/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Munich
 *
 *  All rights reserved.
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EMPIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
 */
#ifdef USE_INTEL_MKL
#include <mkl.h>
#include <mkl_lapacke.h>
#include <mkl_spblas.h>
#endif

#ifndef USE_INTEL_MKL
#include "lapacke.h"
#endif

#ifdef FLANN
#include "flann/flann.hpp"
#endif

#ifdef ANN
#include "ANN/ANN.h"
#endif

#include "MortarMapper.h"
#include "Message.h"
//#include "AuxiliaryFunctions.h"
#include <iostream>
#include <stdlib.h>
#include <set>
#include <assert.h>
#include <omp.h>
#include "Message.h"

using namespace std;
namespace EMPIRE {

// set default value for MKL threads
int MortarMapper::mklSetNumThreads = 1;
// set default value for mapper threads
int MortarMapper::mapperSetNumThreads = 1;

const int MortarMapper::numGPsMassMatrixTri = 6;
const int MortarMapper::numGPsMassMatrixQuad = 4;
const int MortarMapper::numGPsOnClipTri = 6;
const int MortarMapper::numGPsOnClipQuad = 12;

MortarMapper::MortarMapper(int _slaveNumNodes, int _slaveNumElems, const int *_slaveNodesPerElem,
        const double *_slaveNodeCoors, const int *_slaveNodeNumbers, const int *_slaveElemTable,
        int _masterNumNodes, int _masterNumElems, const int *_masterNodesPerElem,
        const double *_masterNodeCoors, const int *_masterNodeNumbers, const int *_masterElemTable,
        bool _oppositeSurfaceNormal, bool _dual, bool _toEnforceConsistency) :
        slaveNumNodes(_slaveNumNodes), slaveNumElems(_slaveNumElems), slaveNodesPerElem(
                _slaveNodesPerElem), slaveNodeCoors(_slaveNodeCoors), slaveNodeNumbers(
                _slaveNodeNumbers), slaveElemTable(_slaveElemTable), masterNumNodes(
                _masterNumNodes), masterNumElems(_masterNumElems), masterNodesPerElem(
                _masterNodesPerElem), masterNodeCoors(_masterNodeCoors), masterNodeNumbers(
                _masterNodeNumbers), masterElemTable(_masterElemTable), oppositeSurfaceNormal(
                _oppositeSurfaceNormal), dual(_dual), toEnforceConsistency(_toEnforceConsistency) {
    // a could-be NULL pointer must be initialized to NULL, otherwise there could be segmentation fault when delete it
    C_BB_A = NULL;
    C_BB_IA = NULL;
    C_BB_JA = NULL;
    C_BB_A_DUAL = NULL;
    C_BA_A = NULL;
    C_BA_IA = NULL;
    C_BA_JA = NULL;
    C_BA_A_DUAL = NULL;
    // check whether the necessary libraries are there
#ifndef USE_INTEL_MKL
    if (!dual) {
        cerr << endl;
        cerr
                << "MortarMapper::MortarMapper: No pardiso library is found, standard mortar mapper cannot be used!"
                << endl << "\t Try dual mortar mapper!" << endl;
        exit(EXIT_FAILURE);
    }
#endif

    // 1. initialize data that could be used later
    initTables();
    initANNTree();

    // 2. compute C_BB
    computeC_BB();
    /*if (!dual) {
        MortarMath::printCSRMatrixSymmetric(C_BB_A, C_BB_IA, C_BB_JA, masterNumNodes);
    } else {
        MortarMath::printDiagonalMatrix(C_BB_A_DUAL, masterNumNodes);
    }*/
    // 3. compute C_BA
    computeC_BA();
    /*if (!dual) {
        MortarMath::printCSRMatrixUnsymmetric(C_BA_A, C_BA_IA, C_BA_JA, masterNumNodes, slaveNumNodes);
    } else {
        MortarMath::printCSRMatrixUnsymmetric(C_BA_A_DUAL, C_BA_IA, C_BA_JA, masterNumNodes, slaveNumNodes);
    }*/
    deleteANNTree();
    deleteTables();

    // 4. use pardiso to do factorization on C_BB
    if (!dual)
        initPardiso();

    // 5. check NULL pointers
    checkNullPointers();
}

MortarMapper::~MortarMapper() {
    if (!dual)
        deletePardiso();
    delete[] C_BB_A;
    delete[] C_BB_IA;
    delete[] C_BB_JA;
    delete[] C_BB_A_DUAL;
    delete[] C_BA_A;
    delete[] C_BA_IA;
    delete[] C_BA_JA;
    delete[] C_BA_A_DUAL;
#ifdef ANN
    for (int i = 0; i < slaveNumNodes; i++) {
        delete[] ANNSlaveNodes[i];
    }
    delete[] ANNSlaveNodes;
#endif
}

void MortarMapper::consistentMapping(const double *slaveField, double *masterField) {
    double *slaveFieldCopy = new double[slaveNumNodes];
    for (int i = 0; i < slaveNumNodes; i++)
        slaveFieldCopy[i] = slaveField[i];

    // 1. matrix vector product (W_tmp = C_BA * W_A)
    int m = masterNumNodes; // number of rows of C_BA
    int n = slaveNumNodes; // number of columns of C_BA
    char noTrans = 'N';
    char descra[] = "G00F"; // general matrix, indexing from 1
    double alpha = 1.0;
    double beta = 0.0;
#ifdef USE_INTEL_MKL
    if (!dual) {
        mkl_set_num_threads(mklSetNumThreads);
        mkl_dcsrmv(&noTrans, &m, &n, &alpha, descra, C_BA_A, C_BA_JA, C_BA_IA, &C_BA_IA[1], slaveFieldCopy,
                &beta, masterField);
    }
    else {
        mkl_set_num_threads(mklSetNumThreads);
        mkl_dcsrmv(&noTrans, &m, &n, &alpha, descra, C_BA_A_DUAL, C_BA_JA, C_BA_IA, &C_BA_IA[1], slaveFieldCopy,
                &beta, masterField);
    }
#else
    if (!dual)
        MortarMath::dcsrmv(noTrans, m, n, C_BA_A, C_BA_JA, C_BA_IA, slaveFieldCopy, masterField);
    else
        MortarMath::dcsrmv(noTrans, m, n, C_BA_A_DUAL, C_BA_JA, C_BA_IA, slaveFieldCopy,
                masterField);
#endif

    delete[] slaveFieldCopy;

    // 2. solve C_BB * W_B = W_tmp
    if (!dual) {
        // pardiso forward and backward substitution
        int phase = 33; // forward and backward substitution
        int idum; // integer dummy
        double *ddum = new double[masterNumNodes]; // dummy but the memory is asked for
        int error = 0;
        iparm[5] = 1; // write solution to b
#ifdef USE_INTEL_MKL
                mkl_set_num_threads(mklSetNumThreads);
                pardiso(pt, &maxfct, &mnum, &mtype, &phase, &neq, C_BB_A, C_BB_IA, C_BB_JA, &idum, &nrhs,
                        iparm, &msglvl, masterField, ddum, &error);
#endif
        if (error != 0) {
            cerr << "Error in MortarMapper: pardiso solver failed!" << error << endl;
            exit(EXIT_FAILURE);
        }
        delete[] ddum;
    } else {
        for (int i = 0; i < masterNumNodes; i++)
            masterField[i] /= C_BB_A_DUAL[i];
    }
}

void MortarMapper::conservativeMapping(const double *masterField, double *slaveField) {
    /*
     * consistent mapping:
     * C_BB * W_B = C_BA * W_A
     * since conservative mapping -- W_B^T * F_B = W_A^T * F_A
     * => F_A = C_BA^T * C_BB^(-1) * F_B
     * F_A --- slaveField
     * F_B --- masterField
     */
    double *masterFieldCopy = new double[masterNumNodes];
    for (int i = 0; i < masterNumNodes; i++)
        masterFieldCopy[i] = masterField[i];
    // 1. solve C_BB * F_tmp = F_B
    if (!dual) {
        // pardiso forward and backward substitution
        int phase = 33; // forward and backward substitution
        int idum; // integer dummy
        double *ddum = new double[masterNumNodes]; // dummy but the memory is asked for
        int error = 0;
        iparm[5] = 1; // write solution to x
#ifdef USE_INTEL_MKL
                mkl_set_num_threads(mklSetNumThreads);
                pardiso(pt, &maxfct, &mnum, &mtype, &phase, &neq, C_BB_A, C_BB_IA, C_BB_JA, &idum, &nrhs,
                        iparm, &msglvl, masterFieldCopy, ddum, &error);
#endif
        if (error != 0) {
            cerr << "Error in MortarMapper: pardiso solver failed!" << error << endl;
            exit(EXIT_FAILURE);
        }
        delete ddum;
    } else {
        for (int i = 0; i < masterNumNodes; i++)
            masterFieldCopy[i] /= C_BB_A_DUAL[i];
    }

    // 2. matrix vector product (F_A = C_BA^T * F_tmp)
    int m = masterNumNodes; // number of rows of C_BA
    int n = slaveNumNodes; // number of columns of C_BA
    char trans = 'T';
    char descra[] = "G00F"; // general matrix, indexing from 1
    double alpha = 1.0;
    double beta = 0.0;
#ifdef USE_INTEL_MKL
    if (!dual) {
        mkl_set_num_threads(mklSetNumThreads);
        mkl_dcsrmv(&trans, &m, &n, &alpha, descra, C_BA_A, C_BA_JA, C_BA_IA, &C_BA_IA[1], masterFieldCopy,
                &beta, slaveField);
    }
    else {
        mkl_set_num_threads(mklSetNumThreads);
        mkl_dcsrmv(&trans, &m, &n, &alpha, descra, C_BA_A_DUAL, C_BA_JA, C_BA_IA, &C_BA_IA[1], masterFieldCopy,
                &beta, slaveField);
    }
#else
    if (!dual)
        MortarMath::dcsrmv(trans, m, n, C_BA_A, C_BA_JA, C_BA_IA, masterFieldCopy, slaveField);
    else
        MortarMath::dcsrmv(trans, m, n, C_BA_A_DUAL, C_BA_JA, C_BA_IA, masterFieldCopy, slaveField);
#endif
    delete[] masterFieldCopy;
}

void MortarMapper::computeC_BB() {
    // 1. compute the sparsity map
    // sparsity map has the information of a, ia, ja in a CSR formated matrix
    map<int, double> **sparsityMapC_BB = NULL; // a "could be" NULL pointer is better to be init to NULL
    if (!dual) {
        sparsityMapC_BB = new map<int, double>*[masterNumNodes];
        for (int i = 0; i < masterNumNodes; i++)
            sparsityMapC_BB[i] = new map<int, double>;
    } else {
        C_BB_A_DUAL = new double[masterNumNodes];
        for (int i = 0; i < masterNumNodes; i++)
            C_BB_A_DUAL[i] = 0.0;
    }

    for (int i = 0; i < masterNumElems; i++) {
        const int numNodesMasterElem = masterNodesPerElem[i];
        double elem[numNodesMasterElem * 3]; // this element
        int pos[numNodesMasterElem]; // the position of the node in the nodeCoors
        for (int j = 0; j < numNodesMasterElem; j++) {
            pos[j] = masterDirectElemTable[i]->at(j);
            for (int k = 0; k < 3; k++)
                elem[j * 3 + k] = masterNodeCoors[pos[j] * 3 + k];
        }
        if (numNodesMasterElem == 4) { // replace the master element by the projection of it on its "element plane"
            double masterElemNormal[3];
            MortarMath::computeNormalOfQuad(elem, true, masterElemNormal);
            double masterQuadCenter[3];
            MortarMath::computePolygonCenter(elem, 4, masterQuadCenter);
            double masterQuadPrj[12];
            MortarMath::projectToPlane(masterQuadCenter, masterElemNormal, elem, 4, masterQuadPrj);
            for (int i = 0; i < 12; i++)
                elem[i] = masterQuadPrj[i];
        }

        // make use of the symmetry
        double massMatrix[numNodesMasterElem * numNodesMasterElem];
        if (numNodesMasterElem == 4)
            MortarMath::computeMassMatrixOfQuad(elem, numGPsMassMatrixQuad, dual, massMatrix);
        else if (numNodesMasterElem == 3)
            MortarMath::computeMassMatrixOfTrianlge(elem, numGPsMassMatrixTri, dual, massMatrix);
        else
            assert(false);
        if (!dual) {
            for (int j = 0; j < numNodesMasterElem; j++) {
                for (int k = j; k < numNodesMasterElem; k++) {
                    int smaller, larger;
                    if (pos[j] > pos[k]) {
                        larger = pos[j];
                        smaller = pos[k];
                    } else {
                        larger = pos[k];
                        smaller = pos[j];
                    }
                    double massMatrixJK = massMatrix[j * numNodesMasterElem + k];
                    bool inserted = sparsityMapC_BB[smaller]->insert(
                            pair<int, double>(larger, massMatrixJK)).second; // *.second is a bool indicating inserted or not
                    if (!inserted)
                        (sparsityMapC_BB[smaller]->at(larger)) += massMatrixJK; // at() returns a reference, so using "+=" is correct
                }
            }
        } else {
            for (int j = 0; j < numNodesMasterElem; j++) {
                C_BB_A_DUAL[pos[j]] += massMatrix[j * numNodesMasterElem + j];
            }
        }
    }

    // 2. according to sparsity map, compute a, ia, ja of CSR formated C_BB
    if (!dual) {
        C_BB_IA = new int[masterNumNodes + 1];
        C_BB_IA[0] = 1; // numbering from 1
        int nnz = 1; // number of non-zero entries
        for (int i = 0; i < masterNumNodes; i++) {
            nnz += sparsityMapC_BB[i]->size();
            C_BB_IA[i + 1] = nnz;
        }
        nnz--;
        C_BB_JA = new int[nnz];
        C_BB_A = new double[nnz];
        int count = 0;
        for (int i = 0; i < masterNumNodes; i++) {
            for (map<int, double>::iterator it = sparsityMapC_BB[i]->begin();
                    it != sparsityMapC_BB[i]->end(); it++) {
                C_BB_JA[count] = it->first + 1; // numbering from 1
                C_BB_A[count] = it->second;
                count++;
            }
        }assert(nnz == count);

        // delete spasity map
        for (int i = 0; i < masterNumNodes; i++)
            delete sparsityMapC_BB[i];
        delete[] sparsityMapC_BB;
    }
}

void MortarMapper::computeC_BA() {
    // 1. initialize sparsity map, which has the information of a, ia, ja in a CSR formated matrix
    map<int, double> **sparsityMapC_BA = NULL;
    map<int, double> **sparsityMapC_BA_DUAL = NULL;
    if (!dual) {
        sparsityMapC_BA = new map<int, double>*[masterNumNodes];
        for (int i = 0; i < masterNumNodes; i++)
            sparsityMapC_BA[i] = new map<int, double>;
    } else {
        sparsityMapC_BA_DUAL = new map<int, double>*[masterNumNodes];
        for (int i = 0; i < masterNumNodes; i++)
            sparsityMapC_BA_DUAL[i] = new map<int, double>;
    }

    // 2. compute entries in the sparsity map by looping over the master elements
#ifdef FLANN
#pragma omp parallel num_threads(mapperSetNumThreads)
#endif
    {
#ifdef FLANN
        //EMPIRE::AuxiliaryFunctions::report_num_threads(2);
#pragma omp for
#endif
        for (int i = 0; i < masterNumElems; i++) {
            // 2.1 compute the searching radius
            int numNodesMasterElem = masterNodesPerElem[i];
            double *masterElem = new double[numNodesMasterElem * 3];
            getElemCoor(i, MortarMapper::MASTER, masterElem);

            double radiusSqr = computeSearchRadiusSquare(masterElem, numNodesMasterElem);
            // if (!(i%1000)) cout << endl << i << "RADIUS: " << radius << endl;

            // 2.2 find all candidates which may overlap the master element
            double masterElemNormal[3];
            if (numNodesMasterElem == 3) {
                MortarMath::computeNormalOfTriangle(masterElem, true, masterElemNormal);
            } else {
                MortarMath::computeNormalOfQuad(masterElem, true, masterElemNormal);
            }
            set<int> *neighborElems = new set<int>;
            findCandidates(masterElem, numNodesMasterElem, masterElemNormal, radiusSqr,
                    neighborElems);

            map<int, double*> *projections = new map<int, double*>; // the projections of all neighboring nodes
            if (numNodesMasterElem == 3)
                projectToElemPlane(masterElem, masterElemNormal, neighborElems, projections);
            if (numNodesMasterElem == 4) { // replace the master element by the projection of it on its "element plane"
                double masterQuadCenter[3];
                MortarMath::computePolygonCenter(masterElem, numNodesMasterElem, masterQuadCenter);
                double masterQuadPrj[12];
                MortarMath::projectToPlane(masterQuadCenter, masterElemNormal, masterElem, 4,
                        masterQuadPrj);
                for (int j = 0; j < 12; j++)
                    masterElem[j] = masterQuadPrj[j];
                projectToElemPlane(masterQuadCenter, masterElemNormal, neighborElems, projections);
            }

            // if (!(i%1000)) cout << "candidates: " << neighborElems->size() << endl;

            int posMasterNodes[numNodesMasterElem]; // the position of the nodes in the master element
            for (int j = 0; j < numNodesMasterElem; j++)
                posMasterNodes[j] = masterDirectElemTable[i]->at(j);

            // 2.3 create the point clipper
            int planeToProject = MortarMath::computePlaneToProject(masterElemNormal);
            MortarMath::PolygonClipper *clipper = new MortarMath::PolygonClipper(masterElem,
                    numNodesMasterElem, planeToProject);

            // 2.4 if dual, compute the coefficient matrix here
            double coeffMatrix[numNodesMasterElem * numNodesMasterElem];
            if (dual)
                computeDualCoeffMatrix(masterElem, numNodesMasterElem, coeffMatrix);
            // 2.5 loop over the candidates, do clipping
            for (set<int>::iterator it = neighborElems->begin(); it != neighborElems->end(); it++) {
                int numNodesSlaveElem = slaveNodesPerElem[*it];
                int posSlaveNodes[numNodesSlaveElem]; // the position of the node in the slave triangle
                double slaveElemPrj[numNodesSlaveElem * 3];

                for (int ii = 0; ii < numNodesSlaveElem; ii++) {
                    posSlaveNodes[ii] = slaveDirectElemTable[*it]->at(ii);
                    for (int jj = 0; jj < 3; jj++) {
                        slaveElemPrj[ii * 3 + jj] = projections->at(posSlaveNodes[ii])[jj];
                    }
                }
                double result[numNodesMasterElem * numNodesSlaveElem];

                vector<double*> *clippedPolygon = new vector<double*>;
                bool overlap = clipper->clip(slaveElemPrj, numNodesSlaveElem, clippedPolygon); // clip
                if (overlap)
                    gaussQuadratureOnClip(masterElem, numNodesMasterElem, slaveElemPrj,
                            numNodesSlaveElem, planeToProject, clippedPolygon, result);
                for (int i = 0; i < clippedPolygon->size(); i++)
                    delete[] clippedPolygon->at(i);
                delete clippedPolygon;
                if (overlap) { // only add the result to sparsity map if overlap happens, that is how C_BA is sparse
                    if (!dual) {
                        for (int ii = 0; ii < numNodesMasterElem; ii++) {
                            for (int jj = 0; jj < numNodesSlaveElem; jj++) {
#ifdef FLANN
#pragma omp critical (map1)
#endif
                                {
                                    bool inserted = sparsityMapC_BA[posMasterNodes[ii]]->insert(
                                            pair<int, double>(posSlaveNodes[jj],
                                                    result[ii * numNodesSlaveElem + jj])).second; // *.second is a bool indicating inserted or not
                                    if (!inserted)
                                        (sparsityMapC_BA[posMasterNodes[ii]]->at(posSlaveNodes[jj])) +=
                                                result[ii * numNodesSlaveElem + jj]; // at() returns a reference, hence using "+=" is correct
                                } //omp critical
                            }
                        }
                    } else {
                        double result_dual[numNodesMasterElem * numNodesSlaveElem];
                        for (int ii = 0; ii < numNodesMasterElem * numNodesSlaveElem; ii++)
                            result_dual[ii] = result[ii]; // copy the value in result
                        MortarMath::computeMatrixProduct(numNodesMasterElem, numNodesSlaveElem,
                                coeffMatrix, result_dual); // now it is dual
                        for (int ii = 0; ii < numNodesMasterElem; ii++) {
                            for (int jj = 0; jj < numNodesSlaveElem; jj++) {
#ifdef FLANN
#pragma omp critical (map2)
#endif
                                {
                                    bool inserted =
                                            sparsityMapC_BA_DUAL[posMasterNodes[ii]]->insert(
                                                    pair<int, double>(
                                                            posSlaveNodes[jj],
                                                            result_dual[ii * numNodesSlaveElem + jj])).second; // *.second is a bool indicating inserted or not
                                    if (!inserted)
                                        (sparsityMapC_BA_DUAL[posMasterNodes[ii]]->at(
                                                posSlaveNodes[jj])) += result_dual[ii
                                                * numNodesSlaveElem + jj]; // at() returns a reference, hence using "+=" is correct
                                } //omp critical
                            }
                        }
                    }
                }
            }
            //delete
            for (map<int, double*>::iterator it = projections->begin(); it != projections->end();
                    it++)
                delete[] it->second;

            delete[] masterElem;
            delete projections;
            delete neighborElems;
            delete clipper;
        }

    } //#pragma omp parallel

    // 3. modify C_BA to enforce consistency
    if (toEnforceConsistency) {
        if (!dual)
            enforceConsistency(sparsityMapC_BA);
        else
            enforceConsistency(sparsityMapC_BA_DUAL);
    }

    // 4. according to sparsity map, compute a, ia, ja of CSR formated C_BA
    if (!dual) {
        C_BA_IA = new int[masterNumNodes + 1];
        C_BA_IA[0] = 1; // numbering from 1
        int nnz = 1; // number of non-zero entries
        for (int i = 0; i < masterNumNodes; i++) {
            nnz += sparsityMapC_BA[i]->size();
            C_BA_IA[i + 1] = nnz;
        }
        nnz--;
        C_BA_JA = new int[nnz];
        C_BA_A = new double[nnz];
        int count = 0;
        for (int i = 0; i < masterNumNodes; i++) {
            for (map<int, double>::iterator it = sparsityMapC_BA[i]->begin();
                    it != sparsityMapC_BA[i]->end(); it++) {
                C_BA_JA[count] = it->first + 1; // numbering from 1
                C_BA_A[count] = it->second;
                count++;
            }
        }
        assert(nnz == count);
    } else {
        C_BA_IA = new int[masterNumNodes + 1];
        C_BA_IA[0] = 1; // numbering from 1
        int nnz = 1; // number of non-zero entries
        for (int i = 0; i < masterNumNodes; i++) {
            nnz += sparsityMapC_BA_DUAL[i]->size();
            C_BA_IA[i + 1] = nnz;
        }
        nnz--;
        C_BA_JA = new int[nnz];
        C_BA_A_DUAL = new double[nnz];
        int count = 0;
        for (int i = 0; i < masterNumNodes; i++) {
            for (map<int, double>::iterator it = sparsityMapC_BA_DUAL[i]->begin();
                    it != sparsityMapC_BA_DUAL[i]->end(); it++) {
                C_BA_JA[count] = it->first + 1; // numbering from 1
                C_BA_A_DUAL[count] = it->second;
                count++;
            }
        }
        assert(nnz == count);
    }

    // delete
    if (!dual) {
        for (int i = 0; i < masterNumNodes; i++)
            delete sparsityMapC_BA[i];
        delete[] sparsityMapC_BA;
    } else {
        for (int i = 0; i < masterNumNodes; i++)
            delete sparsityMapC_BA_DUAL[i];
        delete[] sparsityMapC_BA_DUAL;
    }
}

void MortarMapper::enforceConsistency(map<int, double> **sparsityMapC_BA) {
    // solve the problem that C_BB * 1 != C_BA * 1, this happens due to the projection of A not covering B
    if (!dual) {
        double *factor = new double[masterNumNodes];
        for (int i = 0; i < masterNumNodes; i++) {
            factor[i] = 0.0;
        }
        for (int i = 0; i < masterNumNodes; i++) {
            for (int j = C_BB_IA[i]; j < C_BB_IA[i + 1]; j++) {
                factor[i] += C_BB_A[j - 1];
            }
        }
        for (int i = 0; i < masterNumNodes; i++) {
            for (int j = C_BB_IA[i] + 1; j < C_BB_IA[i + 1]; j++) {
                factor[C_BB_JA[j - 1] - 1] += C_BB_A[j - 1];
            }
        }
        for (int i = 0; i < masterNumNodes; i++) {
            double sum = 0.0;
            for (map<int, double>::iterator it = sparsityMapC_BA[i]->begin();
                    it != sparsityMapC_BA[i]->end(); it++) {
                sum += it->second;
            }

            if (sum < factor[i] * 0.5) { // if the master element is not fully covered by slave elements, use nearest neighbor
                cout << "WARNING(MortarMapper::enforceConsistency): Nearest neighbor is used for node: ";
                MortarMath::printPoint(&masterNodeCoors[i*3]);
                sparsityMapC_BA[i]->clear();
                double dummy;
                int nb;
#ifdef ANN
                slaveNodesTree->annkSearch(&masterNodeCoors[i * 3], 1, &nb, &dummy);
#endif
#ifdef FLANN
                flann::Matrix<double> masterNodeFlann(const_cast<double*>(&masterNodeCoors[i * 3]), 1, 3);
                vector<vector<int> > indexes_tmp;
                vector<vector<double> > dists_tmp;
                FLANNkd_tree->knnSearch(masterNodeFlann, indexes_tmp, dists_tmp, 1, flann::SearchParams(1));
                nb = indexes_tmp[0][0];
#endif
                sparsityMapC_BA[i]->insert(pair<int, double>(nb, factor[i]));
                sum = factor[i];
            }
            // rowsum of C_BA cannot be large than rowsum of C_BB
            assert(sum<factor[i]*(1+1E-2));
            factor[i] /= sum;
        }
        for (int i = 0; i < masterNumNodes; i++) {
            for (map<int, double>::iterator it = sparsityMapC_BA[i]->begin();
                    it != sparsityMapC_BA[i]->end(); it++) {
                it->second *= factor[i];
            }
        }
        delete[] factor;
    } else {
        for (int i = 0; i < masterNumNodes; i++) {
            double factor = C_BB_A_DUAL[i];
            double sum = 0.0;
            for (map<int, double>::iterator it = sparsityMapC_BA[i]->begin();
                    it != sparsityMapC_BA[i]->end(); it++) {
                sum += it->second;
            }

            if ((sum < factor * 0.5) || (sum > factor * 1.5)) { // if the master element is not fully covered by slave elements, use nearest neighbor
                sparsityMapC_BA[i]->clear();
                double dummy;
                int nb;
#ifdef ANN
                slaveNodesTree->annkSearch(&masterNodeCoors[i * 3], 1, &nb, &dummy);
#endif
#ifdef FLANN
                flann::Matrix<double> masterNodeFlann(const_cast<double*>(&masterNodeCoors[i * 3]), 1, 3);
                vector<vector<int> > indexes_tmp;
                vector<vector<double> > dists_tmp;
                FLANNkd_tree->knnSearch(masterNodeFlann, indexes_tmp, dists_tmp, 1, flann::SearchParams(1));
                nb = indexes_tmp[0][0];
#endif
                sparsityMapC_BA[i]->insert(pair<int, double>(nb, factor));
                sum = factor;
            }
            // rowsum of C_BA cannot be large than rowsum of C_BB
            //assert(sum<factor*(1+1E-2)); // This is not the case for dual
            factor /= sum;
            for (map<int, double>::iterator it = sparsityMapC_BA[i]->begin();
                    it != sparsityMapC_BA[i]->end(); it++) {
                it->second *= factor;
            }
        }
    }
}

void MortarMapper::initPardiso() {
#ifdef USE_INTEL_MKL
    mtype = 2; // real symmetric matrix
    // set pardiso default parameters
    pardisoinit(pt, &mtype, iparm);
    iparm[2] = mklSetNumThreads;
    //cout << endl << "iparm[3]: " << iparm[2] << endl;
    maxfct = 1;// max number of factorizations
    mnum = 1;// which factorization to use
    msglvl = 0;// do NOT print statistical information
    neq = masterNumNodes;// number of rows of C_BB
    nrhs = 1;// number of right hand side
    int phase = 12;// analysis and factorization
    double ddum;// double dummy
    int idum;// integer dummy
    int error = 0;
    //cout<<"factorizing"<<endl;
    mkl_set_num_threads(mklSetNumThreads);
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &neq, C_BB_A, C_BB_IA, C_BB_JA, &idum, &nrhs, iparm,
            &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        cerr << "Error in MortarMapper: pardiso factorization failed!" << error << endl;
        exit(EXIT_FAILURE);
    }
#endif
}

void MortarMapper::deletePardiso() {
#ifdef USE_INTEL_MKL
    int phase = -1; // deallocate memory
    double ddum;// double dummy
    int idum;// integer dummy
    int error = 0;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &neq, C_BB_A, C_BB_IA, C_BB_JA, &idum, &nrhs, iparm,
            &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        cerr << "Error in MortarMapper: pardiso factorization failed!" << error << endl;
        exit(EXIT_FAILURE);
    }
#endif
}

void MortarMapper::initTables() {
    // using the map to store the nodeNumbers
    // but here the "key" is the node number, and the value is the position in nodeNumbers
    // the map is sorted automatically, so it is efficient for searching

    { // 1. compute slaveDirectElemTable
        slaveDirectElemTable = new vector<int>*[slaveNumElems];
        for (int i = 0; i < slaveNumElems; i++)
            slaveDirectElemTable[i] = new vector<int>;
        map<int, int> *slaveNodesMap = new map<int, int>();
        for (int i = 0; i < slaveNumNodes; i++)
            slaveNodesMap->insert(slaveNodesMap->end(), pair<int, int>(slaveNodeNumbers[i], i));
        int count = 0;
        for (int i = 0; i < slaveNumElems; i++) {
            const int numNodesSlaveElem = slaveNodesPerElem[i];   
            for (int j = 0; j < numNodesSlaveElem; j++) {
                if(slaveNodesMap->find(slaveElemTable[count + j])== slaveNodesMap->end() ){
                ERROR_OUT()<< "Slave Node Label " << slaveElemTable[count + j] << " is not part of slaveNodesMap." << endl;
                }
                slaveDirectElemTable[i]->push_back(slaveNodesMap->at(slaveElemTable[count + j]));
            }
            count += numNodesSlaveElem;
        }
        delete slaveNodesMap;
    }

    { // 2. compute masterDirectElemTable
        masterDirectElemTable = new vector<int>*[masterNumElems];
        for (int i = 0; i < masterNumElems; i++)
            masterDirectElemTable[i] = new vector<int>;
        map<int, int> *masterNodesMap = new map<int, int>();

        for (int i = 0; i < masterNumNodes; i++)
            masterNodesMap->insert(masterNodesMap->end(), pair<int, int>(masterNodeNumbers[i], i));

        int count = 0;
        for (int i = 0; i < masterNumElems; i++) {
            const int numNodesMasterElem = masterNodesPerElem[i];
            for (int j = 0; j < numNodesMasterElem; j++) {
                masterDirectElemTable[i]->push_back(masterNodesMap->at(masterElemTable[count + j]));
            }
            count += numNodesMasterElem;
        }
        delete masterNodesMap;
    }

    { // 3. compute slaveNodeToElemTable
        slaveNodeToElemTable = new vector<int>*[slaveNumNodes];
        for (int i = 0; i < slaveNumNodes; i++)
            slaveNodeToElemTable[i] = new vector<int>;
        for (int i = 0; i < slaveNumElems; i++) {
            const int numNodesSlaveElem = slaveNodesPerElem[i];
            for (int j = 0; j < numNodesSlaveElem; j++) {
                int nodePos = slaveDirectElemTable[i]->at(j);
                slaveNodeToElemTable[nodePos]->push_back(i);
            }
        }
    }

    // 4. computeSlaveElemNormals
    computeSlaveElemNormals();
}

void MortarMapper::deleteTables() {
    for (int i = 0; i < slaveNumElems; i++)
        delete slaveDirectElemTable[i];
    delete[] slaveDirectElemTable;

    for (int i = 0; i < masterNumElems; i++)
        delete masterDirectElemTable[i];
    delete[] masterDirectElemTable;

    for (int i = 0; i < slaveNumNodes; i++)
        delete slaveNodeToElemTable[i];
    delete[] slaveNodeToElemTable;

    delete[] slaveElemNormals;
}

void MortarMapper::initANNTree() {
#ifdef ANN
    ANNSlaveNodes = new double*[slaveNumNodes]; // ANN uses 2D array
    for (int i = 0; i < slaveNumNodes; i++) {
        ANNSlaveNodes[i] = new double[3];
        for (int j = 0; j<3; j++)
        ANNSlaveNodes[i][j] = slaveNodeCoors[i * 3 + j];
    }
    slaveNodesTree = new ANNkd_tree(ANNSlaveNodes, slaveNumNodes, 3);
#endif
#ifdef FLANN
    FLANNSlaveNodes = new flann::Matrix<double>(const_cast<double*>(slaveNodeCoors), slaveNumNodes, 3);
    FLANNkd_tree = new flann::Index<flann::L2<double> >(*FLANNSlaveNodes, flann::KDTreeSingleIndexParams(1));
    FLANNkd_tree->buildIndex(); // Build binary tree for searching
#endif
}

void MortarMapper::deleteANNTree() {
#ifdef ANN
    delete[] ANNSlaveNodes;
    delete slaveNodesTree;
    annClose();
#endif
#ifdef FLANN
    delete FLANNSlaveNodes;
    delete FLANNkd_tree;
#endif
}

void MortarMapper::gaussQuadratureOnClip(const double *masterElem, int numNodesMasterElem,
        const double *slaveElem, int numNodesSlaveElem, int planeToProject,
        vector<double*> *clippedPolygon, double *result) {
    for (int i = 0; i < numNodesMasterElem * numNodesSlaveElem; i++)
        result[i] = 0.0;

    int numGPs = 0;
    if (numNodesMasterElem == 3 && numNodesSlaveElem == 3)
        numGPs = numGPsOnClipTri;
    else
        numGPs = numGPsOnClipQuad;
    // if the clipped polygon is a triangle, do Gauss quadrature on it directly
    if (clippedPolygon->size() == 3) {
        double clipTriangle[9];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                clipTriangle[i * 3 + j] = clippedPolygon->at(i)[j];

        MortarMath::GaussQuadratureOnTriangle *gaussQuadratureOnTriangle =
                new MortarMath::GaussQuadratureOnTriangle(clipTriangle, numGPs);
        ShapeFunctionProduct *integrand = new ShapeFunctionProduct(masterElem, numNodesMasterElem,
                slaveElem, numNodesSlaveElem, planeToProject);
        integrand->setGaussPoints(gaussQuadratureOnTriangle->gaussPointsGlobal,
                gaussQuadratureOnTriangle->numGaussPoints);
        integrand->computeShapeFunctionProducts();

        for (int i = 0; i < numNodesMasterElem; i++) {
            for (int j = 0; j < numNodesSlaveElem; j++) {
                integrand->setShapeFunctions(i, j);
                gaussQuadratureOnTriangle->setIntegrandFunc(integrand);
                result[i * numNodesSlaveElem + j] = gaussQuadratureOnTriangle->computeIntegral();
            }
        }
        delete gaussQuadratureOnTriangle;
        delete integrand;
    } else { // if not, divide it into triangles, on each of which doing Gauss quadrature
        int size = clippedPolygon->size();
        double tmp[size * 3];
        for (int i = 0; i < size; i++)
            for (int j = 0; j < 3; j++)
                tmp[i * 3 + j] = clippedPolygon->at(i)[j];
        double center[3];
        MortarMath::computePolygonCenter(tmp, size, center);

        for (int i = 0; i < size; i++) {
            double clipTriangle[9];
            MortarMath::buildTrianagle(center, clippedPolygon->at(i),
                    clippedPolygon->at((i + 1) % size), clipTriangle);
            MortarMath::GaussQuadratureOnTriangle *gaussQuadratureOnTriangle =
                    new MortarMath::GaussQuadratureOnTriangle(clipTriangle, numGPs);
            ShapeFunctionProduct *integrand = new ShapeFunctionProduct(masterElem,
                    numNodesMasterElem, slaveElem, numNodesSlaveElem, planeToProject);
            integrand->setGaussPoints(gaussQuadratureOnTriangle->gaussPointsGlobal,
                    gaussQuadratureOnTriangle->numGaussPoints);
            integrand->computeShapeFunctionProducts();
            for (int i = 0; i < numNodesMasterElem; i++) {
                for (int j = 0; j < numNodesSlaveElem; j++) {
                    integrand->setShapeFunctions(i, j);
                    gaussQuadratureOnTriangle->setIntegrandFunc(integrand);
                    result[i * numNodesSlaveElem + j] +=
                            gaussQuadratureOnTriangle->computeIntegral();
                }
            }
            delete gaussQuadratureOnTriangle;
            delete integrand;
        }
    }
}

double MortarMapper::computeSearchRadiusSquare(const double* masterElem, int numNodesMasterElem) {
    // 1. find the longest edge of this element (masterElem)
    double masterElemCopy[3 * numNodesMasterElem];
    MortarMath::copyElem(masterElem, numNodesMasterElem, masterElemCopy);
    double lengthSqr = MortarMath::longestEdgeLengthSquare(masterElemCopy, numNodesMasterElem);

    // 2. find the longest edge of the neighboring elements (the slave side)
    // the neighboring elements are found by
    // first find all neighboring nodes of the nodes of the master element,
    // then all elements on the slave side having these nodes are the neighboring elements.
    // These elements together with the master element itself help to compute a reasonable search radius.
    double dummy;
    int nb;
    vector<int>* elemVec[numNodesMasterElem]; // store nearest elements of each node of the masterElem
    for (int i = 0; i < numNodesMasterElem; i++) {
#ifdef ANN
        slaveNodesTree->annkSearch(&masterElemCopy[i * 3], 1, &nb, &dummy);
#endif
#ifdef FLANN
        flann::Matrix<double> masterElemCopyFlann(&masterElemCopy[i*3], 1, 3);
        vector<vector<int> > indexes_tmp;
        vector<vector<double> > dists_tmp;
        FLANNkd_tree->knnSearch(masterElemCopyFlann, indexes_tmp, dists_tmp, 1, flann::SearchParams(1));
        nb = indexes_tmp[0][0];
#endif

        elemVec[i] = slaveNodeToElemTable[nb]; // do not modify/delete contents in elemVec!!!
    }

    double searchRadiusSqr = lengthSqr; // the goal of the first two steps is to set up this value
    set<int> neighborElemsTmp; // set allows no multiple entries
    for (int i = 0; i < numNodesMasterElem; i++)
        for (unsigned j = 0; j < elemVec[i]->size(); j++)
            neighborElemsTmp.insert(elemVec[i]->at(j));

    for (set<int>::iterator it = neighborElemsTmp.begin(); it != neighborElemsTmp.end(); it++) {
        int numNodesSlaveElem = slaveNodesPerElem[*it];
        double elemTmp[3 * numNodesSlaveElem];
        // for a single element, find its longest edge
        getElemCoor(*it, MortarMapper::SLAVE, elemTmp);
        lengthSqr = MortarMath::longestEdgeLengthSquare(elemTmp, numNodesSlaveElem);
        if (lengthSqr > searchRadiusSqr)
            searchRadiusSqr = lengthSqr;
    }

    // 3. the radius is 1.12 times the longest edge in the first two steps
    searchRadiusSqr *= 1.12;
    searchRadiusSqr *= 1.12;
    return searchRadiusSqr;
}

void MortarMapper::findCandidates(const double* masterElem, int masterElemNumNodes,
        const double* masterElemNormal, double searchRadiusSqr, set<int> *neighborElems) {
    // Use fixed radius search on end points of the element,
    // all elements containing these points are the overlapped candidates
    // 1. find all neighboring elements in the radius
    double masterElemCopy[masterElemNumNodes * 3];
    MortarMath::copyElem(masterElem, masterElemNumNodes, masterElemCopy);
    // OpenMP parallelize this loop
    // Ann is not thread safe
    for (int i = 0; i < masterElemNumNodes; i++) {
        int n_nb = 0;
        int *nbs;
#ifdef ANN
        n_nb = slaveNodesTree->annkFRSearch(&masterElemCopy[i * 3], searchRadiusSqr, 0); // get the number of neighbors in a radius
        nbs = new int[n_nb];
        double *dummy = new double[n_nb];
        slaveNodesTree->annkFRSearch(&masterElemCopy[i * 3], searchRadiusSqr, n_nb, nbs, dummy);// get the real neighbors (ANN uses the square of the radius)
        delete dummy;
#endif

#ifdef FLANN
        flann::Matrix<double> masterElemCopyFlann(&masterElemCopy[i*3], 1, 3);
        vector<vector<int> > indexes_tmp;
        vector<vector<double> > dists_tmp;
        FLANNkd_tree->radiusSearch(masterElemCopyFlann, indexes_tmp, dists_tmp, searchRadiusSqr, flann::SearchParams(1));
        n_nb = indexes_tmp[0].size();
        nbs = new int[indexes_tmp[0].size()];
        for (int j=0; j<indexes_tmp[0].size(); j++)
        nbs[j] = indexes_tmp[0][j];
#endif
        for (int j = 0; j < n_nb; j++) {
            vector<int> * vecTmp = slaveNodeToElemTable[nbs[j]];
            for (unsigned k = 0; k < vecTmp->size(); k++) {
                neighborElems->insert(vecTmp->at(k));
            }
        }
        delete nbs;
    }
    vector<set<int>::iterator> toDelete;
    // 2. kick out the elements that have wrong normal direction
    for (set<int>::iterator it = neighborElems->begin(); it != neighborElems->end(); it++) {
        double *slaveElemNormal = &slaveElemNormals[(*it) * 3];
        if (kickOutCandidate(masterElemNormal, slaveElemNormal, 0.7))
            toDelete.push_back(it); //neighborElems->erase(it) is wrong! because it++ goes to a magic place later
    }
    /*if (toDelete.size() == neighborElems->size()) { // loose the restriction if no neighbors are found
     toDelete.clear();
     for (set<int>::iterator it = neighborElems->begin(); it != neighborElems->end(); it++) {
     double *slaveElemNormal = &slaveElemNormals[(*it) * 3];
     if (kickOutCandidate(masterElemNormal, slaveElemNormal, 0.01))
     toDelete.push_back(it); //neighborElems->erase(it) is wrong! because it++ goes to a magic place later
     }
     }*/
    for (int i = 0; i < toDelete.size(); i++)
        neighborElems->erase(toDelete[i]); // out kicking
    //assert(neighborElems->size()!=0);
}

bool MortarMapper::kickOutCandidate(const double *masterUnitNormal, const double *slaveUnitNormal,
        double bound) {
    //acos(0.7) = pi/4
    //acos(0.01) = pi/2
    if (!oppositeSurfaceNormal)
        return (MortarMath::computeVectorDotProduct(masterUnitNormal, slaveUnitNormal) < bound);
    return (MortarMath::computeVectorDotProduct(masterUnitNormal, slaveUnitNormal) > -bound);
}

void MortarMapper::projectToElemPlane(const double *elem, const double *planeUnitNormal,
        set<int> *neighborElems, map<int, double*> *projections) {
    set<int> neighborNodes;
    for (set<int>::iterator it = neighborElems->begin(); it != neighborElems->end(); it++) {
        for (int i = 0; i < slaveNodesPerElem[*it]; i++)
            neighborNodes.insert(slaveDirectElemTable[*it]->at(i));
    }
    //cout << neighborNodes.size();
    for (set<int>::iterator it = neighborNodes.begin(); it != neighborNodes.end(); it++) {
        double *projectionTmp = new double[3];
        MortarMath::projectToPlane(&elem[0], planeUnitNormal, &slaveNodeCoors[(*it) * 3], 1,
                projectionTmp);
        projections->insert(projections->end(), pair<int, double*>(*it, projectionTmp));
    }
}

void MortarMapper::getElemCoor(int elemIndex, MeshLabel label, double *elem) {
    // compute the coordinates of a triangle by its id
    if (label == MASTER) {
        int numNodesMasterElem = masterNodesPerElem[elemIndex];
        for (int i = 0; i < numNodesMasterElem; i++) {
            int nodePos = masterDirectElemTable[elemIndex]->at(i); // position of the node
            for (int j = 0; j < 3; j++) {
                elem[i * 3 + j] = masterNodeCoors[nodePos * 3 + j];
            }
        }

    } else if (label == SLAVE) {
        int numNodesSlaveElem = slaveNodesPerElem[elemIndex];
        for (int i = 0; i < numNodesSlaveElem; i++) {
            int nodePos = slaveDirectElemTable[elemIndex]->at(i); // position of the node
            for (int j = 0; j < 3; j++) {
                elem[i * 3 + j] = slaveNodeCoors[nodePos * 3 + j];
            }
        }
    } else {
        assert(false);
    }
}

void MortarMapper::computeSlaveElemNormals() {
    slaveElemNormals = new double[slaveNumElems * 3];
    for (int i = 0; i < slaveNumElems; i++) {
        const int numNodesSlaveElem = slaveNodesPerElem[i];
        double slaveElem[numNodesSlaveElem * 3];
        for (int j = 0; j < numNodesSlaveElem; j++) {
            for (int k = 0; k < 3; k++) {
                int pos = slaveDirectElemTable[i]->at(j);
                slaveElem[j * 3 + k] = slaveNodeCoors[pos * 3 + k];
            }
        }
        if (numNodesSlaveElem == 3) {
            MortarMath::computeNormalOfTriangle(slaveElem, true, &slaveElemNormals[i * 3]);
        } else if (numNodesSlaveElem == 4) {
            MortarMath::computeNormalOfQuad(slaveElem, true, &slaveElemNormals[i * 3]);
        } else {
            assert(false);
        }
    }
}

void MortarMapper::checkNullPointers() {
    if (!dual) {
        assert(C_BB_A!=NULL);
        assert(C_BB_IA!=NULL);
        assert(C_BB_JA!=NULL);
        assert(C_BB_A_DUAL==NULL);

        assert(C_BA_A!=NULL);
        assert(C_BA_IA!=NULL);
        assert(C_BA_JA!=NULL);
        assert(C_BA_A_DUAL==NULL);
    } else {
        assert(C_BB_A==NULL);
        assert(C_BB_IA==NULL);
        assert(C_BB_JA==NULL);
        assert(C_BB_A_DUAL!=NULL);

        assert(C_BA_A==NULL);
        assert(C_BA_IA!=NULL);
        assert(C_BA_JA!=NULL);
        assert(C_BA_A_DUAL!=NULL);
    }
}

void MortarMapper::computeDualCoeffMatrix(const double *elem, int numNodesElem,
        double *coeffMatrix) {
    double massMatrix[numNodesElem * numNodesElem];
    if (numNodesElem == 4)
        MortarMath::computeMassMatrixOfQuad(elem, numGPsMassMatrixQuad, false, massMatrix);
    else if (numNodesElem == 3)
        MortarMath::computeMassMatrixOfTrianlge(elem, numGPsMassMatrixTri, false, massMatrix);
    else
        assert(false);

    // store the dual mass matrix in coeffMatrix temporarily
    if (numNodesElem == 4)
        MortarMath::computeMassMatrixOfQuad(elem, numGPsMassMatrixQuad, true, coeffMatrix);
    else if (numNodesElem == 3)
        MortarMath::computeMassMatrixOfTrianlge(elem, numGPsMassMatrixTri, true, coeffMatrix);
    else
        assert(false);

    int dummy[numNodesElem];
    int info = -1;
    info = LAPACKE_dsysv(LAPACK_COL_MAJOR, 'L', numNodesElem, numNodesElem, massMatrix,
            numNodesElem, dummy, coeffMatrix, numNodesElem);
    assert(info == 0);
    // if it fails, then a singular mass matrix is indicated
}

MortarMapper::ShapeFunctionProduct::ShapeFunctionProduct(const double *_masterElem,
        int _numNodesMasterElem, const double *_slaveElem, int _numNodesSlaveElem,
        int _planeToProject) :
        masterElem(_masterElem), numNodesMasterElem(_numNodesMasterElem), slaveElem(_slaveElem), numNodesSlaveElem(
                _numNodesSlaveElem), planeToProject(_planeToProject) {
    gaussPoints = NULL;
    shapeFunctionProducts = NULL;
}

MortarMapper::ShapeFunctionProduct::~ShapeFunctionProduct() {
    assert(gaussPoints!=NULL);
    assert(shapeFunctionProducts!=NULL);
    delete[] gaussPoints;
    for (int i = 0; i < numNodesMasterElem * numNodesSlaveElem; i++)
        delete[] shapeFunctionProducts[i];
    delete[] shapeFunctionProducts;
}

void MortarMapper::ShapeFunctionProduct::computeShapeFunctionProducts() {
    assert(shapeFunctionProducts == NULL);
    assert(gaussPoints != NULL);
    shapeFunctionProducts = new double*[numNodesMasterElem * numNodesSlaveElem];
    for (int i = 0; i < numNodesMasterElem * numNodesSlaveElem; i++)
        shapeFunctionProducts[i] = new double[numGaussPoints];

    for (int i = 0; i < numGaussPoints; i++) {
        double *gaussPoint = &gaussPoints[i * 3];
        // *. compute shape function values on master element (on the Gauss point)
        double shapeFuncValueMasterElem[numNodesMasterElem];
        if (numNodesMasterElem == 3) {
            double localCoor[3];
            bool inside = MortarMath::computeLocalCoorInTriangle(masterElem, planeToProject,
                    gaussPoint, localCoor);
            // debug
            if (!inside){
            	cout<<"Error in computing local coordinates in tria master element"<<endl;
                cout<<"GP coordinates: "<<gaussPoint[0]<<"\t"<<gaussPoint[1]<<"\t"<<gaussPoint[2]<<endl;
                cout<<"Element nodes:"<<endl;
                for(int ctr=0;ctr<numNodesMasterElem;ctr++){
                	cout<<ctr+1<<": "<<masterElem[ctr*3]<<"\t"<<masterElem[ctr*3+1]<<"\t"<<masterElem[ctr*3+2]<<endl;
                }
            }
            // debug end
            assert(inside);
            for (int j = 0; j < 3; j++) {
                shapeFuncValueMasterElem[j] = localCoor[j];
            }
        } else if (numNodesMasterElem == 4) {
            double localCoor[2];
            bool inside = MortarMath::computeLocalCoorInQuad(masterElem, planeToProject, gaussPoint,
                    localCoor);
            // debug
            if (!inside){
            	cout<<"Error in computing local coordinates in quad master element"<<endl;
            	cout<<"GP coordinates: "<<gaussPoint[0]<<"\t"<<gaussPoint[1]<<"\t"<<gaussPoint[2]<<endl;
            	cout<<"Element nodes:"<<endl;
                for(int ctr=0;ctr<numNodesMasterElem;ctr++){
                	cout<<ctr+1<<": "<<masterElem[ctr*3]<<"\t"<<masterElem[ctr*3+1]<<"\t"<<masterElem[ctr*3+2]<<endl;
                }
             }
            // debug end
            assert(inside);
            MortarMath::computeShapeFuncOfQuad(localCoor, shapeFuncValueMasterElem);

        } else {
            assert(false);
        }
        // *. compute shape function values on slave element (on the Gauss point)
        double shapeFuncValueSlaveElem[numNodesSlaveElem];
        if (numNodesSlaveElem == 3) {
            double localCoor[3];
            bool inside = MortarMath::computeLocalCoorInTriangle(slaveElem, planeToProject,
                    gaussPoint, localCoor);
            // debug
            if (!inside){
            	cout<<"Error in computing local coordinates in quad slave element"<<endl;
                cout<<"GP coordinates: "<<gaussPoint[0]<<"\t"<<gaussPoint[1]<<"\t"<<gaussPoint[2]<<endl;
                cout<<"Element nodes:"<<endl;
                for(int ctr=0;ctr<numNodesSlaveElem;ctr++){
                	cout<<ctr+1<<": "<<slaveElem[ctr*3]<<slaveElem[ctr*3+1]<<slaveElem[ctr*3+2]<<endl;
                }
            }
            // debug end
            assert(inside);
            for (int j = 0; j < 3; j++) {
                shapeFuncValueSlaveElem[j] = localCoor[j];
            }
        } else if (numNodesSlaveElem == 4) {
            double localCoor[2];
            bool inside = MortarMath::computeLocalCoorInQuad(slaveElem, planeToProject, gaussPoint,
                    localCoor);
            // debug
            if (!inside){
            	cout<<"Error in computing local coordinates in quad slave element"<<endl;
            	cout<<"GP coordinates: "<<gaussPoint[0]<<"\t"<<gaussPoint[1]<<"\t"<<gaussPoint[2]<<endl;
            	cout<<"Element nodes:"<<endl;
            	for(int ctr=0;ctr<numNodesSlaveElem;ctr++){
            		cout<<ctr+1<<": "<<slaveElem[ctr*3]<<slaveElem[ctr*3+1]<<slaveElem[ctr*3+2]<<endl;
            	}
            }
            // debug end
            assert(inside);

            MortarMath::computeShapeFuncOfQuad(localCoor, shapeFuncValueSlaveElem);
        } else {
            assert(false);
        }
        // *. compute shape function products (on the Gauss point)
        for (int j = 0; j < numNodesMasterElem; j++)
            for (int k = 0; k < numNodesSlaveElem; k++)
                shapeFunctionProducts[j * numNodesSlaveElem + k][i] = shapeFuncValueMasterElem[j]
                        * shapeFuncValueSlaveElem[k];
    }

}

void MortarMapper::ShapeFunctionProduct::setGaussPoints(const double *_gaussPoints,
        int _numGaussPoints) {
    assert(gaussPoints == NULL);
    numGaussPoints = _numGaussPoints;
    gaussPoints = new double[numGaussPoints * 3];
    for (int i = 0; i < numGaussPoints * 3; i++)
        gaussPoints[i] = _gaussPoints[i];
}

void MortarMapper::ShapeFunctionProduct::setShapeFunctions(int _masterShapeFuncID,
        int _slaveShapeFuncID) {
    masterShapeFuncID = _masterShapeFuncID;
    slaveShapeFuncID = _slaveShapeFuncID;
}

double MortarMapper::ShapeFunctionProduct::operator ()(double *gaussPoint) {
    assert(gaussPoints!=NULL);
    assert(shapeFunctionProducts!=NULL);
    for (int i = 0; i < numGaussPoints; i++) {
        if (gaussPoints[i * 3 + 0] == gaussPoint[0] && gaussPoints[i * 3 + 1] == gaussPoint[1]
                && gaussPoints[i * 3 + 2] == gaussPoint[2]) {
            return shapeFunctionProducts[masterShapeFuncID * numNodesSlaveElem + slaveShapeFuncID][i];
        }
    }
    // should not come here
    assert(false);
    return 0.0;
}

} /* namespace EMPIRE */
