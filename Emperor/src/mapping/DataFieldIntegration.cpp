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
#include "DataFieldIntegration.h"
#include "MortarMath.h"
#include <map>
#include <vector>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#ifdef USE_INTEL_MKL
#include <mkl.h>
#include <mkl_lapacke.h>
#include <mkl_spblas.h>
#endif

using namespace std;

namespace EMPIRE {

int DataFieldIntegration::mklSetNumThreads = 1; // set default value for MKL threads
const int DataFieldIntegration::numGPsMassMatrixTri = 6;
const int DataFieldIntegration::numGPsMassMatrixQuad = 4;

DataFieldIntegration::DataFieldIntegration(int _numNodes, int _numElems,
        const int *_numNodesPerElem, const double *_nodes, const int *_nodeIDs, const int *_elems) :
        numNodes(_numNodes), numElems(_numElems), numNodesPerElem(_numNodesPerElem), nodes(_nodes), nodeIDs(
                _nodeIDs), elems(_elems) {
    pardisoInitialized = false;
    // compute directElemTable which link the element to the all its nodes' positions (instead of IDs)
    vector<int> **directElemTable = new vector<int>*[numElems];
    for (int i = 0; i < numElems; i++)
        directElemTable[i] = new vector<int>;
    map<int, int> *nodesMap = new map<int, int>();
    for (int i = 0; i < numNodes; i++)
        nodesMap->insert(nodesMap->end(), pair<int, int>(nodeIDs[i], i));
    int count = 0;
    for (int i = 0; i < numElems; i++) {
        const int numNodesThisElem = numNodesPerElem[i];
        for (int j = 0; j < numNodesThisElem; j++) {
            directElemTable[i]->push_back(nodesMap->at(elems[count + j]));
        }
        count += numNodesThisElem;
    }
    delete nodesMap;

    // compute the sparsity map
    // sparsity map has the information of a, ia, ja in a CSR formated matrix
    map<int, double> **sparsityMap = new map<int, double>*[numNodes];
    for (int i = 0; i < numNodes; i++)
        sparsityMap[i] = new map<int, double>;

    for (int i = 0; i < numElems; i++) {
        const int numNodesThisElem = numNodesPerElem[i];
        double elem[numNodesThisElem * 3]; // this element
        int pos[numNodesThisElem]; // the position of the node in the nodeCoors
        for (int j = 0; j < numNodesThisElem; j++) {
            pos[j] = directElemTable[i]->at(j);
            for (int k = 0; k < 3; k++)
                elem[j * 3 + k] = nodes[pos[j] * 3 + k];
        }

        // replace the master element by the projection of it on its "element plane"
        // we do it here because we have done the same in MortarMapper
        if (numNodesThisElem == 4) {
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
        double massMatrixElem[numNodesThisElem * numNodesThisElem];
        if (numNodesThisElem == 4)
            MortarMath::computeMassMatrixOfQuad(elem, numGPsMassMatrixQuad, false, massMatrixElem);
        else if (numNodesThisElem == 3)
            MortarMath::computeMassMatrixOfTrianlge(elem, numGPsMassMatrixTri, false,
                    massMatrixElem);
        else
            assert(false);
        for (int j = 0; j < numNodesThisElem; j++) {
            for (int k = j; k < numNodesThisElem; k++) {
                int smaller, larger;
                if (pos[j] > pos[k]) {
                    larger = pos[j];
                    smaller = pos[k];
                } else {
                    larger = pos[k];
                    smaller = pos[j];
                }
                double massMatrixJK = massMatrixElem[j * numNodesThisElem + k];
                bool inserted = sparsityMap[smaller]->insert(
                        pair<int, double>(larger, massMatrixJK)).second; // *.second is a bool indicating inserted or not
                if (!inserted)
                    (sparsityMap[smaller]->at(larger)) += massMatrixJK; // at() returns a reference, so using "+=" is correct
            }
        }
    }

    // 2. according to sparsity map, compute a, ia, ja of CSR formated massMatrix
    massMatrix_IA = new int[numNodes + 1];
    massMatrix_IA[0] = 1; // numbering from 1
    int nnz = 1; // number of non-zero entries
    for (int i = 0; i < numNodes; i++) {
        nnz += sparsityMap[i]->size();
        massMatrix_IA[i + 1] = nnz;
    }
    nnz--;
    massMatrix_JA = new int[nnz];
    massMatrix_A = new double[nnz];
    int count2 = 0;
    for (int i = 0; i < numNodes; i++) {
        for (map<int, double>::iterator it = sparsityMap[i]->begin(); it != sparsityMap[i]->end();
                it++) {
            massMatrix_JA[count2] = it->first + 1; // numbering from 1
            massMatrix_A[count2] = it->second;
            count2++;
        }
    }assert(nnz == count2);

    // delete spasity map
    for (int i = 0; i < numNodes; i++)
        delete sparsityMap[i];
    delete[] sparsityMap;

    for (int i = 0; i < numElems; i++)
        delete directElemTable[i];
    delete[] directElemTable;

    //MortarMath::printCSRMatrixSymmetric(massMatrix_A, massMatrix_IA, massMatrix_JA, numNodes);
}

DataFieldIntegration::~DataFieldIntegration() {
    delete[] massMatrix_A;
    delete[] massMatrix_IA;
    delete[] massMatrix_JA;
    if (pardisoInitialized)
        deletePardiso();
}

void DataFieldIntegration::integrate(const double *tractions, double *forces) {
    int n = numNodes;
    char up = 'u';
    // This routine supports only one-based indexing of the input arrays.
#ifdef USE_INTEL_MKL
    mkl_dcsrsymv(&up, &n, massMatrix_A, massMatrix_IA, massMatrix_JA, const_cast<double*>(tractions), forces);
#else
    MortarMath::dcsrsymv(n, massMatrix_A, massMatrix_IA, massMatrix_JA, tractions, forces);
#endif
}

void DataFieldIntegration::deIntegrate(const double *forces, double *tractions) {
    if (!pardisoInitialized) {
        initPardiso();
        pardisoInitialized = true;
    }
    int n = numNodes;
    char up = 'u';
    // use pardiso solver
    // pardiso forward and backward substitution
    int phase = 33; // forward and backward substitution
    int idum; // integer dummy
    double *ddum = new double[n]; // dummy but the memory is asked for
    int error = 0;
    iparm[5] = 1; // write solution to b (i.e. overwrite)
    for (int i = 0; i < n; i++)
        tractions[i] = forces[i];
#ifdef USE_INTEL_MKL
    mkl_set_num_threads(mklSetNumThreads);
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &neq, massMatrix_A, massMatrix_IA, massMatrix_JA, &idum, &nrhs,
            iparm, &msglvl, tractions, ddum, &error);
#endif
    assert(error == 0);
    delete[] ddum;
}

void DataFieldIntegration::initPardiso() {
#ifdef USE_INTEL_MKL
    mtype = 2; // real symmetric matrix
    // set pardiso default parameters
    pardisoinit(pt, &mtype, iparm);
    iparm[2] = mklSetNumThreads;
    //cout << endl << "iparm[3]: " << iparm[2] << endl;
    maxfct = 1;// max number of factorizations
    mnum = 1;// which factorization to use
    msglvl = 0;// do NOT print statistical information
    neq = numNodes;// number of rows of massMatrix
    nrhs = 1;// number of right hand side
    int phase = 12;// analysis and factorization
    double ddum;// double dummy
    int idum;// integer dummy
    int error = 0;
    //cout<<"factorizing"<<endl;
    mkl_set_num_threads(mklSetNumThreads);
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &neq, massMatrix_A, massMatrix_IA, massMatrix_JA, &idum, &nrhs, iparm,
            &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        cerr << "Error in MortarMapper: pardiso factorization failed!" << error << endl;
        exit(EXIT_FAILURE);
    }
#endif
}

void DataFieldIntegration::deletePardiso() {
#ifdef USE_INTEL_MKL
    int phase = -1; // deallocate memory
    double ddum;// double dummy
    int idum;// integer dummy
    int error = 0;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &neq, massMatrix_A, massMatrix_IA, massMatrix_JA, &idum, &nrhs, iparm,
            &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        cerr << "Error in MortarMapper: pardiso factorization failed!" << error << endl;
        exit(EXIT_FAILURE);
    }
#endif
}

} /* namespace EMPIRE */
