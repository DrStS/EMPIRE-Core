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
#include "NearestElementMapper.h"
#include "MortarMath.h"
#ifdef FLANN
#include "flann/flann.hpp"
#endif

#ifdef ANN
#include "ANN/ANN.h"
#endif

#include <assert.h>
#include <math.h>
#include <map>
#include <iostream>
//#include <time.h>

using namespace std;
namespace EMPIRE {
const int NearestElementMapper::MAX_NUM_NEIGHBORS_TO_SEARCH = 10;
// set default value for mapper threads
int NearestElementMapper::mapperSetNumThreads = 1;

NearestElementMapper::NearestElementMapper(int _numNodesA, int _numElemsA,
        const int *_numNodesPerElemA, const double *_nodesA, const int *_nodeIDsA,
        const int *_elemTableA, int _numNodesB, int _numElemsB, const int *_numNodesPerElemB,
        const double *_nodesB, const int *_nodeIDsB, const int *_elemTableB) :
        numNodesA(_numNodesA), numElemsA(_numElemsA), numNodesPerElemA(_numNodesPerElemA), nodesA(
                _nodesA), nodeIDsA(_nodeIDsA), elemTableA(_elemTableA), numNodesB(_numNodesB), numElemsB(
                _numElemsB), numNodesPerElemB(_numNodesPerElemB), nodesB(_nodesB), nodeIDsB(
                _nodeIDsB), elemTableB(_elemTableB) {
    numNodesPerNeighborElem = new int[numNodesB];
    neighborsTable = new vector<int*>(numNodesB);
    weightsTable = new vector<double*>(numNodesB);
    computeNeighborsAndWeights();
}

NearestElementMapper::~NearestElementMapper() {
    delete numNodesPerNeighborElem;
    for (int i = 0; i < numNodesB; i++)
        delete[] neighborsTable->at(i);
    delete neighborsTable;
    for (int i = 0; i < numNodesB; i++)
        delete[] weightsTable->at(i);
    delete weightsTable;
}

void NearestElementMapper::consistentMapping(const double *fieldA, double *fieldB) {
    //time_t timeStart, timeEnd;
    //time(&timeStart);
#pragma omp parallel num_threads(mapperSetNumThreads)
    {
#pragma omp for
        for (int i = 0; i < numNodesB; i++) { // i-th node
            int numNodesMyNeighbor = numNodesPerNeighborElem[i];
            int neighbors[numNodesMyNeighbor];
            for (int j = 0; j < numNodesMyNeighbor; j++) {
                neighbors[j] = neighborsTable->at(i)[j];
            }
            double val = 0.0;
            for (int j = 0; j < numNodesMyNeighbor; j++) { //j-th neighbor
                val += weightsTable->at(i)[j] * fieldA[neighbors[j]];
            }
            fieldB[i] = val;
        }
    } //#pragma omp parallel
    //time(&timeEnd);
    //cout << "It took " << difftime(timeEnd, timeStart) << " seconds for consistentMapping" << endl;
}

void NearestElementMapper::conservativeMapping(const double *fieldB, double *fieldA) {
    for (int i = 0; i < numNodesA; i++) {
        fieldA[i] = 0.0;
    }
    //time_t timeStart, timeEnd;
    //time(&timeStart);
#pragma omp parallel num_threads(mapperSetNumThreads)
    {
#pragma omp for
        for (int i = 0; i < numNodesB; i++) { // i-th node
            int numNodesMyNeighbor = numNodesPerNeighborElem[i];
            int neighbors[numNodesMyNeighbor];
            for (int j = 0; j < numNodesMyNeighbor; j++)
                neighbors[j] = neighborsTable->at(i)[j];
#pragma omp critical
            {
                for (int j = 0; j < numNodesMyNeighbor; j++) { //j-th neighbor
                    fieldA[neighbors[j]] += weightsTable->at(i)[j] * fieldB[i];
                }
            }
        }
    } //#pragma omp parallel
    //time(&timeEnd);
    //cout << "It took " << difftime(timeEnd, timeStart) << " seconds for conservativeMapping"
    //      << endl;
}

void NearestElementMapper::computeNeighborsAndWeights() {
    // compute directElemTableA
    map<int, int> *nodesIDToPosMap = new map<int, int>;
    for (int i = 0; i < numNodesA; i++)
        nodesIDToPosMap->insert(nodesIDToPosMap->end(), pair<int, int>(nodeIDsA[i], i));
    directElemTableA = new vector<int>*[numElemsA];
    for (int i = 0; i < numElemsA; i++)
        directElemTableA[i] = new vector<int>;
    int count = 0;
    for (int i = 0; i < numElemsA; i++) {
        const int numNodesThisElem = numNodesPerElemA[i];
        for (int j = 0; j < numNodesThisElem; j++) {
            directElemTableA[i]->push_back(nodesIDToPosMap->at(elemTableA[count + j]));
        }
        count += numNodesThisElem;
    }
    delete nodesIDToPosMap;
    // compute elementCentroidsA
    double *elementCentroidsA = new double[numElemsA * 3];
    for (int i = 0; i < numElemsA; i++) { // compute the centroid of all elements in A
        int numNodesThisElem = numNodesPerElemA[i];
        double thisElem[numNodesThisElem * 3];
        getElemCoorInA(i, thisElem);
        MortarMath::computePolygonCenter(thisElem, numNodesThisElem, &(elementCentroidsA[i * 3]));
    }

    int NUM_NEIGHBORS_TO_SEARCH = MAX_NUM_NEIGHBORS_TO_SEARCH;
    if (NUM_NEIGHBORS_TO_SEARCH > numElemsA) {
        NUM_NEIGHBORS_TO_SEARCH = numElemsA;
    }
    //time_t timeStart, timeEnd;
    //time(&timeStart);
    {
#ifdef FLANN
        double *nodesBCasted = const_cast<double*>(nodesB);
        flann::Matrix<double> *elementCentroidsA_FLANN = new flann::Matrix<double>(
                const_cast<double*>(elementCentroidsA), numElemsA, 3);
        flann::Index<flann::L2<double> > *ANodesTree = new flann::Index<flann::L2<double> >(
                *elementCentroidsA_FLANN, flann::KDTreeSingleIndexParams(1));
        ANodesTree->buildIndex(); // Build binary tree for searching

#pragma omp parallel num_threads(mapperSetNumThreads)
        {
#pragma omp for
            for (int i = 0; i < numNodesB; i++) {
                flann::Matrix<double> nodeI(&(nodesBCasted[i * 3]), 1, 3);
                vector<vector<int> > indexes_tmp;
                vector<vector<double> > dists_tmp;
                ANodesTree->knnSearch(nodeI, indexes_tmp, dists_tmp, NUM_NEIGHBORS_TO_SEARCH,
                        flann::SearchParams(1));

                for (int j = 0; j < NUM_NEIGHBORS_TO_SEARCH; j++) {
                    int numNodesThisElem = numNodesPerElemA[indexes_tmp[0][j]];
                    if (numNodesThisElem == 3) {
                        double triangle[3 * 3];
                        getElemCoorInA(indexes_tmp[0][j], triangle);
                        double normal[3];
                        MortarMath::computeNormalOfTriangle(triangle, true, normal);
                        int planeToProject = MortarMath::computePlaneToProject(normal);

                        double projection[3];
                        MortarMath::projectToPlane(&triangle[0], normal, &(nodesB[i * 3]), 1,
                                projection);

                        double localCoors[3];
                        MortarMath::computeLocalCoorInTriangle(triangle, planeToProject, projection,
                                localCoors);
                        bool inside = insideElement(3, localCoors);
                        if (inside) {
                            numNodesPerNeighborElem[i] = 3;
                            neighborsTable->at(i) = new int[3];
                            for (int k = 0; k < 3; k++) {
                                neighborsTable->at(i)[k] = directElemTableA[indexes_tmp[0][j]]->at(
                                        k);
                            }
                            weightsTable->at(i) = new double[3];
                            for (int k = 0; k < 3; k++) {
                                weightsTable->at(i)[k] = localCoors[k];
                            }
                            break;
                        }
                    } else if (numNodesThisElem == 4) {
                        double quad[4 * 3];
                        getElemCoorInA(indexes_tmp[0][j], quad);
                        double normal[3];
                        MortarMath::computeNormalOfQuad(quad, true, normal);
                        int planeToProject = MortarMath::computePlaneToProject(normal);

                        { // replace the element by the projection of it on its "element plane"
                            double quadCenter[3];
                            MortarMath::computePolygonCenter(quad, 4, quadCenter);
                            double quadPrj[12];
                            MortarMath::projectToPlane(quadCenter, normal, quad, 4, quadPrj);
                            for (int i = 0; i < 12; i++)
                            quad[i] = quadPrj[i];
                        }
                        double projection[3];
                        MortarMath::projectToPlane(&quad[0], normal, &(nodesB[i * 3]), 1,
                                projection);

                        double localCoors[2];
                        MortarMath::computeLocalCoorInQuad(quad, planeToProject, projection,
                                localCoors);
                        bool inside = insideElement(4, localCoors);
                        if (inside) {
                            numNodesPerNeighborElem[i] = 4;
                            neighborsTable->at(i) = new int[4];
                            for (int k = 0; k < 4; k++) {
                                neighborsTable->at(i)[k] = directElemTableA[indexes_tmp[0][j]]->at(
                                        k);
                            }
                            weightsTable->at(i) = new double[4];
                            MortarMath::computeShapeFuncOfQuad(localCoors, weightsTable->at(i));
                            break;
                        }
                    } else {
                        assert(false);
                    }
                    if (j == NUM_NEIGHBORS_TO_SEARCH - 1) { // projections do not locate inside the neighboring elements, use the nearest one
                        int numNodesThisElem = numNodesPerElemA[indexes_tmp[0][0]];
                        if (numNodesThisElem == 3) {
                            double triangle[3 * 3];
                            getElemCoorInA(indexes_tmp[0][0], triangle);
                            double normal[3];
                            MortarMath::computeNormalOfTriangle(triangle, true, normal);
                            int planeToProject = MortarMath::computePlaneToProject(normal);

                            double projection[3];
                            MortarMath::projectToPlane(&triangle[0], normal, &(nodesB[i * 3]), 1,
                                    projection);

                            double localCoors[3];
                            MortarMath::computeLocalCoorInTriangle(triangle, planeToProject,
                                    projection, localCoors);
                            {
                                numNodesPerNeighborElem[i] = 3;
                                neighborsTable->at(i) = new int[3];
                                for (int k = 0; k < 3; k++) {
                                    neighborsTable->at(i)[k] =
                                    directElemTableA[indexes_tmp[0][0]]->at(k);
                                }
                                weightsTable->at(i) = new double[3];
                                for (int k = 0; k < 3; k++) {
                                    weightsTable->at(i)[k] = localCoors[k];
                                }
                            }
                        } else if (numNodesThisElem == 4) {
                            double quad[4 * 3];
                            getElemCoorInA(indexes_tmp[0][0], quad);
                            double normal[3];
                            MortarMath::computeNormalOfQuad(quad, true, normal);
                            int planeToProject = MortarMath::computePlaneToProject(normal);

                            { // replace the element by the projection of it on its "element plane"
                                double quadCenter[3];
                                MortarMath::computePolygonCenter(quad, 4, quadCenter);
                                double quadPrj[12];
                                MortarMath::projectToPlane(quadCenter, normal, quad, 4, quadPrj);
                                for (int i = 0; i < 12; i++)
                                quad[i] = quadPrj[i];
                            }
                            double projection[3];
                            MortarMath::projectToPlane(&quad[0], normal, &(nodesB[i * 3]), 1,
                                    projection);

                            double localCoors[2];
                            MortarMath::computeLocalCoorInQuad(quad, planeToProject, projection,
                                    localCoors);
                            {
                                numNodesPerNeighborElem[i] = 4;
                                neighborsTable->at(i) = new int[4];
                                for (int k = 0; k < 4; k++) {
                                    neighborsTable->at(i)[k] =
                                    directElemTableA[indexes_tmp[0][0]]->at(k);
                                }
                                weightsTable->at(i) = new double[4];
                                MortarMath::computeShapeFuncOfQuad(localCoors, weightsTable->at(i));
                            }
                        } else {
                            assert(false);
                        }
                    }
                }

            }
        } //#pragma omp parallel
        delete elementCentroidsA_FLANN;
        delete ANodesTree;
#endif
    }
    //time(&timeEnd);
    //cout << "It took " << difftime(timeEnd, timeStart) << " seconds for neighbor search" << endl;
    for (int i = 0; i < numElemsA; i++)
        delete directElemTableA[i];
    delete[] directElemTableA;
    delete[] elementCentroidsA;
}

void NearestElementMapper::getElemCoorInA(int elemIndex, double *elem) {
    // compute the coordinates of an element by its id
    int numNodesThisElem = numNodesPerElemA[elemIndex];
    for (int i = 0; i < numNodesThisElem; i++) {
        int nodePos = directElemTableA[elemIndex]->at(i); // position of the node
        for (int j = 0; j < 3; j++) {
            elem[i * 3 + j] = nodesA[nodePos * 3 + j];
        }
    }
}

bool NearestElementMapper::insideElement(int numNodesThisElem, double *localCoor) {
    const double EPS = 1e-10;
    if (numNodesThisElem == 3) {
        for (int i = 0; i < 3; i++) {
            if (localCoor[i] > 1.0 + EPS)
                return false;
            if (localCoor[i] < 0.0 - EPS)
                return false;
        }
        return true;
    } else if (numNodesThisElem == 4) {
        for (int i = 0; i < 2; i++) {
            if (localCoor[i] > 1.0 + EPS)
                return false;
            if (localCoor[i] < -1.0 - EPS)
                return false;
        }
        return true;
    } else {
        assert(false);
        return false;
    }
}

} /* namespace EMPIRE */
