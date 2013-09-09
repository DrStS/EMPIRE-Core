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
#include "BarycentricInterpolationMapper.h"
#include "MortarMath.h"
#include <assert.h>
#include <math.h>
#include <vector>

#ifdef FLANN
#include "flann/flann.hpp"
#endif

#ifdef ANN
#include "ANN/ANN.h"
#endif

using namespace std;

namespace EMPIRE {

const int BarycentricInterpolationMapper::MAX_NUM_NEIGHBORS_TO_SEARCH = 50;

BarycentricInterpolationMapper::BarycentricInterpolationMapper(int _numNodesA,
        const double *_nodesA, int _numNodesB, const double *_nodesB) :
        numNodesA(_numNodesA), nodesA(_nodesA), numNodesB(_numNodesB), nodesB(_nodesB) {
    neighborsTable = new int[numNodesB * 3];
    weightsTable = new double[numNodesB * 3];

    computeNeighbors();
    computeWeights();
}

BarycentricInterpolationMapper::~BarycentricInterpolationMapper() {
    delete neighborsTable;
    delete weightsTable;
}

void BarycentricInterpolationMapper::consistentMapping(const double *fieldA, double *fieldB) {
    for (int i = 0; i < numNodesB; i++) { // i-th node
        int neighbors[3];
        for (int j = 0; j < 3; j++) {
            neighbors[j] = neighborsTable[i * 3 + j];
        }
        double val = 0.0;
        for (int j = 0; j < 3; j++) { //j-th neighbor
            val += weightsTable[i * 3 + j] * fieldA[neighbors[j]];
        }
        fieldB[i] = val;
    }
}

void BarycentricInterpolationMapper::conservativeMapping(const double *fieldB,
        double *fieldA) {
    for (int i = 0; i < numNodesA; i++) {
        fieldA[i] = 0.0;
    }
    for (int i = 0; i < numNodesB; i++) { // i-th node
        int neighbors[3];
        for (int j = 0; j < 3; j++)
            neighbors[j] = neighborsTable[i * 3 + j];
        for (int j = 0; j < 3; j++) { //j-th neighbor
            fieldA[neighbors[j]] += weightsTable[i * 3 + j] * fieldB[i];
        }
    }
}

void BarycentricInterpolationMapper::computeNeighbors() {
    int NUM_NEIGHBORS_TO_SEARCH = MAX_NUM_NEIGHBORS_TO_SEARCH;
    if (NUM_NEIGHBORS_TO_SEARCH > numNodesA) {
        NUM_NEIGHBORS_TO_SEARCH = numNodesA;
    }
    {
#ifdef ANN
        double **ANodes = new double*[numNodesA]; // ANN uses 2D array
        for (int i = 0; i < numNodesA; i++)
        ANodes[i] = &(nodesA[i * 3]);

        ANNkd_tree *ANodesTree = new ANNkd_tree(ANodes, numNodesA, 3);

        for (int i = 0; i < numNodesB; i++) {
            int nb[NUM_NEIGHBORS_TO_SEARCH];
            double dummy[NUM_NEIGHBORS_TO_SEARCH];
            ANodesTree->annkSearch(&(nodesB[i * 3]), NUM_NEIGHBORS_TO_SEARCH, nb, dummy);
            const double *p1 = &(nodesA[nb[0] * 3]);
            const double *p2 = &(nodesA[nb[1] * 3]);
            int p3Index = 2;
            while (true) {
                const double *p3 = &(nodesA[nb[p3Index] * 3]);
                if (areNotOnTheSameLine(p1, p2, p3)) {
                    break;
                }
                p3Index++;
                assert(p3Index != NUM_NEIGHBORS_TO_SEARCH);
            }

            neighborsTable[i * 3 + 0] = nb[0];
            neighborsTable[i * 3 + 1] = nb[1];
            neighborsTable[i * 3 + 2] = nb[p3Index];
        }

        delete[] ANodes;
        delete ANodesTree;
        annClose();
#endif
    }
    {
#ifdef FLANN
        double *nodesBCasted = const_cast<double*>(nodesB);
        flann::Matrix<double> *ANodes = new flann::Matrix<double>(
                const_cast<double*>(nodesA), numNodesA, 3);
        flann::Index<flann::L2<double> > *ANodesTree = new flann::Index<flann::L2<double> >(*ANodes,
                flann::KDTreeSingleIndexParams(1));
        ANodesTree->buildIndex(); // Build binary tree for searching

        for (int i = 0; i < numNodesB; i++) {
            flann::Matrix<double> nodeI(&(nodesBCasted[i * 3]), 1, 3);
            vector<vector<int> > indexes_tmp;
            vector<vector<double> > dists_tmp;
            ANodesTree->knnSearch(nodeI, indexes_tmp, dists_tmp, NUM_NEIGHBORS_TO_SEARCH,
                    flann::SearchParams(1));

            const double *p1 = &(nodesA[indexes_tmp[0][0] * 3]);
            const double *p2 = &(nodesA[indexes_tmp[0][1] * 3]);
            int p3Index = 2;
            while (true) {
                const double *p3 = &(nodesA[indexes_tmp[0][p3Index] * 3]);
                if (areNotOnTheSameLine(p1, p2, p3)) {
                    break;
                }
                p3Index++;
                assert(p3Index != NUM_NEIGHBORS_TO_SEARCH);
            }

            neighborsTable[i * 3 + 0] = indexes_tmp[0][0];
            neighborsTable[i * 3 + 1] = indexes_tmp[0][1];
            neighborsTable[i * 3 + 2] = indexes_tmp[0][p3Index];
        }

        delete ANodes;
        delete ANodesTree;
#endif
    }
}

void BarycentricInterpolationMapper::computeWeights() {
    for (int i = 0; i < numNodesB; i++) { // i-th node in B
        const double *nodeI = &(nodesB[i * 3]);
        double triangle[9];
        for (int j = 0; j < 3; j++) { //j-th neighbor in A
            int neighbor = neighborsTable[i * 3 + j];
            for (int k = 0; k < 3; k++) { // k-th coordinate
                triangle[j * 3 + k] = nodesA[neighbor * 3 + k];
            }
        }
        double normal[3];
        MortarMath::computeNormalOfTriangle(triangle, false, normal);
        int planeToProject = MortarMath::computePlaneToProject(normal);
        double localCoors[3];
        MortarMath::computeLocalCoorInTriangle(triangle, planeToProject, nodeI, localCoors);
        for (int j = 0; j < 3; j++) {
            weightsTable[i * 3 + j] = localCoors[j];
        }
    }
}

bool BarycentricInterpolationMapper::areNotOnTheSameLine(const double *p1, const double *p2,
        const double *p3) {
    const double TOL = 1E-3; // tolerance for determining whether or not on the same line
    double p1p2Square = MortarMath::distanceSquare(p1, p2);
    double p2p3Square = MortarMath::distanceSquare(p2, p3);
    double p3p1Square = MortarMath::distanceSquare(p3, p1);
    double refSize = 0.0;
    if (p1p2Square <= p2p3Square && p1p2Square <= p3p1Square) {
        refSize = p2p3Square * p3p1Square;
    } else if (p2p3Square <= p1p2Square && p2p3Square <= p3p1Square) {
        refSize = p1p2Square * p3p1Square;
    } else if (p3p1Square <= p1p2Square && p3p1Square <= p2p3Square) {
        refSize = p1p2Square * p2p3Square;
    } else {
        assert(false);
    }
    refSize = sqrt(refSize);
    double triangle[9];
    triangle[0 * 3 + 0] = p1[0];
    triangle[0 * 3 + 1] = p1[1];
    triangle[0 * 3 + 2] = p1[2];
    triangle[1 * 3 + 0] = p2[0];
    triangle[1 * 3 + 1] = p2[1];
    triangle[1 * 3 + 2] = p2[2];
    triangle[2 * 3 + 0] = p3[0];
    triangle[2 * 3 + 1] = p3[1];
    triangle[2 * 3 + 2] = p3[2];
    double area = MortarMath::computeAreaOfTriangle(triangle);

    return area > TOL * refSize;
}

} /* namespace EMPIRE */
