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
#include "NearestNeighborMapper.h"
#include <assert.h>
#include <vector>

#ifdef FLANN
#include "flann/flann.hpp"
#endif

#ifdef ANN
#include "ANN/ANN.h"
#endif

using namespace std;

namespace EMPIRE {

NearestNeighborMapper::NearestNeighborMapper(int _numNodesA, const double *_nodesA, int _numNodesB,
        const double *_nodesB) :
        numNodesA(_numNodesA), nodesA(_nodesA), numNodesB(_numNodesB), nodesB(_nodesB) {
    neighborsTable = new int[numNodesB];
    {
#ifdef ANN
        double **ANodes = new double*[numNodesA]; // ANN uses 2D array
        for (int i = 0; i < numNodesA; i++)
            ANodes[i] = &(nodesA[i * 3]);

        ANNkd_tree *ANodesTree = new ANNkd_tree(ANodes, numNodesA, 3);

        for (int i = 0; i < numNodesB; i++) {
            int nb;
            double dummy;
            ANodesTree->annkSearch(&(nodesB[i * 3]), 1, &nb, &dummy);
            neighborsTable[i] = nb;
        }

        delete[] ANodes;
        delete ANodesTree;
        annClose();
#endif
    }
    {
#ifdef FLANN
        double *nodesBCasted = const_cast<double*>(nodesB);
        flann::Matrix<double> *ANodes = new flann::Matrix<double>(const_cast<double*>(nodesA),
                numNodesA, 3);
        flann::Index<flann::L2<double> > *ANodesTree = new flann::Index<flann::L2<double> >(*ANodes,
                flann::KDTreeSingleIndexParams(1));
        ANodesTree->buildIndex(); // Build binary tree for searching

        for (int i = 0; i < numNodesB; i++) {
            flann::Matrix<double> nodeI(&(nodesBCasted[i * 3]), 1, 3);
            vector<vector<int> > indexes_tmp;
            vector<vector<double> > dists_tmp;
            ANodesTree->knnSearch(nodeI, indexes_tmp, dists_tmp, 1, flann::SearchParams(1));
            int nb = indexes_tmp[0][0];
            neighborsTable[i] = nb;
        }

        delete ANodes;
        delete ANodesTree;
#endif
    }
}

NearestNeighborMapper::~NearestNeighborMapper() {
    delete[] neighborsTable;
}

void NearestNeighborMapper::consistentMapping(const double *DOF_A, double *DOF_B) {
    for (int i = 0; i < numNodesB; i++) {
        int neighbor = neighborsTable[i];
        DOF_B[i] = DOF_A[neighbor];
    }
}

void NearestNeighborMapper::conservativeMapping(const double *DOF_B, double *DOF_A) {
    for (int i = 0; i < numNodesA; i++) {
        DOF_A[i] = 0.0;
    }

    for (int i = 0; i < numNodesB; i++) {
        int neighbor = neighborsTable[i];
        DOF_A[neighbor] += DOF_B[i];
    }

}

} /* namespace EMPIRE */
