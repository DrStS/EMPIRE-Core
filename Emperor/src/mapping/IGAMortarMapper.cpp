/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Munich
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
/*
 * IGAMortarMapper.cpp
 *
 *  Created on: May 8, 2013
 *      Author: chenshen
 */

#include "IGAMortarMapper.h"
#include "IGAMortarMath.h"
#include "MathLibrary.h"
#include "FEMesh.h"
#include "IGAPatchSurface.h"
#include "IGAMesh.h"
#include "DataField.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using namespace std;

namespace EMPIRE {

IGAMortarMapper::IGAMortarMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE,
        double _disTol, int _numGPsTri, int _numGPsQuad) :
        meshIGA(_meshIGA), disTol(_disTol), numGPsTri(_numGPsTri), numGPsQuad(_numGPsQuad) {

    assert(_meshIGA != NULL);
    assert(_meshFE != NULL);
    assert(_meshIGA->type == EMPIRE_Mesh_IGAMesh);
    assert(_meshFE->type == EMPIRE_Mesh_FEMesh);

    if (_meshFE->triangulate() == NULL)
        meshFE = _meshFE;
    else
        meshFE = _meshFE->triangulate();

    projectedCoords = new vector<map<int, double*> >(meshFE->numNodes);

    C_NR = new MathLibrary::SparseMatrix<double>((const size_t) meshFE->numNodes,
            (const size_t) meshIGA->getNumNodes());
    C_NN = new MathLibrary::SparseMatrix<double>((const size_t) meshFE->numNodes, true);

    gaussTriangle = new IGAMortarMath::GaussQuadratureOnTriangle(numGPsTri);
    gaussQuad = new IGAMortarMath::GaussQuadratureOnQuad(numGPsQuad);

    initTables();

    projectPointsToSurface();

    computeCouplingMatrices();

    C_NN->factorize();

}

IGAMortarMapper::~IGAMortarMapper() {

    for (int i = 0; i < meshFE->numElems; i++)
        delete[] meshFEDirectElemTable[i];
    delete[] meshFEDirectElemTable;

    for (int nodeIndex = 0; nodeIndex < meshFE->numNodes; nodeIndex++) {
        for (int patchCount = 0; patchCount < meshIGA->getSurfacePatches().size(); patchCount++) {
            if ((*projectedCoords)[nodeIndex].find(patchCount)
                    != (*projectedCoords)[nodeIndex].end()) {
                delete[] (*projectedCoords)[nodeIndex][patchCount];
            }
        }
    }
    delete projectedCoords;

    delete gaussTriangle;
    delete gaussQuad;

    C_NN->cleanPardiso();
    delete C_NR;
    delete C_NN;
}

void IGAMortarMapper::initTables() {
    /* using the map to store the nodeIDs
     * but here the "key" is the node ID, and the value is the position in nodeIDs
     * the map is sorted automatically, so it is efficient for searching
     */

    // compute direct element table for fluid mesh
    meshFEDirectElemTable = new int*[meshFE->numElems]; // deleted
    for (int i = 0; i < meshFE->numElems; i++)
        meshFEDirectElemTable[i] = new int[meshFE->numNodesPerElem[i]];

    map<int, int> *meshFENodesMap = new map<int, int>(); // deleted
    for (int i = 0; i < meshFE->numNodes; i++)
        meshFENodesMap->insert(meshFENodesMap->end(), pair<int, int>(meshFE->nodeIDs[i], i));
    int count = 0;

    for (int i = 0; i < meshFE->numElems; i++) {
        const int numNodesPerElem = meshFE->numNodesPerElem[i];

        for (int j = 0; j < numNodesPerElem; j++) {
            if (meshFENodesMap->find(meshFE->elems[count + j]) == meshFENodesMap->end()) {
                ERROR_OUT() << "Cannot find node ID " << meshFE->elems[count + j] << endl;
                exit(-1);
            }
            meshFEDirectElemTable[i][j] = meshFENodesMap->at(meshFE->elems[count + j]);
        }
        count += numNodesPerElem;
    }

    delete meshFENodesMap;

}

void IGAMortarMapper::projectPointsToSurface() {
    /* Loop over all patches
     *      Loop over all the elements on the fluid side
     *          1. Loop over all nodes of the current element to check if there exist one node has been projected already
     *          2. Check if there exist one node in the current element has been successfully projected. If not,
     *              2.1 if so, use result of the projected node as the initial guess
     *              2.2 otherwise,  find the nearest knot intersection as initial guess
     *          3. Loop over each node at the current element
     *              3.1 Check if the node has been already projected
     *              3.2 if not, compute the point projection on the IGA patch using the initial guess get from the last step
     * The result of the projected coordinates are stored in the class member: projectedCoords
     */

    int numPatches = meshIGA->getSurfacePatches().size();
    for (int patchCount = 0; patchCount < numPatches; patchCount++) {

        IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchCount];

        bool *isProjected = new bool[meshFE->numNodes]; // deleted
        for (int i = 0; i < meshFE->numNodes; i++)
            isProjected[i] = false;

        /// Loop over all the elements on the fluid side
        for (int i = 0; i < meshFE->numElems; i++) {

            bool isNodeInsideElementProjecteded = false;
            int projectedNode = 0;

            /// 1. Loop over all nodes of the current element to check if there exist one node has been projected already
            for (int j = 0; j < meshFE->numNodesPerElem[i]; j++)
                if (isProjected[meshFEDirectElemTable[i][j]]) {
                    isNodeInsideElementProjecteded = true;
                    projectedNode = meshFEDirectElemTable[i][j];
                    break;
                }

            // Coordinates of the projected nodes on the IGA surface
            double initialU, initialV;

            /// 2. Check if there exist one node in the current element has been successfully projected. If not,
            if (isNodeInsideElementProjecteded) {

                /// 2.1 if so, use result of the projected node as the initial guess
                initialU = (*projectedCoords)[projectedNode][patchCount][0];
                initialV = (*projectedCoords)[projectedNode][patchCount][1];

            } else {

                /// 2.2 otherwise,  find the nearest knot intersection as initial guess
                int nodeIndex = meshFEDirectElemTable[i][0];
                double cartesianCoords[3];
                cartesianCoords[0] = meshFE->nodes[nodeIndex * 3];
                cartesianCoords[1] = meshFE->nodes[nodeIndex * 3 + 1];
                cartesianCoords[2] = meshFE->nodes[nodeIndex * 3 + 2];

                thePatch->findInitialGuess4PointProjection(initialU, initialV, cartesianCoords);
            }

            /// 3. Loop over each node at the current element
            for (int j = 0; j < meshFE->numNodesPerElem[i]; j++) {

                int nodeIndex = meshFEDirectElemTable[i][j];

                /// Check if the node has been already projected
                if (!isProjected[nodeIndex]) {

                    /// if not, compute the point projection on the IGA patch using the initial guess get from the last step
                    /// the result are stored in the class member "projectedCoords"
                    double cartesianCoords[3];
                    cartesianCoords[0] = meshFE->nodes[nodeIndex * 3];
                    cartesianCoords[1] = meshFE->nodes[nodeIndex * 3 + 1];
                    cartesianCoords[2] = meshFE->nodes[nodeIndex * 3 + 2];

                    double projectedU = initialU;
                    double projectedV = initialV;

                    bool converge;
                    bool convergeInside = thePatch->computePointProjectionOnPatch(projectedU,
                            projectedV, cartesianCoords, converge);

                    if (convergeInside
                            && IGAMortarMath::computePointDistance(&meshFE->nodes[nodeIndex * 3],
                                    cartesianCoords) < disTol) {

                        isProjected[nodeIndex] = true;
                        double* coordTmp = new double(2);
                        coordTmp[0] = projectedU;
                        coordTmp[1] = projectedV;

                        (*projectedCoords)[nodeIndex].insert(
                                std::pair<int, double*>(patchCount, coordTmp));

                    }
                }
            }
        }

        delete[] isProjected;
    }

    double U, V;
    double P[3];
    for (int nodeIndex = 0; nodeIndex < meshFE->numNodes; nodeIndex++) {

        bool isProjected = false;
        for (int patchCount = 0; patchCount < numPatches; patchCount++) {

            if ((*projectedCoords)[nodeIndex].find(patchCount)
                    != (*projectedCoords)[nodeIndex].end()) {
                isProjected = true;
                break;
            }
        }

        if (!isProjected) {

            for (int patchCount = 0; patchCount < numPatches; patchCount++) {
                IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchCount];

                P[0] = meshFE->nodes[nodeIndex * 3];
                P[1] = meshFE->nodes[nodeIndex * 3 + 1];
                P[2] = meshFE->nodes[nodeIndex * 3 + 2];

                thePatch->findInitialGuess4PointProjection(U, V, P, REFINED_NUM_PARAMETRIC_LOCATIONS, REFINED_NUM_PARAMETRIC_LOCATIONS);

                bool isConverge;
                bool isConvergeInside = thePatch->computePointProjectionOnPatch(U, V, P,
                        isConverge);

                if (isConvergeInside
                        && IGAMortarMath::computePointDistance(&meshFE->nodes[nodeIndex * 3], P)
                                < disTol) {
                    isProjected = true;
                    double* coordTmp = new double(2);
                    coordTmp[0] = U;
                    coordTmp[1] = V;

                    (*projectedCoords)[nodeIndex].insert(
                            std::pair<int, double*>(patchCount, coordTmp));
                }
            }
            if (!isProjected) {
                ERROR_OUT() << " in IGAMortarMapper::projectPointsToSurface" << endl;
                ERROR_OUT() << "Cannot project node: " << nodeIndex << "  ("
                        << meshFE->nodes[nodeIndex * 3] << ", " << meshFE->nodes[nodeIndex * 3 + 1]
                        << ", " << meshFE->nodes[nodeIndex * 3 + 2] << ")" << endl;
                exit (EXIT_FAILURE);
            } else {
            }

        }
    }
}

void IGAMortarMapper::computeCouplingMatrices() {
    /*
     * Computes the coupling matrices CNR and CNN.
     * Loop over all the elements in the FE side
     * ->
     * 1. Find whether the projected FE element is located on one patch
     *
     * 2. Compute the coupling matrices
     * ->
     * 2i. If the current element can be projected on one patch
     * 2ii. If the current element cannot be projected on one patch
     * <-
     */
    // The vertices of the canonical polygons
    double parentTriangle[6] = { 0, 0, 1, 0, 0, 1 };
    double parentQuadriliteral[8] = { -1, -1, 1, -1, 1, 1, -1, 1 };

    /// Loop over all the elements in the FE side
    for (int elemCount = 0; elemCount < meshFE->numElems; elemCount++) {

        // Compute the number of shape functions. Depending on number of nodes in the current element
        int numNodesElementFE = meshFE->numNodesPerElem[elemCount];

        // The vertices of the FE element in the parent domain
        double* projectedElementFEWZ;

        // Decide whether the parent element is a triangle or a quadrilateral
        if (numNodesElementFE == 3)
            projectedElementFEWZ = parentTriangle;
        else
            projectedElementFEWZ = parentQuadriliteral;

        /// 1. Find whether the projected FE element is located on one patch

        // Initialize the patch that possibly contains that element
        IGAPatchSurface* thePatch = NULL;

        // Initialize the index of this patch
        int patchIndex = 0;

        // Initialize the flag whether the projected FE element is located on one patch
        bool isAllNodesOnPatch = true;

        // Loop over all the patches
        for (int patchCount = 0; patchCount < meshIGA->getSurfacePatches().size(); patchCount++) {
            isAllNodesOnPatch = true;

            // Loop over all the nodes of the unclipped element
            for (int nodeCount = 0; nodeCount < numNodesElementFE; nodeCount++) {
                // Find the index of the node in the FE mesh
                int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];

                // Find whether this index is in the projected nodes array
                bool isNodeOnPatch = (*projectedCoords)[nodeIndex].find(patchCount)
                        != (*projectedCoords)[nodeIndex].end();

                if (!isNodeOnPatch) {
                    isAllNodesOnPatch = false;
                    break;
                }
            }

            // If all nodes are on the patch save this patch
            if (isAllNodesOnPatch) {
                thePatch = meshIGA->getSurfacePatches()[patchCount];
                patchIndex = patchCount;
                break;
            }
        }

        /// 2. Compute the coupling matrices

        /// 2i. If the current element can be projected on one patch
        if (isAllNodesOnPatch) {
            // Get the projected coordinates for the current element
            double clippedByPatchProjElementFEUV[numNodesElementFE * 2];

            for (int nodeCount = 0; nodeCount < numNodesElementFE; nodeCount++) {
                int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];
                clippedByPatchProjElementFEUV[nodeCount * 2] =
                        (*projectedCoords)[nodeIndex][patchIndex][0];
                clippedByPatchProjElementFEUV[nodeCount * 2 + 1] =
                        (*projectedCoords)[nodeIndex][patchIndex][1];
            }
            computeCouplingMatrices4ClippedByPatchProjectedElement(thePatch, numNodesElementFE,
                    clippedByPatchProjElementFEUV, projectedElementFEWZ, elemCount,
                    numNodesElementFE);
        }
        /// 2ii. If the current element cannot be projected on one patch
        else {
            // Loop over all the patches in the IGA mesh
            for (int patchCount = 0; patchCount < meshIGA->getSurfacePatches().size();
                    patchCount++) {
                IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchCount];
                int numNodesClippedByPatchProjElementFE = 0;
                // the parameter coordinates in IGA of the sub-element divided by the patch
                double clippedByPatchProjElementFEUV[16];

                // the parameter coordinates in low order element of the sub-element divided by the patch
                double clippedByPatchProjElementFEWZ[16];

                // if the line segment connecting the node with the next one is a edge in the sub-element
                bool isEdge[8];

                for (int i = 0; i < 8; i++)
                    isEdge[i] = false;

                // Find nodes and edges of low order element inside IGA Patch
                for (int nodeCount = 0; nodeCount < numNodesElementFE; nodeCount++) {
                    int nodeCountNext = (nodeCount + 1) % numNodesElementFE;
                    int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];
                    int nodeIndexNext = meshFEDirectElemTable[elemCount][nodeCountNext];

                    // Flag on whether the node is inside the current NURBS patch
                    bool isInsidePatch = (*projectedCoords)[nodeIndex].find(patchCount)
                            != (*projectedCoords)[nodeIndex].end();
                    bool isNextNodeInsidePatch = (*projectedCoords)[nodeIndexNext].find(patchCount)
                            != (*projectedCoords)[nodeIndexNext].end();

                    if (isInsidePatch) {
                        // if the node is inside the patch, put it into the Clipped By Patch Projected Element
                        clippedByPatchProjElementFEUV[numNodesClippedByPatchProjElementFE * 2] =
                                (*projectedCoords)[nodeIndex][patchCount][0];
                        clippedByPatchProjElementFEUV[numNodesClippedByPatchProjElementFE * 2 + 1] =
                                (*projectedCoords)[nodeIndex][patchCount][1];
                        clippedByPatchProjElementFEWZ[numNodesClippedByPatchProjElementFE * 2] =
                                projectedElementFEWZ[nodeCount * 2];
                        clippedByPatchProjElementFEWZ[numNodesClippedByPatchProjElementFE * 2 + 1] =
                                projectedElementFEWZ[nodeCount * 2 + 1];
                        isEdge[numNodesClippedByPatchProjElementFE] = true;
                        numNodesClippedByPatchProjElementFE++;

                        // if the node is inside and the next node is outside the patch,
                        // find the intersection with patch boundary, and put it into the Clipped By Patch Projected Element
                        if (!isNextNodeInsidePatch) {
                            double u = (*projectedCoords)[nodeIndex][patchCount][0], v =
                                    (*projectedCoords)[nodeIndex][patchCount][1];
                            double div, dis;
                            double* P1 = &(meshFE->nodes[nodeIndex * 3]);
                            double* P2 = &(meshFE->nodes[nodeIndexNext * 3]);

                            bool isProjectedOnPatchBoundary =
                                    thePatch->computePointProjectionOnPatchBoundary(u, v, div, dis,
                                            P1, P2);
                            if (isProjectedOnPatchBoundary && dis <= disTol) {
                                clippedByPatchProjElementFEUV[numNodesClippedByPatchProjElementFE
                                        * 2] = u;
                                clippedByPatchProjElementFEUV[numNodesClippedByPatchProjElementFE
                                        * 2 + 1] = v;
                                double P1x = projectedElementFEWZ[nodeCount * 2];
                                double P1y = projectedElementFEWZ[nodeCount * 2 + 1];
                                double P2x = projectedElementFEWZ[nodeCountNext * 2];
                                double P2y = projectedElementFEWZ[nodeCountNext * 2 + 1];
                                clippedByPatchProjElementFEWZ[numNodesClippedByPatchProjElementFE
                                        * 2] = P1x * (1 - div) + P2x * div;
                                clippedByPatchProjElementFEWZ[numNodesClippedByPatchProjElementFE
                                        * 2 + 1] = P1y * (1 - div) + P2y * div;
                                numNodesClippedByPatchProjElementFE++;
                            } else {
                                ERROR_OUT() << "Error in IGAMortarMapper::computeCouplingMatrices"
                                        << endl;
                                ERROR_OUT() << "Cannot find point projection on patch boundary"
                                        << endl;
                                ERROR_OUT()
                                        << "Cannot find point projection on patch boundary between node ["
                                        << nodeIndex << "]:(" << meshFE->nodes[nodeIndex * 3] << ","
                                        << meshFE->nodes[nodeIndex * 3 + 1] << ","
                                        << meshFE->nodes[nodeIndex * 3 + 2] << ") and node ["
                                        << nodeIndexNext << "]:("
                                        << meshFE->nodes[nodeIndexNext * 3] << ","
                                        << meshFE->nodes[nodeIndexNext * 3 + 1] << ","
                                        << meshFE->nodes[nodeIndexNext * 3 + 2] << ") on patch ["
                                        << patchCount << "] boundary" << endl;
                                exit (EXIT_FAILURE);
                            }
                        }
                    } else if (isNextNodeInsidePatch) {
                        // if this node is outside and the next node is inside, find the intersection with patch boundary
                        // and put it into the Clipped By Patch Projected Element
                        double u = (*projectedCoords)[nodeIndexNext][patchCount][0], v =
                                (*projectedCoords)[nodeIndexNext][patchCount][1];
                        double div, dis;
                        double* P1 = &(meshFE->nodes[nodeIndex * 3]);
                        double* P2 = &(meshFE->nodes[nodeIndexNext * 3]);
                        bool isProjectedOnPatchBoundary =
                                thePatch->computePointProjectionOnPatchBoundary(u, v, div, dis, P1,
                                        P2);
                        if (isProjectedOnPatchBoundary && dis <= disTol) {
                            clippedByPatchProjElementFEUV[numNodesClippedByPatchProjElementFE * 2] =
                                    u;
                            clippedByPatchProjElementFEUV[numNodesClippedByPatchProjElementFE * 2
                                    + 1] = v;
                            double P1x = projectedElementFEWZ[nodeCount * 2];
                            double P1y = projectedElementFEWZ[nodeCount * 2 + 1];
                            double P2x = projectedElementFEWZ[nodeCountNext * 2];
                            double P2y = projectedElementFEWZ[nodeCountNext * 2 + 1];
                            clippedByPatchProjElementFEWZ[numNodesClippedByPatchProjElementFE * 2] =
                                    P1x * (1 - div) + P2x * div;
                            clippedByPatchProjElementFEWZ[numNodesClippedByPatchProjElementFE * 2
                                    + 1] = P1y * (1 - div) + P2y * div;
                            isEdge[numNodesClippedByPatchProjElementFE] = true;
                            clippedByPatchProjElementFEWZ[numNodesClippedByPatchProjElementFE * 2] =
                                    P1x * (1 - div) + P2x * div;
                            clippedByPatchProjElementFEWZ[numNodesClippedByPatchProjElementFE * 2
                                    + 1] = P1y * (1 - div) + P2y * div;
                            numNodesClippedByPatchProjElementFE++;
                        } else {
                            ERROR_OUT() << "Error in IGAMortarMapper::computeCouplingMatrices"
                                    << endl;
                            ERROR_OUT() << "Cannot find point projection on patch boundary" << endl;
                            ERROR_OUT()
                                    << "Cannot find point projection on patch boundary between node ["
                                    << nodeIndex << "]:(" << meshFE->nodes[nodeIndex * 3] << ","
                                    << meshFE->nodes[nodeIndex * 3 + 1] << ","
                                    << meshFE->nodes[nodeIndex * 3 + 2] << ") and node ["
                                    << nodeIndexNext << "]:(" << meshFE->nodes[nodeIndexNext * 3]
                                    << "," << meshFE->nodes[nodeIndexNext * 3 + 1] << ","
                                    << meshFE->nodes[nodeIndexNext * 3 + 2] << ") on patch ["
                                    << patchCount << "] boundary" << endl;
                            exit (EXIT_FAILURE);
                        }
                    }
                }

                // Find the corner node of the IGA patch which is inside the projected FE element
                if (numNodesClippedByPatchProjElementFE >= 3) {
                    // Find direction of FE element(Clock or Counter-Clock)
                    double x1 = clippedByPatchProjElementFEUV[2] - clippedByPatchProjElementFEUV[0];
                    double y1 = clippedByPatchProjElementFEUV[3] - clippedByPatchProjElementFEUV[1];
                    double x2 = clippedByPatchProjElementFEUV[4] - clippedByPatchProjElementFEUV[2];
                    double y2 = clippedByPatchProjElementFEUV[5] - clippedByPatchProjElementFEUV[3];
                    bool isCounterClockWise = IGAMortarMath::computeCrossProduct2D(x1, y1, x2, y2)
                            > 0;

                    // IGA Patch
                    double *knotVectorU =
                            thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
                    double *knotVectorV =
                            thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();
                    int numKnotsU = thePatch->getIGABasis()->getUBSplineBasis1D()->getNoKnots();
                    int numKnotsV = thePatch->getIGABasis()->getVBSplineBasis1D()->getNoKnots();
                    double u0 = knotVectorU[0];
                    double uEnd = knotVectorU[numKnotsU - 1];
                    double v0 = knotVectorV[0];
                    double vEnd = knotVectorV[numKnotsV - 1];
                    // Corner nodes of the IGA Patch
                    double patchCornerNodesUV[8] = { u0, v0, uEnd, v0, uEnd, vEnd, u0, vEnd };
                    // Corner nodes of the IGA Patch which inside the current element
                    int patchCornerNodesInsideElemIndex[4];
                    int numPatchCornerNodesInsideElem = 0;

                    // Loop over the 4 corner node of the IGA Patch to check which of them are in side the current element
                    for (int nodeCount = 0; nodeCount < 4; nodeCount++) {
                        bool isNodeInsideElem = true;
                        for (int edgeCount = 0; edgeCount < numNodesClippedByPatchProjElementFE;
                                edgeCount++) {
                            if (isEdge[edgeCount]) {
                                int node1 = edgeCount;
                                int node2 = (edgeCount + 1) % numNodesClippedByPatchProjElementFE;

                                x1 = clippedByPatchProjElementFEUV[node2 * 2]
                                        - clippedByPatchProjElementFEUV[node1 * 2];
                                y1 = clippedByPatchProjElementFEUV[node2 * 2 + 1]
                                        - clippedByPatchProjElementFEUV[node1 * 2 + 1];
                                x2 = patchCornerNodesUV[nodeCount * 2]
                                        - clippedByPatchProjElementFEUV[node2 * 2];
                                y2 = patchCornerNodesUV[nodeCount * 2 + 1]
                                        - clippedByPatchProjElementFEUV[node2 * 2 + 1];
                                bool isRightHandSide = IGAMortarMath::computeCrossProduct2D(x1, y1,
                                        x2, y2) > 0;
                                if (isRightHandSide != isCounterClockWise) {
                                    isNodeInsideElem = false;
                                    break;
                                }
                            }
                        }
                        if (isNodeInsideElem) {
                            patchCornerNodesInsideElemIndex[numPatchCornerNodesInsideElem] =
                                    nodeCount;
                            numPatchCornerNodesInsideElem++;
                        }
                    }

                    // center of all nodes inside element
                    double centerU = 0.0;
                    double centerV = 0.0;
                    // corner nodes of IGA Patch
                    for (int nodeCount = 0; nodeCount < numPatchCornerNodesInsideElem;
                            nodeCount++) {
                        centerU +=
                                patchCornerNodesUV[patchCornerNodesInsideElemIndex[nodeCount] * 2];
                        centerV += patchCornerNodesUV[patchCornerNodesInsideElemIndex[nodeCount] * 2
                                + 1];
                    }
                    // nodes of FE and their intersection between the edge of FE element and IGA patch
                    for (int nodeCount = 0; nodeCount < numNodesClippedByPatchProjElementFE;
                            nodeCount++) {
                        centerU += clippedByPatchProjElementFEUV[nodeCount * 2];
                        centerV += clippedByPatchProjElementFEUV[nodeCount * 2 + 1];
                    }
                    centerU /= numPatchCornerNodesInsideElem + numNodesClippedByPatchProjElementFE;
                    centerV /= numPatchCornerNodesInsideElem + numNodesClippedByPatchProjElementFE;

                    // angle between the line of node to the center point to the Y-axis
                    double ang[numNodesClippedByPatchProjElementFE];
                    for (int nodeCount = 0; nodeCount < numNodesClippedByPatchProjElementFE;
                            nodeCount++) {
                        double x0 = clippedByPatchProjElementFEUV[nodeCount * 2] - centerU;
                        double y0 = clippedByPatchProjElementFEUV[nodeCount * 2 + 1] - centerV;
                        ang[nodeCount] = atan2(y0, x0);
                    }

                    // Cartesian coordinates of the low order element
                    double elementFEXYZ[12];
                    for (int i = 0; i < numNodesElementFE; i++) {
                        int nodeIndex = meshFEDirectElemTable[elemCount][i];
                        for (int j = 0; j < 3; j++)
                            elementFEXYZ[i * 3 + j] = meshFE->nodes[nodeIndex * 3 + j];
                    }

                    // Loop over all corner nodes of IGA Patch inside the element and put them into the clipped by patch element in a correct order
                    for (int nodeCount = 0; nodeCount < numPatchCornerNodesInsideElem;
                            nodeCount++) {
                        double x0 = patchCornerNodesUV[patchCornerNodesInsideElemIndex[nodeCount]
                                * 2] - centerU;
                        double y0 = patchCornerNodesUV[patchCornerNodesInsideElemIndex[nodeCount]
                                * 2 + 1] - centerV;
                        double ang0 = atan2(y0, x0);

                        // Loop over all edges(nodes) of the clipped-by-patch-element to find the correct position for the corner nodes of the patch
                        for (int edgeCount = numNodesClippedByPatchProjElementFE + nodeCount - 1;
                                edgeCount >= 0; edgeCount--) {
                            int node1 = edgeCount;
                            int node2 = (edgeCount + 2)
                                    % (numNodesClippedByPatchProjElementFE + nodeCount + 1);
                            double x1 = clippedByPatchProjElementFEUV[node1 * 2] - centerU;
                            double y1 = clippedByPatchProjElementFEUV[node1 * 2 + 1] - centerV;
                            double x2 = clippedByPatchProjElementFEUV[node2 * 2] - centerU;
                            double y2 = clippedByPatchProjElementFEUV[node2 * 2 + 1] - centerV;
                            double ang1 = atan2(y1, x1);
                            double ang2 = atan2(y2, x2);

                            // check if the corner nodes are between node1 and node2
                            if (isCounterClockWise
                                    && (ang1 < ang2 && ang1 < ang0 && ang0 < ang2
                                            || ang1 > ang2 && (ang0 > ang1 || ang0 < ang2))
                                    || (!isCounterClockWise
                                            && (ang1 > ang2 && ang1 > ang0 && ang0 > ang2
                                                    || ang1 < ang2 && (ang0 < ang1 || ang0 > ang2)))) {

                                clippedByPatchProjElementFEUV[node1 * 2 + 2] =
                                        patchCornerNodesUV[patchCornerNodesInsideElemIndex[nodeCount]
                                                * 2];
                                clippedByPatchProjElementFEUV[node1 * 2 + 3] =
                                        patchCornerNodesUV[patchCornerNodesInsideElemIndex[nodeCount]
                                                * 2 + 1];
                                double nodeXYZ[3];
                                double normalVec[3];

                                thePatch->computeCartesianCoordinatesAndNormalVector(nodeXYZ,
                                        normalVec, clippedByPatchProjElementFEUV[node1 * 2 + 2],
                                        clippedByPatchProjElementFEUV[node1 * 2 + 3]);

                                double coordFE[2];
                                if (numNodesElementFE == 3)
                                    IGAMortarMath::computeIntersectionBetweenLineAndTriangle(
                                            elementFEXYZ, nodeXYZ, normalVec, coordFE);
                                else
                                    IGAMortarMath::computeIntersectionBetweenLineAndQuad(
                                            elementFEXYZ, nodeXYZ, normalVec, coordFE);

                                clippedByPatchProjElementFEWZ[node1 * 2 + 2] = coordFE[0];
                                clippedByPatchProjElementFEWZ[node1 * 2 + 3] = coordFE[1];

                                break;
                            }

                            clippedByPatchProjElementFEUV[node1 * 2 + 2] =
                                    clippedByPatchProjElementFEUV[node1 * 2];
                            clippedByPatchProjElementFEUV[node1 * 2 + 3] =
                                    clippedByPatchProjElementFEUV[node1 * 2 + 1];
                            clippedByPatchProjElementFEWZ[node1 * 2 + 2] =
                                    clippedByPatchProjElementFEWZ[node1 * 2];
                            clippedByPatchProjElementFEWZ[node1 * 2 + 3] =
                                    clippedByPatchProjElementFEWZ[node1 * 2 + 1];
                        }
                    }
                    numNodesClippedByPatchProjElementFE += numPatchCornerNodesInsideElem;

                    // compute the coupling matrix for the sub-element
                    computeCouplingMatrices4ClippedByPatchProjectedElement(thePatch,
                            numNodesClippedByPatchProjElementFE, clippedByPatchProjElementFEUV,
                            clippedByPatchProjElementFEWZ, elemCount, numNodesElementFE);
                }

            } // end of loop over patches

        } // end of if in same patch

    } // end of loop over all the element

}

void IGAMortarMapper::computeCouplingMatrices4ClippedByPatchProjectedElement(
        IGAPatchSurface* _thePatch, int _numNodesClippedByPatchProjectedElement,
        double* _clippedByPatchProjElementFEUV, double* _clippedByPatchProjElementFEWZ,
        int _elemCount, int _numNodesElementFE) {

    double *knotVectorU = _thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
    double *knotVectorV = _thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();

/// 1.find the knot span which the current element located in.
//      from minSpanu to maxSpanu in U-direction, and from minSpanV to max SpanV in V-direction
    int minSpanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
            _clippedByPatchProjElementFEUV[0]);
    int minSpanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
            _clippedByPatchProjElementFEUV[1]);
    int maxSpanU = minSpanU;
    int maxSpanV = minSpanV;

    for (int nodeCount = 0; nodeCount < _numNodesClippedByPatchProjectedElement; nodeCount++) {

        int spanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
                _clippedByPatchProjElementFEUV[nodeCount * 2]);
        int spanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
                _clippedByPatchProjElementFEUV[nodeCount * 2 + 1]);
        if (spanU < minSpanU)
            minSpanU = spanU;
        if (spanU > maxSpanU)
            maxSpanU = spanU;
        if (spanV < minSpanV)
            minSpanV = spanV;
        if (spanV > maxSpanV)
            maxSpanV = spanV;
    }

    if (minSpanU == maxSpanU & minSpanV == maxSpanV)
        /* 2.1 if the whole element is located in a single knot span, set the coordinates in the linear element as
         the full triangle or full quadriliteral defined at the beginning of this function. */
        integrate(_thePatch, _numNodesClippedByPatchProjectedElement,
                _clippedByPatchProjElementFEUV, minSpanU, minSpanV, _clippedByPatchProjElementFEWZ,
                _elemCount, _numNodesElementFE);
    else {
/// 2.2  if not, cut the element by each knot span.
        for (int spanU = minSpanU; spanU <= maxSpanU; spanU++)
            for (int spanV = minSpanV; spanV <= maxSpanV; spanV++) {

                if (knotVectorU[spanU] != knotVectorU[spanU + 1]
                        && knotVectorV[spanV] != knotVectorV[spanV + 1]) {

                    IGAMortarMath::IGAPolygonClipper clipper(knotVectorU[spanU],
                            knotVectorU[spanU + 1], knotVectorV[spanV], knotVectorV[spanV + 1]);

                    vector<double*> *polygonResult = new vector<double*>; // deleted

                    if (clipper.clip(_clippedByPatchProjElementFEUV,
                            _numNodesClippedByPatchProjectedElement, polygonResult)) {

                        // number of nodes of the clipped polygon
                        int numNodesClippedByKnotSpanProjElementFE = polygonResult->size();

                        // the array to store the coordinates(on IGA patch) of the clipped sub-element
                        double ClippedByKnotSpanProjElementFEUV[numNodesClippedByKnotSpanProjElementFE
                                * 2];
                        for (int nodeCount = 0; nodeCount < numNodesClippedByKnotSpanProjElementFE;
                                nodeCount++) {
                            ClippedByKnotSpanProjElementFEUV[nodeCount * 2] =
                                    (*polygonResult)[nodeCount][0];
                            ClippedByKnotSpanProjElementFEUV[nodeCount * 2 + 1] =
                                    (*polygonResult)[nodeCount][1];
                        }

                        // the array to store the coordinates(Linear Element) of the clipped sub-element
                        double ClippedByKnotSpanProjElementFEWZ[numNodesClippedByKnotSpanProjElementFE
                                * 2];

                        /// find the local coordinates in FE(triangle or quadrilateral) from parameter of IGAPatch
                        double nodeWZ[2];
                        for (int nodeCount = 0; nodeCount < numNodesClippedByKnotSpanProjElementFE;
                                nodeCount++) {
                            if (_numNodesClippedByPatchProjectedElement == 3)
                                IGAMortarMath::computeLocalCoordsInTriangle(
                                        _clippedByPatchProjElementFEUV, (*polygonResult)[nodeCount],
                                        nodeWZ);
                            else
                                IGAMortarMath::computeLocalCoordsInQuad(
                                        _clippedByPatchProjElementFEUV, (*polygonResult)[nodeCount],
                                        nodeWZ);
                            IGAMortarMath::computeLinearCombinationValueFromVerticesValues(
                                    _numNodesClippedByPatchProjectedElement, 2,
                                    _clippedByPatchProjElementFEWZ, nodeWZ,
                                    &ClippedByKnotSpanProjElementFEWZ[nodeCount * 2]);
                        }

                        integrate(_thePatch, numNodesClippedByKnotSpanProjElementFE,
                                ClippedByKnotSpanProjElementFEUV, spanU, spanV,
                                ClippedByKnotSpanProjElementFEWZ, _elemCount, _numNodesElementFE);

                        for (int nodeCount = 0; nodeCount < polygonResult->size(); nodeCount++)
                            delete[] (*polygonResult)[nodeCount];
                        delete polygonResult;
                    }
                }
            }
    } // end of if same span
}

void IGAMortarMapper::integrate(IGAPatchSurface* _thePatch, int _numNodes, double* _polygonUV,
        int _spanU, int _spanV, double* _polygonWZ, int _elementIndex, int _numNodesElementFE) {
    /*
     * 1. Divide the polygon into several quadratures(triangle or quadriliteral) for integration
     * 2. Loop through each quadrature
     *   2.1 Choose a Gauss quadrature (triangle or quadriliteral)
     *   2.2 Loop throught each Gauss point
     *       2.2.1 compute shape functions from Gauss points
     *       2.2.2 evaluate the coordinates in IGA patch from shape functions
     *       2.2.3 evaluate the coordinates in the linear element from shape functions
     *       2.2.4 compute the shape functions in the linear element for the current integration point
     *       2.2.5 Compute the local basis functions(shape functions of IGA) and their derivatives(for Jacobian)
     *       2.2.6 Compute the Jacobian from parameter space on IGA patch to physical
     *       2.2.7 Compute the Jacobian from the canonical space to the parameter space of IGA patch
     *       2.2.8 integrate the shape function product for C_NN(Linear shape function multiply linear shape function)
     *       2.2.9 integrate the shape function product for C_NR(Linear shape function multiply IGA shape function)
     * 3. Assemble the element coupling matrix to the global coupling matrix.
     */

    assert(_polygonUV != NULL);
    assert(_polygonWZ != NULL);

    vector<double*> quadratureVecUV;
    vector<double*> quadratureVecWZ;
    vector<int> numNodesQuadrature;

/// 1. Divide the polygon into several quadratures(triangle or quadriliteral) for integration
    if (_numNodes <= 4) {
        quadratureVecUV.push_back(_polygonUV);
        quadratureVecWZ.push_back(_polygonWZ);
        numNodesQuadrature.push_back(_numNodes);
    } else {
        double centerIGA[2] = { 0, 0 };
        double centerFE[2] = { 0, 0 };
        for (int i = 0; i < _numNodes; i++) {
            centerIGA[0] += _polygonUV[i * 2] / _numNodes;
            centerIGA[1] += _polygonUV[i * 2 + 1] / _numNodes;
            centerFE[0] += _polygonWZ[i * 2] / _numNodes;
            centerFE[1] += _polygonWZ[i * 2 + 1] / _numNodes;
        }

        double *triangleIGA;
        double *triangleFE;
        for (int i = 0; i < _numNodes; i++) {
            int iNext = (i + 1) % _numNodes;
            triangleIGA = new double[6]; // deleted
            triangleFE = new double[6]; // deleted
            triangleIGA[0] = _polygonUV[i * 2];
            triangleIGA[1] = _polygonUV[i * 2 + 1];
            triangleFE[0] = _polygonWZ[i * 2];
            triangleFE[1] = _polygonWZ[i * 2 + 1];

            triangleIGA[2] = _polygonUV[iNext * 2];
            triangleIGA[3] = _polygonUV[iNext * 2 + 1];
            triangleFE[2] = _polygonWZ[iNext * 2];
            triangleFE[3] = _polygonWZ[iNext * 2 + 1];

            triangleIGA[4] = centerIGA[0];
            triangleIGA[5] = centerIGA[1];
            triangleFE[4] = centerFE[0];
            triangleFE[5] = centerFE[1];

            quadratureVecUV.push_back(triangleIGA);
            quadratureVecWZ.push_back(triangleFE);
            numNodesQuadrature.push_back(3);
        }
    }

    int pDegree = _thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = _thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
    int nShapeFuncsIGA = (pDegree + 1) * (qDegree + 1);

    double elementCouplingMatrixNN[_numNodesElementFE * (_numNodesElementFE + 1) / 2];
    double elementCouplingMatrixNR[nShapeFuncsIGA * _numNodesElementFE];

    for (int arrayIndex = 0; arrayIndex < _numNodesElementFE * (_numNodesElementFE + 1) / 2;
            arrayIndex++)
        elementCouplingMatrixNN[arrayIndex] = 0;

    for (int arrayIndex = 0; arrayIndex < nShapeFuncsIGA * _numNodesElementFE; arrayIndex++)
        elementCouplingMatrixNR[arrayIndex] = 0;

    /// 2. Loop through each quadrature
    for (int quadratureCount = 0; quadratureCount < quadratureVecUV.size(); quadratureCount++) {

/// 2.1 Choose gauss triangle or gauss quadriliteral

        IGAMortarMath::GaussQuadrature *theGaussQuadrature;

        int nNodesQuadrature;
        if (numNodesQuadrature[quadratureCount] == 3) {
            theGaussQuadrature = gaussTriangle;
            nNodesQuadrature = 3;
        } else {
            theGaussQuadrature = gaussQuad;
            nNodesQuadrature = 4;
        }

        double *quadratureUV = quadratureVecUV[quadratureCount];
        double *quadratureWZ = quadratureVecWZ[quadratureCount];

/// 2.2 Loop throught each Gauss point
        for (int GPCount = 0; GPCount < theGaussQuadrature->numGaussPoints; GPCount++) {

            /// 2.2.1 compute shape functions from Gauss points(in the quadrature).
            const double *GP = theGaussQuadrature->getGaussPoint(GPCount);

            double shapeFuncs[nNodesQuadrature];
            IGAMortarMath::computeLowOrderShapeFunc(nNodesQuadrature, GP, shapeFuncs);

            /// 2.2.2 evaluate the coordinates in IGA patch from shape functions
            double GPIGA[2];
            IGAMortarMath::computeLinearCombination(nNodesQuadrature, 2, quadratureUV, shapeFuncs,
                    GPIGA);

            /// 2.2.3 evaluate the coordinates in the linear element from shape functions
            double GPFE[2];
            IGAMortarMath::computeLinearCombination(nNodesQuadrature, 2, quadratureWZ, shapeFuncs,
                    GPFE);

            /// 2.2.4 compute the shape function(in the linear element) of the current integration point
            double shapeFuncsFE[_numNodesElementFE];
            IGAMortarMath::computeLowOrderShapeFunc(_numNodesElementFE, GPFE, shapeFuncsFE);
            int derivDegree = 1;

            /// 2.2.5 Compute the local basis functions(shape functions of IGA) and their derivatives(for Jacobian)
            double localBasisFunctionsAndDerivatives[(derivDegree + 1) * (derivDegree + 2)
                    * nShapeFuncsIGA / 2];

            _thePatch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                    localBasisFunctionsAndDerivatives, derivDegree, GPIGA[0], _spanU, GPIGA[1],
                    _spanV);

            /// 2.2.6 Compute the Jacobian from parameter space on IGA patch to physical
            double baseVectors[6];
            _thePatch->computeBaseVectors(baseVectors, localBasisFunctionsAndDerivatives, _spanU,
                    _spanV);

            double JacobianUVToPhysical = IGAMortarMath::computeAreaTriangle(baseVectors[0],
                    baseVectors[1], baseVectors[2], baseVectors[3], baseVectors[4], baseVectors[5])
                    * 2;

            /// 2.2.7 Compute the Jacobian from the canonical space to the parameter space of IGA patch
            double JacobianCanonicalToUV;
            if (nNodesQuadrature == 3) {
                JacobianCanonicalToUV = IGAMortarMath::computeAreaTriangle(
                        quadratureUV[2] - quadratureUV[0], quadratureUV[3] - quadratureUV[1], 0,
                        quadratureUV[4] - quadratureUV[0], quadratureUV[5] - quadratureUV[1], 0);
            } else {
                double dudx = .25
                        * (-(1 - GP[2]) * quadratureUV[0] + (1 - GP[2]) * quadratureUV[2]
                                + (1 + GP[2]) * quadratureUV[4] - (1 + GP[2]) * quadratureUV[6]);
                double dudy = .25
                        * (-(1 - GP[1]) * quadratureUV[0] - (1 + GP[1]) * quadratureUV[2]
                                + (1 + GP[1]) * quadratureUV[4] + (1 - GP[1]) * quadratureUV[6]);
                double dvdx = .25
                        * (-(1 - GP[2]) * quadratureUV[1] + (1 - GP[2]) * quadratureUV[3]
                                + (1 + GP[2]) * quadratureUV[5] - (1 + GP[2]) * quadratureUV[7]);
                double dvdy = .25
                        * (-(1 - GP[1]) * quadratureUV[1] - (1 + GP[1]) * quadratureUV[3]
                                + (1 + GP[1]) * quadratureUV[5] + (1 - GP[1]) * quadratureUV[7]);
                JacobianCanonicalToUV = fabs(dudx * dvdy - dudy * dvdx);
            }
            double Jacobian = JacobianUVToPhysical * JacobianCanonicalToUV;

            /// 2.2.8 integrate the shape function product for C_NN(Linear shape function multiply linear shape function)
            int count = 0;
            for (int i = 0; i < _numNodesElementFE; i++)
                for (int j = i; j < _numNodesElementFE; j++) {
                    elementCouplingMatrixNN[count++] += shapeFuncsFE[i] * shapeFuncsFE[j] * Jacobian
                            * theGaussQuadrature->weights[GPCount];
                }

            /// 2.2.9 integrate the shape function product for C_NR(Linear shape function multiply IGA shape function)
            count = 0;
            for (int i = 0; i < _numNodesElementFE; i++)
                for (int j = 0; j < nShapeFuncsIGA; j++)
                    elementCouplingMatrixNR[count++] +=
                            shapeFuncsFE[i]
                                    * localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                                            1, 0, 0, j)] * Jacobian
                                    * theGaussQuadrature->weights[GPCount];

        }
    }

/// 3.Assemble the element coupling matrix to the global coupling matrix.
    int count = 0;
    for (int i = 0; i < _numNodesElementFE; i++)
        for (int j = i; j < _numNodesElementFE; j++) {
            int dof1 = meshFEDirectElemTable[_elementIndex][i];
            int dof2 = meshFEDirectElemTable[_elementIndex][j];

            if (dof1 < dof2)
                (*C_NN)(dof1, dof2) += elementCouplingMatrixNN[count++];
            else
                (*C_NN)(dof2, dof1) += elementCouplingMatrixNN[count++];
        }

    int dofIGA[nShapeFuncsIGA];
    _thePatch->getIGABasis()->getBasisFunctionsIndex(_spanU, _spanV, dofIGA);

    for (int i = 0; i < nShapeFuncsIGA; i++)
        dofIGA[i] = _thePatch->getControlPointNet()[dofIGA[i]]->getDofIndex();

    count = 0;
    for (int i = 0; i < _numNodesElementFE; i++)
        for (int j = 0; j < nShapeFuncsIGA; j++)
            (*C_NR)(meshFEDirectElemTable[_elementIndex][i], dofIGA[j]) +=
                    elementCouplingMatrixNR[count++];

    if (_numNodes > 4)
        for (int i = 0; i < quadratureVecUV.size(); i++) {
            delete quadratureVecUV[i];
            delete quadratureVecWZ[i];
        }
}

void IGAMortarMapper::consistentMapping(const double* _fieldIGA, double* _fieldFE) {
    /*
     * Mapping of the
     * C_NN * x_FE = C_NR * x_IGA
     */
    double tmpVec[meshFE->numNodes];

// 1. matrix vector product (x_tmp = C_NR * x_IGA)
    C_NR->mulitplyVec(_fieldIGA, tmpVec, meshFE->numNodes);

// 2. solve C_NN * x_FE = x_tmp
    C_NN->solve(_fieldFE, tmpVec);

//    ofstream myfile;
//    myfile.open("D_tmp", ios::app);
//    myfile.precision(14);
//    myfile << std::dec;
//    myfile << _fieldIGA[33] << "\n";
//    myfile.close();

}

void IGAMortarMapper::conservativeMapping(const double* _fieldFE, double* _fieldIGA) {
    /*
     * Mapping of the
     * f_IGA = (C_NN^(-1) * C_NR)^T * f_FE
     */

    double tmpVec[meshFE->numNodes];

// 1. solve C_NN * f_tmp = f_FE;
    C_NN->solve(tmpVec, const_cast<double *>(_fieldFE));

// 2. matrix vector product (f_IGA = C_NR^T * f_tmp)
    C_NR->transposeMulitplyVec(tmpVec, _fieldIGA, meshFE->numNodes);

//    ofstream myfile;
//    myfile.open("F_tmp", ios::app);
//    myfile.precision(14);
//    myfile << std::dec;
//    myfile << _fieldIGA[36] << "\n";
//    myfile.close();
}

void IGAMortarMapper::printCouplingMatrices() {

    ERROR_OUT() << "C_NN" << endl;
    C_NN->printCSR();
    ERROR_OUT() << "C_NR" << endl;
    C_NR->printCSR();
}

}

