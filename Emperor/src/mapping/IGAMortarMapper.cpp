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
#include <stdlib.h>
#include <math.h>

using namespace std;

namespace EMPIRE {

IGAMortarMapper::IGAMortarMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE, double _disTol, int _numGPsTri, int _numGPsQuad):
        meshIGA(_meshIGA), meshFE(_meshFE), disTol(_disTol), numGPsTri(_numGPsTri), numGPsQuad(_numGPsQuad) {

    assert(_meshIGA != NULL);
    assert(_meshFE != NULL);
    assert(_meshIGA->type == EMPIRE_Mesh_IGAMesh);
    assert(_meshFE->type == EMPIRE_Mesh_FEMesh);

    cout << std::dec;

    projectedCoords = new vector<map<int, double*> >(meshFE->numNodes);

    C_NR = new MathLibrary::SparseMatrix<double>(meshFE->numNodes,
            (unsigned long) meshIGA->getNumControlPoints());
    C_NN = new MathLibrary::SparseMatrix<double>(meshFE->numNodes, true);

    gaussTriangle = new IGAMortarMath::GaussQuadratureOnTriangle(numGPsTri);
    gaussQuad = new IGAMortarMath::GaussQuadratureOnQuad(numGPsQuad);

    initTables();
    projectPointsToSurface();
    computeCouplingMatrices();
    C_NN->factorize();
//    printCouplingMatrices();

}

IGAMortarMapper::~IGAMortarMapper() {

//  delete[] projectedCoords;
    for (int i = 0; i < meshFE->numElems; i++)
        delete[] meshFEDirectElemTable[i];
    delete[] meshFEDirectElemTable;
    delete C_NR;
    C_NN->cleanPardiso();
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
        for (int j = 0; j < numNodesPerElem; j++)
            meshFEDirectElemTable[i][j] = meshFENodesMap->at(meshFE->elems[count + j]);
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

                thePatch->findNearestKnotIntersection(initialU, initialV, cartesianCoords);
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
                    bool convergeInside = thePatch->computePointProjectionOnPatch(projectedU, projectedV, cartesianCoords, converge);

                    if (convergeInside && IGAMortarMath::computePointDistance(&meshFE->nodes[nodeIndex * 3], cartesianCoords) < disTol) {

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
}

void IGAMortarMapper::computeCouplingMatrices() {
    /*
     * Loop through the each element on the fluid side(linear element side).
     * 1. find the knot span which the current element located in.
     * 2. Check whether the whole element is located in a single knot span
     *   2.1  if so, integrate it directly.
     *   2.2  if not,  go through each knot span
     *      2.2.1 clip the element by the window defined by the current knot span
     *      2.2.2 calculate coordinates of each node of the clipped sub-element
     *      2.2.3 integrate in the sub element
     */

    // the coordinates of nodes of a complete element.
    double polygonFullTriangle[6] = { 0, 0, 1, 0, 0, 1};
    double polygonFullQuadriliteral[8] = { -1, -1, 1, -1, 1, 1, -1, 1 };

    for (int elemCount = 0; elemCount < meshFE->numElems; elemCount++) {

        // number of shape functions. Depending on number of nodes in the current element
        int nShapeFuncsFE = meshFE->numNodesPerElem[elemCount];
        double* elementFE;
        if (nShapeFuncsFE == 3)
            elementFE = polygonFullTriangle;
        else
            elementFE = polygonFullQuadriliteral;

        IGAPatchSurface* thePatch = NULL;
        int patchIndex = 0;
        bool isAllNodesInPatch = true;

        // Loop over all patches to find out that if there exist a patch that all nodes in the current element can be projected to
        for (int patchCount = 0; patchCount < meshIGA->getSurfacePatches().size(); patchCount++) {

            isAllNodesInPatch = true;
            for (int nodeCount = 0; nodeCount < nShapeFuncsFE; nodeCount++) {
                int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];
                if ((*projectedCoords)[nodeIndex].find(patchCount) == (*projectedCoords)[nodeIndex].end()) {
                    isAllNodesInPatch = false;
                    break;
                }
            }
            if (isAllNodesInPatch) {
                thePatch = meshIGA->getSurfacePatches()[patchCount];
                patchIndex = patchCount;
                break;
            }
        }

        if (isAllNodesInPatch) {
            // All nodes of the current element can be projected one patch
            double elementInPatchIGA[nShapeFuncsFE * 2];
            for (int nodeCount = 0; nodeCount < nShapeFuncsFE; nodeCount ++){
                int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];
                elementInPatchIGA[nodeCount * 2] =  (*projectedCoords)[nodeIndex][patchIndex][0];
                elementInPatchIGA[nodeCount * 2 + 1] =  (*projectedCoords)[nodeIndex][patchIndex][1];
            }
            if (nShapeFuncsFE == 3)
                computeCouplingMatricesInPatch(thePatch, nShapeFuncsFE, elementInPatchIGA, polygonFullTriangle, elemCount, nShapeFuncsFE);
            else
                computeCouplingMatricesInPatch(thePatch, nShapeFuncsFE, elementInPatchIGA, polygonFullQuadriliteral, elemCount, nShapeFuncsFE);

        } else {
            // The element need to be subdivided by patches

            for (int patchCount = 0; patchCount < meshIGA->getSurfacePatches().size(); patchCount++){
                IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchCount];
                int numNodesInPatch = 0;
                double elementInPatchIGA[16];       // the parameter coordinates in IGA of the sub-element divided by the patch
                double elementInPatchFE[16];        // the parameter coordinates in low order element of the sub-element divided by the patch
                bool isEdge[8];             // if the line segment connecting the node with the next one is a edge in the sub-element
                for (int i = 0; i < 8; i++) isEdge[i] = false;

                // Find nodes and edges of low order element inside IGA Patch
                for (int nodeCount = 0; nodeCount < nShapeFuncsFE; nodeCount++) {
                    int nodeCountNext = (nodeCount + 1) % nShapeFuncsFE;
                    int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];
                    int nodeIndexNext = meshFEDirectElemTable[elemCount][nodeCountNext];

                    if ((*projectedCoords)[nodeIndex].find(patchCount) != (*projectedCoords)[nodeIndex].end()) {
                        // if the node is inside the patch, put it into the sub-element
                        elementInPatchIGA[numNodesInPatch * 2] = (*projectedCoords)[nodeIndex][patchCount][0];
                        elementInPatchIGA[numNodesInPatch * 2 + 1] = (*projectedCoords)[nodeIndex][patchCount][1];
                        elementInPatchFE[numNodesInPatch * 2] = elementFE[nodeCount * 2];
                        elementInPatchFE[numNodesInPatch * 2 + 1] = elementFE[nodeCount * 2 + 1];
                        isEdge[numNodesInPatch] = true;
                        numNodesInPatch++ ;

                        if ((*projectedCoords)[nodeIndexNext].find(patchCount) == (*projectedCoords)[nodeIndexNext].end()){
                            // if the node is inside and the next node is outside the patch, find the intersection with patch boundary
                            // and put it into the sub-element
                            double u,v,div,dis;
                            double* P1 = &(meshFE->nodes[nodeIndex * 3]);
                            double* P2 = &(meshFE->nodes[nodeIndexNext * 3]);

                            if (thePatch->computePointProjectionOnPatchBoundary(u, v, div, dis, P1, P2)){
                                if (dis <= disTol) {
                                    elementInPatchIGA[numNodesInPatch * 2] = u;
                                    elementInPatchIGA[numNodesInPatch * 2 + 1] = v;
                                    double P1x = elementFE[nodeCount * 2];
                                    double P1y = elementFE[nodeCount * 2 + 1];
                                    double P2x = elementFE[nodeCountNext * 2];
                                    double P2y = elementFE[nodeCountNext * 2 + 1];
                                    elementInPatchFE[numNodesInPatch * 2] = P1x * (1 - div) + P2x * div;
                                    elementInPatchFE[numNodesInPatch * 2 + 1] = P1y * (1 - div) + P2y * div;

                                    numNodesInPatch ++;
                                } else {
                                    assert(0);
                                }
                            } else
                                assert(0);
                        } else {                }

                    } else if ((*projectedCoords)[nodeIndexNext].find(patchCount) != (*projectedCoords)[nodeIndexNext].end()) {
                        // if this node is outside and the next node is inside, find the intersection with patch boundary
                        // and put it into the sub-element
                        double u,v,div,dis;
                        double* P1 = &(meshFE->nodes[nodeIndex * 3]);
                        double* P2 = &(meshFE->nodes[nodeIndexNext * 3]);
                        if (thePatch->computePointProjectionOnPatchBoundary(u, v, div, dis, P1, P2))
                            if (dis <= disTol) {
                                elementInPatchIGA[numNodesInPatch * 2] = u;
                                elementInPatchIGA[numNodesInPatch * 2 + 1] = v;
                                double P1x = elementFE[nodeCount * 2];
                                double P1y = elementFE[nodeCount * 2 + 1];
                                double P2x = elementFE[nodeCountNext * 2];
                                double P2y = elementFE[nodeCountNext * 2 + 1];
                                elementInPatchFE[numNodesInPatch * 2] = P1x * (1 - div) + P2x * div;
                                elementInPatchFE[numNodesInPatch * 2 + 1] = P1y * (1 - div) + P2y * div;
                                isEdge[numNodesInPatch] = true;
                                elementInPatchFE[numNodesInPatch * 2] = P1x * (1 - div) + P2x * div;
                                elementInPatchFE[numNodesInPatch * 2 + 1] = P1y * (1 - div) + P2y * div;

                                numNodesInPatch ++;
                            } else
                                assert(0);
                        else
                            assert(0);
                    }
                }

                // Find of corner Point of IGA Patch inside low order element
                if (numNodesInPatch > 2) {
                    // Find direction of FE element(Clock or Counter-Clock)
                    double counterClockWise = 0.0;
                    double x1 = elementInPatchIGA[2] - elementInPatchIGA[0];
                    double y1 = elementInPatchIGA[3] - elementInPatchIGA[1];
                    double x2 = elementInPatchIGA[4] - elementInPatchIGA[2];
                    double y2 = elementInPatchIGA[5] - elementInPatchIGA[3];
                    if (IGAMortarMath::computeCrossProduct2D(x1,y1,x2,y2) > 0)
                        counterClockWise = 1.0;
                    else
                        counterClockWise = -1.0;

                    // IGA Patch
                    double *knotVectorU = thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
                    double *knotVectorV = thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();
                    int numKnotsU = thePatch->getIGABasis()->getUBSplineBasis1D()->getNoKnots();
                    int numKnotsV = thePatch->getIGABasis()->getVBSplineBasis1D()->getNoKnots();
                    double u0 = knotVectorU[0];
                    double uEnd = knotVectorU[numKnotsU - 1];
                    double v0 = knotVectorV[0];
                    double vEnd = knotVectorV[numKnotsV - 1];
                    // Corner nodes of the IGA Patch
                    double patchNodesIGA[8] = {u0,v0,uEnd,v0,uEnd,vEnd,u0,vEnd};
                    // Corner nodes of the IGA Patch which inside the current element
                    int patchNodesInsideElem[4];
                    int numPatchNodesInsideElem = 0;

                    // Loop over the 4 corner node of the IGA Patch to check which of them are in side the current element
                    for (int nodeCount = 0; nodeCount < 4; nodeCount ++){
                        bool isNodeInsideElem = true;
                        for (int edgeCount = 0; edgeCount < numNodesInPatch; edgeCount ++){
                            if (isEdge[edgeCount]){
                                int node1 = edgeCount;
                                int node2 = (edgeCount + 1) % numNodesInPatch;

                                x1 = elementInPatchIGA[node2 * 2] - elementInPatchIGA[node1 * 2];
                                y1 = elementInPatchIGA[node2 * 2 + 1] - elementInPatchIGA[node1 * 2 + 1];
                                x2 = patchNodesIGA[nodeCount * 2] - elementInPatchIGA[node2 * 2];
                                y2 = patchNodesIGA[nodeCount * 2 + 1] - elementInPatchIGA[node2 * 2 + 1];
                                if (IGAMortarMath::computeCrossProduct2D(x1,y1,x2,y2) * counterClockWise < 0){
                                    isNodeInsideElem = false;
                                    break;
                                }
                            }
                        }
                        if (isNodeInsideElem){
                            patchNodesInsideElem[numPatchNodesInsideElem] = nodeCount;
                            numPatchNodesInsideElem ++;
                        }
                    }

                    // center of all nodes inside element
                    double centerX = 0.0;
                    double centerY = 0.0;
                    for (int nodeCount = 0; nodeCount < numPatchNodesInsideElem; nodeCount ++) {
                        centerX += patchNodesIGA[patchNodesInsideElem[nodeCount] * 2];
                        centerY += patchNodesIGA[patchNodesInsideElem[nodeCount] * 2 + 1];
                    }
                    for (int nodeCount = 0; nodeCount < numNodesInPatch; nodeCount ++){
                        centerX += elementInPatchIGA[nodeCount * 2];
                        centerY += elementInPatchIGA[nodeCount * 2 + 1];
                    }
                    centerX /= numPatchNodesInsideElem + numNodesInPatch;
                    centerY /= numPatchNodesInsideElem + numNodesInPatch;

                    double ang[numNodesInPatch];
                    for (int nodeCount = 0; nodeCount < numNodesInPatch; nodeCount ++){
                        double x0 = elementInPatchIGA[nodeCount * 2] - centerX;
                        double y0 = elementInPatchIGA[nodeCount * 2 + 1] - centerY;
                        ang[nodeCount] = atan2(y0,x0);
                    }

                    // Cartesian coordinates of the low order element
                    double coordElementFE[12];
                    for (int i = 0; i < nShapeFuncsFE; i++){
                        int nodeIndex = meshFEDirectElemTable[elemCount][i];
                        for (int j = 0; j < 3; j ++)
                            coordElementFE[i * 3 + j] = meshFE->nodes[nodeIndex * 3 + j];
                    }

                    // Loop over all corner nodes of IGA Patch inside the element and put them into the sub-element in a correct order
                    for (int nodeCount = 0; nodeCount < numPatchNodesInsideElem; nodeCount ++){
                        double x0 = patchNodesIGA[patchNodesInsideElem[nodeCount] * 2] - centerX;
                        double y0 = patchNodesIGA[patchNodesInsideElem[nodeCount] * 2 + 1] - centerY;
                        double ang0 = atan2(y0,x0);

                        // Loop over all edges(nodes) of the sub-element to find out the correct position for the corner nodes of the patch
                        for (int edgeCount = numNodesInPatch + nodeCount - 1; edgeCount >= 0 ; edgeCount --){
                            int node1 = edgeCount;
                            int node2 = (edgeCount + 2) % (numNodesInPatch + nodeCount + 1);
                            double x1 = elementInPatchIGA[node1 * 2] - centerX;
                            double y1 = elementInPatchIGA[node1 * 2 + 1] - centerY;
                            double x2 = elementInPatchIGA[node2 * 2] - centerX;
                            double y2 = elementInPatchIGA[node2 * 2 + 1] - centerY;
                            double ang1 = atan2(y1,x1);
                            double ang2 = atan2(y2,x2);

                            // check if the corner nodes are between node1 and node2
                            if ( counterClockWise > 0 && (ang1 < ang2 && ang1 < ang0 && ang0 < ang2
                                ||  ang1 > ang2 && (ang0 > ang1 || ang0 < ang2))
                            || ( counterClockWise < 0 && (ang1 > ang2 && ang1 > ang0 && ang0 > ang2
                                ||  ang1 < ang2 && (ang0 < ang1 || ang0 > ang2)))) {


                                elementInPatchIGA[node1 * 2 + 2] = patchNodesIGA[patchNodesInsideElem[nodeCount] * 2];
                                elementInPatchIGA[node1 * 2 + 3] = patchNodesIGA[patchNodesInsideElem[nodeCount] * 2 + 1];
                                double coord[3];
                                double normal[3];

                                thePatch->computeCartesianCoordinatesAndNormalVector(coord, normal,
                                        elementInPatchIGA[node1 * 2 + 2], elementInPatchIGA[node1 * 2 + 3]);


                                double coordFE[2];
                                if (nShapeFuncsFE == 3)
                                    IGAMortarMath::computeIntersectionBetweenLineAndTriangle(coordElementFE, coord, normal, coordFE);
                                else
                                    IGAMortarMath::computeIntersectionBetweenLineAndQuad(coordElementFE, coord, normal, coordFE);

                                elementInPatchFE[node1 * 2 + 2] = coordFE[0];
                                elementInPatchFE[node1 * 2 + 3] = coordFE[1];

                                break;
                            }

                            elementInPatchIGA[node1 * 2 + 2] = elementInPatchIGA[node1 * 2];
                            elementInPatchIGA[node1 * 2 + 3] = elementInPatchIGA[node1 * 2 + 1];
                            elementInPatchFE[node1 * 2 + 2] = elementInPatchFE[node1 * 2];
                            elementInPatchFE[node1 * 2 + 3] = elementInPatchFE[node1 * 2 + 1];
                        }
                    }
                    numNodesInPatch += numPatchNodesInsideElem;

                    // compute the coupling matrix for the sub-element
                    // if the sub-element has more then 4 nodes, sub divide it into triangles
                    if (numNodesInPatch <= 4)
                        computeCouplingMatricesInPatch(thePatch, numNodesInPatch, elementInPatchIGA,
                                elementInPatchFE, elemCount, nShapeFuncsFE);
                    else {
                        double coord[3];
                        double normal[3];
                        thePatch->computeCartesianCoordinatesAndNormalVector(coord, normal,centerX, centerY);
                        double coordFE[2];
                        if (nShapeFuncsFE == 3)
                            IGAMortarMath::computeIntersectionBetweenLineAndTriangle(coordElementFE, coord, normal, coordFE);
                        else
                            IGAMortarMath::computeIntersectionBetweenLineAndQuad(coordElementFE, coord, normal, coordFE);

                        double elementInPatchIGATri[6];
                        elementInPatchIGATri[0] = centerX;
                        elementInPatchIGATri[1] = centerY;
                        double elementInPatchFETri[6];
                        elementInPatchFETri[0] = coordFE[0];
                        elementInPatchFETri[1] = coordFE[1];
                        for (int i = 0; i < numNodesInPatch; i++){
                            int nextI = (i + 1) % numNodesInPatch;
                            elementInPatchIGATri[2] = elementInPatchIGA[i * 2];
                            elementInPatchIGATri[3] = elementInPatchIGA[i * 2 + 1];
                            elementInPatchIGATri[4] = elementInPatchIGA[nextI * 2];
                            elementInPatchIGATri[5] = elementInPatchIGA[nextI * 2 + 1];

                            elementInPatchFETri[2] = elementInPatchFE[i * 2];
                            elementInPatchFETri[3] = elementInPatchFE[i * 2 + 1];
                            elementInPatchFETri[4] = elementInPatchFE[nextI * 2];
                            elementInPatchFETri[5] = elementInPatchFE[nextI * 2 + 1];
                            computeCouplingMatricesInPatch(thePatch, 3, elementInPatchIGATri,
                                    elementInPatchFETri, elemCount, nShapeFuncsFE);
                        }
                    }
                } else if (numNodesInPatch > 0)
                    assert(0);

            } // end of loop over patches

        } // end of if in same patch

    } // end of loop over all the element

}

void IGAMortarMapper::computeCouplingMatricesInPatch(IGAPatchSurface* _thePatch, int _numNodesInPatch, double* _elementInPatchIGA,
        double* _elementInPatchFE, int _elemCount, int _nShapeFuncsFE){

    double *knotVectorU = _thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
    double *knotVectorV = _thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();

    /// 1.find the knot span which the current element located in.
    //      from minSpanu to maxSpanu in U-direction, and from minSpanV to max SpanV in V-direction
    int minSpanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
            _elementInPatchIGA[0]);
    int minSpanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
            _elementInPatchIGA[1]);
    int maxSpanU = minSpanU;
    int maxSpanV = minSpanV;

    for (int nodeCount = 0; nodeCount < _numNodesInPatch; nodeCount++) {

        int spanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
                _elementInPatchIGA[nodeCount * 2]);
        int spanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
                _elementInPatchIGA[nodeCount * 2 + 1]);
        if (spanU < minSpanU)            minSpanU = spanU;
        if (spanU > maxSpanU)            maxSpanU = spanU;
        if (spanV < minSpanV)            minSpanV = spanV;
        if (spanV > maxSpanV)            maxSpanV = spanV;
    }

    if (minSpanU == maxSpanU & minSpanV == maxSpanV)

        /* 2.1 if the whole element is located in a single knot span, set the coordinates in the linear element as
         the full triangle or full quadriliteral defined at the beginning of this function. */
        integrate(_thePatch, _numNodesInPatch, _elementInPatchIGA, minSpanU, minSpanV, _elementInPatchFE, _elemCount, _nShapeFuncsFE);

    else {
        /// 2.2  if not, cut the element through each knot span.
        for (int spanU = minSpanU; spanU <= maxSpanU; spanU++)
            for (int spanV = minSpanV; spanV <= maxSpanV; spanV++) {

                IGAMortarMath::IGAPolygonClipper clipper(knotVectorU[spanU],
                        knotVectorU[spanU + 1], knotVectorV[spanV], knotVectorV[spanV + 1]);

                vector<double*> *polygonResult = new vector<double*>; // deleted

                if (clipper.clip(_elementInPatchIGA, _numNodesInPatch, polygonResult)) {

                    // number of nodes of the clipped polygon
                    int numNodesInSpan = polygonResult->size();

                    // the array to store the coordinates(on IGA patch) of the clipped sub-element
                    double elementInSpanIGA[numNodesInSpan * 2];
                    for (int nodeCount = 0; nodeCount < numNodesInSpan; nodeCount++) {
                        elementInSpanIGA[nodeCount * 2] = (*polygonResult)[nodeCount][0];
                        elementInSpanIGA[nodeCount * 2 + 1] = (*polygonResult)[nodeCount][1];
                    }

                    // the array to store the coordinates(Linear Element) of the clipped sub-element
                    double elementInSpanFE[numNodesInSpan * 2];

                    /// find the local coordinates in FE(triangle or quadrilateral) from parameter of IGAPatch
                    double localCoordsFE[2];
                    for (int nodeCount = 0; nodeCount < numNodesInSpan; nodeCount++) {
                        if (_numNodesInPatch == 3)
                            IGAMortarMath::computeLocalCoordsInTriangle(_elementInPatchIGA,
                                    (*polygonResult)[nodeCount], localCoordsFE);
                        else
                            IGAMortarMath::computeLocalCoordsInQuad(_elementInPatchIGA,
                                    (*polygonResult)[nodeCount], localCoordsFE);
                        IGAMortarMath::computeLinearCombinationValueFromVerticesValues(_numNodesInPatch, 2, _elementInPatchFE,
                                localCoordsFE, &elementInSpanFE[nodeCount * 2]);
                    }

                    integrate(_thePatch, numNodesInSpan, elementInSpanIGA, spanU, spanV,
                            elementInSpanFE, _elemCount, _nShapeFuncsFE);

                    for (int nodeCount = 0; nodeCount < polygonResult->size(); nodeCount++)
                        delete[] (*polygonResult)[nodeCount];
                    delete polygonResult;
                }
            }
    } // end of if same span
}

void IGAMortarMapper::integrate(IGAPatchSurface* _thePatch, int _numNodes, double* _polygonIGA,
        int _spanU, int _spanV, double* _polygonFE, int _elementIndex, int _nShapeFuncsFE) {
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


    assert(_polygonIGA != NULL);
    assert(_polygonFE != NULL);

    vector<double*> quadratureIGA;
    vector<double*> quadratureFE;
    vector<int> numNodesQuadrature;

    /// 1. Divide the polygon into several quadratures(triangle or quadriliteral) for integration
    if (_numNodes <= 4) {
        quadratureIGA.push_back(_polygonIGA);
        quadratureFE.push_back(_polygonFE);
        numNodesQuadrature.push_back(_numNodes);
    } else {
        double centerIGA[2] = { 0, 0 };
        double centerFE[2] = { 0, 0 };
        for (int i = 0; i < _numNodes; i++) {
            centerIGA[0] += _polygonIGA[i * 2] / _numNodes;
            centerIGA[1] += _polygonIGA[i * 2 + 1] / _numNodes;
            centerFE[0] += _polygonFE[i * 2] / _numNodes;
            centerFE[1] += _polygonFE[i * 2 + 1] / _numNodes;
        }

        double *triangleIGA;
        double *triangleFE;
        for (int i = 0; i < _numNodes; i++) {
            int iNext = (i + 1) % _numNodes;
            triangleIGA = new double[6]; // deleted
            triangleFE = new double[6]; // deleted
            triangleIGA[0] = _polygonIGA[i * 2];
            triangleIGA[1] = _polygonIGA[i * 2 + 1];
            triangleFE[0] = _polygonFE[i * 2];
            triangleFE[1] = _polygonFE[i * 2 + 1];

            triangleIGA[2] = _polygonIGA[iNext * 2];
            triangleIGA[3] = _polygonIGA[iNext * 2 + 1];
            triangleFE[2] = _polygonFE[iNext * 2];
            triangleFE[3] = _polygonFE[iNext * 2 + 1];

            triangleIGA[4] = centerIGA[0];
            triangleIGA[5] = centerIGA[1];
            triangleFE[4] = centerFE[0];
            triangleFE[5] = centerFE[1];

            quadratureIGA.push_back(triangleIGA);
            quadratureFE.push_back(triangleFE);
            numNodesQuadrature.push_back(3);
        }
    }

    int pDegree = _thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = _thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
    int nShapeFuncsIGA = (pDegree + 1) * (qDegree + 1);

    double elementCouplingMatrixNN[_nShapeFuncsFE * (_nShapeFuncsFE + 1) / 2];
    double elementCouplingMatrixNR[nShapeFuncsIGA * _nShapeFuncsFE];

    for (int arrayIndex = 0; arrayIndex < _nShapeFuncsFE * (_nShapeFuncsFE + 1) / 2; arrayIndex++)
        elementCouplingMatrixNN[arrayIndex] = 0;

    for (int arrayIndex = 0; arrayIndex < nShapeFuncsIGA * _nShapeFuncsFE; arrayIndex++)
        elementCouplingMatrixNR[arrayIndex] = 0;

/// 2. Loop through each quadrature
    for (int quadratureCount = 0; quadratureCount < quadratureIGA.size(); quadratureCount++) {

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

        double *coordsIGA = quadratureIGA[quadratureCount];
        double *coordsFE = quadratureFE[quadratureCount];

        /// 2.2 Loop throught each Gauss point
        for (int GPCount = 0; GPCount < theGaussQuadrature->numGaussPoints; GPCount++) {

            /// 2.2.1 compute shape functions from Gauss points(in the quadrature).
            const double *GP = theGaussQuadrature->getGaussPoint(GPCount);

            double shapeFuncs[nNodesQuadrature];
            IGAMortarMath::computeLowOrderShapeFunc(nNodesQuadrature, GP, shapeFuncs);

            /// 2.2.2 evaluate the coordinates in IGA patch from shape functions
            double GPIGA[2];
            IGAMortarMath::computeLinearCombination(nNodesQuadrature, 2, coordsIGA, shapeFuncs, GPIGA);

            /// 2.2.3 evaluate the coordinates in the linear element from shape functions
            double GPFE[2];
            IGAMortarMath::computeLinearCombination(nNodesQuadrature, 2, coordsFE, shapeFuncs, GPFE);

            /// 2.2.4 compute the shape function(in the linear element) of the current integration point
            double shapeFuncsFE[_nShapeFuncsFE];
            IGAMortarMath::computeLowOrderShapeFunc(_nShapeFuncsFE, GPFE, shapeFuncsFE);
            int derivDegree = 1;

            /// 2.2.5 Compute the local basis functions(shape functions of IGA) and their derivatives(for Jacobian)

            double localBasisFunctionsAndDerivatives[(derivDegree + 1) * (derivDegree + 2) * nShapeFuncsIGA / 2];

            _thePatch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                    localBasisFunctionsAndDerivatives, derivDegree, GPIGA[0], _spanU, GPIGA[1], _spanV);

            /// 2.2.6 Compute the Jacobian from parameter space on IGA patch to physical
            double baseVectors[6];
            _thePatch->computeBaseVectors(baseVectors, localBasisFunctionsAndDerivatives, _spanU, _spanV);

            double JacobianUVToPhysical = IGAMortarMath::computeAreaTriangle(baseVectors[0],
                    baseVectors[1], baseVectors[2], baseVectors[3], baseVectors[4], baseVectors[5]) * 2;

            /// 2.2.7 Compute the Jacobian from the canonical space to the parameter space of IGA patch
            double JacobianCanonicalToUV;
            if (nNodesQuadrature == 3) {
                JacobianCanonicalToUV = IGAMortarMath::computeAreaTriangle(coordsIGA[2] - coordsIGA[0],
                        coordsIGA[3] - coordsIGA[1], 0, coordsIGA[4] - coordsIGA[0],
                        coordsIGA[5] - coordsIGA[1], 0);
            } else {
                double dudx = .25
                        * (-(1 - GP[2]) * coordsIGA[0] + (1 - GP[2]) * coordsIGA[2]
                                + (1 + GP[2]) * coordsIGA[4] - (1 + GP[2]) * coordsIGA[6]);
                double dudy = .25
                        * (-(1 - GP[1]) * coordsIGA[0] - (1 + GP[1]) * coordsIGA[2]
                                + (1 + GP[1]) * coordsIGA[4] + (1 - GP[1]) * coordsIGA[6]);
                double dvdx = .25
                        * (-(1 - GP[2]) * coordsIGA[1] + (1 - GP[2]) * coordsIGA[3]
                                + (1 + GP[2]) * coordsIGA[5] - (1 + GP[2]) * coordsIGA[7]);
                double dvdy = .25
                        * (-(1 - GP[1]) * coordsIGA[1] - (1 + GP[1]) * coordsIGA[3]
                                + (1 + GP[1]) * coordsIGA[5] + (1 - GP[1]) * coordsIGA[7]);
                JacobianCanonicalToUV = fabs(dudx * dvdy - dudy * dvdx);
            }
            double Jacobian = JacobianUVToPhysical * JacobianCanonicalToUV;

            /// 2.2.8 integrate the shape function product for C_NN(Linear shape function multiply linear shape function)
            int count = 0;
            for (int i = 0; i < _nShapeFuncsFE; i++)
                for (int j = i; j < _nShapeFuncsFE; j++) {
                    elementCouplingMatrixNN[count++] += shapeFuncsFE[i] * shapeFuncsFE[j] * Jacobian
                            * theGaussQuadrature->weights[GPCount];
                }

            /// 2.2.9 integrate the shape function product for C_NR(Linear shape function multiply IGA shape function)
            count = 0;
            for (int i = 0; i < _nShapeFuncsFE; i++)
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
    for (int i = 0; i < _nShapeFuncsFE; i++)
        for (int j = i; j < _nShapeFuncsFE; j++) {
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
        dofIGA[i] =
                meshIGA->getMapControlPointIDToIndex()[_thePatch->getControlPointNet()[dofIGA[i]]->getId()];

    count = 0;
    for (int i = 0; i < _nShapeFuncsFE; i++)
        for (int j = 0; j < nShapeFuncsIGA; j++)
            (*C_NR)(meshFEDirectElemTable[_elementIndex][i], dofIGA[j]) +=
                    elementCouplingMatrixNR[count++];

    if (_numNodes > 4)
        for (int i = 0; i < quadratureIGA.size(); i++) {
            delete quadratureIGA[i];
            delete quadratureFE[i];
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

}

void IGAMortarMapper::printCouplingMatrices() {

    cout << "C_NN" << endl;
    C_NN->printCSR();
    cout << "C_NR" << endl;
    C_NR->printCSR();
}

}

