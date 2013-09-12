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

IGAMortarMapper::IGAMortarMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE) :
        meshIGA(_meshIGA), meshFE(_meshFE) {

    assert(_meshIGA != NULL);
    assert(_meshFE != NULL);
    assert(_meshIGA->type == EMPIRE_Mesh_IGAMesh);
    assert(_meshFE->type == EMPIRE_Mesh_FEMesh);

    projectedCoords = new vector<map<int, double*> >(meshFE->numNodes);

    C_NR = new MathLibrary::SparseMatrix<double>(meshFE->numNodes,
            (unsigned long) meshIGA->getNumControlPoints());
    C_NN = new MathLibrary::SparseMatrix<double>(meshFE->numNodes, true);

    gaussTriangle = new IGAMortarMath::GaussQuadratureOnTriangle(numGPsTri);
    gaussQuad = new IGAMortarMath::GaussQuadratureOnQuad(numGPsQuad);

    initTables();
    projectPointsToSurface();
    computeCouplingMatrices();

}

IGAMortarMapper::~IGAMortarMapper() {

//	delete[] projectedCoords;
    for (int i = 0; i < meshFE->numElems; i++)
        delete[] meshFEDirectElemTable[i];
    delete[] meshFEDirectElemTable;
    delete C_NR;
    delete C_NN;

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

    // the coordinates of nodes of a whole element. see comment 2.1
    double polygonFullTriangle[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    double polygonFullQuadriliteral[8] = { -1, -1, 1, -1, 1, 1, -1, 1 };

    for (int elemCount = 0; elemCount < meshFE->numElems; elemCount++) {

        // number of shape functions. Depending on number of nodes in the current element
        int nShapeFuncsFE = meshFE->numNodesPerElem[elemCount];

        IGAPatchSurface* thePatch = NULL;
        int patchIndex = 0;

        for (int patchCount = 0; patchCount < meshIGA->getSurfacePatches().size(); patchCount++) {
            bool isInPatch = true;
            for (int nodeCount = 0; nodeCount < meshFE->numNodesPerElem[elemCount]; nodeCount++) {
                int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];
                if ((*projectedCoords)[nodeIndex].find(patchCount)
                        == (*projectedCoords)[nodeIndex].end()) {
                    isInPatch = false;
                    break;
                }
            }
            if (isInPatch) {
                thePatch = meshIGA->getSurfacePatches()[patchCount];
                patchIndex = patchCount;
                break;
            }
        }

        if (thePatch != NULL) {

            double *knotVectorU = thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
            double *knotVectorV = thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();

            /// 1.find the knot span which the current element located in.
            // 		from minSpanu to maxSpanu in U-direction, and from minSpanV to max SpanV in V-direction
            int minSpanU = thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
                    (*projectedCoords)[meshFEDirectElemTable[elemCount][0]][patchIndex][0]);
            int minSpanV = thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
                    (*projectedCoords)[meshFEDirectElemTable[elemCount][0]][patchIndex][1]);
            int maxSpanU = minSpanU;
            int maxSpanV = minSpanV;

            for (int nodeCount = 0; nodeCount < meshFE->numNodesPerElem[elemCount]; nodeCount++) {
                int spanU =
                        thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
                                (*projectedCoords)[meshFEDirectElemTable[elemCount][nodeCount]][patchIndex][0]);
                int spanV =
                        thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
                                (*projectedCoords)[meshFEDirectElemTable[elemCount][nodeCount]][patchIndex][1]);
                if (spanU < minSpanU)
                    minSpanU = spanU;
                if (spanU > maxSpanU)
                    maxSpanU = spanU;
                if (spanV < minSpanV)
                    minSpanV = spanV;
                if (spanV > maxSpanV)
                    maxSpanV = spanV;
            }

            // the array to store the coordinates(on IGA patch) of the current element
            double elementPolygonIGA[nShapeFuncsFE * 2];
            for (int nodeCount = 0; nodeCount < nShapeFuncsFE; nodeCount++) {
                elementPolygonIGA[nodeCount * 2] =
                        (*projectedCoords)[meshFEDirectElemTable[elemCount][nodeCount]][patchIndex][0];
                elementPolygonIGA[nodeCount * 2 + 1] =
                        (*projectedCoords)[meshFEDirectElemTable[elemCount][nodeCount]][patchIndex][1];
            }

            if (minSpanU == maxSpanU & minSpanV == maxSpanV) {

                /* 2.1 if the whole element is located in a single knot span, set the coordinates in the linear element as
                 the full triangle or full quadriliteral defined at the beginning of this function. */
                if (nShapeFuncsFE == 3) {
                    integrate(thePatch, nShapeFuncsFE, elementPolygonIGA, minSpanU, minSpanV,
                            polygonFullTriangle, elemCount, nShapeFuncsFE);
                } else {
                    integrate(thePatch, nShapeFuncsFE, elementPolygonIGA, minSpanU, minSpanV,
                            polygonFullQuadriliteral, elemCount, nShapeFuncsFE);
                }

            } else {
                /// 2.2  if not, cut the element through each knot span.
                for (int spanU = minSpanU; spanU <= maxSpanU; spanU++)
                    for (int spanV = minSpanV; spanV <= maxSpanV; spanV++) {

                        IGAMortarMath::IGAPolygonClipper clipper(knotVectorU[spanU],
                                knotVectorU[spanU + 1], knotVectorV[spanV], knotVectorV[spanV + 1]);

                        vector<double*> *polygonResult = new vector<double*>; // deleted

                        if (clipper.clip(elementPolygonIGA, nShapeFuncsFE, polygonResult)) {

                            // number of nodes of the clipped polygon
                            int numNodesClipped = polygonResult->size();

                            // the array to store the coordinates(on IGA patch) of the clipped sub-element
                            double polygonIGA[numNodesClipped * 2];
                            for (int nodeCount = 0; nodeCount < numNodesClipped; nodeCount++) {
                                polygonIGA[nodeCount * 2] = (*polygonResult)[nodeCount][0];
                                polygonIGA[nodeCount * 2 + 1] = (*polygonResult)[nodeCount][1];
                            }

                            // the array to store the coordinates(Linear Element) of the clipped sub-element
                            double polygonFE[numNodesClipped * 2];

                            /// find the local coordinates in FE(triangle or quadrilateral) from parameter of IGAPatch
                            if (nShapeFuncsFE == 3) {
                                /// for triangle element
                                double localCoordsFE[2];
                                for (int nodeCount = 0; nodeCount < numNodesClipped; nodeCount++) {
                                    IGAMortarMath::computeLocalCoordsInTriangle(elementPolygonIGA,
                                            (*polygonResult)[nodeCount], localCoordsFE);
                                    polygonFE[nodeCount * 2] = localCoordsFE[0];
                                    polygonFE[nodeCount * 2 + 1] = localCoordsFE[1];
                                }
                            } else {
                                /// for quadrilateral element
                                double localCoordsFE[2];
                                for (int nodeCount = 0; nodeCount < numNodesClipped; nodeCount++) {
                                    IGAMortarMath::computeLocalCoordsInQuad(elementPolygonIGA,
                                            (*polygonResult)[nodeCount], localCoordsFE);
                                    polygonFE[nodeCount * 2] = localCoordsFE[0];
                                    polygonFE[nodeCount * 2 + 1] = localCoordsFE[1];
                                }
                            }

                            integrate(thePatch, numNodesClipped, polygonIGA, spanU, spanV,
                                    polygonFE, elemCount, nShapeFuncsFE);

                            for (int nodeCount = 0; nodeCount < polygonResult->size(); nodeCount++)
                                delete[] (*polygonResult)[nodeCount];
                            delete polygonResult;
                        }
                    }
            } // end of if same span

        } // end of if same patch

    } // end of loop over all the element

}

void IGAMortarMapper::integrate(IGAPatchSurface* _thePatch, int _numNodes, double* _polygonIGA,
        int _spanU, int _spanV, double* _polygonFE, int _elementIndex, int _nShapeFuncsFE) {
    /*
     * 1. Divide the polygon into several quadratures(triangle or quadriliteral) for integration
     * 2. Loop through each quadrature
     *   2.1 Choose a Gauss quadrature (triangle or quadriliteral)
     *   2.2 Loop throught each Gauss point
     *   	 2.2.1 compute shape functions from Gauss points
     *   	 2.2.2 evaluate the coordinates in IGA patch from shape functions
     *   	 2.2.3 evaluate the coordinates in the linear element from shape functions
     *   	 2.2.4 compute the shape functions in the linear element for the current integration point
     *   	 2.2.5 Compute the local basis functions(shape functions of IGA) and their derivatives(for Jacobian)
     *   	 2.2.6 Compute the Jacobian from parameter space on IGA patch to physical
     *   	 2.2.7 Compute the Jacobian from the canonical space to the parameter space of IGA patch
     *   	 2.2.8 integrate the shape function product for C_NN(Linear shape function multiply linear shape function)
     *   	 2.2.9 integrate the shape function product for C_NR(Linear shape function multiply IGA shape function)
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
            triangleIGA = new double[6]; // deleted
            triangleFE = new double[6]; // deleted
            triangleIGA[0] = _polygonIGA[i * 2];
            triangleIGA[1] = _polygonIGA[i * 2 + 1];
            triangleFE[0] = _polygonFE[i * 2];
            triangleFE[1] = _polygonFE[i * 2 + 1];

            triangleIGA[2] = _polygonIGA[i * 2 + 2];
            triangleIGA[3] = _polygonIGA[i * 2 + 3];
            triangleFE[2] = _polygonFE[i * 2 + 2];
            triangleFE[3] = _polygonFE[i * 2 + 3];

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

        /// 2.2 Loop throught each Gauss point
        for (int GPCount = 0; GPCount < theGaussQuadrature->numGaussPoints; GPCount++) {

            /// 2.2.1 compute shape functions from Gauss points(in the quadrature).
            const double *GP = theGaussQuadrature->getGaussPoint(GPCount);

            double shapeFuncs[nNodesQuadrature];
            IGAMortarMath::computeLowOrderShapeFunc(nNodesQuadrature, GP, shapeFuncs);

            /// 2.2.2 evaluate the coordinates in IGA patch from shape functions
            double GPIGA[2];
            IGAMortarMath::computeLinearCombination(nNodesQuadrature, 2,
                    quadratureIGA[quadratureCount], shapeFuncs, GPIGA);

            /// 2.2.3 evaluate the coordinates in the linear element from shape functions
            double GPFE[2];
            IGAMortarMath::computeLinearCombination(nNodesQuadrature, 2,
                    quadratureFE[quadratureCount], shapeFuncs, GPFE);

            /// 2.2.4 compute the shape function(in the linear element) of the current integration point
            double shapeFuncsFE[_nShapeFuncsFE];
            IGAMortarMath::computeLowOrderShapeFunc(_nShapeFuncsFE, GPFE, shapeFuncsFE);
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
            if (_nShapeFuncsFE == 3) {
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
            int node1 = meshFEDirectElemTable[_elementIndex][i];
            int node2 = meshFEDirectElemTable[_elementIndex][j];

            if (node1 < node2)
                (*C_NN)(node1, node2) += elementCouplingMatrixNN[count++];
            else
                (*C_NN)(node2, node1) += elementCouplingMatrixNN[count++];
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
// f_IGA = (C_NN^(-1) * C_NR)^T * f_FE

    double tmpVec[meshFE->numNodes];

// 1. solve C_NN * f_tmp = f_FE;
    C_NN->solve(tmpVec, const_cast<double *>(_fieldFE));

// 2. matrix vector product (f_IGA = C_NR^T * f_tmp)
    C_NR->transposeMulitplyVec(tmpVec, _fieldIGA, meshFE->numNodes);

}

void IGAMortarMapper::projectPointsToSurface() {
    /* Loop over all the elements on the fluid side
     *   1. Loop over all nodes of the current element to check if there exist one node has been projected already
     *	 2. Check if there exist one node in the current element has been successfully projected. If not,
     *	 	2.1 if so, use result of the projected node as the initial guess
     *	 	2.2 otherwise,  find the nearest knot intersection as initial guess
     *	 3. Loop over each node at the current element
     *	      Check if the node has been already projected
     *	      if not, compute the point projection on the IGA patch using the initial guess get from the last step
     *	 The result of the projected coordinates are stored in the class member: projectedCoords
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
                    bool convergeInside = thePatch->computePointProjectionOnPatch(projectedU,
                            projectedV, cartesianCoords, converge);
                    if (IGAMortarMath::computePointDistance(&meshFE->nodes[nodeIndex * 3], cartesianCoords)
                            < disTol) {
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

void IGAMortarMapper::initTables() {
    /* using the map to store the nodeIDs
     * but here the "key" is the node ID, and the value is the position in nodeIDs
     * the map is sorted automatically, so it is efficient for searching
     */

// 1. compute direct element table for fluid mesh
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

void IGAMortarMapper::printCouplingMatrices() {
    cout << "C_NN" << endl;
    C_NN->printCSR();
    cout << "C_NR" << endl;
    C_NR->printCSR();
}

}

