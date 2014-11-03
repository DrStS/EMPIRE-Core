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
#include "MortarMath.h"
#include "MathLibrary.h"
#include "FEMesh.h"
#include "IGAPatchSurface.h"
#include "IGAMesh.h"
#include "DataField.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <set>
#include <algorithm>
#include "clipper.hpp"
using namespace std;

namespace EMPIRE {

/// Declaration statement
static const string HEADER_DECLARATION = "Author: Andreas Apostolatos";

IGAMortarMapper::IGAMortarMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE,
        double _disTol, int _numGPsTri, int _numGPsQuad, bool _isMappingIGA2FEM) :
        name(_name), meshIGA(_meshIGA), disTol(_disTol), numGPsTri(_numGPsTri), numGPsQuad(
                _numGPsQuad), isMappingIGA2FEM(_isMappingIGA2FEM) {

    assert(_meshIGA != NULL);
    assert(_meshFE != NULL);
    assert(_meshIGA->type == EMPIRE_Mesh_IGAMesh);
    assert(_meshFE->type == EMPIRE_Mesh_FEMesh);

    DEBUG_OUT()<<"------------------------------"<<endl;
    DEBUG_OUT()<<"DEBUG for mapper "<<_name<<endl;
    DEBUG_OUT()<<"------------------------------"<<endl;

    if (_meshFE->triangulate() == NULL)
        meshFE = _meshFE;
    else
        meshFE = _meshFE->triangulate();

    projectedCoords = new vector<map<int, double*> >(meshFE->numNodes);

    if (isMappingIGA2FEM) {
        numNodesSlave = meshIGA->getNumNodes();
        numNodesMaster = meshFE->numNodes;
    } else {
        numNodesSlave = meshFE->numNodes;
        numNodesMaster = meshIGA->getNumNodes();
    }

    C_NR = new MathLibrary::SparseMatrix<double>(numNodesMaster, numNodesSlave);
    C_NN = new MathLibrary::SparseMatrix<double>(numNodesMaster, true);

    gaussTriangle = new IGAMortarMath::GaussQuadratureOnTriangle(numGPsTri);
    gaussQuad = new IGAMortarMath::GaussQuadratureOnQuad(numGPsQuad);

    initTables();

    projectPointsToSurface();

    // Write the projected points on to a file
    writeProjectedNodesOntoIGAMesh();

    computeCouplingMatrices();


    printCouplingMatricesToFile();

    C_NN->factorize();

    checkConsistency();
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
    /*  Projects all nodes of the FE side onto the IGA mesh.
     *
     *  0. Read input
     *
     *  1. Loop over all the patches in the IGA mesh
     *     1i. Get the IGA patch
     *    1ii. Initialize all projection flags to false
     *   1iii. Loop over all the elements on the fluid side
     *         1iii.1. Initialize the flag to false and the node id to zero
     *         1iii.2. Loop over all nodes of the current element to check if there exist one node has been projected already
     *                 1iii.2i. If the node has already been projected set projection flag to true
     *                1iii.2ii. Get the global ID of the projected node
     *               1iii.2iii. Break the loop
     *         1iii.3. Check if there exist one node in the current element has been successfully projected
     *                 1iii.3i. If so, use result of the projected node as the initial guess for the projection step
     *                1iii.3ii. Otherwise, find the nearest knot intersection as initial guess for the projection step
     *         1iii.4. Loop over each node at the current element in the FE side
     *                 1iii.4i. Get the node ID from the EFT
     *                1iii.4ii. Check if the node has been already projected. If not, compute the point projection on the IGA patch using the initial guess get from the last step the results are stored in the class member "projectedCoords"
     *                          1iii.4ii.1. get the Cartesian coordinates of the node in the FE side
     *                          1iii.4ii.2. Get an initial guess for the parametric location of the projected node of the FE side on the NURBS patch
     *                          1iii.4ii.3. Compute point projection on the NURBS patch using the Newton-Rapshon iteration method
     *                          1iii.4ii.4. Check if the Newton-Rapshon iterations have converged and if the points are coinciding
     *                                      1iii.4ii.4i. Set projection flag to true
     *                                     1iii.4ii.4ii. Insert the parametric locations of the projected FE nodes into the map projectedCoords
     *
     *   2. Loop over all the nodes in the FE side
     *      2i. Initialize projection flag to false
     *     2ii. Loop over all the patches in the IGA mesh and check if the node is has been projected to any patch
     *    2iii. If the node has not been projected to any patch of the IGA mesh, loop over all the patches in the mesh again and try to project it with better initial guess
     *          2iii.1. Get the NURBS patch
     *          2iii.2. Get the Cartesian coordinates of the node in the FE side
     *          2iii.3. Get a better initial guess for the Newton-Rapshon iterations using the variable REFINED_NUM_PARAMETRIC_LOCATIONS
     *          2iii.4. Compute point projection on the NURBS patch using the Newton-Rapshon iteration method
     *          2iii.5. Check if the Newton-Rapshon iterations have converged and if the points are coinciding
     *                  2iii.5i. Set projection flag to true
     *                 2iii.5ii. Insert the parametric locations of the projected FE nodes into the map projectedCoords
     *           2iv. If the node can still not be projected assert an error in the projection phase
     *
     *  3. Deallocate dynamically allocated memory
     */

    /// 0. Read input
    // Initialization of variables
    // Array of booleans containing flags on the projection of the FE nodes onto the NURBS patch
    bool *isProjected = new bool[meshFE->numNodes]; // deleted

    // Flag on whether a node is projected on the IGA mesh
    bool isNodeProjected = false;

    // Flag (??)
    bool isNodeInsideElementProjecteded;

    // id of the projected node
    int projectedNode;

    // Coordinates of the projected nodes on the NURBS patch
    double initialU, initialV;

    // Initialize the parametric coordinates of the FE node on the NURBS patch
    double projectedU;
    double projectedV;

    // Initialize flag on the convergence of the Newton-Rapshon iterations for the projection of a node on the NURBS patch
    bool isConverged, isConvergedInside;

    // Initialize the array of the Cartesian coordinates of a node in the FE side
    double cartesianCoords[3];

    // Initialize the node ID
    int nodeIndex;

    /// Initialize parametric coordinates of the projected node onto the NURBS patch
    double U, V;
    double P[3];

    // Get the number of patches in the IGA mesh
    int numPatches = meshIGA->getSurfacePatches().size();

    /// 1. Loop over all the patches in the IGA mesh
    for (int patchCount = 0; patchCount < numPatches; patchCount++) {
        /// 1i. Get the IGA patch
        IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchCount];

        /// 1ii. Initialize all projection flags to false
        for (int i = 0; i < meshFE->numNodes; i++)
            isProjected[i] = false;

        /// 1iii. Loop over all the elements on the fluid side
        for (int i = 0; i < meshFE->numElems; i++) {

            /// 1iii.1. Initialize the flag to false and the node id to zero
            isNodeInsideElementProjecteded = false;
            projectedNode = 0;

            /// 1iii.2. Loop over all nodes of the current element to check if there exist one node has been projected already
            for (int j = 0; j < meshFE->numNodesPerElem[i]; j++)
                if (isProjected[meshFEDirectElemTable[i][j]]) {
                    /// 1iii.2i. If the node has already been projected set projection flag to true
                    isNodeInsideElementProjecteded = true;

                    /// 1iii.2ii. Get the global ID of the projected node
                    projectedNode = meshFEDirectElemTable[i][j];

                    /// 1iii.2iii. Break the loop
                    break;
                }

            /// 1iii.3. Check if there exist one node in the current element has been successfully projected
            if (isNodeInsideElementProjecteded) {
                /// 1iii.3i. If so, use result of the projected node as the initial guess for the projection step
                initialU = (*projectedCoords)[projectedNode][patchCount][0];
                initialV = (*projectedCoords)[projectedNode][patchCount][1];

            } else {
                /// 1iii.3ii. Otherwise, find the nearest knot intersection as initial guess for the projection step

                // Get the node ID from the EFT
                nodeIndex = meshFEDirectElemTable[i][0];

                // Get the Cartesian coordinates of that node
                cartesianCoords[0] = meshFE->nodes[nodeIndex * 3];
                cartesianCoords[1] = meshFE->nodes[nodeIndex * 3 + 1];
                cartesianCoords[2] = meshFE->nodes[nodeIndex * 3 + 2];

                // Get accordingly an initial guess for the projection onto the NURBS patch
                thePatch->findInitialGuess4PointProjection(initialU, initialV, cartesianCoords);
            }

            /// 1iii.4. Loop over each node at the current element in the FE side
            for (int j = 0; j < meshFE->numNodesPerElem[i]; j++) {

                // 1iii.4i. Get the node ID from the EFT
                nodeIndex = meshFEDirectElemTable[i][j];

                /// 1iii.4ii. Check if the node has been already projected. If not, compute the point projection on the IGA patch using the initial guess get from the last step the results are stored in the class member "projectedCoords"
                if (!isProjected[nodeIndex]) {

                    /// 1iii.4ii.1. get the Cartesian coordinates of the node in the FE side
                    cartesianCoords[0] = meshFE->nodes[nodeIndex * 3];
                    cartesianCoords[1] = meshFE->nodes[nodeIndex * 3 + 1];
                    cartesianCoords[2] = meshFE->nodes[nodeIndex * 3 + 2];

                    /// 1iii.4ii.2. Get an initial guess for the parametric location of the projected node of the FE side on the NURBS patch
                    projectedU = initialU;
                    projectedV = initialV;

                    /// 1iii.4ii.3. Compute point projection on the NURBS patch using the Newton-Rapshon iteration method

                    // isConverged --> projected onto patch
                    isConvergedInside = thePatch->computePointProjectionOnPatch(projectedU,
                            projectedV, cartesianCoords, isConverged);

                    /// 1iii.4ii.4. Check if the Newton-Rapshon iterations have converged
                    if (isConvergedInside
                            && IGAMortarMath::computePointDistance(&meshFE->nodes[nodeIndex * 3],
                                    cartesianCoords) < disTol) {

                        /// 1iii.4ii.4i. Set projection flag to true
                        isProjected[nodeIndex] = true;

                        /// 1iii.4ii.4ii. Insert the parametric locations of the projected FE nodes into the map projectedCoords
                        double* coordTmp = new double(2);
                        coordTmp[0] = projectedU;
                        coordTmp[1] = projectedV;
                        (*projectedCoords)[nodeIndex].insert(
                                std::pair<int, double*>(patchCount, coordTmp));

                    }
                }
            }
        }
    }

    /// 2. Loop over all the nodes in the FE side
    for (int nodeIndex = 0; nodeIndex < meshFE->numNodes; nodeIndex++) {
        /// 2i.1 Initialize projection flag to false
        isNodeProjected = false;
        /// 2ii. Loop over all the patches in the IGA mesh and check if the node is has been projected to any patch
        for (int patchCount = 0; patchCount < numPatches; patchCount++) {
            if ((*projectedCoords)[nodeIndex].find(patchCount)
                    != (*projectedCoords)[nodeIndex].end()) {
                /// 2ii.1. If the node has been projected to a patch in the IGA mesh, set the flag to true
                isNodeProjected = true;

                /// 2ii.2. Break the loop
                break;
            }
        }

        /// 2iii. If the node has not been projected to any patch of the IGA mesh, loop over all the patches in the mesh again and try to project it with better initial guess
        if (!isNodeProjected) {
            /// 2iii.0 Get the Cartesian coordinates of the node in the FE side
            P[0] = meshFE->nodes[nodeIndex * 3];
            P[1] = meshFE->nodes[nodeIndex * 3 + 1];
            P[2] = meshFE->nodes[nodeIndex * 3 + 2];
            for (int patchCount = 0; patchCount < numPatches; patchCount++) {
                /// 2iii.1. Get the NURBS patch
                IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchCount];

                /// 2iii.2. Get the Cartesian coordinates of the node in the FE side
                double P_out[3];
                P_out[0] = P[0];
                P_out[1] = P[1];
				P_out[2] = P[2];

                /// 2iii.3. Get a better initial guess for the Newton-Rapshon iterations using the variable REFINED_NUM_PARAMETRIC_LOCATIONS
                thePatch->findInitialGuess4PointProjection(U, V, P,
                        REFINED_NUM_PARAMETRIC_LOCATIONS, REFINED_NUM_PARAMETRIC_LOCATIONS);

                /// 2iii.4. Compute point projection on the NURBS patch using the Newton-Rapshon iteration method
                isConvergedInside = thePatch->computePointProjectionOnPatch(U, V, P_out, isConverged);

                /// 2iii.5. Check if the Newton-Rapshon iterations have converged and if the points are coinciding
                bool hasConverged=isConvergedInside;
                if (hasConverged && IGAMortarMath::computePointDistance(P, P_out) < disTol) {
                    /// 2iii.5i. Set projection flag to true
                    isNodeProjected = true;

                    /// 2iii.5ii. Insert the parametric locations of the projected FE nodes into the map projectedCoords
                    double* coordTmp = new double(2);
                    coordTmp[0] = U;
                    coordTmp[1] = V;
                    (*projectedCoords)[nodeIndex].insert(make_pair(patchCount, coordTmp));
                }
            }

            /// 2iv. If the node can still not be projected assert an error in the projection phase
            if (!isNodeProjected) {
                ERROR_OUT() << " in IGAMortarMapper::projectPointsToSurface" << endl;
                ERROR_OUT() << "Cannot project node: " << nodeIndex << "  ("
                        << P[0] << ", " << P[1] << ", " << P[2] << ")" << endl;
                ERROR_OUT() << "Projection failed in IGA mapper " << name << endl;
                exit (EXIT_FAILURE);
            }
        }
    }

    /// 3. Deallocate dynamically allocated memory
    delete[] isProjected;

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

	//List of integrated element
	set<int> elementIntegrated;
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
        // Initialize the flag whether the projected FE element is not located at all on the patch
        bool isAllNodesOut= true;

        set<int> patchWithFullElt;
        set<int> patchWithSplitElt;
        // Loop over all the patches
        for (int patchCount = 0; patchCount < meshIGA->getSurfacePatches().size(); patchCount++) {
            isAllNodesOnPatch = true;
            isAllNodesOut= true;
			// Loop over all the nodes of the unclipped element
            for (int nodeCount = 0; nodeCount < numNodesElementFE; nodeCount++) {
                // Find the index of the node in the FE mesh
                int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];

                // Find whether this index is in the projected nodes array
                bool isNodeOnPatch = (*projectedCoords)[nodeIndex].find(patchCount)
                        != (*projectedCoords)[nodeIndex].end();

                if (!isNodeOnPatch) {
                    isAllNodesOnPatch = false;
                } else {
                	isAllNodesOut=false;
                }
            }
            // If all nodes are on the patch "patchCount", save this patch and go for next patch
            if (isAllNodesOnPatch) {
                patchWithFullElt.insert(patchCount);
                continue;
            }
            // If element is splitted for patch "patchCount", save this patch and go for next patch
            if(!isAllNodesOut) {
                patchWithSplitElt.insert(patchCount);
                continue;
            }
        }
        DEBUG_OUT()<<"Number of patch in set FULL "<<patchWithFullElt.size()<<endl;
        DEBUG_OUT()<<"Number of patch in set SPLIT "<<patchWithSplitElt.size()<<endl;

        /// 2. Compute the coupling matrices
        /// 2i. If the current element can be projected on one patch
		for (set<int>::iterator it = patchWithFullElt.begin();
				it != patchWithFullElt.end(); ++it) {
			int patchIndex=*it;
			IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchIndex];
			// Get the projected coordinates for the current element
			double clippedByPatchProjElementFEUV[numNodesElementFE * 2];

			for (int nodeCount = 0; nodeCount < numNodesElementFE;
					nodeCount++) {
				int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];
				clippedByPatchProjElementFEUV[nodeCount * 2] =
						(*projectedCoords)[nodeIndex][patchIndex][0];
				clippedByPatchProjElementFEUV[nodeCount * 2 + 1] =
						(*projectedCoords)[nodeIndex][patchIndex][1];
			}
			bool isIntegrated = computeCouplingMatrices4ClippedByPatchProjectedElement(thePatch,
					numNodesElementFE, clippedByPatchProjElementFEUV,
					projectedElementFEWZ, elemCount, numNodesElementFE);
			if(isIntegrated)
				elementIntegrated.insert(elemCount);
		}
        /// 2ii. If the current element is split in some patches
            // Loop over all the patches in the IGA Mesh having a piece of the FE element
		for (set<int>::iterator it = patchWithSplitElt.begin(); it != patchWithSplitElt.end(); it++) {
			int patchIndex=*it;
			IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchIndex];
			int numNodesClippedByPatchProjElementFE = 0;
			// the parameter coordinates in IGA of the sub-element divided by the patch
			double clippedByPatchProjElementFEUV[16];

			// the parameter coordinates in low order element of the sub-element divided by the patch
			double clippedByPatchProjElementFEWZ[16];

			// Cartesian coordinates of the low order element
			double elementFEXYZ[12];
			for (int i = 0; i < numNodesElementFE; i++) {
				int nodeIndex = meshFEDirectElemTable[elemCount][i];
				for (int j = 0; j < 3; j++)
					elementFEXYZ[i * 3 + j] = meshFE->nodes[nodeIndex * 3 + j];
			}
			// Stores points of the polygon clipped by the nurbs patch
			vector<pair<double,double> > polygonUV, polygonWZ;
			// Find nodes and edges of low order element inside IGA Patch
			for (int nodeCount = 0; nodeCount < numNodesElementFE; nodeCount++) {
				int nodeCountNext = (nodeCount + 1) % numNodesElementFE;
				int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];
				int nodeIndexNext = meshFEDirectElemTable[elemCount][nodeCountNext];

				// Flag on whether the node is inside the current NURBS patch
				bool isInsidePatch = (*projectedCoords)[nodeIndex].find(patchIndex)
							!= (*projectedCoords)[nodeIndex].end();
				bool isNextNodeInsidePatch = (*projectedCoords)[nodeIndexNext].find(patchIndex)
							!= (*projectedCoords)[nodeIndexNext].end();

/////////////////////////////////////////////////////////////////////////
////////////////////// MY METHOD TO CLIP BY PATCH ///////////////////////
/////////////////////////////////////////////////////////////////////////
			bool isProjectedOnPatchBoundary=1;
			double u, v;
			double w, z;
			double div, dis;
			double P0[3] = {0, 0, 0};
			int edge[2] = {-1, -1};
			double* P1 = &(meshFE->nodes[nodeIndex * 3]);
			double* P2 = &(meshFE->nodes[nodeIndexNext * 3]);
			DEBUG_OUT()<<"inside "<<isInsidePatch<<" / "<<isNextNodeInsidePatch<<endl;
			double u0=thePatch->getIGABasis()->getUBSplineBasis1D()->getFirstKnot();
			double uN=thePatch->getIGABasis()->getUBSplineBasis1D()->getLastKnot();
			double v0=thePatch->getIGABasis()->getVBSplineBasis1D()->getFirstKnot();
			double vN=thePatch->getIGABasis()->getVBSplineBasis1D()->getLastKnot();
			bool corner;
			if(!isInsidePatch && !isNextNodeInsidePatch) {
				u = 0;
				v = 0;
				thePatch->computePointMinimumDistanceToPatchBoundary(u, v, dis, P1, edge);
				corner=(u==u0 || u==uN) && (v==v0 || v==vN);
				if(corner) {
					polygonUV.push_back(make_pair(u,v));
					double nodeXYZ[3];
					double normalVec[3];
					thePatch->computeCartesianCoordinatesAndNormalVector(nodeXYZ, normalVec, u, v);
					double coordFE[2];
					if (numNodesElementFE == 3)
						IGAMortarMath::computeIntersectionBetweenLineAndTriangle(
								elementFEXYZ, nodeXYZ, normalVec, coordFE);
					else
						IGAMortarMath::computeIntersectionBetweenLineAndQuad(
								elementFEXYZ, nodeXYZ, normalVec, coordFE);
					w = coordFE[0];
					z = coordFE[1];
					polygonWZ.push_back(make_pair(w,z));
				}
				thePatch->computeLineMinimumDistanceToPatchBoundary(u, v, dis, P1, P2);
				corner=(u==u0 || u==uN) && (v==v0 || v==vN);
				if(corner) {
					polygonUV.push_back(make_pair(u,v));
					double nodeXYZ[3];
					double normalVec[3];
					thePatch->computeCartesianCoordinatesAndNormalVector(nodeXYZ, normalVec, u, v);
					double coordFE[2];
					if (numNodesElementFE == 3)
						IGAMortarMath::computeIntersectionBetweenLineAndTriangle(
								elementFEXYZ, nodeXYZ, normalVec, coordFE);
					else
						IGAMortarMath::computeIntersectionBetweenLineAndQuad(
								elementFEXYZ, nodeXYZ, normalVec, coordFE);
					w = coordFE[0];
					z = coordFE[1];
					polygonWZ.push_back(make_pair(w,z));
				}
				thePatch->computePointMinimumDistanceToPatchBoundary(u, v, dis, P2, edge);
				corner=(u==u0 || u==uN) && (v==v0 || v==vN);
				if(corner) {
					polygonUV.push_back(make_pair(u,v));
					double nodeXYZ[3];
					double normalVec[3];
					thePatch->computeCartesianCoordinatesAndNormalVector(nodeXYZ, normalVec, u, v);
					double coordFE[2];
					if (numNodesElementFE == 3)
						IGAMortarMath::computeIntersectionBetweenLineAndTriangle(
								elementFEXYZ, nodeXYZ, normalVec, coordFE);
					else
						IGAMortarMath::computeIntersectionBetweenLineAndQuad(
								elementFEXYZ, nodeXYZ, normalVec, coordFE);
					w = coordFE[0];
					z = coordFE[1];
					polygonWZ.push_back(make_pair(w,z));
				}

			} else if(isInsidePatch && isNextNodeInsidePatch) {
				// First point
				u = (*projectedCoords)[nodeIndex][patchIndex][0];
				v = (*projectedCoords)[nodeIndex][patchIndex][1];
				polygonUV.push_back(make_pair(u,v));
				w=projectedElementFEWZ[nodeCount * 2];
				z=projectedElementFEWZ[nodeCount * 2 + 1];
				polygonWZ.push_back(make_pair(w,z));
				// Second point
				u = (*projectedCoords)[nodeIndexNext][patchIndex][0];
				v = (*projectedCoords)[nodeIndexNext][patchIndex][1];
				polygonUV.push_back(make_pair(u,v));
				w=projectedElementFEWZ[nodeCountNext * 2];
				z=projectedElementFEWZ[nodeCountNext * 2 + 1];
				polygonWZ.push_back(make_pair(w,z));

			} else if(isInsidePatch && !isNextNodeInsidePatch) {
				// First point
				u = (*projectedCoords)[nodeIndex][patchIndex][0];
				v = (*projectedCoords)[nodeIndex][patchIndex][1];
				polygonUV.push_back(make_pair(u,v));
				w=projectedElementFEWZ[nodeCount * 2];
				z=projectedElementFEWZ[nodeCount * 2 + 1];
				polygonWZ.push_back(make_pair(w,z));
				// Second point
				isProjectedOnPatchBoundary=thePatch->computePointProjectionOnPatchBoundary_Brute(u, v, div, dis, P1, P2);
				if (isProjectedOnPatchBoundary && dis <= disTol) {
					polygonUV.push_back(make_pair(u,v));
					double P1x = projectedElementFEWZ[nodeCount * 2];
					double P1y = projectedElementFEWZ[nodeCount * 2 + 1];
					double P2x = projectedElementFEWZ[nodeCountNext * 2];
					double P2y = projectedElementFEWZ[nodeCountNext * 2 + 1];
					w = P1x * (1 - div) + P2x * div;
					z = P1y * (1 - div) + P2y * div;
					polygonWZ.push_back(make_pair(w,z));
				}
				thePatch->computePointMinimumDistanceToPatchBoundary(u, v, dis, P2, edge);
				corner=(u==u0 || u==uN) && (v==v0 || v==vN);
				if(corner) {
					polygonUV.push_back(make_pair(u,v));
					double nodeXYZ[3];
					double normalVec[3];
					thePatch->computeCartesianCoordinatesAndNormalVector(nodeXYZ, normalVec, u, v);
					double coordFE[2];
					if (numNodesElementFE == 3)
						IGAMortarMath::computeIntersectionBetweenLineAndTriangle(
								elementFEXYZ, nodeXYZ, normalVec, coordFE);
					else
						IGAMortarMath::computeIntersectionBetweenLineAndQuad(
								elementFEXYZ, nodeXYZ, normalVec, coordFE);
					w = coordFE[0];
					z = coordFE[1];
					polygonWZ.push_back(make_pair(w,z));
				}


			} else if(!isInsidePatch && isNextNodeInsidePatch) {
				thePatch->computePointMinimumDistanceToPatchBoundary(u, v, dis, P2, edge);
				corner=(u==u0 || u==uN) && (v==v0 || v==vN);
				if(corner) {
					polygonUV.push_back(make_pair(u,v));
					double nodeXYZ[3];
					double normalVec[3];
					thePatch->computeCartesianCoordinatesAndNormalVector(nodeXYZ, normalVec, u, v);
					double coordFE[2];
					if (numNodesElementFE == 3)
						IGAMortarMath::computeIntersectionBetweenLineAndTriangle(
								elementFEXYZ, nodeXYZ, normalVec, coordFE);
					else
						IGAMortarMath::computeIntersectionBetweenLineAndQuad(
								elementFEXYZ, nodeXYZ, normalVec, coordFE);
					w = coordFE[0];
					z = coordFE[1];
					polygonWZ.push_back(make_pair(w,z));
				}
				// First point
				u = (*projectedCoords)[nodeIndexNext][patchIndex][0];
				v = (*projectedCoords)[nodeIndexNext][patchIndex][1];
				isProjectedOnPatchBoundary=thePatch->computePointProjectionOnPatchBoundary_Brute(u, v, div, dis, P2, P1);
				if (isProjectedOnPatchBoundary && dis <= disTol) {
					polygonUV.push_back(make_pair(u,v));
					double P1x = projectedElementFEWZ[nodeCountNext * 2];
					double P1y = projectedElementFEWZ[nodeCountNext * 2 + 1];
					double P2x = projectedElementFEWZ[nodeCount * 2];
					double P2y = projectedElementFEWZ[nodeCount * 2 + 1];
					w = P1x * (1 - div) + P2x * div;
					z = P1y * (1 - div) + P2y * div;
					polygonWZ.push_back(make_pair(w,z));
				}
				// Second point
				u = (*projectedCoords)[nodeIndexNext][patchIndex][0];
				v = (*projectedCoords)[nodeIndexNext][patchIndex][1];
				polygonUV.push_back(make_pair(u,v));
				w=projectedElementFEWZ[nodeCountNext * 2];
				z=projectedElementFEWZ[nodeCountNext * 2 + 1];
				polygonWZ.push_back(make_pair(w,z));
			}
			if(!isProjectedOnPatchBoundary) {
				if(thePatch->isTrimmed()) {
					WARNING_OUT() << "Warning in IGAMortarMapper::computeCouplingMatrices"
							<< endl;
					WARNING_OUT() << "Cannot find point projection on patch boundary" << endl;
					WARNING_OUT() << "But proceed anyway !" << endl;
					break;//break loop over node
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
							<< patchIndex << "] boundary" << endl;
					ERROR_OUT() << "Projection failed in IGA mapper " << name << endl;
					exit (EXIT_FAILURE);
				}
			}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
			}
			DEBUG_OUT()<<"Display coordinates of FE element in parametric domain"<<endl;
			for(int nodeCount=0;nodeCount<polygonUV.size();nodeCount++) {
				DEBUG_OUT()<<polygonUV[nodeCount].first<<" / "<<polygonUV[nodeCount].second<<endl;
			}
			DEBUG_OUT()<<"Display coordinates of normalized FE element"<<endl;
			for(int nodeCount=0;nodeCount<polygonWZ.size();nodeCount++) {
				DEBUG_OUT()<<polygonWZ[nodeCount].first<<" / "<<polygonWZ[nodeCount].second<<endl;
			}
			IGAMortarMath::cleanPolygon(polygonUV,polygonWZ);
			//IGAMortarMath::cleanPolygon(polygonWZ);

			DEBUG_OUT()<<"["<<elemCount<<"] Element of size "<<polygonUV.size()<<endl;
			numNodesClippedByPatchProjElementFE=polygonUV.size();

			//Proceed toward integration if the polygon is valid, i.e. at least a triangle
			if (numNodesClippedByPatchProjElementFE >= 3) {
				// Check orientation of the polygon and reverse it if clockwise
				double v1[2],v2[2];
				v1[0]=polygonUV[1].first-polygonUV[0].first;
				v1[1]=polygonUV[1].second-polygonUV[0].second;
				v2[0]=polygonUV[2].first-polygonUV[1].first;
				v2[1]=polygonUV[2].second-polygonUV[1].second;
				if(IGAMortarMath::computeCrossProduct2D(v1[0],v1[1],v2[0],v2[1])<0){
					reverse(polygonUV.begin(),polygonUV.end());
					reverse(polygonWZ.begin(),polygonWZ.end());
				}
				///Copy the output polygon into output single pointer array
				DEBUG_OUT()<<"Display coordinates of FE element in parametric domain"<<endl;
				for(int nodeCount=0;nodeCount<polygonUV.size();nodeCount++) {
					DEBUG_OUT()<<polygonUV[nodeCount].first<<" / "<<polygonUV[nodeCount].second<<endl;
					clippedByPatchProjElementFEUV[2*nodeCount]=polygonUV[nodeCount].first;
					clippedByPatchProjElementFEUV[2*nodeCount+1]=polygonUV[nodeCount].second;

				}
				DEBUG_OUT()<<"Display coordinates of normalized FE element"<<endl;
				for(int nodeCount=0;nodeCount<polygonWZ.size();nodeCount++) {
					DEBUG_OUT()<<polygonWZ[nodeCount].first<<" / "<<polygonWZ[nodeCount].second<<endl;
					clippedByPatchProjElementFEWZ[2*nodeCount]  =polygonWZ[nodeCount].first;
					clippedByPatchProjElementFEWZ[2*nodeCount+1]=polygonWZ[nodeCount].second;

				}
				// compute the coupling matrix for the sub-element
				bool isIntegrated=computeCouplingMatrices4ClippedByPatchProjectedElement(thePatch,
						numNodesClippedByPatchProjElementFE, clippedByPatchProjElementFEUV,
						clippedByPatchProjElementFEWZ, elemCount, numNodesElementFE);
				if(isIntegrated)
					elementIntegrated.insert(elemCount);

			} //end of if polygon has more than 3 edges

		} // end of loop over set of split patch

    } // end of loop over all the element
    if(elementIntegrated.size()!=meshFE->numElems) {
    	ERROR_OUT()<<"Number of FE mesh integrated is "<<elementIntegrated.size()<<" over "<<meshFE->numElems<<endl;
    	for(int i=0;i<meshFE->numElems;i++) {
    		if(!elementIntegrated.count(i))
    			ERROR_OUT()<<"Missing element number "<<i<<endl;
    	}
    	ERROR_BLOCK_OUT("IGAMortarMapper","ComputeCouplingMatrices","Not all element in FE mesh integrated ! Coupling matrices invalid");
    	exit(-1);
    }
}

bool IGAMortarMapper::computeCouplingMatrices4ClippedByPatchProjectedElement(
        IGAPatchSurface* _thePatch, int _numNodesClippedByPatchProjectedElement,
        double* _clippedByPatchProjElementFEUV, double* _clippedByPatchProjElementFEWZ,
        int _elemCount, int _numNodesElementFE) {
    /*
     * Compute the coupling matrices for the input polygon which has been clipped by patch boundary previously
     * 2. Find the knot span which the current element located in.
     *   2.1 If the whole element is located in a single knot span
     *   	2.1.1  If patch is trimmed
     *   		2.1.1.1	Clip the input polygon by the trimming polygon
     *   		2.1.1.2 Integrate if the polygon is valid (at least a triangle)
     *		2.1.2 Else integrate the input polygon
     *   2.2 If the polygon is living inside several knot spans
     *       2.2.1 Loop over the knot span
     *       	2.2.1.1 Clip by Knot span
     *       	2.2.1.2	Try to integrate following same procedure as 2.1.1 and 2.1.2
     */

	// Flag indicating if current polygon has been integrated
	bool isIntegrated=false;

    double *knotVectorU = _thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
    double *knotVectorV = _thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();

/// 1.find the knot span which the current element located in.
//      from minSpanu to maxSpanu in U-direction, and from minSpanV to max SpanV in V-direction
    int span[4];
    computeKnotSpanOfProjElement(_thePatch, _numNodesClippedByPatchProjectedElement, _clippedByPatchProjElementFEUV,span);

    int minSpanU = span[0];
    int maxSpanU = span[1];
    int minSpanV = span[2];
    int maxSpanV = span[3];

	/* 2.1 if the whole element is located in a single knot span, set the coordinates in the linear element as
     the full triangle or full quadrilateral defined at the beginning of this function. */
    if (minSpanU == maxSpanU && minSpanV == maxSpanV) {
    	if(_thePatch->isTrimmed()) {
			int numNodesClippedByTrim;
			double* ClippedByTrim;
        	///Clip to trimming window
			computeTrimmedPolygon2(_thePatch,_numNodesClippedByPatchProjectedElement,_clippedByPatchProjElementFEUV,numNodesClippedByTrim,ClippedByTrim);
			// If polygon is not degenerated (line or point)
			if(numNodesClippedByTrim >2) {
				double* ClippedByPatchProjElementFEWZ;
				/// Update local coord and integrate
				computeLocalElementCoord(numNodesClippedByTrim,ClippedByTrim,_numNodesClippedByPatchProjectedElement,_clippedByPatchProjElementFEUV,_clippedByPatchProjElementFEWZ,ClippedByPatchProjElementFEWZ);
				integrate(_thePatch, numNodesClippedByTrim,
						ClippedByTrim, minSpanU, minSpanV, ClippedByPatchProjElementFEWZ,
						_elemCount, _numNodesElementFE);
				delete ClippedByPatchProjElementFEWZ;
				isIntegrated=true;
			} else {
//				DEBUG_OUT()<<"Elem not integrated is "<<_elemCount<<endl;
//				for(int i=0;i<_numNodesClippedByPatchProjectedElement;i++)
//				{
//					DEBUG_OUT()<<"U: "<<_clippedByPatchProjElementFEUV[i*2]<<"/ V: "<<_clippedByPatchProjElementFEUV[i*2+1]<<endl;
//				}
			}
			delete ClippedByTrim;
    	} else {
			integrate(_thePatch, _numNodesClippedByPatchProjectedElement,
					_clippedByPatchProjElementFEUV, minSpanU, minSpanV, _clippedByPatchProjElementFEWZ,
					_elemCount, _numNodesElementFE);
			isIntegrated=true;
    	}
    } else {
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
                        for (int nodeCount = 0; nodeCount < polygonResult->size(); nodeCount++)
                            delete[] (*polygonResult)[nodeCount];
                        delete polygonResult;

                        /// Proceed to integration
                        if(_thePatch->isTrimmed()) {
                        	///Clip to trimming window
							int numNodesClippedByTrim;
							double* ClippedByTrim;
							computeTrimmedPolygon2(_thePatch,numNodesClippedByKnotSpanProjElementFE,
									ClippedByKnotSpanProjElementFEUV,numNodesClippedByTrim,ClippedByTrim);
							// If polygon is not degenerated (line or point)
							if(numNodesClippedByTrim>2) {
								double* ClippedByTrimProjElementFEWZ;
								/// Update local coord and integrate
								computeLocalElementCoord(numNodesClippedByTrim,ClippedByTrim,
										_numNodesClippedByPatchProjectedElement,_clippedByPatchProjElementFEUV,
										_clippedByPatchProjElementFEWZ,ClippedByTrimProjElementFEWZ);
								integrate(_thePatch, numNodesClippedByTrim,
										ClippedByTrim, spanU, spanV, ClippedByTrimProjElementFEWZ,
										_elemCount, _numNodesElementFE);
								delete ClippedByTrimProjElementFEWZ;
								isIntegrated=true;
							} else {
//								DEBUG_OUT()<<"Elem not integrated is "<<_elemCount<<endl;
//								for(int i=0;i<_numNodesClippedByPatchProjectedElement;i++)
//								{
//									INFO_OUT()	<<"U: "	<<_clippedByPatchProjElementFEUV[i*2]
//												<<"/ V: "<<_clippedByPatchProjElementFEUV[i*2+1]<<endl;
//								}
							}
							delete ClippedByTrim;
                        } else {
                            // the array to store the coordinates(Linear Element) of the clipped sub-element
                            double* ClippedByKnotSpanProjElementFEWZ;
                            /// find the local coordinates in FE(triangle or quadrilateral) from parameter of IGAPatch
                            computeLocalElementCoord(numNodesClippedByKnotSpanProjElementFE,ClippedByKnotSpanProjElementFEUV,
                            		_numNodesClippedByPatchProjectedElement,_clippedByPatchProjElementFEUV,
                            		_clippedByPatchProjElementFEWZ,ClippedByKnotSpanProjElementFEWZ);
                            /// integrate
							integrate(_thePatch, numNodesClippedByKnotSpanProjElementFE,
									ClippedByKnotSpanProjElementFEUV, spanU, spanV,
									ClippedByKnotSpanProjElementFEWZ, _elemCount, _numNodesElementFE);
							delete ClippedByKnotSpanProjElementFEWZ;
							isIntegrated=true;
                        } // end if/else trimmed

                    } // end if clip by knot span successful

                } // end if different knot span

            } // end for knot span of element

    } // end if/else same span
    return isIntegrated;
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

    // Definitions
    double IGABasisFctsI = 0;
    double IGABasisFctsJ = 0;
    double basisFctsMaster = 0;
    double basisFctsSlave = 0;
    int numNodesElMaster = 0;
    int numNodesElSlave = 0;

    int pDegree = _thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = _thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
    int nShapeFuncsIGA = (pDegree + 1) * (qDegree + 1);

    if (isMappingIGA2FEM) {
        numNodesElMaster = _numNodesElementFE;
        numNodesElSlave = nShapeFuncsIGA;
    } else {
        numNodesElMaster = nShapeFuncsIGA;
        numNodesElSlave = _numNodesElementFE;
    }

    double elementCouplingMatrixNN[numNodesElMaster * (numNodesElMaster + 1) / 2];
    double elementCouplingMatrixNR[numNodesElSlave * numNodesElMaster];

    for (int arrayIndex = 0; arrayIndex < numNodesElMaster * (numNodesElMaster + 1) / 2;
            arrayIndex++)
        elementCouplingMatrixNN[arrayIndex] = 0;

    for (int arrayIndex = 0; arrayIndex < numNodesElSlave * numNodesElMaster; arrayIndex++)
        elementCouplingMatrixNR[arrayIndex] = 0;

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

            for (int i = 0; i < numNodesElMaster; i++)
                for (int j = i; j < numNodesElMaster; j++) {
                    if (isMappingIGA2FEM)
                        elementCouplingMatrixNN[count++] += shapeFuncsFE[i] * shapeFuncsFE[j]
                                * Jacobian * theGaussQuadrature->weights[GPCount];
                    else {
                        IGABasisFctsI =
                                localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                                        1, 0, 0, i)];
                        IGABasisFctsJ =
                                localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                                        1, 0, 0, j)];
                        elementCouplingMatrixNN[count++] += IGABasisFctsI * IGABasisFctsJ * Jacobian
                                * theGaussQuadrature->weights[GPCount];
                    }

                }

            /// 2.2.9 integrate the shape function product for C_NR(Linear shape function multiply IGA shape function)
            count = 0;

            for (int i = 0; i < numNodesElMaster; i++)
                for (int j = 0; j < numNodesElSlave; j++) {
                    if (isMappingIGA2FEM) {
                        basisFctsMaster = shapeFuncsFE[i];
                        basisFctsSlave =
                                localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                                        1, 0, 0, j)];
                    } else {
                        basisFctsMaster =
                                localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                                        1, 0, 0, i)];
                        basisFctsSlave = shapeFuncsFE[j];
                    }
                    elementCouplingMatrixNR[count++] += basisFctsMaster * basisFctsSlave * Jacobian
                            * theGaussQuadrature->weights[GPCount];
                }

        }
    }

/// 3.Assemble the element coupling matrix to the global coupling matrix.
    int dofIGA[nShapeFuncsIGA];
    _thePatch->getIGABasis()->getBasisFunctionsIndex(_spanU, _spanV, dofIGA);

    for (int i = 0; i < nShapeFuncsIGA; i++)
        dofIGA[i] = _thePatch->getControlPointNet()[dofIGA[i]]->getDofIndex();

    int count = 0;
    int dof1, dof2;
    for (int i = 0; i < numNodesElMaster; i++)
        for (int j = i; j < numNodesElMaster; j++) {
            if (isMappingIGA2FEM) {
                dof1 = meshFEDirectElemTable[_elementIndex][i];
                dof2 = meshFEDirectElemTable[_elementIndex][j];
            } else {
                dof1 = dofIGA[i];
                dof2 = dofIGA[j];
            }
            if (dof1 < dof2)
                (*C_NN)(dof1, dof2) += elementCouplingMatrixNN[count++];
            else
                (*C_NN)(dof2, dof1) += elementCouplingMatrixNN[count++];
        }

    count = 0;
    for (int i = 0; i < numNodesElMaster; i++)
        for (int j = 0; j < numNodesElSlave; j++) {
            if (isMappingIGA2FEM) {
                dof1 = meshFEDirectElemTable[_elementIndex][i];
                dof2 = dofIGA[j];
            } else {
                dof1 = dofIGA[i];
                dof2 = meshFEDirectElemTable[_elementIndex][j];
            }
            (*C_NR)(dof1, dof2) += elementCouplingMatrixNR[count++];
        }

    if (_numNodes > 4)
        for (int i = 0; i < quadratureVecUV.size(); i++) {
            delete quadratureVecUV[i];
            delete quadratureVecWZ[i];
        }
}

void IGAMortarMapper::computeKnotSpanOfProjElement(IGAPatchSurface* _thePatch, int _numNodesClippedByPatchProjectedElement,
        double* _clippedByPatchProjElementFEUV, int* _span) {
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
    _span[0]=minSpanU;
    _span[1]=maxSpanU;
    _span[2]=minSpanV;
    _span[3]=maxSpanV;
}

void IGAMortarMapper::computeTrimmedPolygon(IGAPatchSurface* _thePatch,int _numNodesPolygonToClip,
        double* _nodesPolygonToClip, int& _numNodesClippedByTrimming, double*& _clippedByTrimming) {

    vector<double*> *polygonResult = new vector<double*>; // deleted
	/// Copy input in working structure
    int numNodesClippedByKnotSpanProjElementFE=_numNodesPolygonToClip;
    double* ClippedByKnotSpanProjElementFEUV = new double[_numNodesPolygonToClip*3];
    for (int nodeCount = 0; nodeCount < numNodesClippedByKnotSpanProjElementFE;
            nodeCount++) {
        ClippedByKnotSpanProjElementFEUV[nodeCount * 3] =
                _nodesPolygonToClip[nodeCount * 2];
        ClippedByKnotSpanProjElementFEUV[nodeCount * 3 + 1] =
                _nodesPolygonToClip[nodeCount*2 +1 ];
        ClippedByKnotSpanProjElementFEUV[nodeCount * 3 + 2] = 0;
    }
    /// Go over all trim loops of the patch to clip
	for(int loop=0;loop<_thePatch->getTrimming().getNumOfLoops();loop++) {
		const std::vector<double> polygonWindow=_thePatch->getTrimming().getLoop(loop).getPolylines();
		int nPointsClipWindow=polygonWindow.size()/2;
		/// Restructure the data from the loop such that it is clipper compliant
		double fullPolygon[nPointsClipWindow*3];
		for(int p=0;p<nPointsClipWindow;p++) {
			fullPolygon[3*p]=polygonWindow[2*p];
			fullPolygon[3*p+1]=polygonWindow[2*p+1];
			fullPolygon[3*p+2]=0;
		}
		///Create the clipper based on trimming polygon
		MortarMath::PolygonClipper clipper=MortarMath::PolygonClipper(fullPolygon,nPointsClipWindow,2);
		/// Apply clipping on the input
		clipper.clip(ClippedByKnotSpanProjElementFEUV,numNodesClippedByKnotSpanProjElementFE,polygonResult);
		delete[] ClippedByKnotSpanProjElementFEUV;
        // number of nodes of the clipped polygon
        numNodesClippedByKnotSpanProjElementFE = polygonResult->size();
        // the array to store the coordinates(on IGA patch) of the clipped sub-element
        ClippedByKnotSpanProjElementFEUV=new double[numNodesClippedByKnotSpanProjElementFE * 3];
        for (int nodeCount = 0; nodeCount < numNodesClippedByKnotSpanProjElementFE;
                nodeCount++) {
            ClippedByKnotSpanProjElementFEUV[nodeCount * 3] =
                    (*polygonResult)[nodeCount][0];
            ClippedByKnotSpanProjElementFEUV[nodeCount * 3 + 1] =
                    (*polygonResult)[nodeCount][1];
            ClippedByKnotSpanProjElementFEUV[nodeCount * 3 + 2] =
                    (*polygonResult)[nodeCount][2];
        }
        /// Delete Polygon nodes
        for (int nodeCount = 0; nodeCount < polygonResult->size(); nodeCount++)
            delete[] (*polygonResult)[nodeCount];
        polygonResult->resize(0);
	}
	/// Copy result to output
	_numNodesClippedByTrimming=numNodesClippedByKnotSpanProjElementFE;
	_clippedByTrimming=new double[numNodesClippedByKnotSpanProjElementFE*2];
    for (int nodeCount = 0; nodeCount < numNodesClippedByKnotSpanProjElementFE;
            nodeCount++) {
    	_clippedByTrimming[nodeCount * 2] =
    			ClippedByKnotSpanProjElementFEUV[nodeCount * 3];
    	_clippedByTrimming[nodeCount * 2 + 1] =
    			ClippedByKnotSpanProjElementFEUV[nodeCount * 3 +1];
    }
	/// Delete result
	delete[] ClippedByKnotSpanProjElementFEUV;
    /// Delete Polygon
    delete polygonResult;
}

void IGAMortarMapper::computeTrimmedPolygon2(IGAPatchSurface* _thePatch,int _numNodesPolygonToClip,
        double* _nodesPolygonToClip, int& _numNodesClippedByTrimming, double*& _clippedByTrimming) {

	using namespace ClipperLib;

	Path subj;
	Paths clip(_thePatch->getTrimming().getNumOfLoops()), solution;

    /// Go over all trim loops of the patch to clip
	for(int loop=0;loop<_thePatch->getTrimming().getNumOfLoops();loop++) {
		const std::vector<double> clippingWindow=_thePatch->getTrimming().getLoop(loop).getPolylines();
		int nPointsClipWindow=clippingWindow.size()/2;
		for(int p=0;p<nPointsClipWindow;p++) {
			clip[loop]<<IntPoint((cInt)(clippingWindow[2*p]*1e16),(cInt)(clippingWindow[2*p+1]*1e16));
		}
	}
	for(int p=0;p<_numNodesPolygonToClip;p++) {
		subj<<IntPoint((cInt)(_nodesPolygonToClip[2*p]*1e16),(cInt)(_nodesPolygonToClip[2*p+1]*1e16));
	}
	bool isCounterClockwise=1;
	if(Orientation(subj)!=isCounterClockwise) ReversePath(subj);
	Clipper c;
	c.AddPath(subj, ptSubject, true);
	c.AddPaths(clip, ptClip, true);
	c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);
	if(solution.empty()) {
		_numNodesClippedByTrimming=0;
		_clippedByTrimming=new double[2*_numNodesClippedByTrimming];
		return;
	}
	_numNodesClippedByTrimming=solution[0].size();
	_clippedByTrimming=new double[2*_numNodesClippedByTrimming];
	for(int p=0;p<_numNodesClippedByTrimming;p++) {
		_clippedByTrimming[2*p]=1e-16*solution[0][p].X;
		_clippedByTrimming[2*p+1]=1e-16*solution[0][p].Y;
	}
}


void IGAMortarMapper::computeLocalElementCoord(int _numNodesClippedByKnotSpanProjElementFE, double* _clippedByKnotSpanProjElementFEUV,
		int _numNodesClippedByPatchProjectedElement, double* _clippedByPatchProjElementFEUV,
		double* _clippedByPatchProjElementFEWZ, double*&  ClippedByKnotSpanProjElementFEWZ) {
    ClippedByKnotSpanProjElementFEWZ=new double[_numNodesClippedByKnotSpanProjElementFE * 2];
    /// find the local coordinates in FE(triangle or quadrilateral) from parameter of IGAPatch
    double nodeWZ[2];
    for (int nodeCount = 0; nodeCount < _numNodesClippedByKnotSpanProjElementFE;
            nodeCount++) {
        if (_numNodesClippedByPatchProjectedElement == 3)
            IGAMortarMath::computeLocalCoordsInTriangle(
                    _clippedByPatchProjElementFEUV, &_clippedByKnotSpanProjElementFEUV[nodeCount * 2],
                    nodeWZ);
        else
            IGAMortarMath::computeLocalCoordsInQuad(
                    _clippedByPatchProjElementFEUV, &_clippedByKnotSpanProjElementFEUV[nodeCount * 2],
                    nodeWZ);
        IGAMortarMath::computeLinearCombinationValueFromVerticesValues(
                _numNodesClippedByPatchProjectedElement, 2,
                _clippedByPatchProjElementFEWZ, nodeWZ,
                &ClippedByKnotSpanProjElementFEWZ[nodeCount * 2]);
    }
}

void IGAMortarMapper::consistentMapping(const double* _slaveField, double *_masterField) {
    /*
     * Mapping of the
     * C_NN * x_master = C_NR * x_slave
     */
    double* tmpVec = new double[numNodesMaster];

// 1. matrix vector product (x_tmp = C_NR * x_slave)
    C_NR->mulitplyVec(_slaveField, tmpVec, numNodesMaster);

// 2. solve C_NN * x_master = x_tmp
    C_NN->solve(_masterField, tmpVec);

    delete[] tmpVec;

}

void IGAMortarMapper::conservativeMapping(const double* _masterField, double *_slaveField) {
    /*
     * Mapping of the
     * f_slave = (C_NN^(-1) * C_NR)^T * f_master
     */

    double* tmpVec = new double[numNodesMaster];

// 1. solve C_NN * f_tmp = f_master;
    C_NN->solve(tmpVec, const_cast<double *>(_masterField));

// 2. matrix vector product (f_slave = C_NR^T * f_tmp)
    C_NR->transposeMulitplyVec(tmpVec, _slaveField, numNodesMaster);

    delete[] tmpVec;

}

void IGAMortarMapper::writeProjectedNodesOntoIGAMesh() {
    // Initializations
    IGAPatchSurface* IGAPatch;
    int numXiKnots, numEtaKnots, numXiControlPoints, numEtaControlPoints;
    double xiKnot, etaKnot;

    // Open file for writing the projected nodes
    ofstream projectedNodesFile;
    const string UNDERSCORE = "_";
    string projectedNodesFileName = "projectedNodesOntoNURBSSurface" + UNDERSCORE + name + ".m";
    projectedNodesFile.open(projectedNodesFileName.c_str());
    projectedNodesFile.precision(14);
    projectedNodesFile << std::dec;

    projectedNodesFile << HEADER_DECLARATION << endl << endl;

    // Loop over all the patches
    for (int patchCounter = 0; patchCounter < meshIGA->getSurfacePatches().size(); patchCounter++) {
        // Get the IGA patch
        IGAPatch = meshIGA->getSurfacePatches()[patchCounter];

        // Get number of knots in each parametric direction
        numXiKnots = IGAPatch->getIGABasis()->getUBSplineBasis1D()->getNoKnots();
        numEtaKnots = IGAPatch->getIGABasis()->getVBSplineBasis1D()->getNoKnots();

        // Write out the knot vector in u-direction
        projectedNodesFile << "Patch" << patchCounter << endl << endl;
        projectedNodesFile << "xiKnotVector" << endl;

        for (int xiCounter = 0; xiCounter < numXiKnots; xiCounter++) {
            xiKnot = IGAPatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[xiCounter];
            projectedNodesFile << xiKnot << " ";
        }
        projectedNodesFile << endl << endl;

        // Write out the knot vector in v-direction
        projectedNodesFile << "etaKnotVector" << endl;
        for (int etaCounter = 0; etaCounter < numEtaKnots; etaCounter++) {
            etaKnot = IGAPatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[etaCounter];
            projectedNodesFile << etaKnot << " ";
        }
        projectedNodesFile << endl << endl;

        // Loop over all the nodes
        for (int nodeIndex = 0; nodeIndex < meshFE->numNodes; nodeIndex++) {
            // Loop over all the patches on which this node has been successfully projected
            for (map<int, double*>::iterator it = (*projectedCoords)[nodeIndex].begin();
                    it != (*projectedCoords)[nodeIndex].end(); it++)
                if (it->first == patchCounter) {
                    projectedNodesFile << nodeIndex << "\t" << it->first << "\t" << it->second[0]
                            << "\t" << it->second[1] << endl;
                }
        }
        projectedNodesFile << endl;
    }

    // Close file
    projectedNodesFile.close();
}

void IGAMortarMapper::printCouplingMatrices() {

    ERROR_OUT() << "C_NN" << endl;
    C_NN->printCSR();
    ERROR_OUT() << "C_NR" << endl;
    C_NR->printCSR();
}

void IGAMortarMapper::printCouplingMatricesToFile() {
	DEBUG_OUT()<<"Size of C_NR is "<<numNodesMaster<<" by "<<numNodesSlave<<endl;
    if(Message::userSetOutputLevel==Message::DEBUG) {
		C_NR->printToFile("C_NR.dat");
		C_NN->printToFile("C_NN.dat");
	}
}

void IGAMortarMapper::checkConsistency() {
    double ones[numNodesSlave];
    for(int i=0;i<numNodesSlave;i++) {
    	ones[i]=1.0;
    }
    double output[numNodesMaster];
    this->consistentMapping(ones,output);
    double norm=0;
    for(int i=0;i<numNodesMaster;i++) {
    	norm+=output[i]*output[i];
    }
    norm=sqrt(norm/numNodesMaster);
    DEBUG_OUT()<<"Norm of output field = "<<norm<<endl;
    if(fabs(norm-1.0)>1e-6) {
    	ERROR_OUT()<<"Coupling not consistent !"<<endl;
    	ERROR_OUT()<<"Coupling of unit field deviating from 1 of "<<fabs(norm-1.0)<<endl;
    	exit(-1);
    }
}
}

