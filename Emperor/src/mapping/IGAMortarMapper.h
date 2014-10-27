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
/******************************************************************************//**
 * \file IGAMortarMapper.h
 * The header file of class IGAMortarMapper.
 * \date 3/6/2013
 *********************************************************************************/

#ifndef IGAMORTARMAPPER_H_
#define IGAMORTARMAPPER_H_

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <assert.h>
#include "AbstractMapper.h"

namespace EMPIRE {

namespace MathLibrary {
template<class T> class SparseMatrix;
}

class IGAPatchSurface;
class IGAMesh;
class FEMesh;
class DataField;

namespace IGAMortarMath {

class GaussQuadratureOnTriangle;
class GaussQuadratureOnQuad;

}

/***********************************************************************************************
 * \brief This class is related to IGA mortar mapping
 * ***********/
class IGAMortarMapper: public AbstractMapper {

private:
    /// Name of the mapper
    std::string name;

    /// IGA Mesh
    IGAMesh *meshIGA;

    /// Fluid Mesh
    FEMesh *meshFE;

    /// The element freedom table for the fluid mesh
    int **meshFEDirectElemTable;

    /// The mass-like matrix
    MathLibrary::SparseMatrix<double> *C_NN;

    /// The right-hand side matrix
    MathLibrary::SparseMatrix<double> *C_NR;

    /// Quadrature rule over the triangulated subdomains
    IGAMortarMath::GaussQuadratureOnTriangle *gaussTriangle;

    /// Quadrature rule over the non-triangulated subdomains
    IGAMortarMath::GaussQuadratureOnQuad *gaussQuad;

    /// The parametric coordinates of the projected nodes on the surface
    std::vector<std::map<int, double*> > *projectedCoords;
//    std::vector<std::map<int, bool> > *isProjectionOrthogonal;

/// Tolerance up to which projection is trusted
    double disTol;

    /// number of Gauss points used for computing triangle element
    int numGPsTri;

    /// number of Gauss points used for computing quad element
    int numGPsQuad;

    /// Flag on the mapping direction
    bool isMappingIGA2FEM;

    size_t numNodesSlave;
    size_t numNodesMaster;

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name The name of the mapper
     * \param[in] _meshIGA The IGAMesh
     * \param[in] _meshFE The FEMesh
     * \param[in] _disTol Tolerance up to which projection is trusted
     * \param[in] _numGPsTri The number of Gauss points used for computing triangle element
     * \param[in] _numGPsQuad The number of Gauss points used for computing quad element
     * \author Chenshen Wu
     ***********/
    IGAMortarMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE, double _disTol,
            int _numGPsTri, int _numGPsQuad, bool _isMappingIGA2FEM);

    /***********************************************************************************************
     * \brief Destructor Chenshen Wu
     ***********/
    virtual ~IGAMortarMapper();

    /***********************************************************************************************
     * \brief Initialization of the element freedom tables
     * \author Chenshen Wu
     ***********/
    void initTables();

    /***********************************************************************************************
     * \brief Fills up the array projectedCoords by performing closest point projection
     * \author Chenshen Wu
     ***********/
    void projectPointsToSurface();

    /***********************************************************************************************
     * \brief Compute matrices C_NN and C_NR
     * \author Chenshen Wu
     ***********/
    void computeCouplingMatrices();

    /***********************************************************************************************
     * \brief Computes coupling matrices in the given patch for the element which is split into more than one patches
     * \author Chenshen Wu
     ***********/
    bool computeCouplingMatrices4ClippedByPatchProjectedElement(IGAPatchSurface* _thePatch,
            int _numNodesInPatch, double* _elementInPatchIGA, double* _elementInPatchFE,
            int _elemCount, int _nShapeFuncsFE);

    /***********************************************************************************************
     * \brief Integrate the element coupling matrices and assemble them to the global one
     * \param[in] _igaPatchSurface The patch to compute the coupling matrices for
     * \param[in] _numNodes The number of nodes of the clipped polygon
     * \param[in] _polygonIGA The resulting from the clipping polygon at each knot span in the NURBS space
     * \param[in] _spanU The knot span index in the u-direction
     * \param[in] _spanV The knot span index in the v-direction
     * \param[in] _polygonFE The resulting from the clipping polygon at each knot span in the bilinear/linear space
     * \param[in] _elementIndex The global numbering of the element from the FE mesh
     * \param[in] _nShapeFuncsFE The number of the shape functions in the clipped element of the FE side
     * \author Chenshen Wu
     ***********/
    void integrate(IGAPatchSurface* _igaPatchSurface, int _numNodes, double* _polygonIGA,
            int _spanU, int _spanV, double* _polygonFE, int _elementIndex, int _nShapeFuncsFE);
    /***********************************************************************************************
     * \brief Perform consistent mapping from IGA to FE (map displacements)
     * \param[in] fieldIGA is the input data
     * \param[out] fieldFE is the output data
     * \author Chenshen Wu
     ***********/
    void consistentMapping(const double *fieldIGA, double *fieldFE);

    /***********************************************************************************************
     * \brief Perform conservative mapping from FE to IGA (map forces)
     * \param[in] fieldFE is the input data
     * \param[out] fieldIGA is the output data
     * \author Chenshen Wu
     ***********/
    void conservativeMapping(const double *fieldFE, double *fieldIGA);

    /// intern function used for mapping
private:
    /***********************************************************************************************
     * \brief Compute the span of the projected element living in _thePatch
     * \param[in] _igaPatchSurface The patch to compute the coupling matrices for
     * \param[in] _numNodes The number of nodes of the clipped polygon
     * \param[in] _polygonIGA The resulting from the clipping polygon at each knot span in the NURBS space
     * \param[out] _span An array size 4 containing [minSpanU maxSpanU minSpanV maxSpanV]
     * \author Fabien Pean
     ***********/
    void computeKnotSpanOfProjElement(IGAPatchSurface* _thePatch, int _numNodesClippedByPatchProjectedElement,
            double* _clippedByPatchProjElementFEUV, int* _span);
    /***********************************************************************************************
     * \brief Clip the input polygon by the trimming polygons of the patch
     * \param[in] _igaPatchSurface The patch to compute the coupling matrices for
     * \param[in] _numNodes The number of nodes of the polygon clipped by the knot span in input
     * \param[in] _polygonIGA The polygon clipped at each knot span as input.
     * \param[out] _numNodesClippedByTrimming The number of nodes of the output polygon after clipping by trimming
     * \param[out] _clippedByTrimming The polygon after clipping by the trimming loops of the patch
     * \author Fabien Pean
     ***********/
    void computeTrimmedPolygon(IGAPatchSurface* _thePatch,int _numNodesClippedByPatchProjectedElement,
            double* _clippedByPatchProjElementFEUV, int& _numNodesClippedByTrimming, double*& _clippedByTrimming);
    void computeTrimmedPolygon2(IGAPatchSurface* _thePatch,int _numNodesPolygonToClip,
            double* _nodesPolygonToClip, int& _numNodesClippedByTrimming, double*& _clippedByTrimming);
    /***********************************************************************************************
     * \brief Compute the local generalized coordinates for the input polygon
     * \param[in] _numNodesClippedByKnotSpanProjElementFE The number of nodes in the input polygon
     * \param[in] _clippedByKnotSpanProjElementFEUV The points of the polygon in the parametric space clipped
     * \param[in] _numNodesClippedByPatchProjectedElement The number of initial points in the polygon
     * \param[in] _clippedByPatchProjElementFEUV The initial points of the polygon in the parametric space
     * \param[in] _clippedByPatchProjElementFEWZ The initial points of the polygon in the generalized coordinates
     * \param[out] ClippedByKnotSpanProjElementFEWZ The output points of the polygon in the generalized coordinates
     * \author Fabien Pean
     ***********/
    void computeLocalElementCoord(int _numNodesClippedByKnotSpanProjElementFE, double* _clippedByKnotSpanProjElementFEUV,
    		int _numNodesClippedByPatchProjectedElement, double* _clippedByPatchProjElementFEUV,
    		double* _clippedByPatchProjElementFEWZ, double*&  ClippedByKnotSpanProjElementFEWZ);

    /// Writing output functions
public:
    /***********************************************************************************************
     * \brief Writes the projected nodes of the FE mesh onto the IGA surface into a file
     * \author Andreas Apostolatos
     ***********/
    void writeProjectedNodesOntoIGAMesh();

    /// Debugging functions
public:
    /***********************************************************************************************
     * \brief Print both coupling matrices C_NN and C_NR
     * \author Chenshen Wu
     ***********/
    void printCouplingMatrices();
    /***********************************************************************************************
     * \brief Print both coupling matrices C_NN and C_NR in file in csv format with space delimiter
     * \author Fabien Pean
     ***********/
    void printCouplingMatricesToFile();
    /***********************************************************************************************
     * \brief Check consistency of the coupling, constant field gives constant field
     * \author Fabien Pean
     ***********/
    void checkConsistency();
    /// unit test class
    friend class TestIGAMortarMapperTube;
    friend class TestIGAMortarMapperMultiPatchPlanarSurface;
    friend class TestIGAMortarMapperCylinder;

    /// Number of refined parametric locations where to find the candidate closest points for the projection
    static const int REFINED_NUM_PARAMETRIC_LOCATIONS = 10;
};
}

#endif /* IGAMORTARMAPPER_H_ */
