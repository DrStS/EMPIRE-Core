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

    /// Tolerance up to which projection is trusted
    double disTol;

    /// number of Gauss points used for computing triangle element
    int numGPsTri;

    /// number of Gauss points used for computing quad element
    int numGPsQuad;

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
            int _numGPsTri, int _numGPsQuad);

    /***********************************************************************************************
     * \brief Destructor
     * \author Chenshen Wu
     ***********/
    virtual ~IGAMortarMapper();

    /***********************************************************************************************
     * \brief Initializes of the element freedom tables
     * \author Chenshen Wu
     ***********/
    void initTables();

    /***********************************************************************************************
     * \brief Fills up the array projectedCoords by performing closest point projection
     * \author Chenshen Wu
     ***********/
    void projectPointsToSurface();

    /***********************************************************************************************
     * \brief Computes matrices C_NN and C_NR
     * \author Chenshen Wu
     ***********/
    void computeCouplingMatrices();

    /***********************************************************************************************
     * \brief Computes the coupling matrices in the given patch for the element
     * \param[in] _thePatch The NURBS patch
     * \param[in] _numNodesInPatch No. of nodes of the clipped by patch projected element
     * \param[in] _elementInPatchIGA The projected surface parameters on the NURBS patch of the clipped by patch projected element
     * \param[in] _elementInPatchFE The vertices of the clipped by patch projected element
     * \param[in] _elemCount The id of the unclipped element
     * \param[in] _nShapeFuncsFE The number of shape functions for the unclipped element
     * \author Chenshen Wu
     ***********/
    void computeCouplingMatrices4ClippedByPatchProjectedElement(IGAPatchSurface* _thePatch,
            int _numNodesInPatch, double* _elementInPatchIGA, double* _elementInPatchFE,
            int _elemCount, int _nShapeFuncsFE);

    /***********************************************************************************************
     * \brief Integrates the element coupling matrices and assemble them to the global one
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
     * \brief Performs consistent mapping from IGA to FE (map displacements)
     * \param[in] _fieldIGA is the input data
     * \param[in/out] _fieldFE is the output data
     * \author Chenshen Wu
     ***********/
    void consistentMapping(const double* _fieldIGA, double* _fieldFE);

    /***********************************************************************************************
     * \brief Performs conservative mapping from FE to IGA (map forces)
     * \param[in] _fieldFE is the input data
     * \param[out] _fieldIGA is the output data
     * \author Chenshen Wu
     ***********/
    void conservativeMapping(const double* _fieldFE, double* _fieldIGA);

    /***********************************************************************************************
     * \brief Prints both coupling matrices C_NN and C_NR
     * \author Chenshen Wu
     ***********/
    void printCouplingMatrices();

    /// Unit test class
    friend class TestIGAMortarMapperTube;
    friend class TestIGAMortarMapperMultiPatch;
};
}

#endif /* IGAMORTARMAPPER_H_ */
