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

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name The name of the mapper
     * \param[in] _meshIGA The IGAMesh
     * \param[in] _meshFE The FEMesh
     * \author Chenshen Wu
     ***********/
    IGAMortarMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE);

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
     * \brief Do consistent mapping from IGA to FE (map displacements)
     * \param[in] fieldIGA is the input data
     * \param[out] fieldFE is the output data
     * \author Chenshen Wu
     ***********/
    void consistentMapping(const double *fieldIGA, double *fieldFE);

    /***********************************************************************************************
     * \brief Do conservative mapping from FE to IGA (map forces)
     * \param[in] fieldFE is the input data
     * \param[out] fieldIGA is the output data
     * \author Chenshen Wu
     ***********/
    void conservativeMapping(const double *fieldFE, double *fieldIGA);

    /***********************************************************************************************
     * \brief Print both coupling matrices C_NN and C_NR
     * \author Chenshen Wu
     ***********/
    void printCouplingMatrices();

    /// Tolerance up to which projection is trusted
    const static double disTol = 1e-6;

    /// number of Gauss points used for computing triangle element
    static const int numGPsTri = 6;

    /// number of Gauss points used for computing quad element
    static const int numGPsQuad = 25;

    /// unit test class
    friend class TestIGAMortarMapper;
};
}

#endif /* IGAMORTARMAPPER_H_ */
