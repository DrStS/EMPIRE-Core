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
 * \brief This class performs mortar mapping for IGA
 * ***********/
class IGAMortarMapper: public AbstractMapper {
public:
	/***********************************************************************************************
	 * \brief Constructor
	 * \param[in] _name name of the mapper
	 * \param[in] _meshS meshS, IGAMesh
	 * \param[in] _meshF meshF, FEMesh
	 * \author Chenshen Wu
	 ***********/
	IGAMortarMapper(std::string _name, IGAMesh *_meshS, FEMesh *_meshF);

	/***********************************************************************************************
	 * \brief Destructor Chenshen Wu
	 ***********/
	virtual ~IGAMortarMapper();
	/***********************************************************************************************
	 * \brief Do consistent mapping from A to B (map displacements)
	 * \param[in] fieldA is the input data
	 * \param[out] fieldB is the output data
	 * \author Chenshen Wu
	 ***********/
	void consistentMapping(const double *fieldS, double *fieldF);
	/***********************************************************************************************
	 * \brief Do conservative mapping from B to A (map forces)
	 * \param[in] fieldB is the input data
	 * \param[out] fieldA is the output data
	 * \author Chenshen Wu
	 ***********/
	void conservativeMapping(const double *fieldF, double *fieldS);

	void printCouplingMatrix();

private:
	/// IGA Mesh
	IGAMesh *meshS;

	/// Fluid Mesh
	FEMesh *meshF;

	MathLibrary::SparseMatrix<double> *C_NN;

	MathLibrary::SparseMatrix<double> *C_NR;

	IGAMortarMath::GaussQuadratureOnTriangle *gaussTriangle;
	IGAMortarMath::GaussQuadratureOnQuad *gaussQuad;

	/// coordinates of parameter space of the projected nodes on the surface
	std::vector<std::map<int, double*> > *projectedCoords;

	void projectPointsToSurface();
//    /// number of Gauss points used for computing triangle element mass matrix
//	static const int numGPsMassMatrixTri;
//	/// number of Gauss points used for computing quad element mass matrix
//	static const int numGPsMassMatrixQuad;

//	// number of Gauss points used for computing shape function (tri) products on a clip
//	static const int numGPsOnClipTri;
//	// number of Gauss points used for computing shape function (quad) products on a clip
//	static const int numGPsOnClipQuad;

/// unit test class
	friend class TestIGAMortarMapper;

	/***********************************************************************************************
	 * \brief Compute matrix C_NN and C_NR
	 * \author Chenshen Wu
	 ***********/
	void computeCouplingMatrix();
	void initTables();

	void integrate(IGAPatchSurface* thePatch, double *polygonIGA, int spanU, int spanV, double *polygonFE,
			int elementIndex, int numNodes, int nShapeFuncsFE);

	int **meshFDirectElemTable;

	const static double disTol = 1e-6;

};
}

#endif /* IGAMORTARMAPPER_H_ */
