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
/******************************************************************************//**
 * \file MortarMapper.h
 * The header file of class MortarMapper.
 * \date 9/21/2011
 *********************************************************************************/
#ifndef MORTARMAPPER_H_
#define MORTARMAPPER_H_

#include <vector>
#include <set>
#include <map>
#include "MortarMath.h"
#include "AbstractMapper.h"

class ANNkd_tree;

namespace flann {
template<typename Distance> class Index;
template<class T> struct L2;
template<typename T> class Matrix;
}

namespace EMPIRE {

/***********************************************************************************************
 * \brief This class performs mortar mapping
 * ***********/
class MortarMapper : public AbstractMapper {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _slaveNumNodes number of nodes on slave side
     * \param[in] _slaveNumElems number of elements on slave side
     * \param[in] _slaveNodesPerElem number of nodes for each element (only 3 or 4, enable hybrid mesh)
     * \param[in] _slaveNodeCoors x,y,z coordinates of all slave nodes
     * \param[in] _slaveNodeNumbers the id of each slave node
     * \param[in] _slaveElemTable the information of how the elements are constructed by the nodes
     *
     * \param[in] _masterNumNodes number of nodes on master side
     * \param[in] _masterNumElems number of elements on master side
     * \param[in] _masterNodesPerElem number of nodes for each element (only 3 or 4, enable hybrid mesh)
     * \param[in] _masterNodeCoors x,y,z coordinates of all master nodes
     * \param[in] _masterNodeNumbers the id of each master node
     * \param[in] _masterElemTable the information of how the elements are constructed by the nodes
     *
     * \param[in] _oppositeSurfaceNormal whether the interface of master side and of master side have opposite normals  or not (true or false)
     * \param[in] _dual whether or not to use dual mortar (true or false)
     * \param[in] _toEnforceConsistency whether or not to enforce consistency
     * \author Tianyang Wang
     ***********/
    MortarMapper(int _slaveNumNodes, int _slaveNumElems, const int *_slaveNodesPerElem,
            const double *_slaveNodeCoors, const int *_slaveNodeNumbers, const int *_slaveElemTable,
            int _masterNumNodes, int _masterNumElems, const int *_masterNodesPerElem,
            const double *_masterNodeCoors, const int *_masterNodeNumbers,
            const int *_masterElemTable, bool _oppositeSurfaceNormal, bool _dual,
            bool _toEnforceConsistency);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~MortarMapper();
    /***********************************************************************************************
     * \brief Do consistent mapping on fields (e.g. displacements or tractions) --- C_BB * masterField = C_BA * slaveField
     * \param[in] slaveField the field of the slave side (e.g. x-displacements on all structure nodes)
     * \param[out] masterField the field of the master side (e.g. x-displacements on all fluid nodes)
     * \author Tianyang Wang
     ***********/
    void consistentMapping(const double *slaveField, double *masterField);
    /***********************************************************************************************
     * \brief Do conservative mapping on integral fields (e.g. forces)--- slaveField = C_BA^T * C_BB^(-1) * masterField
     * \param[in] masterField the field of the master side (e.g. x-forces on all fluid nodes)
     * \param[out] slaveField the field of the slave side (e.g. x-forces on all structure nodes)
     * \author Tianyang Wang
     * ***********/
    void conservativeMapping(const double *masterField, double *slaveField);
    /// defines number of threads used for MKL routines
    static int mklSetNumThreads;
    /// defines number of threads used for mapper routines
    static int mapperSetNumThreads;

private:
    /// number of nodes on slave side
    int slaveNumNodes;
    /// number of elements on slave side
    int slaveNumElems;
    /// number of nodes for each element
    const int *slaveNodesPerElem;
    /// x,y,z coordinates of all slave nodes
    const double *slaveNodeCoors;
    /// the id of each slave node
    const int *slaveNodeNumbers;
    /// the connectivity table
    const int *slaveElemTable;

    /// number of nodes on master side
    int masterNumNodes;
    /// number of elements on master side
    int masterNumElems;
    /// number of nodes for each element
    const int *masterNodesPerElem;
    /// x,y,z coordinates of all master nodes
    const double *masterNodeCoors;
    /// the id of each master node
    const int *masterNodeNumbers;
    /// the connectivity table
    const int *masterElemTable;

    /// whether the slave mesh and the master mesh have the opposite normal directions or not
    bool oppositeSurfaceNormal;
    /// whether to use dual mortar or not
    bool dual;
    /// whether or not to enforce consistency
    bool toEnforceConsistency;

    /// nearest neighbors searching tree of FLAN libraray
    flann::Index<flann::L2<double> > *FLANNkd_tree;
    /// nodes constructing the searching tree
    flann::Matrix<double> *FLANNSlaveNodes;

    /// nearest neighbors searching tree of ANN libraray
    ANNkd_tree *slaveNodesTree;
    /// nodes constructing the searching tree
    double **ANNSlaveNodes;

    /// directElemTable means the entries is not the node number, but the position in nodeCoors
    std::vector<int> **slaveDirectElemTable;
    /// directElemTable means the entries is not the node number, but the position in nodeCoors
    std::vector<int> **masterDirectElemTable;
    /// given a node, all the elements containing it are got
    std::vector<int> **slaveNodeToElemTable;
    /// compute normals of all slave elements
    double *slaveElemNormals;

    /// C_BB csr format
    double *C_BB_A;
    /// C_BB csr format
    int *C_BB_IA;
    /// C_BB csr format
    int *C_BB_JA;
    /// dual version of C_BB_A, which is diagonal
    double *C_BB_A_DUAL;

    /// C_BA csr format
    double *C_BA_A;
    /// C_BA csr format
    int *C_BA_IA;
    /// C_BA csr format
    int *C_BA_JA;
    /// dual version of C_BA_A, it shares the same IA and JA with C_BA_A
    double *C_BA_A_DUAL;

    /// number of Gauss points used for computing triangle element mass matrix
    static const int numGPsMassMatrixTri;
    /// number of Gauss points used for computing quad element mass matrix
    static const int numGPsMassMatrixQuad;
    /// number of Gauss points used for computing shape function (tri) products on a clip
    static const int numGPsOnClipTri;
    /// number of Gauss points used for computing shape function (quad) products on a clip
    static const int numGPsOnClipQuad;

    /// pardiso variable
    void *pt[64]; // this is related to internal memory management, see PARDISO manual
    /// pardiso variable
    int iparm[64];
    /// pardiso variable
    int mtype;
    /// pardiso variable
    int maxfct;
    /// pardiso variable
    int mnum;
    /// pardiso variable
    int msglvl;
    /// pardiso variable
    int neq;
    /// pardiso variable
    int nrhs;

    /// slave or master
    enum MeshLabel {
        SLAVE, MASTER
    };

    /// unit test class
    friend class TestMortarMapper;

    /***********************************************************************************************
     * \brief Compute matrix C_BB
     * \author Tianyang Wang
     ***********/
    void computeC_BB();
    /***********************************************************************************************
     * \brief Compute matrix C_BA
     * \author Tianyang Wang
     ***********/
    void computeC_BA();
    /***********************************************************************************************
     * \brief force C_BB * 1 == C_BA * 1 by modifying C_BA
     * \param[in] sparsity map of C_BA
     * \author Tianyang Wang
     ***********/
    void enforceConsistency(std::map<int, double> **sparsityMapC_BA);
    /***********************************************************************************************
     * \brief Initialize pardiso to factorize C_BB
     * \author Tianyang Wang
     ***********/
    void initPardiso();
    /***********************************************************************************************
     * \brief Deallocate the memory of pardiso
     * \author Tianyang Wang
     ***********/
    void deletePardiso();
    /***********************************************************************************************
     * \brief Initialize all tables that help referring from an element to its nodes or vice versa
     * \author Tianyang Wang
     ***********/
    void initTables();
    /***********************************************************************************************
     * \brief Deallocate the memory of all tables
     * \author Tianyang Wang
     ***********/
    void deleteTables();
    /***********************************************************************************************
     * \brief Initialize the ANN nearest neighbor searching tree
     * \author Tianyang Wang
     ***********/
    void initANNTree();
    /***********************************************************************************************
     * \brief Deallocate the memory of the searching tree
     * \author Tianyang Wang
     ***********/
    void deleteANNTree();
    /***********************************************************************************************
     * \brief Compute the Gauss quadrature of the shape function product.
     * \param[in] masterElem the master element
     * \param[in] numNodesMasterElem number of nodes of the master element
     * \param[in] slaveElem the slave element
     * \param[in] numNodesSlaveElem number of nodes of the slave element
     * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
     * \param[in] polygon the clipped polygon
     * \param[out] result the Gauss quadrature of the shape function product. 9 entries if triangle, 16 entries if quadrilateral.
     * \author Tianyang Wang
     ***********/
    void gaussQuadratureOnClip(const double *masterElem, int numNodesMasterElem,
            const double *slaveElem, int numNodesSlaveElem, int planeToProject,
            std::vector<double*> *clippedPolygon, double *result);
    /***********************************************************************************************
     * \brief Compute the search radius to find the overlapping candidates of the master element.
     * \param[in] masterElem master element
     * \param[in] numNodesMasterElem number of nodes of the master element
     * \return the square of the search radius
     * \author Tianyang Wang
     ***********/
    double computeSearchRadiusSquare(const double* masterElem, int numNodesMasterElem);
    /***********************************************************************************************
     * \brief Find the overlapping candidates of the master element with the searching radius.
     * \param[in] masterElem the master element
     * \param[in] numNodesMasterElem number of nodes of the master element
     * \param[in] masterElemNormal the normal direction of the master element
     * \param[in] searchRadiusSqr the square of the search radius
     * \param[out] neighborElems the neighboring elements which are the overlapping candidates
     * \author Tianyang Wang
     ***********/
    void findCandidates(const double* masterElem, int numNodesMasterElem,
            const double *masterElemNormal, double searchRadiusSqr, std::set<int> *neighborElems);
    /***********************************************************************************************
     * \brief Kick out the candidates who have wrong normal direction. Only here the field "oppositeSurfaceNormal" is used.
     * \brief This helps to avoid getting wrong overlap on the meshes like Turek benchmark 3D mesh.
     * \param[in] masterUnitNormal the unit normal of the master element
     * \param[in] slaveUnitNormal the unit normal of the slave element
     * \param[in] bound the bound determining the upper limit of the normal direction difference
     * \return should be kicked out or not (true or false).
     * \author Tianyang Wang
     ***********/
    bool kickOutCandidate(const double *masterUnitNormal, const double *slaveUnitNormal,
            double bound);
    /***********************************************************************************************
     * \brief Project all nodes in the elements to the plane of the element.
     * \param[in] elem the element (this can also be a point on the element plane)
     * \param[in] planeUnitNormal the unit normal of the element plane
     * \param[in] neighborElems the neighboring elements which are the overlapping candidates
     * \param[out] projections the projections of all nodes in the elements
     * \author Tianyang Wang
     ***********/
    void projectToElemPlane(const double *elem, const double *planeUnitNormal,
            std::set<int> *neighborElems, std::map<int, double*> *projections);
    /***********************************************************************************************
     * \brief Given the element index/id, return the element
     * \param[in] elemIndex the element index/id
     * \param[in] label saying it is a slave element or a master element
     * \param[out] elem the element corresponds to the element index
     * \author Tianyang Wang
     ***********/
    void getElemCoor(int elemIndex, MeshLabel label, double *elem);
    /***********************************************************************************************
     * \brief Compute the table of normals of all slave elements. Must be implemented by concrete child classes.
     * \author Tianyang Wang
     ***********/
    void computeSlaveElemNormals();
    /***********************************************************************************************
     * \brief Check C_BB and C_BA related pointers, in order to verify the logic
     * \author Tianyang Wang
     ***********/
    void checkNullPointers();
    /***********************************************************************************************
     * \brief Compute the coefficients to transfer the shape functions to dual shape functions
     * \param[in] elem the element
     * \param[in] numNodesElem number of nodes of the element
     * \param[out] coeffMatrix coefficient matrix
     * \author Tianyang Wang
     ***********/
    void computeDualCoeffMatrix(const double *elem, int numNodesElem, double *coeffMatrix);
    /********//**
     * \brief Class shapeFunctionProduct computes the shape function products of two elements
     ***********/
    class ShapeFunctionProduct: public MortarMath::IntegrandFunction {
    public:
        /***********************************************************************************************
         * \brief Constructor
         * \param[in] _masterElem the master element
         * \param[in] _numNodesMasterElem number of nodes of the master element
         * \param[in] _slaveElem the slave element
         * \param[in] _numNodesSlaveElem number of nodes of the slave element
         * \param[in] _planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
         * \author Tianyang Wang
         ***********/
        ShapeFunctionProduct(const double *_masterElem, int _numNodesMasterElem,
                const double *_slaveElem, int _numNodesSlaveElem, int _planeToProject);
        /***********************************************************************************************
         * \brief Destructor
         * \author Tianyang Wang
         ***********/
        virtual ~ShapeFunctionProduct();
        /***********************************************************************************************
         * \brief Compute shape function products on all Gauss points. These are computed here to avoid
         *        computing the local coordinates of the same Gauss point for multiple times.
         * \author Tianyang Wang
         ***********/
        void computeShapeFunctionProducts();
        /***********************************************************************************************
         * \brief set all Gauss points
         * \param[in] _gaussPoints x,y,z coordinates of all Gauss points
         * \param[in] _numGaussPoints x,y,z coordinates of all Gauss points
         * \author Tianyang Wang
         ***********/
        void setGaussPoints(const double *_gaussPoints, int _numGaussPoints);
        /***********************************************************************************************
         * \brief Set the shape functions
         * \param[in] _masterShapeFuncID shape function ID of the master element
         * \param[in] _slaveShapeFuncID shape function ID of the slave element
         * \param[in] gaussPoint x,y,z coordinates of the Gauss point
         * \author Tianyang Wang
         ***********/
        void setShapeFunctions(int _masterShapeFuncID, int _slaveShapeFuncID);
        /***********************************************************************************************
         * \brief Compute the shape function product on the Gauss point. In fact, the shape function
         *        products have been computed by computeAllShapeFunctionProducts().
         * \param[in] gaussPoint x,y,z coordinates of the Gauss point
         * \return the function value
         * \author Tianyang Wang
         ***********/
        double operator()(double *gaussPoint);
    private:
        /// the master element
        const double *masterElem;
        /// number of nodes of the master element
        int numNodesMasterElem;
        /// the slave element
        const double *slaveElem;
        /// number of nodes of the slave element
        int numNodesSlaveElem;
        /// project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
        int planeToProject;
        /// number of Gauss points
        int numGaussPoints;
        /// coordinates of all gauss points
        double * gaussPoints;
        /// different shape function products on all Gauss points
        double **shapeFunctionProducts;
        /// shape function ID of the master element
        int masterShapeFuncID;
        /// shape function ID of the slave element
        int slaveShapeFuncID;
    };
};

} /* namespace EMPIRE */
#endif /* MORTARMAPPER_H_ */
