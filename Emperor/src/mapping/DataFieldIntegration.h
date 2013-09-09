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
/***********************************************************************************************//**
 * \file DataFieldIntegration.h
 * This file holds the class DataFieldIntegration
 * \date 7/17/2013
 **************************************************************************************************/
#ifndef DATAFIELDINTEGRATION_H_
#define DATAFIELDINTEGRATION_H_

namespace EMPIRE {
/********//**
 * \brief Class DataFieldIntegration is an operator from traction to force or vice versa
 * \author Tianyang Wang
 ***********/
class DataFieldIntegration {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _numNodes number of nodes
     * \param[in] _numElems number of elements
     * \param[in] _numNodesPerElem number of nodes per element
     * \param[in] _nodes coordinates of nodes
     * \param[in] _nodeIDs IDs of nodes
     * \param[in] _elems element table
     * \author Tianyang Wang
     ***********/
    DataFieldIntegration(int _numNodes, int _numElems, const int *_numNodesPerElem,
            const double *_nodes, const int *_nodeIDs, const int *_elems);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~DataFieldIntegration();
    /***********************************************************************************************
     * \brief Do integration (massMatrix*tractions=forces).
     * \param[in] tractions tractions (in one direction)
     * \param[out] forces forces (in one direction)
     * \author Tianyang Wang
     ***********/
    void integrate(const double *tractions, double *forces);
    /***********************************************************************************************
     * \brief Do deintegration (massMatrix^(-1)*forces=tractions).
     * \param[in] forces forces (in one direction)
     * \param[out] tractions tractions (in one direction)
     * \author Tianyang Wang
     ***********/
    void deIntegrate(const double *forces, double *tractions);
    /// defines number of threads used for MKL routines
    static int mklSetNumThreads;
private:
    /// number of nodes
    int numNodes;
    /// number of elements
    int numElems;
    /// number of nodes per element
    const int *numNodesPerElem;
    /// coordinates of nodes
    const double *nodes;
    /// IDs of nodes
    const int *nodeIDs;
    /// element table
    const int *elems;

    /// massMatrix csr format
    double *massMatrix_A;
    /// massMatrix csr format
    int *massMatrix_IA;
    /// massMatrix csr format
    int *massMatrix_JA;

    /// number of Gauss points used for computing triangle element mass matrix
    static const int numGPsMassMatrixTri;
    /// number of Gauss points used for computing quad element mass matrix
    static const int numGPsMassMatrixQuad;

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
    /// whether pardiso is initialized
    bool pardisoInitialized;
    /***********************************************************************************************
     * \brief Initialize pardiso to factorize massMatrix
     * \author Tianyang Wang
     ***********/
    void initPardiso();
    /***********************************************************************************************
     * \brief Deallocate the memory of pardiso
     * \author Tianyang Wang
     ***********/
    void deletePardiso();
};

} /* namespace EMPIRE */
#endif /* DATAFIELDINTEGRATION_H_ */
