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
 * \file Empire.h
 * This file holds the class of Empire for the API
 * \date 3/18/2012
 **************************************************************************************************/

#ifndef EMPIRE_H_
#define EMPIRE_H_

#include <string>

namespace EMPIRE {
/********//**
 * \brief This class is the main class for EMPIRE (Client = API)
 *
 * \author Stefan Sicklinger
 ***********/
class Empire {
public:
    /***********************************************************************************************
     * \brief Constructor
     *
     * \author Stefan Sicklinger
     ***********/
    Empire();
    /***********************************************************************************************
     * \brief Destructor
     *
     * \author Stefan Sicklinger
     ***********/
    virtual ~Empire();
    /***********************************************************************************************
     * \brief connects API with Emperor
     *
     * \author Stefan Sicklinger
     ***********/
    void connect();
    /***********************************************************************************************
     * \brief disconnects API from Emperor
     *
     * \author Stefan Sicklinger
     ***********/
    void disconnect();
    /***********************************************************************************************
     * \brief Initializes Meta-database (parsing done) and ClientCommunication
     *
     * \author Stefan Sicklinger
     ***********/
    void initEnvironment(char *inputFileName);
    /***********************************************************************************************
     * \brief Get user defined text by the element name in the XML input file
     * \param[in] elementName name of the XML element
     * \return user defined text
     * \author Tianyang Wang
     ***********/
    std::string getUserDefinedText(std::string elementName);
    /***********************************************************************************************
     * \brief Send the mesh to the server
     * \param[in] numNodes number of nodes
     * \param[in] numElems number of elements
     * \param[in] nodes coordinates of all nodes
     * \param[in] nodeIDs IDs of all nodes
     * \param[in] numNodesPerElem number of nodes per element
     * \param[in] elems connectivity table of all elements
     * \author Tianyang Wang
     ***********/
    void sendMesh(int numNodes, int numElems, double *nodes, int *nodeIDs, int *numNodesPerElem,
            int *elems);
    /***********************************************************************************************
     * \brief Send the IGA patch to the server
     * \param[in] _pDegree The polynomial degree of the IGA 2D patch in the u-direction
     * \param[in] _uNumKnots The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVector The underlying knot vector of the IGA 2D patch in the u-direction
     * \param[in] _qDegree The polynomial degree of the IGA 2D patch in the v-direction
     * \param[in] _vNumKnots The number of knots for the knot vector in the v-direction
     * \param[in] _vKnotVector The underlying knot vector of the IGA 2D patch in the v-direction
     * \param[in] _uNumControlPoints The number of the Control Points for the 2D NURBS patch in the u-direction
     * \param[in] _vNumControlPoints The number of the Control Points for the 2D NURBS patch in the v-direction
     * \param[in] _cpNet The set of the Control Points related to the 2D NURBS patch
     * \param[in] _nodeNet The set of the dof index Control Points related to the 2D NURBS patch
     * \author Chenshen Wu 
     ***********/
    void sendIGAPatch(int _pDegree, int _uNumKnots, double* _uKnotVector, int _qDegree, int _vNumKnots,
            double* _vKnotVector, int _uNumControlPoints, int _vNumControlPoints, double* _cpNet, int* _nodeNet);
    /***********************************************************************************************
     * \brief Send the IGA mesh to the server
     * \param[in] _numPatches The number of the patches out of which the IGA mesh consists
     * \param[in] _numNodes The number of the Control Points which are needed for the computation of the coupling matrices
     * \author Chenshen Wu
     ***********/
    void sendIGAMesh(int _numPatches, int _numNodes);
    /***********************************************************************************************
     * \brief Send the IGA trimming information to the server
     * \param[in] _isTrimmed Whether the current considered patch is trimmed
     * \param[in] _numLoops The number of loops defining boundary
     * \author Fabien Pean
     ***********/
    void sendIGATrimmingInfo(int _isTrimmed, int _numLoops);
    /***********************************************************************************************
     * \brief Send the IGA trimming information about patch to the server
     * \param[in] _uNumKnots The number of knots in U direction
     * \param[in] _vNumKnots The number of knots in V direction
     * \param[in] _knotSpanBelonging The array indicating the knots state, inside,trimmed,outside
     * \author Fabien Pean
     ***********/
    void sendIGATrimmingPatchInfo(int _uNumKnots, int _vNumKnots, int* _knotSpanBelonging);
    /***********************************************************************************************
     * \brief Send the IGA trimming information about the loop to the server
     * \param[in] _inner whether loop is outter boundary loop or inner
     * \param[in] _numCurves The number of curves defining the loop
     * \author Fabien Pean
     ***********/
    void sendIGATrimmingLoopInfo(int _inner, int _numCurves);
    /***********************************************************************************************
     * \brief Send a IGA trimming curve to the server
     * \param[in] direction The direction of the curve if is following standard or not
     * \param[in] _pDegree The polynomial degree of the IGA 1D curve in the u-direction
     * \param[in] _uNumKnots The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVector The underlying knot vector of the IGA 1D curve in the u-direction
     * \param[in] _uNoControlPoints The number of the Control Points for the 1D NURBS patch in the u-direction
     * \param[in] _controlPointNet The set of the Control Points related to the 1D NURBS patch
     * \author Fabien Pean
     ***********/
    void sendIGATrimmingCurve(int _direction, int _pDegree, int _uNumKnots, double* _uKnotVector, int _uNumControlPoints, double* _cpNet);
    /***********************************************************************************************
     * \brief Send data field to the server
     * \param[in] sizeOfArray size of the array (data field)
     * \param[in] dataField the data field to be sent
     * \author Tianyang Wang
     ***********/
    void sendDataField(int sizeOfArray, double *dataField);
    /***********************************************************************************************
     * \brief Receive data field from the server
     * \param[in] sizeOfArray size of the array (data field)
     * \param[out] dataField the data field to be received
     * \author Tianyang Wang
     ***********/
    void recvDataField(int sizeOfArray, double *dataField);
    /***********************************************************************************************
     * \brief Send signal to the server
     * \param[in] name name of the signal
     * \param[in] sizeOfArray size of the array (signal)
     * \param[in] signal the signal
     ***********/
    void sendSignal_double(char *name, int sizeOfArray, double *signal);
    /***********************************************************************************************
     * \brief Receive signal from the server
     * \param[in] name name of the signal
     * \param[in] sizeOfArray size of the array (signal)
     * \param[in] signal the signal
     ***********/
    void recvSignal_double(char *name, int sizeOfArray, double *signal);
    /***********************************************************************************************
     * \brief Send the convergence signal of an loop
     * \param[in] signal 1 means convergence, 0 means non-convergence
     ***********/
    void sendConvergenceSignal(int signal);
    /***********************************************************************************************
     * \brief Receive the convergence signal of an loop
     * \return 1 means convergence, 0 means non-convergence
     * \author Tianyang Wang
     ***********/
    int recvConvergenceSignal();
    /***********************************************************************************************
     * \brief A simple debug function showing the content of the data field
     * \param[in] name name of the data field
     * \param[in] sizeOfArray size of the array (data field)
     * \param[in] dataField the data field to be printed
     * \author Tianyang Wang
     ***********/
    void printDataField(char *name, int sizeOfArray, double *dataField);
};

}/* namespace EMPIRE */

#endif /* EMPIRE_H_ */
