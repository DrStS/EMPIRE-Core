/*           This file has been prepared for Doxygen automatic documentation generation.          */
/***********************************************************************************************//**
 * \mainpage
 *
* \section LICENSE
 *  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis, Stefan Sicklinger, Tianyang Wang, Munich \n
 *  All rights reserved. \n
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify \n
 *  it under the terms of the GNU General Public License as published by \n
 *  the Free Software Foundation, either version 3 of the License, or \n
 *  (at your option) any later version. \n
 *
 *  EMPIRE is distributed in the hope that it will be useful, \n
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of \n
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n
 *  GNU General Public License for more details. \n
 *
 *  You should have received a copy of the GNU General Public License \n
 *  along with EMPIRE.  If not, see <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses/</a>.
 *
 *
 * \section DESCRIPTION
 * This is the API of EMPIRE. It consists of dynamic link library and a header file called EMPIRE_API.h
 *
 *
 * \section COMPILATION
 *  export CC=icc
 *  export CXX=icpc
 *  cd build
 *  cmake ..
 * There are the following make targets available:
 * - make (compilation and linking)
 * - make clean (remove object files and executable including all folders)
 * - make doc (generates documentation) html main file is  /EMPEROR/doc/html/index.html
 * - make cleandoc (removes documentation)
 *
 *
 * \section HOWTO
 * Please find all further information on
 * <a href="http://empire.st.bv.tum.de">EMPIRE Project</a>
 *
 *
 * <EM> Note: The Makefile suppresses per default all compile and linking command output to the terminal.
 *       You may enable this information by make VEREBOSE=1</EM>
 *
 **************************************************************************************************/

/***********************************************************************************************//**
 * \file EMPIRE_API.h
 * This file defines the EMPIRE API for Co-Simulation
 * \date 2/22/2012
 **************************************************************************************************/

#ifndef EMPIRE_API_H_
#define EMPIRE_API_H_

#ifdef __cplusplus
extern "C" { ///Define extern C if C++ compiler is used
#endif

/***********************************************************************************************
 * \brief Establishes the necessary connection with the Emperor
 ***********/
void EMPIRE_API_Connect(char* inputFileName);

/***********************************************************************************************
 * \brief Get user defined text by the element name in the XML input file
 * \param[in] elementName name of the XML element
 * \return user defined text
 ***********/
char *EMPIRE_API_getUserDefinedText(char *elementName);

/***********************************************************************************************
 * \brief Send the mesh to the server
 * \param[in] name name of the mesh
 * \param[in] numNodes number of nodes
 * \param[in] numElems number of elements
 * \param[in] nodes coordinates of all nodes
 * \param[in] nodeIDs IDs of all nodes
 * \param[in] numNodesPerElem number of nodes per element
 * \param[in] elems connectivity table of all elements
 ***********/
void EMPIRE_API_sendMesh(char *name, int numNodes, int numElems, double *nodes, int *nodeIDs,
        int *numNodesPerElem, int *elems);

/***********************************************************************************************
 * \brief Send the IGA patch to the server
 * \param[in] _name name of the field
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
 ***********/
void EMPIRE_API_sendIGAPatch(int _pDegree,  int _uNumKnots, double* _uKnotVector, int _qDegree, int _vNumKnots,
        double* _vKnotVector, int _uNumControlPoints, int _vNumControlPoints, double* _cpNet, int* _nodeNet);

/***********************************************************************************************
 * \brief Send the IGA patch to the server
 * \param[in] _name name of the field
 * \param[in] _numPatches The number of the patches contained in the IGA mesh
 * \param[in] _numNodes The number of nodes of the analysis model, i.e. merged Control Points are seen as one node
 ***********/
void EMPIRE_API_sendIGAMesh(char *_name, int _numPatches, int _numNodes);

/***********************************************************************************************
 * \brief Send the IGA trimming information to the server
 * \param[in] _isTrimmed Whether the current considered patch is trimmed
 * \param[in] _numLoops The number of loops defining boundary
 * \author Fabien Pean
 ***********/
void EMPIRE_API_sendIGATrimmingInfo(int _isTrimmed, int _numLoops);
/***********************************************************************************************
 * \brief Send the IGA trimming information about patch to the server
 * \param[in] _uNumKnots The number of knots in U direction
 * \param[in] _vNumKnots The number of knots in V direction
 * \param[in] _knotSpanBelonging The array indicating the knots state, inside,trimmed,outside
 * \author Fabien Pean
 ***********/
void EMPIRE_API_sendIGATrimmingPatchInfo(int _uNumKnots, int _vNumKnots, int* _knotSpanBelonging);
/***********************************************************************************************
 * \brief Send the IGA trimming information about the loop to the server
 * \param[in] _inner whether loop is outter boundary loop or inner
 * \param[in] _numCurves The number of curves defining the loop
 * \author Fabien Pean
 ***********/
void EMPIRE_API_sendIGATrimmingLoopInfo(int _inner, int _numCurves);
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
void EMPIRE_API_sendIGATrimmingCurve(int _direction, int _pDegree, int _uNumKnots, double* _uKnotVector, int _uNumControlPoints, double* _cpNet);

/***********************************************************************************************
 * \brief Send data field to the server
 * \param[in] name name of the field
 * \param[in] sizeOfArray size of the array (data field)
 * \param[in] dataField the data field to be sent
 ***********/

void EMPIRE_API_sendDataField(char *name, int sizeOfArray, double *dataField);
/***********************************************************************************************
 * \brief Receive data field from the server
 * \param[in] name name of the field
 * \param[in] sizeOfArray size of the array (data field)
 * \param[out] dataField the data field to be received
 ***********/
void EMPIRE_API_recvDataField(char *name, int sizeOfArray, double *dataField);
/***********************************************************************************************
 * \brief Send signal to the server
 * \param[in] name name of the signal
 * \param[in] sizeOfArray size of the array (signal)
 * \param[in] signal the signal
 ***********/
void EMPIRE_API_sendSignal_double(char *name, int sizeOfArray, double *signal);
/***********************************************************************************************
 * \brief Receive signal from the server
 * \param[in] name name of the signal
 * \param[in] sizeOfArray size of the array (signal)
 * \param[in] signal the signal
 ***********/
void EMPIRE_API_recvSignal_double(char *name, int sizeOfArray, double *signal);
/***********************************************************************************************
 * \brief Receive the convergence signal of an loop
 * \return 1 means convergence, 0 means non-convergence
 ***********/
int EMPIRE_API_recvConvergenceSignal();
/***********************************************************************************************
 * \brief Send the convergence signal of an loop
 * \param[in] signal 1 means convergence, 0 means non-convergence
 ***********/
void EMPIRE_API_sendConvergenceSignal(int signal);
/***********************************************************************************************
 * \brief A simple debug function showing the content of the data field
 * \param[in] name name of the data field
 * \param[in] sizeOfArray size of the array (data field)
 * \param[in] dataField the data field to be printed
 ***********/
void EMPIRE_API_printDataField(char *name, int sizeOfArray, double *dataField);
/***********************************************************************************************
 * \brief Performs disconnection and finalization operations to the Emperor
 ***********/
void EMPIRE_API_Disconnect(void);

/// maximum length of a name string
static const int EMPIRE_API_NAME_STRING_LENGTH = 80;

#ifdef __cplusplus
}
#endif

#endif /* EMPIRE_API_H_ */
