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
 * \file ClientCode.h
 * This file holds the class ClientCode
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef CLIENTCODE_H_
#define CLIENTCODE_H_

//#include <mpi.h> --> -DMPICH_IGNORE_CXX_SEEK solves problem
#include <string>
#include <map>
#include <vector>
#include "EMPEROR_Enum.h"

namespace EMPIRE {

class ServerCommunication;
class AbstractMesh;
class Signal;

/********//**
 * \brief Class ClientCode performs communication to an EMPIRE client code (e.g. OpenFOAM, Carat++ ...)
 * \author Tianyang Wang
 ***********/
class ClientCode {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name name of the client code, which is specified in the input file
     * \author Tianyang Wang
     ***********/
    ClientCode(std::string _name);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~ClientCode();
    /***********************************************************************************************
     * \brief Return the name of the client code
     * \author Tianyang Wang
     ***********/
    inline std::string getName() {
        return name;
    }
    /***********************************************************************************************
     * \brief Add the reference to the ServerCommunication instance
     * \param[in] _serverComm pointer to the ServerCommunication instance
     * \author Tianyang Wang
     ***********/
    void setServerCommunication(ServerCommunication *_serverComm);
    /***********************************************************************************************
     * \brief Receive the mesh from a real client
     * \param[in] meshName name of the mesh to be received
     * \param[in] _triangulateAll triangulate all elements
     * \author Tianyang Wang
     ***********/
    void recvFEMesh(std::string meshName, bool triangulateAll);
    /***********************************************************************************************
     * \brief Receive the mesh from a real client
     * \param[in] meshName name of the mesh to be received
     * \author Fabien Pean, Chenshen Wu
     ***********/
    void recvIGAMesh(std::string meshName);
    /***********************************************************************************************
     * \brief Receive the data field from a real client
     * \param[in] meshName name of the mesh which owns the data field
     * \param[in] dataFieldName name of the data field
     * \author Tianyang Wang
     ***********/
    void recvDataField(std::string meshName, std::string dataFieldName);
    /***********************************************************************************************
     * \brief Send the data field to a real client
     * \param[in] meshName name of the mesh which owns the data field
     * \param[in] dataFieldName name of the data field
     * \author Tianyang Wang
     ***********/
    void sendDataField(std::string meshName, std::string dataFieldName);
    /***********************************************************************************************
     * \brief Send the convergence signal to the client
     * \param[in] convergent true if convergent, false otherwise
     * \author Tianyang Wang
     ***********/
    void sendConvergenceSignal(bool convergent);
    /***********************************************************************************************
     * \brief Receive the convergence signal from the client
     * \return true if convergent, false otherwise
     * \author Tianyang Wang
     ***********/
    bool recvConvergenceSignal();
    /***********************************************************************************************
     * \brief Get mesh by its name
     * \return a pointer to the mesh
     * \author Tianyang Wang
     ***********/
    AbstractMesh *getMeshByName(std::string meshName);
    /***********************************************************************************************
     * \brief Initialize the array
     * \param[in] signalName name of the signal
     * \param[in] size0 size0
     * \param[in] size1 size1
     * \param[in] size2 size2
     * \author Tianyang Wang
     ***********/
    void addSignal(std::string signalName, int size0, int size1, int size2);
    /***********************************************************************************************
     * \brief Receive the array from a real client
     * \param[in] signalName name of the signal
     * \author Tianyang Wang
     ***********/
    void recvSignal(std::string signalName);
    /***********************************************************************************************
     * \brief Send the array to a real client
     * \param[in] signalName name of the signal
     * \author Tianyang Wang
     ***********/
    void sendSignal(std::string signalName);
    /***********************************************************************************************
     * \brief Get array by its name
     * \return a pointer to the signal
     * \author Tianyang Wang
     ***********/
    Signal *getSignalByName(std::string signalName);


private:
    /// the name of the client
    std::string name;
    /// a pointer to the ServerCommunication instance
    ServerCommunication *serverComm;
    /// name to mesh map
    std::map<std::string, AbstractMesh*> nameToMeshMap;
    /// name to array map
    std::map<std::string, Signal*> nameToSignalMap;
    /// the unit test class
    friend class TestClientCode;
    /// the unit test class
    friend class TestLoops;
    /// the unit test class
    friend class TestEmperor;
};

} /* namespace EMPIRE */

#endif /* CLIENTCODE_H_ */
