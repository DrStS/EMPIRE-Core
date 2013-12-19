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
 * \file ConnectionIO.h
 * This file holds the class ConnectionIO
 * \date 12/19/2013
 **************************************************************************************************/
#ifndef CONNECTIONIO_H_
#define CONNECTIONIO_H_

#include "EMPEROR_Enum.h"
#include <map>
#include <string>

namespace EMPIRE {

class Signal;
class DataField;
class AbstractMesh;
class ClientCode;
/********//**
 * \brief Class ConnectionIO is the input/output of connection, filter, coupling algorithm ...
 *              It contains a reference to a signal or a dataField.
 ***********/
class ConnectionIO {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] nameToClientCodeMap name to client code map, which is constructed by Emperor
     * \param[in] clientCodeName name of the client code
     * \param[in] meshName name of the mesh
     * \param[in] dataFieldName name of the data field
     * \author Tianyang Wang
     ***********/
    ConnectionIO(std::map<std::string, ClientCode*> &nameToClientCodeMap
            , std::string clientCodeName, std::string meshName, std::string dataFieldName);
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] nameToClientCodeMap name to client code map, which is constructed by Emperor
     * \param[in] clientCodeName name of the client code
     * \param[in] signalName name of the signal
     * \author Tianyang Wang
     ***********/
    ConnectionIO(std::map<std::string, ClientCode*> &nameToClientCodeMap
            , std::string clientCodeName, std::string signalName);
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _clientCode reference to the client code
     * \param[in] _mesh reference to the mesh
     * \param[in] _dataField reference to the data field
     * \author Tianyang Wang
     ***********/
    ConnectionIO(ClientCode *_clientCode, AbstractMesh *_mesh, DataField *_dataField);
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _clientCode reference to the client code
     * \param[in] _signal reference to the signal
     * \author Tianyang Wang
     ***********/
    ConnectionIO(ClientCode *_clientCode, Signal *_signal);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~ConnectionIO();
    /***********************************************************************************************
     * \brief Send the data field to the external client
     * \author Tianyang Wang
     ***********/
    void send();
    /***********************************************************************************************
     * \brief Receive the data field from the external client
     * \author Tianyang Wang
     ***********/
    void receive();
    /// type of the connectionIO (dataField or signal)
    EMPIRE_ConnectionIO_Type type;
    /// reference to the clientCode
    ClientCode *clientCode;
    /// reference to the mesh
    AbstractMesh *mesh;
    /// reference to the dataField
    DataField *dataField;
    /// reference to the signal
    Signal *signal;
private:
    /***********************************************************************************************
     * \brief Constructor, for unit test only
     * \author Tianyang Wang
     ***********/
    ConnectionIO();
    /// unit test class
    friend class ConnectionIOSetup;
};

} /* namespace EMPIRE */
#endif /* CONNECTIONIO_H_ */
