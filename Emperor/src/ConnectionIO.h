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

class ConnectionIO {
public:
    ConnectionIO(std::map<std::string, ClientCode*> &nameToClientCodeMap
            , std::string clientCodeName, std::string meshName, std::string dataFieldName);
    ConnectionIO(std::map<std::string, ClientCode*> &nameToClientCodeMap
            , std::string clientCodeName, std::string signalName);
    ConnectionIO(ClientCode *_clientCode, AbstractMesh *_mesh, DataField *_dataField);
    ConnectionIO(ClientCode *_clientCode, Signal *_signal);
    virtual ~ConnectionIO();
    void send();
    void receive();

    EMPIRE_ConnectionIO_Type type;
    ClientCode *clientCode;
    /*
     * defined when type == EMPIRE_ConnectionIO_DataField
     */
    AbstractMesh *mesh;
    DataField *dataField;
    /*
     * defined when type == EMPIRE_ConnectionIO_Signal
     */
    Signal *signal;
private:
    ConnectionIO(); // for unit test only
    friend class ConnectionIOSetup;
};

} /* namespace EMPIRE */
#endif /* CONNECTIONIO_H_ */
