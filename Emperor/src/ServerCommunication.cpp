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
#include "ServerCommunication.h"
#include "Message.h"
#include "MetaDatabase.h"
#include "MPIErrorHandling.h"

#include <assert.h>
using namespace std;

namespace EMPIRE {

ServerCommunication *ServerCommunication::serverComm = NULL;

void ServerCommunication::init(int *argc, char ***argv) {
    assert(serverComm == NULL);
    serverComm = new ServerCommunication(argc, argv);
}

ServerCommunication *ServerCommunication::getSingleton() {
    assert(serverComm != NULL);
    return serverComm;
}

ServerCommunication::ServerCommunication(int *argc, char ***argv) {
    terminateListening = false;
    clientNameCommMap = new map<string, MPI_Comm>;
    int providedThreadSupport;
    assert(MetaDatabase::getSingleton() != NULL);
    totalNumClients = MetaDatabase::getSingleton()->settingClientCodeVec.size();
    MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &providedThreadSupport);
    if (MPI_THREAD_MULTIPLE != providedThreadSupport) {
        WARNING_OUT() << "Requested MPI thread support is not guaranteed." << endl;
    }
    portFile.open((const char*) (MetaDatabase::getSingleton()->serverPortFile.c_str()));
    if (portFile.is_open()) {
        char portNameChar[MPI_MAX_PORT_NAME];
        MPI_Open_port(MPI_INFO_NULL, portNameChar);
        portName = portNameChar;
        portFile << portName << endl;
        portFile.close();
    } else {
        ERROR_OUT() << "Unable to open port file " << MetaDatabase::getSingleton()->serverPortFile
                << endl;
    }
}

ServerCommunication::~ServerCommunication() {
    delete clientNameCommMap;
}

void ServerCommunication::startListening() {
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, *getOwnMPIHandler()); // set own error handler to MPI_COMM_WORLD
    INFO_OUT() << "intra-communicator: " << "MPI_COMM_WORLD" << "---" << hex << MPI_COMM_WORLD << endl;
    int connectionSuccessful = 0;
    while (1) {
        /// set up inter-communicator for a connecting client
        MPI_Status status;
        MPI_Comm interClient = MPI_COMM_NULL;
        MPI_Comm world = MPI_COMM_NULL; // a duplication of MPI_COMM_WORLD
        MPI_Comm_dup(MPI_COMM_WORLD, &world);
        MPI_Comm_set_errhandler(world, MPI_ERRORS_RETURN); // we set the error handler of it to MPI_ERRORS_RETURN
        MPI_Comm_accept((char*) (portName.c_str()), MPI_INFO_NULL, 0, world, &interClient);
        if (terminateListening == true) {
            INFO_OUT() << "terminateListening" << endl;
            break;
        }
        MPI_Comm_set_errhandler(interClient, *getOwnMPIHandler()); // set own error handler to this inter-communicator
        MPI_Comm_test_inter(interClient, &isInterComm);
        if (isInterComm == false) {
            ERROR_OUT() << "&server is not an inter-communicator, it is a intra-communicator" << endl;
        }

        /// Receive client name
        char clientName[NAME_STRING_LENGTH];
        MPI_Recv(clientName, NAME_STRING_LENGTH, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, interClient,
                &status);
        string clientNameS(clientName);

        /// Add client into clientNameCommMap
        if (MetaDatabase::getSingleton()->checkForClientCodeName(clientNameS)) { /// Check if client code is defined in input file
            bool inserted = clientNameCommMap->insert(
                    pair<string, MPI_Comm>(clientNameS, interClient)).second;
            if (!inserted) { // this means the client with the same name has been connected
                WARNING_OUT() << "Client " << clientNameS << " is connected again!!!" << " The old connection is removed!!!" << endl;
                int MPIerror = MPI_Comm_free(&(clientNameCommMap->at(clientNameS)));
                if (MPIerror != MPI_SUCCESS) {
                    WARNING_OUT() << "MPI_Comm_free fails!" << endl;
                } else {
                    WARNING_OUT() << "MPI_Comm_free done!" << endl;
                }
                clientNameCommMap->erase(clientNameS);
                inserted = clientNameCommMap->insert(
                        pair<string, MPI_Comm>(clientNameS, interClient)).second;
                assert(inserted == true);
            }
            /// Tell client that everything is okay
            connectionSuccessful = 1;
            MPI_Send(&connectionSuccessful, 1, MPI_INT, 0, 0, interClient);
            INFO_OUT() << "Client " << clientNameS << " connected" << ", inter-communicator: " << hex << interClient << ", intra-communicator: " << hex << world << endl;
        } else {
            /// Tell client to disconnect
            connectionSuccessful = 0;
            MPI_Send(&connectionSuccessful, 1, MPI_INT, 0, 0, interClient);
            /// Disconnect client as it is not defined in the input file
            MPI_Comm_disconnect(&interClient);
            WARNING_OUT()
                    << "Client "
                    << clientNameS
                    << "is not defined in the input file. Hence, the client will not be connected to the emperor!"
                    << endl;
        }
        INFO_OUT() << dec;
    }
}

void ServerCommunication::disconnectAllClients() {
    for (map<string, MPI_Comm>::iterator it = clientNameCommMap->begin();
            it != clientNameCommMap->end(); it++) {
        MPI_Comm_disconnect(&(it->second));
    }
    terminateListening = true;
    /// Connect to your self in order to terminate listening
    MPI_Comm dummy;
    MPI_Comm_connect((char*) (portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_WORLD, &dummy);
    MPI_Comm_disconnect(&dummy);
    MPI_Close_port((char*) (portName.c_str()));
}

bool ServerCommunication::allClientsConnected() {
    if (totalNumClients == clientNameCommMap->size()) {
        return (true);
    }
    return (false);
}

void ServerCommunication::getClientNames(set<string> *names) {
    for (map<string, MPI_Comm>::iterator it = clientNameCommMap->begin();
            it != clientNameCommMap->end(); it++) {
        names->insert(it->first);
    }
}

} /* namespace EMPIRE */
