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
#include "ClientCommunication.h"
#include "ClientMetaDatabase.h"
#include <assert.h>
#include <string.h>
#include "EMPIRE_API.h"
using namespace std;

namespace EMPIRE {

ClientCommunication *ClientCommunication::clientComm = NULL;

ClientCommunication *ClientCommunication::getSingleton() {
    if (clientComm == NULL)
        clientComm = new ClientCommunication();
    return (clientComm);
}

void ClientCommunication::free() {
    delete clientComm;
    clientComm = NULL;
}

ClientCommunication::ClientCommunication() {

    portFile.open((const char*) (ClientMetaDatabase::getSingleton()->getServerPortFile().c_str()));
    MPI_Initialized(&isMpiInitCalledByClient);
	myRank = -1;
    if (!isMpiInitCalledByClient) {
        MPI_Init(NULL, NULL);
		MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    } else {
	    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	    if (myRank == 0) {
			cout << "EMPIRE_INFO: Client has already called MPI_Init." << endl;
		}
    }
    isMpiCallLegal = 1;
}

ClientCommunication::~ClientCommunication() {
}

void ClientCommunication::connect() {
    MPI_Status status;
    int connectionSuccessful = 0;
    if (portFile.is_open()) {
        getline(portFile, portName);
        portFile.close();
        MPI_Comm_connect((char*) (portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_WORLD, &server);
        cout << "EMPIRE_INFO: process with rank " << myRank << " is connected." << endl;
        MPI_Comm_test_inter(server, &isInterComm);
        if (isInterComm == false) {
            cout << "EMPIRE_ERROR: &server is not an intercommunicator, it is a intracommunicator" << endl;
        }

        if (myRank == 0) {
            char clientNameBuffer[EMPIRE_API_NAME_STRING_LENGTH];
            strcpy(clientNameBuffer, ClientMetaDatabase::getSingleton()->getClientName().c_str());
            MPI_Send(clientNameBuffer, EMPIRE_API_NAME_STRING_LENGTH, MPI_CHAR, 0, 0, server);
            MPI_Recv(&connectionSuccessful, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, server,
                    &status);

            if (!connectionSuccessful) {
                cout << "EMPIRE_ERROR: Emperor did not acknowledge connection" << endl;
                MPI_Comm_disconnect(&server);
                if (!isMpiInitCalledByClient) {
                    MPI_Finalize();
                }
                isMpiCallLegal = 0;
            } else {
                cout << "EMPIRE_INFO: Emperor acknowledged connection" << endl;
            }
        }
    } else
        cout << "EMPIRE_ERROR: Unable to open port file" << endl;
}

void ClientCommunication::disconnect() {

	if (isMpiCallLegal) {
		if (myRank == 0) {
			int size;
			cout << "EMPIRE_INFO: Disconnect from Emperor" << endl;
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			cout << "EMPIRE_INFO: COMM size MPI_COMM_WORLD: " << size << endl;
			MPI_Comm_size(server, &size);
			cout << "EMPIRE_INFO: COMM size server: " << size << endl;
			MPI_Comm_remote_size(server, &size);
			cout << "EMPIRE_INFO: COMM size remote server: " << size << endl;
			MPI_Comm_disconnect(&server);
		}
		if (!isMpiInitCalledByClient) {
			MPI_Finalize();
		}
	}

}

} /* namespace EMPIRE */
