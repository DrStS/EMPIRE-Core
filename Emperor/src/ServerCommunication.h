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
 * \file ServerCommunication.cpp
 * This file holds the class implementation of the EMPIRE communication layer
 * \date 2/22/2012
 **************************************************************************************************/
#ifndef SERVERCOMMUNICATION_H_
#define SERVERCOMMUNICATION_H_
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>
#include <typeinfo>

namespace EMPIRE {
/********//**
 * \brief This class is the communication layer class for the EMPIRE (Server = Emperor)
 * \author Stefan Sicklinger
 ***********/
class MetaDatabase;
class ServerCommunication {
public:
    /***********************************************************************************************
     * \brief Initialize the serverCommunication
     * \author Tianyang Wang
     ***********/
    static void init(int *argc, char ***argv);
    /***********************************************************************************************
     * \brief return the singleton
     * \return the singleton
     * \author Tianyang Wang
     ***********/
    static ServerCommunication *getSingleton();
    /***********************************************************************************************
     * \brief Destructor
     * \author Stefan Sicklinger
     ***********/
    virtual ~ServerCommunication();
    /***********************************************************************************************
     * \brief This starts the listening thread
     * \author Stefan Sicklinger
     ***********/
    void startListening();
    /***********************************************************************************************
     * \brief Disconnect the all clients and terminate the listening thread
     * \author Stefan Sicklinger
     ***********/
    void disconnectAllClients();
    /***********************************************************************************************
     * \brief Return status of clients connections
     * \return True if all clients are connected
     * \author Stefan Sicklinger
     ***********/
    bool allClientsConnected();
    /***********************************************************************************************
     * \brief Template function for sending data to client
     * \param[in] clientName is a string which holds the name of the receiver client
     * \param[in] message is a pointer to the message which is going to be send
     * \param[in] size is the length of the message
     * \author Stefan Sicklinger
     ***********/
    template<class T>
    void sendToClientBlocking(const std::string &clientName, int size, T* message) {
        MPI_Comm client = clientNameCommMap->at(clientName);
        /*
         #define MPI_BYTE           ...
         #define MPI_PACKED         ...
         */

        if (typeid(T) == typeid(int)) {
            MPI_Ssend(message, size, MPI_INT, 0, 0, client);
        } else if (typeid(T) == typeid(double)) {
            MPI_Ssend(message, size, MPI_DOUBLE, 0, 0, client);
        } else if (typeid(T) == typeid(float)) {
            MPI_Ssend(message, size, MPI_FLOAT, 0, 0, client);
        } else if (typeid(T) == typeid(unsigned char)) {
            MPI_Ssend(message, size, MPI_UNSIGNED_CHAR, 0, 0, client);
        } else if (typeid(T) == typeid(long double)) {
            MPI_Ssend(message, size, MPI_LONG_DOUBLE, 0, 0, client);
        } else if (typeid(T) == typeid(short)) {
            MPI_Ssend(message, size, MPI_SHORT, 0, 0, client);
        } else if (typeid(T) == typeid(char)) {
            MPI_Ssend(message, size, MPI_CHAR, 0, 0, client);
        }
    }
    /***********************************************************************************************
     * \brief Template function for sending data to client
     * \param[in] clientName is a string which holds the name of the receiver client
     * \param[in] size is the length of the message
     * \param[out] message is a pointer to the message which is going to be received
     * \author Tianyang Wang
     ***********/
    template<class T>
    void receiveFromClientBlocking(const std::string &clientName, int size, T* message) {
        /*
         #define MPI_BYTE           ...
         #define MPI_PACKED         ...
         */
        MPI_Comm client = clientNameCommMap->at(clientName);
        if (typeid(T) == typeid(int)) {
            MPI_Recv(message, size, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, client, &status);
        } else if (typeid(T) == typeid(double)) {
            MPI_Recv(message, size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, client, &status);
        } else if (typeid(T) == typeid(float)) {
            MPI_Recv(message, size, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, client, &status);
        } else if (typeid(T) == typeid(unsigned char)) {
            MPI_Recv(message, size, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, client,
                    &status);
        } else if (typeid(T) == typeid(long double)) {
            MPI_Recv(message, size, MPI_LONG_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, client, &status);
        } else if (typeid(T) == typeid(short)) {
            MPI_Recv(message, size, MPI_SHORT, MPI_ANY_SOURCE, MPI_ANY_TAG, client, &status);
        } else if (typeid(T) == typeid(char)) {
            MPI_Recv(message, size, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, client, &status);
        }
    }
    /***********************************************************************************************
     * \brief Get names of all clients
     * \param[out] names is the container of names of all client codes.
     * \author Tianyang Wang
     ***********/
    void getClientNames(std::set<std::string> *names);
    /// length of the name string
    static const int NAME_STRING_LENGTH = 80;
private:
    /// The singleton of this class
    static ServerCommunication* serverComm;
    /// This holds a map of ClientName <=> intercommunicators  to the clients
    std::map<std::string, MPI_Comm>* clientNameCommMap;
    /// MPI status for the MPI calls
    MPI_Status status;
    /// The handle to the port file which is generated by the Emperor
    std::ofstream portFile;
    /// Stores the port name written to the port file
    std::string portName;
    /// Flag which stores the result if a communicator is a intercommunicator
    int isInterComm;
    /// Total number of clients
    unsigned int totalNumClients;
    /// Listening termination signal
    bool terminateListening;
    /***********************************************************************************************
     * \brief Disallow public constructor
     *        This does also the communication initialization
     * \param[in] argc pointer to an integer with holds the number of command line arguments
     * \param[in] argv pointer to an double char array with holds the command line arguments
     * \param[in] Pointer to meta database.
     * \author Stefan Sicklinger
     ***********/
    ServerCommunication(int *argc, char ***argv);
    /***********************************************************************************************
     * \brief disallow public copy constructor
     ***********/
    ServerCommunication(const ServerCommunication&);
    /***********************************************************************************************
     * \brief disallow assignment operator
     ***********/
    ServerCommunication& operator=(const ServerCommunication&);
};

} /* namespace EMPIRE */

#endif /* SERVERCOMMUNICATION_H_*/

