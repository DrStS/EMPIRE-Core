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
 * \file ClientCommunication.cpp
 * This file holds the class implementation of the EMPIRE API communication layer
 * The C++ API is wrapped afterwards in C see EMPIRE_API_WrapperToC.cpp
 * \date 2/22/2012
 **************************************************************************************************/
#ifndef CLIENTCOMMUNICATION_H_
#define CLIENTCOMMUNICATION_H_
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>

namespace EMPIRE {
class ClientMetaDatabase;
/********//**
 * \brief This class is the communication layer class for the EMPIRE API
 ***********/
class ClientCommunication {
public:
    /***********************************************************************************************
     * \brief return the singleton of ClientCommunication
     * \return the singleton
     * \author Stefan Sicklinger
     ***********/
    static ClientCommunication *getSingleton();
    /***********************************************************************************************
     * \brief Delete the singleton
     * \author Tianyang Wang
     ***********/
    static void free();
	/***********************************************************************************************
	 * \brief This does the communication initiation
	 *
	 * \author Stefan Sicklinger
	 ***********/
	void connect();
    /***********************************************************************************************
     * \brief This disconnects the client form the server
     *
     * \author Stefan Sicklinger
     ***********/
    void disconnect();
	/***********************************************************************************************
	 * \brief Template function for receiving data from client
	 *
	 * \param[in] size is the length of the message
	 * \param[out] message is a pointer to the message which is going to be received
	 * \author Stefan Sicklinger
	 ***********/
	template<class T> void receiveFromServerBlocking(int size, T* message) {
		/*
         #define MPI_BYTE           ...
         #define MPI_PACKED         ...
		 */
		if (isMpiCallLegal) {

			if (typeid(T) == typeid(int)) {
				MPI_Recv(message, size, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, server, &status);
			} else if (typeid(T) == typeid(double)) {
				MPI_Recv(message, size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, server, &status);
			} else if (typeid(T) == typeid(float)) {
				MPI_Recv(message, size, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, server, &status);
			} else if (typeid(T) == typeid(unsigned char)) {
				MPI_Recv(message, size, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, server,
						&status);
			} else if (typeid(T) == typeid(long double)) {
				MPI_Recv(message, size, MPI_LONG_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, server,
						&status);
			} else if (typeid(T) == typeid(short)) {
				MPI_Recv(message, size, MPI_SHORT, MPI_ANY_SOURCE, MPI_ANY_TAG, server, &status);
			} else if (typeid(T) == typeid(char)) {
				MPI_Recv(message, size, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, server, &status);
			}
		}
	}

	/***********************************************************************************************
	 * \brief Template function for sending data to client
	 *
	 * \param[in] clientName is a string which holds the name of the receiver client
	 * \param[in] message is a pointer to the message which is going to be sent
	 * \author Tianyang Wang
	 ***********/
	template<class T> void sendToServerBlocking(int size, T* message) {
		/*
         #define MPI_BYTE           ...
         #define MPI_PACKED         ...
		 */
		if (isMpiCallLegal) {
			if (typeid(T) == typeid(int)) {
				MPI_Ssend(message, size, MPI_INT, 0, 0, server);
			} else if (typeid(T) == typeid(double)) {
				MPI_Ssend(message, size, MPI_DOUBLE, 0, 0, server);
			} else if (typeid(T) == typeid(float)) {
				MPI_Ssend(message, size, MPI_FLOAT, 0, 0, server);
			} else if (typeid(T) == typeid(unsigned char)) {
				MPI_Ssend(message, size, MPI_UNSIGNED_CHAR, 0, 0, server);
			} else if (typeid(T) == typeid(long double)) {
				MPI_Ssend(message, size, MPI_LONG_DOUBLE, 0, 0, server);
			} else if (typeid(T) == typeid(short)) {
				MPI_Ssend(message, size, MPI_SHORT, 0, 0, server);
			} else if (typeid(T) == typeid(char)) {
				MPI_Ssend(message, size, MPI_CHAR, 0, 0, server);
			}
		}
	}

private:
	/// This holds a intercommunicator to the server
	MPI_Comm server;
	/// MPI status for the MPI calls
	MPI_Status status;
	/// The handle to the port file which is generated by the Emperor
	std::ifstream portFile;
	/// Stores the port name read from the port file
	std::string portName;
	/// Flag which stores the result if a communicator is a intercommunicator
	int isInterComm;
	/// Flag is MPI calls are allowed or if MPI_Finalize was already called
	bool isMpiCallLegal;
	/// Flag if client called already MPI_Init
	int isMpiInitCalledByClient;
	/// The singleton of this class
	static ClientCommunication* clientComm;
	/// The Rank of the current process
    int myRank;

	/***********************************************************************************************
	 * \brief disallow public constructor
	 * \param[in] Name of port file generated by emperor at runtime
	 * \param[in] Pointer to meta database.
	 *
	 * \author Stefan Sicklinger
	 ***********/
	ClientCommunication();
	/***********************************************************************************************
	 * \brief disallow public copy constructor
	 *
	 * \author Stefan Sicklinger
	 ***********/
	ClientCommunication(const ClientCommunication&);
	/***********************************************************************************************
	 * \brief disallow assignment operator
	 *
	 * \author Stefan Sicklinger
	 ***********/
	ClientCommunication& operator=(const ClientCommunication&);
	/***********************************************************************************************
	 * \brief Destructor
	 *
	 * \author Stefan Sicklinger
	 ***********/
	virtual ~ClientCommunication();

};

} /* namespace EMPIRE */

#endif /* CLIENTCOMMUNICATION_H_*/
