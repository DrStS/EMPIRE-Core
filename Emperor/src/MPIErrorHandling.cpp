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
#include "MPIErrorHandling.h"
#include <assert.h>
#include <sstream>

namespace EMPIRE {

std::string dummyLine = "\n================================================================\n";

//ToDO add message functions here

MPIException::MPIException(MPI_Comm *_comm, int MPIErrorCode) :
        comm(_comm), MPIErrorString(dummyLine) {

    assert(MPIErrorCode!=MPI_SUCCESS);
    char tmp[MPI_MAX_ERROR_STRING];
    int length;
    MPI_Error_string(MPIErrorCode, tmp, &length);
    MPIErrorString.append(tmp, length);
    MPIErrorString.append(dummyLine);
    std::stringstream ss;
    ss << "Error in communicator: " << std::hex << *comm;
    MPIErrorString.append(ss.str());
    MPIErrorString.append(dummyLine);
}

MPIException::~MPIException() {
}

std::string MPIException::what() {
    return MPIErrorString;
}

MPI_Comm *MPIException::getFailedCommunicator() {
    return comm;
}

static void ownMPIErrorHandlingFunction(MPI_Comm *comm, int *errorCode, ...) {
    //if (*comm != MPI_COMM_WORLD)
    throw MPIException(comm, *errorCode);
}

static MPI_Errhandler *ownMPIHandler = NULL;

MPI_Errhandler *getOwnMPIHandler() {
    if (ownMPIHandler != NULL)
        return ownMPIHandler;
    ownMPIHandler = new MPI_Errhandler;
    int tmp = MPI_Comm_create_errhandler(&ownMPIErrorHandlingFunction, ownMPIHandler);
    assert(tmp==MPI_SUCCESS);
    return ownMPIHandler;
}

} /* namespace EMPIRE */
