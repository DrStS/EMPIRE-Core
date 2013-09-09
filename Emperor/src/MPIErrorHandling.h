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
#ifndef MPIERRORHANDLING_H_
#define MPIERRORHANDLING_H_

#include <exception>
#include <mpi.h>
#include <string>

namespace EMPIRE {

class MPIException {
public:
    MPIException(MPI_Comm *_comm, int MPIErrorCode);
    virtual ~MPIException();
    std::string what();
    MPI_Comm *getFailedCommunicator();
private:
    std::string MPIErrorString;
    MPI_Comm *comm;
};

MPI_Errhandler *getOwnMPIHandler();

} /* namespace EMPIRE */
#endif /* MPIERRORHANDLING_H_ */
