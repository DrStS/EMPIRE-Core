/*           This file has been prepared for Doxygen automatic documentation generation.          */
/***********************************************************************************************//**
 * \mainpage
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
 *  This is the server application of EMPIRE called Emperor
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
 * \section HOWTO
 * Please find all further information on
 * <a href="http://empire.st.bv.tum.de">EMPIRE Project</a>
 *
 *
 * <EM> Note: The Makefile suppresses per default all compile and linking command output to the terminal.
 *       You may enable this information by make VEREBOSE=1</EM>
 *
 *
 *
 **************************************************************************************************/
/***********************************************************************************************//**
 * \file main.cpp
 * This file holds the main function of Emperor.
 * \author Stefan Sicklinger, Tianyang Wang
 * \date 3/5/2012
 * \version beta
 **************************************************************************************************/

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include "Emperor.h"
#include "Message.h"
#include "AuxiliaryFunctions.h"
#include "MPIErrorHandling.h"

using namespace std;
using namespace EMPIRE;

int main(int argc, char **argv) {
    Emperor *emperor = new Emperor();

    emperor->initEnvironment(&argc, &argv);
    omp_set_nested(1); /// Enable nested parallelism
    omp_set_max_active_levels(2); /// Max two levels of nested parallelism
    omp_set_dynamic(0); /// I take full control

#pragma omp parallel num_threads(2) shared(emperor)
    {
        AuxiliaryFunctions::report_num_threads(1);
        /// Use OpemMP section construct for function parallelism (nowait = no barrier is called after the section ends)
#pragma omp sections
        {
#pragma omp section
            {
                try {
                    INFO_OUT("Open the listening thread by using openmp!");
                    emperor->startServerListening();
                } catch (MPIException &e) {
                    ERROR_OUT("Caught MPI exception!");
                    ERROR_OUT() << e.what() << endl;
                    ERROR_OUT("Error in the LISTENING thread, EMPEROR terminates!");
                    exit(EXIT_FAILURE);
                }
            }
#pragma omp section
            {
                try {
                    INFO_OUT("Open the master thread by using openmp!");
                    emperor->startServerCoupling();
                } catch (MPIException &e) {
                    ERROR_OUT("Caught MPI exception!");
                    ERROR_OUT() << e.what() << endl;
                    ERROR_OUT("Error in COUPLING thread, EMPEROR terminates!");
                    exit(EXIT_FAILURE);
                }
            }
        } /// End of sections
    } /// End of parallel section

    cout << "Stop openmp threading" << endl;
    delete emperor;
    MPI_Finalize();

    return (0);
}
