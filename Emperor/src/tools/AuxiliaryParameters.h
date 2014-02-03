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
 * \file AuxiliaryParameters.h
 * This file holds the class of AuxiliaryParameters.
 * \date 1/4/2013
 **************************************************************************************************/
#ifndef AUXILIARYPARAMETERS_H_
#define AUXILIARYPARAMETERS_H_
#include <string>

namespace EMPIRE {
/********//**
 * \brief Class AuxiliaryParameters provides a central place for EMPEROR wide parameters
 ***********/
class AuxiliaryParameters {
public:
    /// How many threads are used for MKL thread parallel routines
    static const int mklSetNumThreads;

    /// How many threads are used for mortar mapper thread parallel routines
    static const int mapperSetNumThreads;

    /// Machine epsilon (the difference between 1 and the least value greater than 1 that is representable).
    static const double machineEpsilon;

    /// Git hash is determined during configure by cmake
    static const std::string gitSHA1;

    /// Git tag is determined during configure by cmake
    static const std::string gitTAG;
};

} /* namespace EMPIRE */
#endif /* AUXILIARYFUNCTIONS_H_ */
