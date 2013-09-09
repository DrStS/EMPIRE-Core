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
#include "AuxiliaryFunctions.h"

#include <iostream>
#include <algorithm>
#include <omp.h>
#include "Message.h"

namespace EMPIRE {


bool AuxiliaryFunctions::CompareStringInsensitive(std::string strFirst, std::string strSecond)
{
    /// Convert both strings to upper case by transfrom() before compare.
    transform(strFirst.begin(), strFirst.end(), strFirst.begin(), toupper);
    transform(strSecond.begin(), strSecond.end(), strSecond.begin(), toupper);
    if(strFirst == strSecond) return true; else return false;
}

void AuxiliaryFunctions::report_num_threads(int level)

{
#pragma omp single
    {
         infoOut()<<"Thread level :"<<level<<" number of threads is: "<<omp_get_num_threads()<<std::endl;
    }
}

} /* namespace EMPIRE */
