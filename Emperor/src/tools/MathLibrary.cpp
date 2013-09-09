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
#ifdef USE_INTEL_MKL
#include <mkl.h>
#endif

#ifndef USE_INTEL_MKL
#include "cblas.h"
#endif

#include "MathLibrary.h"
using namespace std;

namespace EMPIRE {
namespace MathLibrary {

double computeDenseDotProduct(const double *vec1, const double *vec2, const int elements) {
    return cblas_ddot(elements, vec1, 1, vec2, 1);
}

double computeDenseDotProduct(const std::vector<double> &vec1, const std::vector<double> &vec2) {
    return cblas_ddot(vec1.size(), &vec1[0], 1, &vec2[0], 1);
}

} /* namespace Math */
} /* namespace EMPIRE */
