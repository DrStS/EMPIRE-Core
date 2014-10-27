/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Munich
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
 * \file MatlabIGAFileIO.h
 * This file holds functions for files read by MATLAB Code
 * \date 3/15/2013
 **************************************************************************************************/
#ifndef MATLABIGAFILEIO_H_
#define MATLABIGAFILEIO_H_

#include <string>

namespace EMPIRE {
class IGAMesh;
class DataField;
namespace MatlabIGAFileIO {
/***********************************************************************************************
 * \brief Write file
 * \param[in] igaMesh Information on the IGA mesh (polynomial orders, knot vectors, Control Points etc.)
 * \author Chenshen Wu
 ***********/
void writeIGAMesh(IGAMesh* igaMesh);

/***********************************************************************************************
 * \brief Write file
 * \param[in] dataFieldName Name of the data field
 * \param[in] _step Number of time steps
 * \param[in] _dataField The array containing the Control Point displacement field
 * \author Chenshen Wu
 ***********/
void writeVectorFieldOnCPs(std::string dataFieldName, int _step, DataField* _dataField);

} /* namespace MatlabIGAFileIO */
} /* namespace EMPIRE */

#endif /* MATLABIGAFILEIO_H_ */
