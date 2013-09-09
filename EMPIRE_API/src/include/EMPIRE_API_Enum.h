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
#ifndef EMPIRE_API_ENUM_H_
#define EMPIRE_API_ENUM_H_

/*
 * The rule of naming:
 * EMPIRE_[className]_[variableName] and
 * EMPIRE_[className]_[variableValue].
 * For example, in class "DataField", there is a member variable "dataType", it could
 * have value "scalar" or "vector", the enum below is named according to this rule.
 */
enum EMPIRE_DataField_dimension {
    EMPIRE_DataField_scalar = 1,
    EMPIRE_DataField_vector = 3
};






#endif /* EMPIRE_API_ENUM_H_ */
