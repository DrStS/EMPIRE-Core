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
#include "DataField.h"

namespace EMPIRE {
using namespace std;

DataField::DataField(std::string _name, EMPIRE_DataField_location _location, int _numLocations,
        EMPIRE_DataField_dimension _dimension, EMPIRE_DataField_typeOfQuantity _typeOfQuantity) :
        name(_name), location(_location), numLocations(_numLocations), dimension(_dimension), typeOfQuantity(
                _typeOfQuantity), data(new double[_numLocations * _dimension]) {
}

DataField::~DataField() {
    delete[] data;
}

Message &operator<<(Message &message, DataField &dataField) {
    message << "\t+" << "DataField name: " << dataField.name << endl;
    for (int i = 0; i < dataField.numLocations; i++) {
        message << "\t\t+" << "\t";
        for (int j = 0; j < dataField.dimension; j++) {
            message << dataField.data[i * (dataField.dimension) + j] << "\t";
        }
        message << endl;
    }
    message() << "\t+" << "---------------------------------" << endl;

    return message;
}

} /* namespace EMPIRE */
