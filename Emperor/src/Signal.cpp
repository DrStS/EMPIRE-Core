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
#include "Signal.h"
#include <assert.h>
#include "Message.h"

using namespace std;

namespace EMPIRE {
Signal::Signal(std::string _name, int size0, int size1, int size2) :
        name(_name) {
    assert(size0 >= 1);
    assert(size1 >= 1);
    assert(size2 >= 1);
    size3D = new int[3];
    size3D[0] = size0;
    size3D[1] = size1;
    size3D[2] = size2;
    size = size0 * size1 * size2;
    array = new double[size];
    if (size0 != 1) {
        dimension = EMPIRE_Signal_3D;
    } else {
        if (size1 != 1) {
            dimension = EMPIRE_Signal_2D;
        } else {
            if (size2 != 1) {
                dimension = EMPIRE_Signal_1D;
            } else {
                dimension = EMPIRE_Signal_0D;
            }
        }
    }
}
Signal::~Signal() {
    delete[] size3D;
    delete[] array;
}

double &Signal::entry(int i, int j, int k) {
    return array[i * size3D[1] * size3D[2] + j * size3D[2] + k];
}

Message &operator<<(Message &message, Signal &signal) {
    message << "\t+" << "Signal name: " << signal.name << endl;
    for (int i = 0; i < signal.size3D[0]; i++) {
        //message << "\t\t+" << "\t" << signal.name << "[" << i << "]" << endl;
        for (int j = 0; j < signal.size3D[1]; j++) {
            message << "\t\t+" << "\t";
            for (int k = 0; k < signal.size3D[2]; k++) {
                message << signal.entry(i, j, k) << "\t";
            }
            message << endl;
        }
        message << "\t\t+" << endl;
    }
    message() << "\t+" << "---------------------------------" << endl;

    return message;
}

} /* namespace EMPIRE */
