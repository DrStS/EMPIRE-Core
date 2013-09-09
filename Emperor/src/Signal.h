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
 * \file Signal.h
 * This file holds the class Signal
 * \date 11/15/2012
 **************************************************************************************************/

#ifndef ARRAY3D_H_
#define ARRAY3D_H_

#include <string>
#include "EMPEROR_Enum.h"

namespace EMPIRE {
/********//**
 * \brief Class Signal maintains a double array and its size
 * \author Tianyang Wang
 ***********/
class Message;
class Signal {
public:
    /***********************************************************************************************
     * \brief Constructor, allocating the storage of the array
     * \param[in] _name name of the array
     * \param[in] _size0 size of the signal
     * \param[in] _size1 size of the signal[i]
     * \param[in] _size2 size of the signal[i][j]
     * \author Tianyang Wang
     ***********/
    Signal(std::string _name, int size0, int size1, int size2);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~Signal();
    /***********************************************************************************************
     * \brief Return reference to the (i,j,k) entry (index starts from 0 in C language)
     *        TODO test the efficiency of this function!
     * \param[in] i i
     * \param[in] j j
     * \param[in] k k
     * \author Tianyang Wang
     ***********/
    double &entry(int i, int j, int k);
    /// name
    std::string name;
    /// 3D array
    double *array;
    /// the size of the array
    int *size3D;
    int size;
    /// dimension of the array
    EMPIRE_Signal_dimension dimension;
};

/***********************************************************************************************
 * \brief Allows for nice debug output later
 * \param[in] message the stream to do output
 * \param[in] signal the signal to output
 * \author Tianyang Wang
 ***********/
Message &operator<<(Message &message, Signal &signal);

} /* namespace EMPIRE */
#endif /* ARRAY3D_H_ */
