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
 * \file AbstractPoint.h
 * This file holds the class AbstractPoint.h
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef ABSTRACTPOINT_H_
#define ABSTRACTPOINT_H_

namespace EMPIRE {

/********//**
 * \brief class AbstractPoint is an abstraction of a point in 3D space
 ***********/

class AbstractPoint {

protected:
    /// The id of the point
    int ID;

    /// The constructor and the destructor
public:
    /***********************************************************************************************
     * \brief Default Constructor
     * \author Andreas Apostolatos
     ***********/
    AbstractPoint() :
            ID(0) {
    }

    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _ID The id of the point
     * \author Andreas Apostolatos
     ***********/
    AbstractPoint(int _ID) :
            ID(_ID) {
    }

    /***********************************************************************************************
     * \brief Virtual Destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~AbstractPoint() {
    }

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Get the id of the point
     * \author Andreas Apostolatos
     ***********/
    inline int getId() {
        return ID;
    }

    /***********************************************************************************************
     * \brief Get the X-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
    virtual double getX() = 0;

    /***********************************************************************************************
     * \brief Get the Y-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
    virtual double getY() = 0;

    /***********************************************************************************************
     * \brief Get the Z-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
    virtual double getZ() = 0;

    /***********************************************************************************************
     * \brief Set the id of the point
     * \author Andreas Apostolatos
     ***********/
    inline void setId(int _ID) {
        ID = _ID;
    }

    /***********************************************************************************************
     * \brief Set the X-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
    virtual void setX(double) = 0;

    /***********************************************************************************************
     * \brief Set the Y-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
    virtual void setY(double) = 0;

    /***********************************************************************************************
     * \brief Set the Z-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
    virtual void setZ(double) = 0;
};

}/* namespace EMPIRE */

#endif /* ABSTRACTPOINT_H_ */
