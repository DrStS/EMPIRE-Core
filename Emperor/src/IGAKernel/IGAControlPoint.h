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
 * \file IGAControlPoint.h
 * This file holds the class IGAControlPoint.h
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef IGACONTROLPOINT_H_
#define IGACONTROLPOINT_H_

// Inclusion of user defined libraries
#include "AbstractPoint.h"

namespace EMPIRE {

/********//**
 * \brief class IGAControlPoint holds information for a Control Point in Isogeometric Analysis
 ***********/

class IGAControlPoint: public AbstractPoint {

protected:
	/// x-Cartesian coordinate of the Control Point
	double X;

	/// y-Cartesian coordinate of the Control Point
	double Y;

	/// z-Cartesian coordinate of the Control Point
	double Z;

	/// Weight of the Control Point
	double W;

	/// The constructor and the desctructor
public:
    /***********************************************************************************************
     * \brief Default Constructor
     * \author Andreas Apostolatos
     ***********/
	IGAControlPoint():AbstractPoint(0),X(0),Y(0),Z(0),W(0) { }

    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _ID The id of the point
     * \param[in] _X The x-Cartesian coordinate of the Control Point
     * \param[in] _Y The y-Cartesian coordinate of the Control Point
     * \param[in] _Z The z-Cartesian coordinate of the Control Point
     * \param[in] _W The weight of the Control Point
     * \author Andreas Apostolatos
     ***********/
	IGAControlPoint(int _ID,double _X, double _Y, double _Z, double _W):AbstractPoint(_ID),X(_X),Y(_Y),Z(_Z),W(_W) { }

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
	~IGAControlPoint() { }

	/// Get and set function
public:
    /***********************************************************************************************
     * \brief Get the x-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
	double getX() { return X;}

    /***********************************************************************************************
     * \brief Get the y-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
	double getY() { return Y;}

    /***********************************************************************************************
     * \brief Get the z-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
	double getZ() { return Z;}

    /***********************************************************************************************
     * \brief Get the weight of the point
     * \author Andreas Apostolatos
     ***********/
	double getW() { return W;}

    /***********************************************************************************************
     * \brief Set the x-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
	void setX(double _X) { X = _X; };

    /***********************************************************************************************
     * \brief Set the y-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
	void setY(double _Y) { Y = _Y; };

    /***********************************************************************************************
     * \brief Set the z-Cartesian coordinate of the point
     * \author Andreas Apostolatos
     ***********/
	void setZ(double _Z) { Z = _Z; };

    /***********************************************************************************************
     * \brief Set the weight of the point
     * \author Andreas Apostolatos
     ***********/
	void setW(double _W) { W = _W; };

};

}

#endif /* IGACONTROLPOINT_H_ */
