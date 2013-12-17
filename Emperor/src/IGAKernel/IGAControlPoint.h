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
 * \file IGAControlPoint.h
 * This file holds the class IGAControlPoint.h
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef IGACONTROLPOINT_H_
#define IGACONTROLPOINT_H_

// Inclusion of user defined libraries

namespace EMPIRE {

/********//**
 * \brief class IGAControlPoint holds information for a Control Point in Isogeometric Analysis
 ***********/

class IGAControlPoint {

protected:

	/// The id of the point
	int ID;

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
//	IGAControlPoint():ID(0),X(0),Y(0),Z(0),W(0) { }
	/***********************************************************************************************
	 * \brief Constructor
	 * \param[in] _ID The id of the point
	 * \param[in] _X The x-Cartesian coordinate of the Control Point
	 * \param[in] _Y The y-Cartesian coordinate of the Control Point
	 * \param[in] _Z The z-Cartesian coordinate of the Control Point
	 * \param[in] _W The weight of the Control Point
	 * \author Andreas Apostolatos
	 ***********/
	IGAControlPoint(int _ID, double _X, double _Y, double _Z, double _W) :
			ID(_ID), X(_X), Y(_Y), Z(_Z), W(_W) {
	}

	/***********************************************************************************************
	 * \brief Constructor
	 * \param[in] _ID The id of the point
	 * \param[in] _X The Cartesian coordinate of the Control Point and its weight
	 * \author Andreas Apostolatos
	 ***********/
	IGAControlPoint(int _ID, double *_X):ID(_ID),X(_X[0]),Y(_X[1]),Z(_X[2]),W(_X[3]) {}

	/***********************************************************************************************
	 * \brief Copy Constructor
	 * \param[in] _CP Pointer to a Control Point
	 * \author Chenshen Wu
	 ***********/
	IGAControlPoint(IGAControlPoint *_CP) :
			ID(_CP->getId()), X(_CP->getX()), Y(_CP->getY()), Z(_CP->getZ()), W(
					_CP->getW()) {
	}

	/***********************************************************************************************
	 * \brief Destructor
	 * \author Andreas Apostolatos
	 ***********/
	~IGAControlPoint() {
	}

	/// Get and set function
public:
	/***********************************************************************************************
	 * \brief Get the id of the point
	 * \author Andreas Apostolatos
	 ***********/
	inline int getId() {
		return ID;
	}

	/***********************************************************************************************
	 * \brief Get the x-Cartesian coordinate of the point
	 * \author Andreas Apostolatos
	 ***********/
	inline double getX() {
		return X;
	}

	/***********************************************************************************************
	 * \brief Get the y-Cartesian coordinate of the point
	 * \author Andreas Apostolatos
	 ***********/
	inline double getY() {
		return Y;
	}

	/***********************************************************************************************
	 * \brief Get the z-Cartesian coordinate of the point
	 * \author Andreas Apostolatos
	 ***********/
	inline double getZ() {
		return Z;
	}

	/***********************************************************************************************
	 * \brief Get the weight of the point
	 * \author Andreas Apostolatos
	 ***********/
	inline double getW() {
		return W;
	}

	/***********************************************************************************************
	 * \brief Set the id of the point
	 * \author Andreas Apostolatos
	 ***********/
	inline void setId(int _ID) {
		ID = _ID;
	}

	/***********************************************************************************************
	 * \brief Set the x-Cartesian coordinate of the point
	 * \author Andreas Apostolatos
	 ***********/
	inline void setX(double _X) {
		X = _X;
	}
	;

	/***********************************************************************************************
	 * \brief Set the y-Cartesian coordinate of the point
	 * \author Andreas Apostolatos
	 ***********/
	inline void setY(double _Y) {
		Y = _Y;
	}
	;

	/***********************************************************************************************
	 * \brief Set the z-Cartesian coordinate of the point
	 * \author Andreas Apostolatos
	 ***********/
	inline void setZ(double _Z) {
		Z = _Z;
	}
	;

	/***********************************************************************************************
	 * \brief Set the weight of the point
	 * \author Andreas Apostolatos
	 * WARNING!!!!!!!!!!!!!!!!!!!!!! The weight is also stored in NurbsBasis!!!!!!!!!!!!!!!!!!!
	 ***********/
	inline void setW(double _W) {
		W = _W;
	}
	;

};

}

#endif /* IGACONTROLPOINT_H_ */
