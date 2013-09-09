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
 * \file IGAPatch1D.h
 * This file holds the class IGAPatch1D.h
 * \date 11/3/2013
 **************************************************************************************************/

#ifndef IGAPATCH1D_H_
#define IGAPATCH1D_H_

// Inclusion of user defined libraries
#include "AbstractIGAPatch.h"
#include "NurbsBasis1D.h"
#include "IGAControlPoint.h"

namespace EMPIRE {

/********//**
 * \brief class IGAPatch1D is a specialization of the class AbstractIGAPatch used for 1D NURBS or B-Spline patches
 ***********/

class IGAPatch1D: public AbstractIGAPatch {

protected:
    /// The basis functions of the 1D NURBS patch
    BSplineBasis1D* IGABasis;

    /// Number of Control Points
    int NoControlPoints;

    /// The set of the Control Points of the patch
    IGAControlPoint* ControlPointNet;

    /// The constructor and the destructor and the copy constructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _ID The id of the IGA patch
     * \param[in] _IDBasis The id of the underlying basis to the IGA 1D patch
     * \param[in] _pDegree The polynomial degree of the IGA 1D patch
     * \param[in] _noKnots The number of knots for the knot vector
     * \param[in] _knotVector The underlying knot vector of the IGA 1D patch
     * \param[in] _noControlPoints The number of the Control Points for the 1D NURBS patch
     * \param[in] _controlPointNet The set of the Control Points related to the 1D NURBS patch
     * \author Andreas Apostolatos
     ***********/
    IGAPatch1D(int, int, int, int, double*, int, IGAControlPoint*);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    ~IGAPatch1D();

    /***********************************************************************************************
     * \brief The copy constructor
     * \param[in] _igaPatch1D Constant reference to an object of class IGAPatch1D
     * \author Andreas Apostolatos
     ***********/
    IGAPatch1D(const IGAPatch1D&);

    /// Basis related functions
public:
    /***********************************************************************************************
     * \brief Returns the Cartesian components of a point on the NURBS surface given its parametric coordinate
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the Patch whose parameter is _uPrm
     * \param[in] _uPrm The parameter on the IGA 1D patch
     * \param[in] _knotSpanIndex The index of the knot span where the parametric coordinates lives in
     * \author Andreas Apostolatos
     ***********/
    void computeCartesianCoordinates(double*, double, int);

    /***********************************************************************************************
     * \brief Returns the Cartesian components of a point on the NURBS surface given its parametric coordinate and the basis functions
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the Patch whose parameter is _uPrm
     * \param[in] _localBasisFunctions The basis functions at _uPrm
     * \param[in] _uPrm The parameter on the IGA 1D patch
     * \param[in] _knotSpanIndex The index of the knot span where the parametric coordinates lives in
     * \author Andreas Apostolatos
     ***********/
    void computeCartesianCoordinates(double*, double*, double, int);

    /***********************************************************************************************
     * \brief Returns the base vector of the NURBS curve in the Cartesian Coordinate system given the parametric coordinate only
     * \param[in/out] _baseVector The Cartesian coordinates of the base vector at the parametric coordinate _uPrm
     * \param[in] _uPrm The parameter on the IGA 1D patch
     * \param[in] _knotSpanIndex The index of the knot span where the parametric coordinates lives in
     * \author Andreas Apostolatos
     ***********/
    void computeBaseVector(double*, double, int);

    /***********************************************************************************************
     * \brief Returns the base vector of the NURBS curve in the Cartesian Coordinate system given the parametric coordinate, the basis functions and their derivatives
     * \param[in/out] _baseVector The Cartesian coordinates of the base vector at the parametric coordinate _uPrm
     * \param[in/out] _localBasisFunctionsAndDerivatives The basis functions and their derivatives at the parametric coordinate _uPrm
     * \param[in] _uPrm The parameter on the IGA 1D patch
     * \param[in] _knotSpanIndex The index of the knot span where the parametric coordinates lives in
     * \author Andreas Apostolatos
     ***********/
    void computeBaseVector(double*, double*, double, int);

    /***********************************************************************************************
     * \brief Returns the base vector and its derivatives up to requested order given the parametric coordinate only
     * \param[in/out] _baseVectorAndDerivatives The Cartesian coordinates of the base and the acceleration vector at the parametric coordinate _uPrm
     * \param[in] _derivDegree The maximal derivative order of the base vector
     * \param[in] _uPrm The parameter on the IGA 1D patch
     * \param[in] _knotSpanIndex The index of the knot span where the parametric coordinates lives in
     * \author Andreas Apostolatos
     ***********/
    void computeBaseVectorAndDerivatives(double*, int, double, int);

    /***********************************************************************************************
     * \brief Returns the base vector and its derivatives up to requested order given the parametric coordinate, the basis functions and their derivatives
     * \param[in/out] _baseVectorAndDerivatives The Cartesian coordinates of the base and the acceleration vector at the parametric coordinate _uPrm
     * \param[in] _localBasisFunctionsAndDerivatives The basis functions and their derivatives at the parametric coordinate _uPrm
     * \param[in] _derivDegree The maximal derivative order of the base vector
     * \param[in] _uPrm The parameter on the IGA 1D patch
     * \param[in] _knotSpanIndex The index of the knot span where the parametric coordinates lives in
     * \author Andreas Apostolatos
     ***********/
    void computeBaseVectorAndDerivatives(double*, double*, int, double, int);

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Get the underlying IsoGeometric basis of the patch
     * \author Andreas Apostolatos
     ***********/
    inline BSplineBasis1D* getIGABasis() {
        return IGABasis;
    }

    /***********************************************************************************************
     * \brief Get the number of the Control Points of the patch
     * \author Andreas Apostolatos
     ***********/
    inline int getNoControlPoints() {
        return NoControlPoints;
    }

    /***********************************************************************************************
     * \brief Get the Control Points of the patch
     * \author Andreas Apostolatos
     ***********/
    inline IGAControlPoint* getControlPointNet() {
        return ControlPointNet;
    }

    /// DEBUGGING functions
public:
    /***********************************************************************************************
     * \brief Prints the Control Point net for the 1D NURBS patch
     * \author Andreas Apostolatos
     ***********/
    void printControlPointNet();
};

}/* namespace EMPIRE */

#endif /* IGAPATCH1D_H_ */
