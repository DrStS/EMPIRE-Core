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
 * \file AbstractIGABasis.h
 * This file holds the class AbstractIGABasis
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef ABSTRACTIGABASIS_H_
#define ABSTRACTIGABASIS_H_

namespace EMPIRE {

class AbstractIGABasis {

/********//**
* \brief class AbstractIGABasis is used as an abstraction of all Isogeometric Analysis bases
***********/

protected:
    /// The ID of the basis
    int ID;

    /// The constructor, the destructor and the copy constructor
public:
    /***********************************************************************************************
     * \brief Default constructor
     * \author Andreas Apostolatos
     ***********/
    AbstractIGABasis() {
        ID = 0;
    }

    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _ID The id of the basis
     * \author Andreas Apostolatos
     ***********/
    AbstractIGABasis(int _ID) {
        ID = _ID;
    }

    /***********************************************************************************************
     * \brief Virtual destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~AbstractIGABasis() {
    }

    /***********************************************************************************************
     * \brief The copy constructor of the abstract class functions as a virtual member
     * \param[in] _abstractIGABasis Constant reference to an object of class AbstractIGABasis
     * \author Andreas Apostolatos
     ***********/
    AbstractIGABasis(const AbstractIGABasis& _abstractIGABasis) :
            ID(_abstractIGABasis.ID) {
    }

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Returns the id of the basis
     * \author Andreas Apostolatos
     ***********/
    int getId() {
        return ID;
    }

    /***********************************************************************************************
     * \brief Modifies the id of the basis
     * \param[in] _ID The new ID of the basis
     * \author Andreas Apostolatos
     ***********/
    void setId(int _ID) {
        ID = _ID;
    }

};

}/* namespace EMPIRE */

#endif /* ABSTRACTIGABASIS_H_ */
