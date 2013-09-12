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
