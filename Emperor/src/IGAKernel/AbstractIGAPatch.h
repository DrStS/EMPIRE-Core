/***********************************************************************************************//**
 * \file AbstractIGAPatch.h
 * This file holds the class AbstractIGAPatch.h
 * \date 11/3/2013
 **************************************************************************************************/

#ifndef ABSTRACTIGAPATCH_H_
#define ABSTRACTIGAPATCH_H_

namespace EMPIRE {

/********//**
 * \brief class AbstractIGAPatch is an abstraction of the classes IGAPatch1D and IGAPatch2D
 ***********/

class AbstractIGAPatch {

protected:
    /// The ID of the IGA patch
    int ID;

    /// The constructor, the virtual destructor and the copy constructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _ID The id of the IGA patch
     * \author Andreas Apostolatos
     ***********/
    AbstractIGAPatch(int _ID = 0) :
            ID(_ID) {
    }

    /***********************************************************************************************
     * \brief Virtual destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~AbstractIGAPatch() {
    }

    /***********************************************************************************************
     * \brief The copy constructor of the abstract class functions as a virtual member
     * \param[in] _abstractIGAPatch Constant reference to an object of class AbstractIGAPatch
     * \author Andreas Apostolatos
     ***********/
    AbstractIGAPatch(const AbstractIGAPatch& _abstractIGAPatch) :
            ID(_abstractIGAPatch.ID) {
    }

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Get the ID of the patch
     * \author Andreas Apostolatos
     ***********/
    inline int getId() {
        return ID;
    }

    /***********************************************************************************************
     * \brief Get the ID of the patch
     * \param[in] _ID The id of the IGA patch
     * \author Andreas Apostolatos
     ***********/
    inline void setId(int _ID) {
        ID = _ID;
    }
};

}/* namespace EMPIRE */

#endif /* ABSTRACTIGAPATCH_H_ */
