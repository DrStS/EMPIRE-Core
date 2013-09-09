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
 * \file CouplingLogicSequence.h
 * This file holds the class CouplingLogicSequence
 * \date 5/15/2012
 **************************************************************************************************/
#ifndef COUPLINGLOGICSEQUENCE_H_
#define COUPLINGLOGICSEQUENCE_H_

#include "AbstractCouplingLogic.h"

namespace EMPIRE {

/********//**
 * \brief Class CouplingLogicSequence simply calls all coupling logics in a sequence, one instance of
 *              it could be the global coupling logic setting ("coSimulation" block in the XML file)
 ***********/
class CouplingLogicSequence: public AbstractCouplingLogic {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \author Tianyang Wang
     ***********/
    CouplingLogicSequence() : AbstractCouplingLogic() {
    }
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~CouplingLogicSequence() {
    }
    /***********************************************************************************************
     * \brief Do the coupling by calling all coupling logics in a sequence
     * \author Tianyang Wang
     ***********/
    void doCoupling() {
        for (int i=0; i<couplingLogicSequence.size(); i++)
            couplingLogicSequence[i]->doCoupling();
    }
};

} /* namespace EMPIRE */
#endif /* COUPLINGLOGICSEQUENCE_H_ */
