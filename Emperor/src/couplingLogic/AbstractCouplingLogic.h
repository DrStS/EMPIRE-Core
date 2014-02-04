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
 * \file AbstractCouplingLogic.h
 * This file holds the class AbstractCouplingLogic
 * \date 5/15/2012
 **************************************************************************************************/
#ifndef ABSTRACTCOUPLINGLOGIC_H_
#define ABSTRACTCOUPLINGLOGIC_H_

#include <vector>

namespace EMPIRE {
class DataOutput;
/***********************************************************************************************
 * \brief Class AbstractCouplingLogic is the mother class of all coupling logics. A coupling logic is
 *              the topmost control unit of the total coupling process. It may contain several sub
 *              coupling logics and call them in a special way (loop, etc.)
 ***********/
class AbstractCouplingLogic {
public:
    /***********************************************************************************************
     * \brief Constructor, initialize the memory of couplingLogicSequence
     * \author Tianyang Wang
     ***********/
    AbstractCouplingLogic();
    /***********************************************************************************************
     * \brief Destructor. All objects of coupling logics are destructed in the destructor of Emperor.
     * \author Tianyang Wang
     ***********/
    virtual ~AbstractCouplingLogic();
    /***********************************************************************************************
     * \brief Do the coupling
     * \author Tianyang Wang
     ***********/
    virtual void doCoupling() = 0;
    /***********************************************************************************************
     * \brief Add a coupling logic to the end of couplingLogicSequence
     * \param[in] couplingLogic the coupling logic to be added
     * \author Tianyang Wang
     ***********/
    void addCouplingLogic(AbstractCouplingLogic *couplingLogic);
    /***********************************************************************************************
     * \brief Return the size of couplingLogicSequence
     * \author Tianyang Wang
     ***********/
    int size();
    /***********************************************************************************************
     * \brief Add a dataOutput which will write data of current iteration
     * \param[in] dataOutput the data output writer
     * \author Tianyang Wang
     ***********/
    void addDataOutput(DataOutput *dataOutput);
protected:
    /// a sequence of coupling logics
    std::vector<AbstractCouplingLogic*> couplingLogicSequence;
    /// dataOutputs
    std::vector<DataOutput*> dataOutputVec;
    /// the unit test class
    friend class TestEmperor;
};

} /* namespace EMPIRE */
#endif /* ABSTRACTCOUPLINGLOGIC_H_ */
