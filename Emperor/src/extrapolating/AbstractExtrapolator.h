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
 * \file AbstractExtrapolator.h
 * This file holds the class AbstractExtrapolator
 * \date 3/5/2012
 **************************************************************************************************/
#ifndef ABSTRACTEXTRAPOLATOR_H_
#define ABSTRACTEXTRAPOLATOR_H_

#include "EMPEROR_Enum.h"
#include <vector>
#include <string>

namespace EMPIRE {

class ConnectionIO;

/********//**
 * \brief Class AbstractExtrapolator is the superclass of all extrapolators
 ***********/
class AbstractExtrapolator {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name the name
     * \author Tianyang Wang
     ***********/
    AbstractExtrapolator(std::string _name);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~AbstractExtrapolator();
    /***********************************************************************************************
     * \brief Set the input and output
     * \param[in] _input
     * \param[in] _output
     * \author Tianyang Wang
     ***********/
    virtual void setInputAndOutput(const ConnectionIO *_input, ConnectionIO *_output);
    /***********************************************************************************************
     * \brief Set the input
     * \param[in] _input
     * \author Michael Andre
     ***********/
    virtual void addInput(const ConnectionIO *_input);
    /***********************************************************************************************
     * \brief Set the output
     * \param[in] _output
     * \author Michael Andre
     ***********/
    virtual void addOutput(ConnectionIO *_output);
    /***********************************************************************************************
     * \brief Do extrapolation
     * \author Tianyang Wang
     ***********/
    virtual void extrapolate() = 0;
    /***********************************************************************************************
     * \brief Set to do extrapolation (by instance of TimeStepLoop)
     * \author Tianyang Wang
     ***********/
    void setDoExtrapolate();

protected:
    /// name
    std::string name;
    /// input
    std::vector<const ConnectionIO*> inputVec;
    /// output
    std::vector<ConnectionIO*> outputVec;
    /// a sign of doing Extrapolation or not
    bool doExtrapolate;
    // the unit test class
    friend class TestExtrapolator;
};

} /* namespace EMPIRE */
#endif /* ABSTRACTEXTRAPOLATOR_H_ */
