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
 * \file Aitken.h
 * This file holds the class Aitken
 * \date 5/15/2012
 **************************************************************************************************/

#ifndef AITKEN_H_
#define AITKEN_H_

#include "AbstractCouplingAlgorithm.h"

namespace EMPIRE {

/********//**
 * \brief Class Aitken uses Aitken algorithms to perform relaxation in an implicit coupling loop
 ***********/
class Aitken: public AbstractCouplingAlgorithm {
public:
    /***********************************************************************************************
     * \brief Constructor, set input and output (the input and output could be the same memory)
     * \param[in] _name name of the coupling algorithm
     * \param[in] _initialAitkenFactor the initial aitken factor
     * \author Tianyang Wang
     ***********/
    Aitken(std::string _name, double _initialAitkenFactor);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~Aitken();
    /***********************************************************************************************
     * \brief Set the input and output (they can be the same!)
     * \param[in] _input
     * \param[in] _output
     * \author Tianyang Wang
     ***********/
    virtual void setInputAndOutput(const ConnectionIO *_input, ConnectionIO *_output);
    /***********************************************************************************************
     * \brief Calculate the output by relaxation on the input
     * \author Tianyang Wang
     ***********/
    void calcNewValue();
private:
    /// the initial Aitken factor
    const double INIT_AITKEN_FACTOR;
    /// at the beginning, it is the output of last step. At the end, it is the computed at this step
    double *X_out;
    /// residual of last step
    double *R_0;
    /// aikten factor of last step
    double w_0;
    /// size of the array
    int SIZE;
    /// limit of aikten factor amplification
    static const double LIMIT;
    /// whether output numbers or not
    bool debugMe;
    /// the unit test class
    friend class TestAitken;
    friend class TestRelaxationMethods;
    friend class TestConvergenceChecker;
    /***********************************************************************************************
     * \brief Sets all the Atiken vectors to zero (needs to be done at the beginning of every
     * timestep)
     * \author Stefan Sicklinger
     ***********/
    void setZeroState();

};

} /* namespace EMPIRE */
#endif /* AITKEN_H_ */
