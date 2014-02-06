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
 * This file holds the class IJCSA
 * \date 6/2/2014
 **************************************************************************************************/

#ifndef IJCSA_H_
#define IJCSA_H_

#include "AbstractCouplingAlgorithm.h"

namespace EMPIRE {
/********//**
 * \brief Class IJCSA does Interface-Jacobian based Co-Simulation
 ***********/
class IJCSA: public AbstractCouplingAlgorithm {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name the name of the coupling algorithm
     * \author Stefan Sicklinger
     ***********/
	IJCSA(std::string _name);
    /***********************************************************************************************
     * \brief Destructor
     * \author Stefan Sicklinger
     ***********/
    virtual ~IJCSA();
    /***********************************************************************************************
     * \brief Calculate the new value of the output (improved value for next iteration/time step)
     * \author Stefan Sicklinger
     ***********/
    void calcNewValue();
    /***********************************************************************************************
     * \brief Init aitken relaxation
     * \author Stefan Sicklinger
     ***********/
    void init();
private:
    /***********************************************************************************************
     * \brief Assemble interface system
     * \author Stefan Sicklinger
     ***********/
    void assembleInterfaceJSystem();
    /// whether output numbers or not
    bool debugMe;
    /// size of global residual vector
    int globalResidualSize;
    /// current global residual vector
    double *globalResidual;
    /// temp vector of size globalResidualSize
    double *tmpVec;
    /// friend class in unit test
    friend class TestIJCSA;
};
}/* namespace EMPIRE */
#endif /* IJCSA_H_ */
