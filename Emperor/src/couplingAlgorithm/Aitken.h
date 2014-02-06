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
 * \brief Class Aitken does a dynamic relaxation
 ***********/
class Aitken: public AbstractCouplingAlgorithm {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name the name of the coupling algorithm
     * \param[in] _initRelaxationFactor the initial relaxation factor
     * \author Stefan Sicklinger
     ***********/
    Aitken(std::string _name, double _initRelaxationFactor);
    /***********************************************************************************************
     * \brief Destructor
     * \author Stefan Sicklinger
     ***********/
    virtual ~Aitken();
    /***********************************************************************************************
     * \brief Calculate the new value of the output
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
     * \brief Reset in order to start new time step
     * \author Stefan Sicklinger
     ***********/
    void startNewTimeStep();
    /***********************************************************************************************
     * \brief Compute new optimal (Aitken) relaxation factor
     * \author Stefan Sicklinger
     ***********/
    void computeRelaxationFactor();
    /// initial relaxation factor
    const double INIT_RELAXATION_FACTOR;
    /// current relaxation factor
    double relaxationFactor;
    /// old relaxation factor
    double relaxationFactorOld;
    /// whether output numbers or not
    bool debugMe;
    /// size of global residual vector
    int globalResidualSize;
    /// current global residual vector
    double *globalResidual   ;
    /// old global residual vector
    double *globalResidualOld;
    /// temp vector of size globalResidualSize
    double *tmpVec;
    /// friend class in unit test
    friend class TestAitken;
};
}/* namespace EMPIRE */
#endif /* AITKEN_H_ */
