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
 * \file OptimizationLoop.h
 * This file holds the class OptimizationLoop
 * \date 2/4/2014
 **************************************************************************************************/
#ifndef OPTIMIZATIONLOOP_H_
#define OPTIMIZATIONLOOP_H_

#include "AbstractCouplingLogic.h"

namespace EMPIRE {
class ClientCode;

/********//**
 * \brief Class OptimizationLoop performs optimization loop
 ***********/
class OptimizationLoop: public AbstractCouplingLogic {
public:
    /***********************************************************************************************
     * \brief Constructor initialize memories of all containers
     * \param[in] _maxNumOfIterations maximum number of iterations
     * \author Tianyang Wang
     ***********/
    OptimizationLoop(int _maxNumOfIterations);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~OptimizationLoop();
    /***********************************************************************************************
     * \brief Do iterative coupling
     * \author Tianyang Wang
     ***********/
    void doCoupling();
    /***********************************************************************************************
     * \brief Set the convergence signal sender
     * \param[in] _convergenceSignalSender the convergence signal sender
     * \author Tianyang Wang
     ***********/
    void setConvergenceSignalSender(ClientCode *_convergenceSignalSender);
    /***********************************************************************************************
     * \brief Set the convergence signal receiver
     * \param[in] convergenceSignalReceiver the convergence signal receiver
     * \author Tianyang Wang
     ***********/
    void addConvergenceSignalReceiver(ClientCode *convergenceSignalReceiver);

private:
    /// convergence signal sender
    ClientCode *convergenceSignalSender;
    /// vector of convergence signal receivers
    std::vector<ClientCode*> convergenceSignalReceivers;
    /// current number of iterations
    int currentNumOfIterations;
    /// maximun number of iterations
    int maxNumOfIterations;

    /// unit test class
    friend class TestEmperor;
};

} /* namespace EMPIRE */
#endif /* OPTIMIZATIONLOOP_H_ */
