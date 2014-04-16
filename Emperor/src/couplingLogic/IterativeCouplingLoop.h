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
 * \file IterativeCouplingLoop.h
 * This file holds the class IterativeCouplingLoop
 * \date 5/15/2012
 **************************************************************************************************/
#ifndef ITERATIVECOUPLINGLOOP_H_
#define ITERATIVECOUPLINGLOOP_H_

#include <string>
#include <vector>

#include "AbstractCouplingLogic.h"

namespace EMPIRE {

class ConvergenceChecker;
class ClientCode;
class AbstractCouplingAlgorithm;

/********//**
 * \brief Class IterativeCouplingLoop performs iterative coupling of a certain time step
 ***********/
class IterativeCouplingLoop: public AbstractCouplingLogic {
public:
    /***********************************************************************************************
     * \brief Constructor initialize memories of all containers
     * \author Tianyang Wang
     ***********/
    IterativeCouplingLoop();
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~IterativeCouplingLoop();
    /***********************************************************************************************
     * \brief Do iterative coupling
     * \author Tianyang Wang
     ***********/
    void doCoupling();
    /***********************************************************************************************
     * \brief Set the convergence checker of the iterative coupling
     * \param[in] convergenceChecker the convergence checker which is initialized outside
     * \author Tianyang Wang
     ***********/
    void setConvergenceChecker(ConvergenceChecker *_convergenceChecker);
    /***********************************************************************************************
     * \brief Add the convergence observer (the loop tells the observers whether it is convergent)
     * \param[in] clientCode the observer to be added
     * \author Tianyang Wang
     ***********/
    void addConvergenceObserver(ClientCode *clientCode);
    /***********************************************************************************************
     * \brief Add the coupling algorithm (the loop tells the coupling algorithm whether to
     *        start new loop)
     * \param[in] _couplingAlgorithm the coupling algorithm to be set
     * \author Tianyang Wang
     ***********/
    void addCouplingAlgorithm(AbstractCouplingAlgorithm *_couplingAlgorithm);

private:
    /// convergence checker
    ConvergenceChecker *convergenceChecker;
    /// vector of coupling algorithms
    std::vector<AbstractCouplingAlgorithm*>  couplingAlgorithmVec;
    /// vector of convergence observers
    std::vector<ClientCode*> convergenceObserverVec;
    /// output counter
    int outputCounter;

    /***********************************************************************************************
     * \brief Broadcast the convergence signal to all clients/observers
     * \param[in] convergent the convergence signal
     * \author Tianyang Wang
     ***********/
    void broadcastConvergenceToClients(bool convergent);
    /// the unit test class
    friend class TestLoops;
    /// the unit test class
    friend class TestEmperor;
};

} /* namespace EMPIRE */
#endif /* ITERATIVECOUPLINGLOOP_H_ */
