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
 * \file ConvergenceChecker.h
 * This file holds the class ConvergenceChecker
 * \date 5/15/2012
 **************************************************************************************************/
#ifndef CONVERGENCECHECKER_H_
#define CONVERGENCECHECKER_H_

#include <string>
#include <vector>
#include "EMPEROR_Enum.h"

namespace EMPIRE {

class DataField;
class Signal;
class AbstractCouplingAlgorithm;

/********//**
 * \brief Class ConvergenceChecker checks whether convergence is got or not
 ***********/
class ConvergenceChecker {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] maxNumOfIters maximum number of iterations
     * \author Tianyang Wang
     ***********/
    ConvergenceChecker(double maxNumOfIters);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~ConvergenceChecker();
    /***********************************************************************************************
     * \brief Calculate whether there is convergence
     * \return convergence signal
     * \author Tianyang Wang
     ***********/
    bool isConvergent();
    /***********************************************************************************************
     * \brief Add a checkResidual object to checkResiduals
     * \param[in] _absoluteTolerance the absolute residual
     * \param[in] _relativeTolerance the relative residual
     * \param[in] _couplingAlgorithm the coupling algorithm which contains the residual
     * \param[in] _residualIndex the index of the residual
     * \author Tianyang Wang
     ***********/
    void addCheckResidual(double _absoluteTolerance, double _relativeTolerance,
            AbstractCouplingAlgorithm *_couplingAlgorithm, int _residualIndex);
    /***********************************************************************************************
     * \brief Get the current number of iterations
     * \return the current number of iterations
     * \author Tianyang Wang
     ***********/
    int getCurrentNumOfIterations();

private:
    /********//**
     * \brief Class ConvergenceChecker checks whether convergence is got or not
     ***********/
    class CheckResidual {
    public:
        /***********************************************************************************************
         * \brief Constructor
         * \param[in] _absoluteTolerance the absolute residual
         * \param[in] _relativeTolerance the relative residual
         * \param[in] _couplingAlgorithm the coupling algorithm which contains the residual
         * \param[in] _residualIndex the index of the residual
         * \author Tianyang Wang
         ***********/
        CheckResidual(double _absoluteTolerance, double _relativeTolerance,
                AbstractCouplingAlgorithm *_couplingAlgorithm, int _residualIndex);
        /***********************************************************************************************
         * \brief Destructor
         * \author Tianyang Wang
         ***********/
        virtual ~CheckResidual();
        /***********************************************************************************************
         * \brief update the initial residual
         * \author Tianyang Wang
         ***********/
        void updateInitialResidual();
        /***********************************************************************************************
         * \brief get the absolute residual
         * \return the absolute residual
         * \author Tianyang Wang
         ***********/
        double getAbsoluteResidual();
        /***********************************************************************************************
         * \brief get the relative residual
         * \return the relative residual
         * \author Tianyang Wang
         ***********/
        double getRelativeResidual();
        /***********************************************************************************************
         * \brief is convergent?
         * \return true if convergent, otherwise false
         * \author Tianyang Wang
         ***********/
        bool isConvergent();
        /***********************************************************************************************
         * \brief write residual to the shell
         * \author Tianyang Wang
         ***********/
        void writeResidualToShell();
    private:
        /// reference to the coupling algorithm
        AbstractCouplingAlgorithm *couplingAlgorithm;
        /// index of the residual
        int residualIndex;
        /// absolute tolerance
        const double ABS_TOL;
        /// relative tolerance
        const double REL_TOL;
        /// initial residual
        double initialResidual;
        /// the unit test class
        friend class TestEmperor;
    };
    /// maximum number of iterations
    const double MAX_NUM_ITERATIONS;
    /// vector of checkResiduals
    std::vector<CheckResidual*> checkResiduals;
    /// current number of loops
    int currentNumOfIterations;
    /// whether always output the residual or not
    bool debugResidual;
    /// the file where the residual at each time step are written
    std::string residualFileName;
    /// time step number, only used when writing the residual in a file
    int timeStepNumber;

    /// the unit test class
    friend class TestConvergenceChecker;
    /// the unit test class
    friend class TestEmperor;
};

} /* namespace EMPIRE */
#endif /* CONVERGENCECHECKER_H_ */
