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
     * \param[in] absTol absolute tolerance of convergence
     * \param[in] relTol relative tolerance of convergence
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
    void addCheckResidual(double _absoluteTolerance, double _relativeTolerance,
            AbstractCouplingAlgorithm *_couplingAlgorithm, int _residualIndex);
    int getcurrentNumOfIterations();

private:
    class CheckResidual {
    public:
        CheckResidual(double _absoluteTolerance, double _relativeTolerance,
                AbstractCouplingAlgorithm *_couplingAlgorithm, int _residualIndex);
        virtual ~CheckResidual();
        void updateInitialResidual();
        double getAbsoluteResidual();
        double getRelativeResidual();
        bool isConvergent();
        void writeResidualToShell();

    private:
        AbstractCouplingAlgorithm *couplingAlgorithm;
        int residualIndex;
        const double ABS_TOL;
        const double REL_TOL;
        double initialResidual;
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

    /// the unit test classes
    friend class TestConvergenceChecker;
    friend class TestEmperor;
};

} /* namespace EMPIRE */
#endif /* CONVERGENCECHECKER_H_ */
