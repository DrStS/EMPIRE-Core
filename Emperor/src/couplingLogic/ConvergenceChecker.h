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
    ConvergenceChecker(double absTol, double relTol, double maxNumOfIters);
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
     * \brief return currentNumOfIterations, this is only used for debugging, that this value is
     *        the same as that inside iterativeCouplingLogic
     * \return currentNumOfIterations
     * \author Tianyang Wang
     ***********/
    int getcurrentNumOfIterations();
    /***********************************************************************************************
     * \brief Set the data field
     * \param[in] _dataField the data field on which convergence is checked
     * \author Tianyang Wang
     ***********/
    void setDataField(DataField *_dataField);
    /***********************************************************************************************
     * \brief Set the signal
     * \param[in] _signal the signal on which convergence is checked
     * \author Tianyang Wang
     ***********/
    void setSignal(Signal *_signal);
    /***********************************************************************************************
     * \brief Set the coupling algorithm
     * \param[in] _couplingAlgorithm the coupling algorithm which computes the residual
     * \author Tianyang Wang
     ***********/
    void setCouplingAlgorithm(AbstractCouplingAlgorithm *_couplingAlgorithm);
    /// default absolute tolerance of convergence
    static const double DEFAULT_ABS_TOL;
    /// default relative tolerance of convergence
    static const double DEFAULT_REL_TOL;
    /// default maximum number of iterations
    static const double DEFAULT_MAX_NUM_ITERATIONS;

private:
    /// has reference to dataField or signal or couplingAlgorithm
    EMPIRE_ConvergenceChecker_whichRef whichRef;
    /// data field on which convergence is checked
    DataField *dataField;
    /// signal on which convergence is checked
    Signal *signal;
    /// coupling algorithm which computes the residual
    AbstractCouplingAlgorithm *couplingAlgorithm;
    /// absolute tolerance of convergence
    const double ABS_TOL;
    /// relative tolerance of convergence
    const double REL_TOL;
    /// maximum number of iterations
    const double MAX_NUM_ITERATIONS;
    /// current number of loops
    int currentNumOfIterations;
    /// data of last inner loop step
    double *dataLastInnerLoopStep;
    /// initial residual in fixed point iteration
    double initialResidual;
    /// current residual in fixed point iteration
    double currentResidual;
    /// whether always output the residual or not
    bool debugResidual;
    /// the file where the residual at each time step are written
    std::string residualFileName;
    /// time step number, only used when writing the residual in a file
    int timeStepNumber;
    /***********************************************************************************************
     * \brief Calculate the L2 norm of the difference
     * \return the L2 norm of the difference
     * \author Tianyang Wang
     ***********/
    static double calcDifferenceL2Norm(const double *array1, const double *array2, int size);

    /// the unit test classes
    friend class TestConvergenceChecker;
    friend class TestEmperor;
};

} /* namespace EMPIRE */
#endif /* CONVERGENCECHECKER_H_ */
