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
#include <string>
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

#include "DataField.h"
#include "Signal.h"
#include "Aitken.h"
#include "ConvergenceChecker.h"
#include "Residual.h"
#include "ConstantRelaxation.h"

#include <math.h>
#include <iostream>

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test the class ConvergenceChecker
 ***********/
class TestConvergenceChecker: public CppUnit::TestFixture {
private:

public:
    void setUp() {
    }
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief Test case: Test whether convergence happens at the expected step. We have 5 time
     *        steps, if at time step 3, the value passed to the checker at each iterative step
     *        is 3.1, 3.01, 3.001, etc.
     ***********/
    void testConvergenceCheck() {
        const int NUM_TIME_STEP = 5;
        { // relative residual
            const double REL_TOL = 2e-6;
            ConvergenceChecker *checker = new ConvergenceChecker(100);
            Residual *residual = new Residual(1);
            ConstantRelaxation *constantRelaxation = new ConstantRelaxation("", 1.0);
            constantRelaxation->addResidual(residual, 1);
            checker->addCheckResidual(-1.0, REL_TOL, constantRelaxation, 1);

            checker->debugResidual = false;
            for (int i = 1; i <= NUM_TIME_STEP; i++) {
                int count = 0;
                do {
                    count++;
                    residual->residualVectorL2Norm = pow(10, (double) (-count));
                } while (!checker->isConvergent());
                /*std::cout << std::endl;
                 std::cout << "Time step: " << i << ", no. of inner loop steps: " << count
                 << std::endl;*/
                CPPUNIT_ASSERT(count==7);
            }

            delete checker;
            delete constantRelaxation;
        }

        { // absolute residual
            const double ABS_TOL = 2e-6;
            ConvergenceChecker *checker = new ConvergenceChecker(100);
            Residual *residual = new Residual(1);
            ConstantRelaxation *constantRelaxation = new ConstantRelaxation("", 1.0);
            constantRelaxation->addResidual(residual, 1);
            checker->addCheckResidual(ABS_TOL, -1.0, constantRelaxation, 1);

            checker->debugResidual = false;
            for (int i = 1; i <= NUM_TIME_STEP; i++) {
                int count = 0;
                do {
                    count++;
                    residual->residualVectorL2Norm = pow(10, (double) (-count));
                } while (!checker->isConvergent());
                /*std::cout << std::endl;
                 std::cout << "Time step: " << i << ", no. of inner loop steps: " << count
                 << std::endl;*/
                CPPUNIT_ASSERT(count==6);
            }

            delete checker;
            delete constantRelaxation;
        }

        { // maximum number of iterations
            ConvergenceChecker *checker = new ConvergenceChecker(100.0);
            Residual *residual = new Residual(1);
            ConstantRelaxation *constantRelaxation = new ConstantRelaxation("", 1.0);
            constantRelaxation->addResidual(residual, 1);
            checker->addCheckResidual(-1.0, -1.0, constantRelaxation, 1);

            checker->debugResidual = false;
            for (int i = 1; i <= NUM_TIME_STEP; i++) {
                int count = 0;
                do {
                    count++;
                    residual->residualVectorL2Norm = pow(10, (double) (-count));
                } while (!checker->isConvergent());
                /*std::cout << std::endl;
                 std::cout << "Time step: " << i << ", no. of inner loop steps: " << count
                 << std::endl;*/
                CPPUNIT_ASSERT(count==100.0);
            }

            delete checker;
            delete constantRelaxation;
        }
    }

CPPUNIT_TEST_SUITE( TestConvergenceChecker );
        CPPUNIT_TEST( testConvergenceCheck);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestConvergenceChecker);
