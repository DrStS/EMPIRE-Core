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
        const double REL_TOL = 2e-6;
        const int NUM_TIME_STEP = 5;
        { // data field with relative residual
            DataField *df = new DataField("", EMPIRE_DataField_atNode, 10, EMPIRE_DataField_vector,
                    EMPIRE_DataField_field);

            ConvergenceChecker *checker = new ConvergenceChecker(-1.0, REL_TOL, 100);
            checker->debugResidual = false;
            checker->setDataField(df);
            CPPUNIT_ASSERT(checker->whichRef == EMPIRE_ConvergenceChecker_dataFieldRef);
            int size = df->numLocations * df->dimension;
            for (int i = 1; i <= NUM_TIME_STEP; i++) {
                double solution = (double) i;
                for (int j = 0; j < size; j++)
                    df->data[j] = solution;

                int count = 0;
                do {
                    count++;
                    for (int j = 0; j < size; j++)
                        df->data[j] = solution + pow(10, (double) (-count));
                } while (!checker->isConvergent());
                /*std::cout << std::endl;
                 std::cout << "Time step: " << i << ", no. of inner loop steps: " << count
                 << std::endl;*/
                CPPUNIT_ASSERT(count==8);
            }

            delete checker;
            delete df;
        }
        { // signal with relative residual
            Signal *signal = new Signal("", 2, 1, 1);
            ConvergenceChecker *checker = new ConvergenceChecker(-1.0, REL_TOL, 100);
            checker->debugResidual = false;
            checker->setSignal(signal);
            CPPUNIT_ASSERT(checker->whichRef == EMPIRE_ConvergenceChecker_signalRef);

            int size = signal->size;
            for (int i = 1; i <= NUM_TIME_STEP; i++) {
                double solution = (double) i;
                for (int j = 0; j < size; j++)
                    signal->array[j] = solution;

                int count = 0;
                do {
                    count++;
                    for (int j = 0; j < size; j++)
                        signal->array[j] = solution + pow(10, (double) (-count));
                } while (!checker->isConvergent());
                /*std::cout << std::endl;
                 std::cout << "Time step: " << i << ", no. of inner loop steps: " << count
                 << std::endl;*/
                CPPUNIT_ASSERT(count==8);
            }
            delete signal;
            delete checker;
        }
        { // coupling algorithm with relative residual
            Aitken *aitken = new Aitken("", 0.9);
            ConvergenceChecker *checker = new ConvergenceChecker(-1.0, REL_TOL, 100);
            checker->debugResidual = false;
            checker->setCouplingAlgorithm(aitken);
            CPPUNIT_ASSERT(checker->whichRef == EMPIRE_ConvergenceChecker_couplingAlgorithmRef);

            for (int i = 1; i <= NUM_TIME_STEP; i++) {
                aitken->initialResidual = 0.01;
                int count = 0;
                do {
                    count++;
                    aitken->currentResidual = pow(10, (double) (-count));
                } while (!checker->isConvergent());
                /*std::cout << std::endl;
                 std::cout << "Time step: " << i << ", no. of inner loop steps: " << count
                 << std::endl;*/
                CPPUNIT_ASSERT(count==8);
            }
            delete aitken;
            delete checker;
        }
        { // signal with absolute tolerance
            const double ABS_TOL = 2e-6;
            Signal *signal = new Signal("", 2, 1, 1);
            ConvergenceChecker *checker = new ConvergenceChecker(ABS_TOL, -1.0, 100);
            checker->debugResidual = false;
            checker->setSignal(signal);
            CPPUNIT_ASSERT(checker->whichRef == EMPIRE_ConvergenceChecker_signalRef);

            int size = signal->size;
            for (int i = 1; i <= NUM_TIME_STEP; i++) {
                double solution = (double) i;
                for (int j = 0; j < size; j++)
                    signal->array[j] = solution;

                int count = 0;
                do {
                    count++;
                    for (int j = 0; j < size; j++)
                        signal->array[j] = solution + pow(10, (double) (-count));
                } while (!checker->isConvergent());
                /*std::cout << std::endl;
                 std::cout << "Time step: " << i << ", no. of inner loop steps: " << count
                 << std::endl;*/
                CPPUNIT_ASSERT(count==7);
            }
            delete signal;
            delete checker;
        }
        { // signal with maximum number of iterations
            const double ABS_TOL = 2e-6;
            Signal *signal = new Signal("", 2, 1, 1);
            ConvergenceChecker *checker = new ConvergenceChecker(ABS_TOL, REL_TOL, 3);
            checker->debugResidual = false;
            checker->setSignal(signal);
            CPPUNIT_ASSERT(checker->whichRef == EMPIRE_ConvergenceChecker_signalRef);

            int size = signal->size;
            for (int i = 1; i <= 1; i++) {
                double solution = (double) i;
                for (int j = 0; j < size; j++)
                    signal->array[j] = solution;

                int count = 0;
                do {
                    count++;
                    for (int j = 0; j < size; j++)
                        signal->array[j] = solution + pow(10, (double) (-count));
                } while (!checker->isConvergent());
                /*std::cout << std::endl;
                 std::cout << "Time step: " << i << ", no. of inner loop steps: " << count
                 << std::endl;*/
                CPPUNIT_ASSERT(count==3);
            }
            delete signal;
            delete checker;
        }
    }
    /***********************************************************************************************
     * \brief Test case: Test the function ConvergenceChecker::calcDifferenceL2Norm().
     *        Compare the result with hand calculated result.
     ***********/
    void testCalcDifferenceL2Norm() {
        double array1[] = { 4, 4, 4, 4 };
        double array2[] = { 2, 2, 2, 2 };
        int size = 4;
        double result = ConvergenceChecker::calcDifferenceL2Norm(array1, array2, size);
        double resultByHand = 2.0;
        CPPUNIT_ASSERT(fabs(result-resultByHand)<1e-12);
    }

CPPUNIT_TEST_SUITE( TestConvergenceChecker );
        CPPUNIT_TEST( testConvergenceCheck);
        CPPUNIT_TEST( testCalcDifferenceL2Norm);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestConvergenceChecker);
