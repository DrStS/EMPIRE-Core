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
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

#include <string>
#include <math.h>
#include <iostream>

#include "ConstantRelaxation.h"
#include "DataField.h"
#include "Signal.h"
#include "EMPEROR_Enum.h"
#include "ComaAitken.h"
#include "Message.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"
#include "Residual.h"

namespace EMPIRE {

using namespace std;

class TestConstantRelaxation: public CppUnit::TestFixture {
private:
    DataField *dfIn;
    DataField *dfOut;
    ConnectionIO *in1;
    ConnectionIO *in2;
    ConnectionIO *out;
    double coeffIn;
    double coeffOut;
    Residual *residual;
    ConstantRelaxation *constantRelaxation;
    double relaxationFactor;

public:
    void setUp() {
        relaxationFactor = 0.1;

        dfIn = new DataField("input", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        dfOut = new DataField("output", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        in1 = ConnectionIOSetup::constructDummyConnectionIO(dfIn);
        in2 = ConnectionIOSetup::constructDummyConnectionIO(dfIn);
        out = ConnectionIOSetup::constructDummyConnectionIO(dfOut);
        coeffIn = -1.0;
        coeffOut = 1.0;

        residual = new Residual(1);
        residual->addComponent(coeffIn, "iterationBeginning", in1);
        residual->addComponent(coeffOut, "iterationEnd", out);
        residual->init();

        constantRelaxation = new ConstantRelaxation("", relaxationFactor);
        constantRelaxation->addResidual(residual, 1);
        constantRelaxation->addOutput(in2, 1); // using in1 will cause segmentation fault

        constantRelaxation->debugMe = false;
    }
    void tearDown() {
        delete dfIn;
        delete dfOut;
        delete constantRelaxation;
    }
    /*
     * Test case:
     * Test check whether the result is as expected
     * In and out are: (10 20) (11 21) (12 22) (13 23) (14 24)
     * output from constant relaxtion should be: 11 12 13 14 15
     */
    void testRelaxation() {
        const int SIZE = dfIn->dimension * dfIn->numLocations;
        int NUM_INNER_LOOPS = 5;
        const double EPS = 1E-10;
        int NUM_OUTER_LOOPS = 2; // can be used to test memory leak
        // test case 1, all entries in an array are equal (e.g. 1,1,1, ...)
        for (int time = 0; time < NUM_OUTER_LOOPS; time++) {
            // initial value
            for (int j = 0; j < SIZE; j++) {
                dfIn->data[j] = 10.0;
            }
            // iterative loop ( dfIn = dfIn + omega * (dfOut - dfIn) )
            for (int i = 0; i < NUM_INNER_LOOPS; i++) {
                constantRelaxation->updateAtIterationBeginning();
                for (int j = 0; j < SIZE; j++) {
                    dfOut->data[j] = 20.0 + double(i);
                }
                constantRelaxation->updateAtIterationEnd();
                //debugOut << *dfIn;
                if (i == 0)
                    constantRelaxation->setNewLoop(); // useful only for aitken
                constantRelaxation->calcNewValue();
                //debugOut << *dfOut;
                double valueRef = 10.0 + double(i) + 1.0;

                for (int j = 0; j < SIZE; j++) {
                    CPPUNIT_ASSERT( fabs(dfIn->data[j] - valueRef) < EPS);
                }
            }
        }
    }

CPPUNIT_TEST_SUITE( class TestConstantRelaxation );
        CPPUNIT_TEST( testRelaxation);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestConstantRelaxation);
