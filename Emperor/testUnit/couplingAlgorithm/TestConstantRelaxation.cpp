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

namespace EMPIRE {

using namespace std;

class TestConstantRelaxation: public CppUnit::TestFixture {
private:
    DataField *dfIn;
    DataField *dfOut;
    ConstantRelaxation *constantRelaxation;
    double relaxationFactor;

public:
    void setUp() {
        relaxationFactor = 0.1;

        dfIn = new DataField("input", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        dfOut = new DataField("output", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        constantRelaxation = new ConstantRelaxation("", relaxationFactor);
        ConnectionIOSetup::setupIOForCouplingAlgorithm(constantRelaxation, NULL, dfIn, NULL, dfOut);
        constantRelaxation->debugMe = false;
    }
    void tearDown() {
        delete dfIn;
        delete dfOut;
        delete constantRelaxation;
    }
    /*
     * Test case:
     * Test check whether the step number is correct
     */
    void testNewTimeStep() {
        const int SIZE = dfIn->dimension * dfIn->numLocations;
        int NUM_INNER_LOOPS = 5;
        const double EPS = 1E-10;
        int NUM_OUTER_LOOPS = 2;

        for (int k = 0; k < NUM_OUTER_LOOPS; k++) {
            for (int i = 1; i <= NUM_INNER_LOOPS; i++) {
                for (int j = 0; j < SIZE; j++) {
                    dfIn->data[j] = 0.0;
                }
                if (i == 1)
                    constantRelaxation->setNewLoop();
                CPPUNIT_ASSERT(constantRelaxation->step == i);
                constantRelaxation->calcNewValue();
                constantRelaxation->initialResidual = 1.0; // dummy
            }
        }
    }
    /*
     * Test case:
     * Test check whether the result is as expected
     */
    void testRelaxation() {
        const int SIZE = dfIn->dimension * dfIn->numLocations;
        int NUM_INNER_LOOPS = 5;
        const double EPS = 1E-10;
        int NUM_OUTER_LOOPS = 2;
        // test case 1, all entries in an array are equal (e.g. 1,1,1, ...)
        for (int k = 0; k < NUM_OUTER_LOOPS; k++) {
            for (int i = 0; i < NUM_INNER_LOOPS; i++) {
                double valueIn = 0.0;
                if (i == 0)
                    valueIn = 10.0;
                else
                    valueIn = 20.0 + double(i) - 1;
                for (int j = 0; j < SIZE; j++) {
                    dfIn->data[j] = valueIn;
                }
                //debugOut << *dfIn;
                if (i == 0)
                    constantRelaxation->setNewLoop();
                constantRelaxation->calcNewValue();
                //debugOut << *dfOut;
                double valueRef = 0.0;
                if (i == 0)
                    valueRef = 10.0;
                else
                    valueRef = 10.0 + double(i);

                for (int j = 0; j < SIZE; j++) {
                    CPPUNIT_ASSERT( fabs(dfOut->data[j] - valueRef) < EPS);
                }
            }
        }
    }

CPPUNIT_TEST_SUITE( class TestConstantRelaxation );
        CPPUNIT_TEST( testNewTimeStep);
        CPPUNIT_TEST( testRelaxation);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestConstantRelaxation);
