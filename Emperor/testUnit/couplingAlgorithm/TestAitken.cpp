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

#include "Aitken.h"
#include "DataField.h"
#include "Signal.h"
#include "EMPEROR_Enum.h"
#include "ComaAitken.h"
#include "Message.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"

namespace EMPIRE {

using namespace std;

class TestAitken: public CppUnit::TestFixture {
private:
    DataField *dfIn;
    DataField *dfOut;
    Aitken *aitken;
    double aitkenFactor;
    double VALUE;

    DataField *dfCurrent;
    DataField *dfOld;
    ComaAitken *comaAitken;

    static const double EPS = 1E-10;

public:
    void setUp() {
        /*aitkenFactor = 0.5;

        dfIn = new DataField("input", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        dfOut = new DataField("output", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        aitken = new Aitken("aitken", aitkenFactor);
        aitken->debugMe = false;
        ConnectionIOSetup::setupIOForCouplingAlgorithm(aitken, NULL, dfIn, NULL,dfOut);

        dfCurrent = new DataField("current", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        dfOld = new DataField("old", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        comaAitken = new ComaAitken(dfCurrent, dfOld, aitkenFactor);

        VALUE = 100.0;*/
    }
    void tearDown() {
        /*delete dfIn;
        delete dfOut;
        delete aitken;
        delete dfCurrent;
        delete dfOld;*/
    }
    /*
     * Test case:
     * Test Aitken in a nested loop, check whether the step number is correct
     */
    void testNewAitken() {
        /*const int SIZE = dfIn->dimension * dfIn->numLocations;
        int NUM_INNER_LOOPS = 5;
        const double EPS = 1E-10;
        int NUM_OUTER_LOOPS = 2;

        for (int k = 0; k < NUM_OUTER_LOOPS; k++) {
            for (int i = 1; i <= NUM_INNER_LOOPS; i++) {
                for (int j = 0; j < SIZE; j++) {
                    dfIn->data[j] = VALUE + 1.0 / (double) i;
                }
                if (i == 1)
                    aitken->setNewLoop();
                CPPUNIT_ASSERT(aitken->step == i);
                aitken->calcNewValue();
                //debugOut << *dfIn;
                //debugOut << *dfOut;
            }
        }*/
    }
    /*
     * Test case:
     * Test Aitken in a nested loop, check whether the result from own Aitken and CoMA's
     * Aitken is the same or not
     */
    void compareWithComaAitken() {
        /*const int SIZE = dfIn->dimension * dfIn->numLocations;
        int NUM_INNER_LOOPS = 5;
        const double EPS = 1E-10;
        int NUM_OUTER_LOOPS = 2;
        // test case 1, all entries in an array are equal (e.g. 1,1,1, ...)
        for (int k = 0; k < NUM_OUTER_LOOPS; k++) {
            for (int i = 1; i <= NUM_INNER_LOOPS; i++) {
                for (int j = 0; j < SIZE; j++) {
                    dfIn->data[j] = VALUE + 1.0 / (double) i;
                    dfCurrent->data[j] = dfIn->data[j];
                }
                //debugOut << *dfCurrent;
                // own aitken
                if (i == 1)
                    aitken->setNewLoop();
                aitken->calcNewValue();
                // coma aitken
                if (i != 1)
                    comaAitken->calcNewValues(i);
                for (int j = 0; j < SIZE; j++) {
                    dfOld->data[j] = dfCurrent->data[j];
                }
                // assert the results from two aitken implementations are equal
                //debugOut << *dfCurrent;
                for (int j = 0; j < SIZE; j++) {
                    CPPUNIT_ASSERT( fabs(dfOut->data[j] - dfCurrent->data[j]) < EPS);
                }
            }
        }

        // test case 2, all entries in an array are unequal (e.g. 1,2,3, ...)
        for (int k = 0; k < NUM_OUTER_LOOPS; k++) {
            for (int i = 1; i <= NUM_INNER_LOOPS; i++) {
                for (int j = 0; j < SIZE; j++) {
                    dfIn->data[j] = VALUE * j + 1.0 / (double) i;
                    dfCurrent->data[j] = dfIn->data[j];
                }
                // own aitken
                if (i == 1)
                    aitken->setNewLoop();
                aitken->calcNewValue();
                // coma aitken
                if (i != 1)
                    comaAitken->calcNewValues(i);
                for (int j = 0; j < SIZE; j++) {
                    dfOld->data[j] = dfCurrent->data[j];
                }
                // assert the results from two aitken implementations are equal
                //debugOut << *dfCurrent;
                //debugOut << *dfOut;
                for (int j = 0; j < SIZE; j++) {
                    CPPUNIT_ASSERT( fabs(dfOut->data[j] - dfCurrent->data[j]) < EPS);
                }
            }
        }*/
    }

CPPUNIT_TEST_SUITE( TestAitken );
        CPPUNIT_TEST( testNewAitken);
        CPPUNIT_TEST( compareWithComaAitken);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestAitken);
