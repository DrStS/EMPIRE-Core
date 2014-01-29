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

#include "AbstractExtrapolator.h"
#include "LinearExtrapolator.h"
#include "DataField.h"
#include "Signal.h"
#include "ConnectionIOSetup.h"
#include "Message.h"

#include <cmath>
using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test the class AbstractExtrapolator and its sub-classes
 ***********/
class TestExtrapolator: public CppUnit::TestFixture {
private:
public:
    void setUp() {
    }
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief Test case: Test whether linear extarpolator works or not
     ***********/
    void testLinearExtrapolator() {

        {
            /*
             * The value of df should be sequentially:
             * 1. set df=1
             * 2. set df=2
             * 3. extrapolate df=2*2-1=3
             * 4. extrapolate df=2*3-2=4
             * 5. extrapolate df=2*4-3=5
             */
            DataField *df = new DataField("dummy", EMPIRE_DataField_atNode, 10,
                    EMPIRE_DataField_vector, EMPIRE_DataField_field);

            LinearExtrapolator *extrapolator = new LinearExtrapolator("");
            extrapolator->addConnectionIO(ConnectionIOSetup::constructDummyConnectionIO(df));
            extrapolator->init();

            int size = df->dimension * df->numLocations;
            for (int i = 1; i <= 5; i++) {
                extrapolator->extrapolate();
                if (i == 1) {
                    for (int j = 0; j < size; j++)
                        df->data[j] = 1.0;
                } else if (i == 2) {
                    for (int j = 0; j < size; j++)
                        df->data[j] = 2.0;
                }

                for (int j = 0; j < size; j++)
                    CPPUNIT_ASSERT(df->data[j] == double(i));
            }
            delete df;
            delete extrapolator;
        }
        {
            /*
             * The value of df should be sequentially:
             * 1. set df=1
             * 2. set df=2
             * 3. extrapolate df=2*2-1=3, set df=6
             * 4. extrapolate df=2*6-2=10, set df=53
             * 5. extrapolate df=2*53-6=100
             */
            DataField *df = new DataField("dummy", EMPIRE_DataField_atNode, 10,
                    EMPIRE_DataField_vector, EMPIRE_DataField_field);

            LinearExtrapolator *extrapolator = new LinearExtrapolator("");
            extrapolator->addConnectionIO(ConnectionIOSetup::constructDummyConnectionIO(df));
            extrapolator->init();

            int size = df->dimension * df->numLocations;
            for (int i = 1; i <= 5; i++) {
                extrapolator->extrapolate();
                if (i == 3) {
                    for (int j = 0; j < size; j++)
                        CPPUNIT_ASSERT(df->data[j] == 3.0);
                    ;
                } else if (i == 4) {
                    for (int j = 0; j < size; j++)
                        CPPUNIT_ASSERT(df->data[j] == 10.0);
                    ;
                } else if (i == 5) {
                    for (int j = 0; j < size; j++)
                        CPPUNIT_ASSERT(df->data[j] == 100.0);
                    ;
                }

                if (i == 1) {
                    for (int j = 0; j < size; j++)
                        df->data[j] = 1.0;
                } else if (i == 2) {
                    for (int j = 0; j < size; j++)
                        df->data[j] = 2.0;
                } else if (i == 3) {
                    for (int j = 0; j < size; j++)
                        df->data[j] = 6.0;
                } else if (i == 4) {
                    for (int j = 0; j < size; j++)
                        df->data[j] = 53.0;
                }
            }
            delete df;
            delete extrapolator;
        }
    }

CPPUNIT_TEST_SUITE( TestExtrapolator );
        CPPUNIT_TEST( testLinearExtrapolator);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestExtrapolator);
