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
#include <iostream>

#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

#include "AdditionFilter.h"
#include "DataField.h"
#include "Signal.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"

using namespace std;

namespace EMPIRE {
/********//**
 ***********************************************************************************************
 * \brief Test the class CopyFilter
 ***********/
class TestAdditionFilter: public CppUnit::TestFixture {
private:

public:
    void setUp() {
    }
    void tearDown() {
    }
    /********//**
     ***********************************************************************************************
     * \brief Test case: Test y = a*x + b*y
     ***********/
    void testFilter() {
        { // on dataField
            double a = 2.0;
            double b = -1.0;
            AdditionFilter *filter = new AdditionFilter(a, b);
            DataField *x, *y;
            x = new DataField("", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                    EMPIRE_DataField_field);
            y = new DataField("", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                    EMPIRE_DataField_field);

            ConnectionIO *in0, *in1, *out0;

            in0 = ConnectionIOSetup::constructDummyConnectionIO(x);
            in1 = ConnectionIOSetup::constructDummyConnectionIO(y);
            out0 = ConnectionIOSetup::constructDummyConnectionIO(y);

            filter->addInput(in0);
            filter->addInput(in1);
            filter->addOutput(out0);
            filter->init();
            int n = x->numLocations * EMPIRE_DataField_vector;
            for (int i = 0; i < n; i++)
                x->data[i] = 1.0;
            for (int i = 0; i < n; i++)
                y->data[i] = -1.0;

            filter->filtering();

            for (int i = 0; i < n; i++) {
                CPPUNIT_ASSERT(y->data[i] == 3.0);
            }
            delete filter;
            delete x, y;
        }
        { // on signal
            double a = 2.0;
            double b = -1.0;
            AdditionFilter *filter = new AdditionFilter(a, b);
            Signal *x, *y;
            x = new Signal("", 1, 1, 5);
            y = new Signal("", 1, 1, 5);

            ConnectionIO *in0, *in1, *out0;

            in0 = ConnectionIOSetup::constructDummyConnectionIO(x);
            in1 = ConnectionIOSetup::constructDummyConnectionIO(y);
            out0 = ConnectionIOSetup::constructDummyConnectionIO(y);

            filter->addInput(in0);
            filter->addInput(in1);
            filter->addOutput(out0);
            filter->init();

            int n = x->size;
            for (int i = 0; i < n; i++)
                x->array[i] = 1.0;
            for (int i = 0; i < n; i++)
                y->array[i] = -1.0;

            filter->filtering();

            for (int i = 0; i < n; i++) {
                CPPUNIT_ASSERT(y->array[i] == 3.0);
            }
            delete filter;
            delete x, y;
        }
    }

CPPUNIT_TEST_SUITE( TestAdditionFilter );
        CPPUNIT_TEST( testFilter);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestAdditionFilter);
