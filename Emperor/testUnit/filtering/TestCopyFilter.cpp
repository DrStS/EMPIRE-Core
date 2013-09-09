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

#include "CopyFilter.h"
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
class TestCopyFilter: public CppUnit::TestFixture {
private:

public:
    void setUp() {
    }
    void tearDown() {
    }
    /********//**
     ***********************************************************************************************
     * \brief Test case: Test dummy filter when size of the input array is larger
     ***********/
    void testFilter1() {
        {
            CopyFilter *filter = new CopyFilter();
            DataField *in;
            DataField *out;
            in = new DataField("", EMPIRE_DataField_atNode, 2, EMPIRE_DataField_vector,
                    EMPIRE_DataField_field);
            out = new DataField("", EMPIRE_DataField_atNode, 1, EMPIRE_DataField_vector,
                    EMPIRE_DataField_field);
            ConnectionIOSetup::setupIOForFilter(filter, NULL, in, NULL, out);

            for (int i = 0; i < in->numLocations * EMPIRE_DataField_vector; i++)
                in->data[i] = 0.1 * i;

            filter->filtering();

            for (int i = 0; i < out->numLocations * EMPIRE_DataField_vector; i++) {
                CPPUNIT_ASSERT(out->data[i] == in->data[i]);
                CPPUNIT_ASSERT(out->data[i] == 0.1 * i);
            }
            delete filter;
            delete in;
            delete out;
        }
        {
            CopyFilter *filter = new CopyFilter();
            Signal *in;
            Signal *out;
            in = new Signal("", 1, 1, 10);
            out = new Signal("", 1, 1, 5);
            ConnectionIOSetup::setupIOForFilter(filter, in, out);

            for (int i = 0; i < in->size; i++)
                in->array[i] = 0.1 * i;

            filter->filtering();

            for (int i = 0; i < out->size; i++) {
                CPPUNIT_ASSERT(out->array[i] == in->array[i]);
                CPPUNIT_ASSERT(out->array[i] == 0.1 * i);
            }
            delete filter;
            delete in;
            delete out;
        }
    }
    /********//**
     ***********************************************************************************************
     * \brief Test case: Test dummy filter when size of the output array is larger
     ***********/
    void testFilter2() {
        {
            CopyFilter *filter = new CopyFilter();
            DataField *in;
            DataField *out;
            in = new DataField("", EMPIRE_DataField_atNode, 1, EMPIRE_DataField_scalar,
                    EMPIRE_DataField_field);
            out = new DataField("", EMPIRE_DataField_atNode, 2, EMPIRE_DataField_scalar,
                    EMPIRE_DataField_field);
            ConnectionIOSetup::setupIOForFilter(filter, NULL, in, NULL, out);

            for (int i = 0; i < in->numLocations * EMPIRE_DataField_scalar; i++)
                in->data[i] = 0.01 * i;

            filter->filtering();
            for (int i = 0; i < in->numLocations * EMPIRE_DataField_scalar; i++) {
                CPPUNIT_ASSERT(out->data[i] == in->data[i]);
                CPPUNIT_ASSERT(out->data[i] == 0.01 * i);
            }
            for (int i = in->numLocations * EMPIRE_DataField_scalar;
                    i < out->numLocations * EMPIRE_DataField_scalar; i++) {
                CPPUNIT_ASSERT(out->data[i] == 0.0);
            }
            delete filter;
            delete in;
            delete out;
        }
        {
            CopyFilter *filter = new CopyFilter();
            Signal *in;
            Signal *out;
            in = new Signal("", 1, 1, 5);
            out = new Signal("", 1, 1, 10);
            ConnectionIOSetup::setupIOForFilter(filter, in, out);

            for (int i = 0; i < in->size; i++)
                in->array[i] = 0.01 * i;

            filter->filtering();
            for (int i = 0; i < in->size; i++) {
                CPPUNIT_ASSERT(out->array[i] == in->array[i]);
                CPPUNIT_ASSERT(out->array[i] == 0.01 * i);
            }
            for (int i = in->size; i < out->size; i++) {
                CPPUNIT_ASSERT(out->array[i] == 0.0);
            }
            delete filter;
            delete in;
            delete out;
        }
    }
CPPUNIT_TEST_SUITE( TestCopyFilter );
        CPPUNIT_TEST( testFilter1);
        CPPUNIT_TEST( testFilter2);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestCopyFilter);
