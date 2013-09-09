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
#include "Signal.h"
#include "EMPEROR_Enum.h"
#include <string>

namespace EMPIRE {
/********//**
 * \brief Test the class Signal
 ***********/
class TestSignal: public CppUnit::TestFixture {
private:
public:
    void setUp() {
    }
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/
    void testConstructor() {
        {
            /*
             * |1,2| |5,6|
             * |3,4| |7,8|
             */
            Signal *signal = new Signal("array", 2, 2, 2);
            CPPUNIT_ASSERT(signal->name == "array");
            CPPUNIT_ASSERT(signal->dimension == EMPIRE_Signal_3D);
            CPPUNIT_ASSERT(signal->size == 8);
            CPPUNIT_ASSERT(signal->size3D[0] == 2);
            CPPUNIT_ASSERT(signal->size3D[1] == 2);
            CPPUNIT_ASSERT(signal->size3D[2] == 2);
            for (int i = 0; i < signal->size; i++)
                signal->array[i] = (double) (i + 1);
            CPPUNIT_ASSERT(signal->entry(0,0,0) == 1.0);
            CPPUNIT_ASSERT(signal->entry(0,0,1) == 2.0);
            CPPUNIT_ASSERT(signal->entry(0,1,0) == 3.0);
            CPPUNIT_ASSERT(signal->entry(0,1,1) == 4.0);
            CPPUNIT_ASSERT(signal->entry(1,0,0) == 5.0);
            CPPUNIT_ASSERT(signal->entry(1,0,1) == 6.0);
            CPPUNIT_ASSERT(signal->entry(1,1,0) == 7.0);
            CPPUNIT_ASSERT(signal->entry(1,1,1) == 8.0);
            delete signal;
        }
        {
            /*
             * |1,2|
             * |3,4|
             */
            Signal *signal = new Signal("array", 1, 2, 2);
            CPPUNIT_ASSERT(signal->name == "array");
            CPPUNIT_ASSERT(signal->dimension == EMPIRE_Signal_2D);
            CPPUNIT_ASSERT(signal->size == 4);
            CPPUNIT_ASSERT(signal->size3D[0] == 1);
            CPPUNIT_ASSERT(signal->size3D[1] == 2);
            CPPUNIT_ASSERT(signal->size3D[2] == 2);
            for (int i = 0; i < signal->size; i++)
                signal->array[i] = (double) (i + 1);
            CPPUNIT_ASSERT(signal->entry(0,0,0) == 1.0);
            CPPUNIT_ASSERT(signal->entry(0,0,1) == 2.0);
            CPPUNIT_ASSERT(signal->entry(0,1,0) == 3.0);
            CPPUNIT_ASSERT(signal->entry(0,1,1) == 4.0);
            delete signal;
        }
        {
            /*
             * |1,2|
             */
            Signal *signal = new Signal("array", 1, 1, 2);
            CPPUNIT_ASSERT(signal->name == "array");
            CPPUNIT_ASSERT(signal->dimension == EMPIRE_Signal_1D);
            CPPUNIT_ASSERT(signal->size == 2);
            CPPUNIT_ASSERT(signal->size3D[0] == 1);
            CPPUNIT_ASSERT(signal->size3D[1] == 1);
            CPPUNIT_ASSERT(signal->size3D[2] == 2);
            for (int i = 0; i < signal->size; i++)
                signal->array[i] = (double) (i + 1);
            CPPUNIT_ASSERT(signal->entry(0,0,0) == 1.0);
            CPPUNIT_ASSERT(signal->entry(0,0,1) == 2.0);
            delete signal;
        }
        {
            /*
             * |1|
             */
            Signal *signal = new Signal("array", 1, 1, 1);
            CPPUNIT_ASSERT(signal->name == "array");
            CPPUNIT_ASSERT(signal->dimension == EMPIRE_Signal_0D);
            CPPUNIT_ASSERT(signal->size == 1);
            CPPUNIT_ASSERT(signal->size3D[0] == 1);
            CPPUNIT_ASSERT(signal->size3D[1] == 1);
            CPPUNIT_ASSERT(signal->size3D[2] == 1);
            for (int i = 0; i < signal->size; i++)
                signal->array[i] = (double) (i + 1);
            CPPUNIT_ASSERT(signal->entry(0,0,0) == 1.0);
            delete signal;
        }
    }
CPPUNIT_TEST_SUITE( TestSignal );
        CPPUNIT_TEST( testConstructor);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestSignal);
