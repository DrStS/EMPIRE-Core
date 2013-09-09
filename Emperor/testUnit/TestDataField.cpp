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
#include "DataField.h"
#include <string>

namespace EMPIRE {
/********//**
 * \brief Test the class DataField
 ***********/
class TestDataField: public CppUnit::TestFixture {
private:
    DataField *df;
public:
    void setUp() {
        std::string name("testData");
        int numLocations = 3;
        EMPIRE_DataField_dimension dimension = EMPIRE_DataField_vector;
        EMPIRE_DataField_location location = EMPIRE_DataField_atNode;
        EMPIRE_DataField_typeOfQuantity typeOfQuantity = EMPIRE_DataField_field;
        df = new DataField(name, location, numLocations, dimension, typeOfQuantity);
    }
    void tearDown() {
        delete df;
    }
    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/
    void testConstructor() {
        CPPUNIT_ASSERT(df->name == "testData");
        CPPUNIT_ASSERT(df->dimension == EMPIRE_DataField_vector);
        CPPUNIT_ASSERT(df->dimension == 3);
        CPPUNIT_ASSERT(df->typeOfQuantity == EMPIRE_DataField_field);
        CPPUNIT_ASSERT(df->data != NULL);
    }
    /***********************************************************************************************
     * \brief Test case: Test the data
     ***********/
    void testData() {
        for (int i = 0; i < df->numLocations; i++)
            for (int j = 0; j < df->dimension; j++)
                df->data[i * 3 + j] = i / 10.0;
        for (int i = 0; i < df->numLocations; i++)
            for (int j = 0; j < df->dimension; j++)
                CPPUNIT_ASSERT(df->data[i*3+j] == i / 10.0);
    }
CPPUNIT_TEST_SUITE( TestDataField );
        CPPUNIT_TEST( testConstructor);
        CPPUNIT_TEST( testData);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestDataField);
