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
#include "ClientCode.h"
#include "DataField.h"
#include "Connection.h"
#include "CopyFilter.h"
#include "Signal.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"

namespace EMPIRE {
/********//**
 * \brief Test the class Connection
 ***********/
class TestConnection: public CppUnit::TestFixture {
private:
    Connection *connection;
    AbstractFilter *filter1;
    AbstractFilter *filter2;
    AbstractFilter *filter3;

public:
    /***********************************************************************************************
     * \brief Set up dummy objects of ClientCode, Connection, DataField and CopyFilter.
     ***********/
    void setUp() {
        connection = new Connection("connection");
        filter1 = new CopyFilter();
        filter2 = new CopyFilter();
        filter3 = new CopyFilter();
    }
    void tearDown() {
        delete connection;
        //delete filter1, filter2, filter3; // filters are already deleted by connection
    }
    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/
    void testInitialization() {
        CPPUNIT_ASSERT(connection->name=="connection");
    }
    /***********************************************************************************************
     * \brief Test case: Test the filter sequence
     ***********/
    void testFilterSequence() {
        CPPUNIT_ASSERT(connection->filterVec.size() == 0);
        connection->addFilter(filter1);
        connection->addFilter(filter2);
        connection->addFilter(filter3);
        CPPUNIT_ASSERT(connection->filterVec.size() == 3);
        CPPUNIT_ASSERT(connection->filterVec.at(0) == filter1);
        CPPUNIT_ASSERT(connection->filterVec.at(1) == filter2);
        CPPUNIT_ASSERT(connection->filterVec.at(2) == filter3);
    }

CPPUNIT_TEST_SUITE( TestConnection );
        CPPUNIT_TEST( testInitialization);
        CPPUNIT_TEST( testFilterSequence);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestConnection);
