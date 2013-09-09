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

#include "Message.h"
#include <iostream>
#include <sstream>
#include <string>

namespace EMPIRE {
using namespace std;

class Integer {
public:
    Integer(int _i) :
            i(_i) {

    }
    virtual ~Integer() {
    }
    friend Message &operator<<(Message &m, Integer &integer);
private:
    int i;
};

Message &operator<<(Message &m, Integer &integer) {
    return m << "Test integer text "<< integer.i << " " ;
}


/********//**
 * \brief Test the class Message
 ***********/
class TestMessage: public CppUnit::TestFixture {
private:
public:
    void setUp() {
    }
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief Test case 1
     ***********/
    void testMessage() {
        Integer interger(10);

        //Redirect buffer to stringstream buffer
        stringstream buffer;
        streambuf * old = infoOut.rdbuf(buffer.rdbuf());

        infoOut() << "TestString1" << endl;
        string testString = buffer.str();

        //Restore old buffer
        infoOut.rdbuf(old);

        //cout << endl;
        //cout << testString << endl;
        CPPUNIT_ASSERT((testString.compare("EMPIRE_INFO: TestString1\n") == 0 ));

        //Redirect buffer to stringstream buffer
        stringstream buffer2;
        old = infoOut.rdbuf(buffer2.rdbuf());

        infoOut() << interger << "TestString2" << endl;
        string testString2 = buffer2.str();

        //Restore old buffer
        infoOut.rdbuf(old);

        //cout << testString2 << endl;
        CPPUNIT_ASSERT((testString2.compare("EMPIRE_INFO: Test integer text 10 TestString2\n") == 0 ) );
    }
    CPPUNIT_TEST_SUITE( TestMessage );
    CPPUNIT_TEST(testMessage);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestMessage);
