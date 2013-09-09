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
#include "Connection.h"
#include "IterativeCouplingLoop.h"
#include "TimeStepLoop.h"
#include "DataField.h"
#include "ConvergenceChecker.h"
#include "EMPEROR_Enum.h"
#include "ConnectionIO.h"

namespace EMPIRE {
/********//**
 * \brief Test the class IterativeCouplingLoop and TimeStepLoop. The coupling procedure can only
 *        be fully tested with outside clients. So we can only test whether the members are set up
 *        correctly.
 ***********/
class TestLoops: public CppUnit::TestFixture {
private:
    Connection *connection1;
    Connection *connection2;
    Connection *connection3;
    IterativeCouplingLoop *iterativeCouplingLoop;
    TimeStepLoop *timeStepLoop;

public:
    void setUp() {

        connection1 = new Connection("connection1");
        connection2 = new Connection("connection2");
        connection3 = new Connection("connection3");
        iterativeCouplingLoop = new IterativeCouplingLoop();
        timeStepLoop = new TimeStepLoop(10);
    }
    void tearDown() {
        delete connection1, connection2, connection3;
        delete iterativeCouplingLoop, timeStepLoop;
    }
    /***********************************************************************************************
     * \brief Test case: Test adding child coupling logics
     ***********/
    void testCouplingLogicSequence() {
        iterativeCouplingLoop->addCouplingLogic(connection1);
        iterativeCouplingLoop->addCouplingLogic(connection2);
        timeStepLoop->addCouplingLogic(iterativeCouplingLoop);
        timeStepLoop->addCouplingLogic(connection3);

        CPPUNIT_ASSERT(iterativeCouplingLoop->couplingLogicSequence.size() == 2);
        CPPUNIT_ASSERT(iterativeCouplingLoop->couplingLogicSequence[0] == connection1);
        CPPUNIT_ASSERT(iterativeCouplingLoop->couplingLogicSequence[1] == connection2);

        CPPUNIT_ASSERT(timeStepLoop->couplingLogicSequence.size() == 2);
        CPPUNIT_ASSERT(timeStepLoop->couplingLogicSequence[0] == iterativeCouplingLoop);
        CPPUNIT_ASSERT(timeStepLoop->couplingLogicSequence[1] == connection3);
    }

CPPUNIT_TEST_SUITE( TestLoops );
        CPPUNIT_TEST( testCouplingLogicSequence);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestLoops);
