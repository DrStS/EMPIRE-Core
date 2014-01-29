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

#include "DataField.h"
#include "Signal.h"
#include "EMPEROR_Enum.h"
#include "Message.h"
#include "Residual.h"
#include "ConnectionIOSetup.h"

namespace EMPIRE {

using namespace std;

class TestResidual: public CppUnit::TestFixture {
private:
    DataField *df01;
    DataField *df02;
    DataField *df11;
    ConnectionIO *io01;
    ConnectionIO *io02;
    ConnectionIO *io11;

    double coeff01;
    double coeff02;
    double coeff11;
    std::string timeToUpdate01;
    std::string timeToUpdate02;
    std::string timeToUpdate11;
    Residual *dataFieldResidual;

public:
    void setUp() {
        df01 = new DataField("", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        df02 = new DataField("", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        df11 = new DataField("", EMPIRE_DataField_atNode, 5, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        for (int i = 0; i < 5 * EMPIRE_DataField_vector; i++) {
            df01->data[i] = 1.0;
            df02->data[i] = 2.0;
            df11->data[i] = 3.0;
        }
        coeff01 = 1.0;
        coeff02 = 1.0;
        coeff11 = -1.0;
        timeToUpdate01 = "iterationBeginning";
        timeToUpdate02 = "iterationBeginning";
        timeToUpdate11 = "iterationEnd";
        io01 = ConnectionIOSetup::constructDummyConnectionIO(df01);
        io02 = ConnectionIOSetup::constructDummyConnectionIO(df02);
        io11 = ConnectionIOSetup::constructDummyConnectionIO(df11);
        dataFieldResidual = new Residual(1);
        dataFieldResidual->addComponent(coeff01, timeToUpdate01, io01);
        dataFieldResidual->addComponent(coeff02, timeToUpdate02, io02);
        dataFieldResidual->addComponent(coeff11, timeToUpdate11, io11);
        dataFieldResidual->init();
        dataFieldResidual->updateAtIterationBeginning();
        dataFieldResidual->updateAtIterationEnd();
    }
    void tearDown() {
        delete df01;
        delete df02;
        delete df11;
        delete io01;
        delete io02;
        delete io11;
    }
    /*
     * Test case:
     * Test the residual computation (1+2-3==0)
     */
    void testResidualComputation() {
        CPPUNIT_ASSERT(dataFieldResidual->components.size() == 3);
        CPPUNIT_ASSERT(dataFieldResidual->index == 1);
        CPPUNIT_ASSERT(dataFieldResidual->size == 5 * EMPIRE_DataField_vector);

        dataFieldResidual->computeCurrentResidual();
        for (int i = 0; i < 5 * EMPIRE_DataField_vector; i++) {
            CPPUNIT_ASSERT(fabs(dataFieldResidual->residualVector[i] - 0.0) < 1e-10);
        }
    }
    /*
     * Test case:
     * Test update
     */
    void testUpdate() {
        const double EPS = 1e-10;
        // 3+3-3 == 3
        for (int i = 0; i < 5 * EMPIRE_DataField_vector; i++) {
            df01->data[i] = 3.0;
            df02->data[i] = 3.0;
            df11->data[i] = 0.0; // should not affect the result
        }
        dataFieldResidual->updateAtIterationBeginning();
        dataFieldResidual->computeCurrentResidual();

        for (int i = 0; i < 5 * EMPIRE_DataField_vector; i++) {
            CPPUNIT_ASSERT(fabs(dataFieldResidual->residualVector[i] - 3.0) < EPS);
        }
        CPPUNIT_ASSERT(fabs(dataFieldResidual->residualVectorL2Norm - 3.0) < EPS);
        // 3+3-4 == 2
        for (int i = 0; i < 5 * EMPIRE_DataField_vector; i++) {
            df01->data[i] = 0.0; // should not affect the result
            df02->data[i] = 0.0; // should not affect the result
            df11->data[i] = 4.0;
        }
        dataFieldResidual->updateAtIterationEnd();
        dataFieldResidual->computeCurrentResidual();
        for (int i = 0; i < 5 * EMPIRE_DataField_vector; i++) {
            CPPUNIT_ASSERT(fabs(dataFieldResidual->residualVector[i] - 2.0) < EPS);
        }
        CPPUNIT_ASSERT(fabs(dataFieldResidual->residualVectorL2Norm - 2.0) < EPS);
    }

CPPUNIT_TEST_SUITE( class TestResidual );
        CPPUNIT_TEST( testResidualComputation);
        CPPUNIT_TEST( testUpdate);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestResidual);
