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
#include <iomanip>
#include <math.h>

#include "Aitken.h"
#include "ConstantRelaxation.h"
#include "DataField.h"
#include "Signal.h"
#include "EMPEROR_Enum.h"
#include "ComaAitken.h"
#include "Message.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"

namespace EMPIRE {

using namespace std;

class TestRelaxationMethods: public CppUnit::TestFixture {
private:
    bool debugMe;
public:
    void setUp() {
        debugMe = false;
    }
    void tearDown() {
    }
    /*
     * solution of x=cos(x) is x=0.73908513321516
     */
    void testRelaxationMethods() {
        DataField *dfIn1 = new DataField("", EMPIRE_DataField_atNode, 1, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        DataField *dfOut1 = new DataField("", EMPIRE_DataField_atNode, 1, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        DataField *dfIn2 = new DataField("", EMPIRE_DataField_atNode, 1, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        DataField *dfOut2 = new DataField("", EMPIRE_DataField_atNode, 1, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        DataField *dfCurrent = new DataField("", EMPIRE_DataField_atNode, 1,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        DataField *dfOld = new DataField("", EMPIRE_DataField_atNode, 1, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);

        double &xIn1 = dfIn1->data[0];
        double &xOut1 = dfOut1->data[0];
        double &xIn2 = dfIn2->data[0];
        double &xOut2 = dfOut2->data[0];
        double &xCurrent = dfCurrent->data[0];
        double &xOld = dfOld->data[0];

        Aitken *aitken = new Aitken("", 1.0);
        aitken->debugMe = false;
        ConnectionIOSetup::setupIOForCouplingAlgorithm(aitken, NULL, dfIn1, NULL, dfOut1);
        ComaAitken *comaAitken = new ComaAitken(dfCurrent, dfOld, 1.0);
        ConstantRelaxation *cr = new ConstantRelaxation("", 0.8);
        ConnectionIOSetup::setupIOForCouplingAlgorithm(cr, NULL, dfIn2, NULL, dfOut2);

        int NUM_LOOPS = 100;
        double x = 0.0;
        xIn1 = x;
        xIn2 = x;
        xCurrent = x;
        if (debugMe) {
            cout << endl;
            cout << setprecision(14);
            cout << setw(20) << "n" << setw(20) << "fixed-point" << setw(20) << "empire_aitken"
                    << setw(20) << "coma_aitken" << setw(20) << "constant_relaxation" << endl;
        }
        for (int i = 1; i <= NUM_LOOPS; i++) {
            if (i == 1) {
                x = 0.0;
                xIn1 = x;
                xIn2 = x;
                xCurrent = x;
            } else {
                x = cos(x);
                xIn1 = cos(xOut1);
                xIn2 = cos(xOut2);
                xCurrent = cos(xOld);
            }
            // own aitken
            if (i == 1)
                aitken->setNewLoop();
            aitken->calcNewValue();
            // coma aitken
            if (i != 1)
                comaAitken->calcNewValues(i);
            xOld = xCurrent;
            if (i == 1)
                cr->setNewLoop();
            cr->calcNewValue();
            if (debugMe) {
                cout << setw(20) << i - 1 << setw(20) << x << setw(20) << xOut1 << setw(20) << xOld
                        << setw(20) << xOut2 << endl;
            }
        }

        const double EPS = 1E-10;
        CPPUNIT_ASSERT(fabs(xOut1-0.73908513321516) < EPS);
        CPPUNIT_ASSERT(fabs(xOut2-0.73908513321516) < EPS);
        CPPUNIT_ASSERT(fabs(xCurrent-0.73908513321516) < EPS);
        delete aitken;
        delete comaAitken;
        delete cr;
        delete dfCurrent, dfIn1, dfIn2, dfOld, dfOut1, dfOut2;
    }

CPPUNIT_TEST_SUITE( TestRelaxationMethods );
        CPPUNIT_TEST( testRelaxationMethods);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestRelaxationMethods);
