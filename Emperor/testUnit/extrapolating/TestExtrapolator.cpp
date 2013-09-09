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
#include "SimpleExtrapolator.h"
#include "GenMSExtrapolator.h"
#include "DataField.h"
#include "Signal.h"
#include "ConnectionIOSetup.h"

#include <cmath>
using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test the class AbstractExtrapolator and its sub-classes
 ***********/
class TestExtrapolator: public CppUnit::TestFixture {
private:
    AbstractExtrapolator *extrapolator;
public:
    void setUp() {
    }
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief Test case: Test whether setDoPredict() works
     ***********/
    void testSimpleExtrapolate() {
        DataField *df = new DataField("dummy", EMPIRE_DataField_atNode, 10, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        extrapolator = new SimpleExtrapolator("");
        ConnectionIOSetup::setupIOForExtrapolator(extrapolator, NULL, df, NULL, df);
        int size = df->dimension * df->numLocations;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < size; j++)
                df->data[j] = static_cast<double>(i);
            CPPUNIT_ASSERT(extrapolator->doExtrapolate==false);
            extrapolator->setDoExtrapolate();
            CPPUNIT_ASSERT(extrapolator->doExtrapolate==true);
            extrapolator->extrapolate();
            for (int j = 0; j < size; j++)
                CPPUNIT_ASSERT(df->data[j] == static_cast<double>(i));
        }
        delete df;
        clearExtrapolator(extrapolator);
    }

    void clearExtrapolator(AbstractExtrapolator *e) {
        for (int i = 0; i < e->inputVec.size(); i++) {
            delete e->inputVec[i];
        }
        for (int i = 0; i < e->outputVec.size(); i++) {
            delete e->outputVec[i];
        }
        delete e;
    }

    void fill(double* arr, double val, int n) {
        for (int i = 0; i < n; i++) {
            arr[i] = val;
        }
    }

    /***********************************************************************************************
     * \brief GenMSExtrapolator test case: Test DataField extrapolation with different input and output.
     * \author Michael Andre
     ***********/
    void testDataFieldGenMSExtrapolateDiffIO() {
        int numInput = 3;
        int seqLen = 6;
        bool sumOutput = true;
        double deltaTime = 0.1;
        int dataSize = 4;
        double eps = 1.e-20;
        vector<double> coefDot0(6, 0.);
        coefDot0[0] = 1.;
        vector<double> coefDot1(6, 0.);
        coefDot1[1] = 0.5;
        vector<double> coefDot2(6, 0.);
        coefDot2[2] = 3.;
        vector<double> coefOut(6, 0.);
        coefOut[2] = 1.;
        coefOut[3] = 2.;
        coefOut[4] = 4.;

        DataField *dfInDot0 = new DataField("dummyIn0", EMPIRE_DataField_atNode, 4,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        DataField *dfInDot1 = new DataField("dummyIn1", EMPIRE_DataField_atNode, 4,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        DataField *dfInDot2 = new DataField("dummyIn2", EMPIRE_DataField_atNode, 4,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        DataField *dfOut = new DataField("dummyOut", EMPIRE_DataField_atNode, 4,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);

        extrapolator = new GenMSExtrapolator("dummyExtrapolator", numInput, sumOutput, seqLen,
                deltaTime, &coefDot0, &coefDot1, &coefDot2, &coefOut);

        ConnectionIOSetup::addInputForExtrapolator(extrapolator, NULL, dfInDot0);
        ConnectionIOSetup::addInputForExtrapolator(extrapolator, NULL, dfInDot1);
        ConnectionIOSetup::addInputForExtrapolator(extrapolator, NULL, dfInDot2);
        ConnectionIOSetup::addOutputForExtrapolator(extrapolator, NULL, dfOut);

        // step 1
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInDot0->data, 1., dataSize);
        fill(dfInDot1->data, (2. / deltaTime), dataSize);
        fill(dfInDot2->data, (2. / (3 * deltaTime * deltaTime)), dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfOut->data[i] - 1.) < eps);
        }
        // step 2
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInDot0->data, 1., dataSize);
        fill(dfInDot1->data, (2. / deltaTime), dataSize);
        fill(dfInDot2->data, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfOut->data[i] - 2.) < eps);
        }
        // step 3
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInDot0->data, 1., dataSize);
        fill(dfInDot1->data, 0., dataSize);
        fill(dfInDot2->data, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfOut->data[i] - 4.) < eps);
        }
        // step 4
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInDot0->data, 0., dataSize);
        fill(dfInDot1->data, 0., dataSize);
        fill(dfInDot2->data, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfOut->data[i] - 1.) < eps);
        }
        // step 5
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInDot0->data, 0., dataSize);
        fill(dfInDot1->data, 0., dataSize);
        fill(dfInDot2->data, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfOut->data[i] - 4.) < eps);
        }
        // step 6
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInDot0->data, 0., dataSize);
        fill(dfInDot1->data, 0., dataSize);
        fill(dfInDot2->data, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfOut->data[i] - 12.) < eps);
        }

        delete dfInDot0;
        delete dfInDot1;
        delete dfInDot2;
        delete dfOut;
        clearExtrapolator(extrapolator);
    }

    /***********************************************************************************************
     * \brief GenMSExtrapolator test case: Test Signal extrapolation with different input and output.
     * \author Michael Andre
     ***********/
    void testSignalGenMSExtrapolateDiffIO() {
        int numInput = 3;
        int seqLen = 6;
        bool sumOutput = true;
        double deltaTime = 0.1;
        int dataSize = 4;
        double eps = 1.e-20;
        vector<double> coefDot0(6, 0.);
        coefDot0[0] = 1.;
        vector<double> coefDot1(6, 0.);
        coefDot1[1] = 0.5;
        vector<double> coefDot2(6, 0.);
        coefDot2[2] = 3.;
        vector<double> coefOut(6, 0.);
        coefOut[2] = 1.;
        coefOut[3] = 2.;
        coefOut[4] = 4.;

        Signal *sigInDot0 = new Signal("dummyIn0", dataSize, 1, 1);
        Signal *sigInDot1 = new Signal("dummyIn1", dataSize, 1, 1);
        Signal *sigInDot2 = new Signal("dummyIn2", dataSize, 1, 1);
        Signal *sigOut = new Signal("dummyOut", dataSize, 1, 1);

        extrapolator = new GenMSExtrapolator("dummyExtrapolator", numInput, sumOutput, seqLen,
                deltaTime, &coefDot0, &coefDot1, &coefDot2, &coefOut);

        ConnectionIOSetup::addInputForExtrapolator(extrapolator, sigInDot0);
        ConnectionIOSetup::addInputForExtrapolator(extrapolator, sigInDot1);
        ConnectionIOSetup::addInputForExtrapolator(extrapolator, sigInDot2);
        ConnectionIOSetup::addOutputForExtrapolator(extrapolator, sigOut);

        // step 1
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInDot0->array, 1., dataSize);
        fill(sigInDot1->array, (2. / deltaTime), dataSize);
        fill(sigInDot2->array, (2. / (3 * deltaTime * deltaTime)), dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigOut->array[i] - 1.) < eps);
        }
        // step 2
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInDot0->array, 1., dataSize);
        fill(sigInDot1->array, (2. / deltaTime), dataSize);
        fill(sigInDot2->array, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigOut->array[i] - 2.) < eps);
        }
        // step 3
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInDot0->array, 1., dataSize);
        fill(sigInDot1->array, 0., dataSize);
        fill(sigInDot2->array, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigOut->array[i] - 4.) < eps);
        }
        // step 4
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInDot0->array, 0., dataSize);
        fill(sigInDot1->array, 0., dataSize);
        fill(sigInDot2->array, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigOut->array[i] - 1.) < eps);
        }
        // step 5
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInDot0->array, 0., dataSize);
        fill(sigInDot1->array, 0., dataSize);
        fill(sigInDot2->array, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigOut->array[i] - 4.) < eps);
        }
        // step 6
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInDot0->array, 0., dataSize);
        fill(sigInDot1->array, 0., dataSize);
        fill(sigInDot2->array, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigOut->array[i] - 12.) < eps);
        }

        delete sigInDot0;
        delete sigInDot1;
        delete sigInDot2;
        delete sigOut;
        clearExtrapolator(extrapolator);
    }

    /***********************************************************************************************
     * \brief GenMSExtrapolator test case: Test DataField extrapolation with same input and output.
     * \author Michael Andre
     ***********/
    void testDataFieldGenMSExtrapolateSameIO() {
        int numInput = 3;
        int seqLen = 6;
        bool sumOutput = true;
        double deltaTime = 0.1;
        int dataSize = 4;
        double eps = 1.e-20;
        vector<double> coefDot0(6, 0.);
        coefDot0[0] = 1.;
        vector<double> coefDot1(6, 0.);
        coefDot1[1] = 0.5;
        vector<double> coefDot2(6, 0.);
        coefDot2[2] = 3.;
        vector<double> coefOut(6, 0.);
        coefOut[2] = 1.;
        coefOut[3] = 2.;
        coefOut[4] = 4.;

        DataField *dfInOutDot0 = new DataField("dummyInOut0", EMPIRE_DataField_atNode, 4,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        DataField *dfInDot1 = new DataField("dummyIn1", EMPIRE_DataField_atNode, 4,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        DataField *dfInDot2 = new DataField("dummyIn2", EMPIRE_DataField_atNode, 4,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);

        extrapolator = new GenMSExtrapolator("dummyExtrapolator", numInput, sumOutput, seqLen,
                deltaTime, &coefDot0, &coefDot1, &coefDot2, &coefOut);

        ConnectionIOSetup::addInputForExtrapolator(extrapolator, NULL, dfInOutDot0);
        ConnectionIOSetup::addInputForExtrapolator(extrapolator, NULL, dfInDot1);
        ConnectionIOSetup::addInputForExtrapolator(extrapolator, NULL, dfInDot2);
        ConnectionIOSetup::addOutputForExtrapolator(extrapolator, NULL, dfInOutDot0);

        // step 1
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInOutDot0->data, 1., dataSize);
        fill(dfInDot1->data, (2. / deltaTime), dataSize);
        fill(dfInDot2->data, (2. / (3 * deltaTime * deltaTime)), dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfInOutDot0->data[i] - 1.) < eps);
        }
        // step 2
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInOutDot0->data, 1., dataSize);
        fill(dfInDot1->data, (2. / deltaTime), dataSize);
        fill(dfInDot2->data, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfInOutDot0->data[i] - 2.) < eps);
        }
        // step 3
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInOutDot0->data, 1., dataSize);
        fill(dfInDot1->data, 0., dataSize);
        fill(dfInDot2->data, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfInOutDot0->data[i] - 4.) < eps);
        }
        // step 4
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInOutDot0->data, 0., dataSize);
        fill(dfInDot1->data, 0., dataSize);
        fill(dfInDot2->data, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfInOutDot0->data[i] - 1.) < eps);
        }
        // step 5
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInOutDot0->data, 0., dataSize);
        fill(dfInDot1->data, 0., dataSize);
        fill(dfInDot2->data, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfInOutDot0->data[i] - 4.) < eps);
        }
        // step 6
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(dfInOutDot0->data, 0., dataSize);
        fill(dfInDot1->data, 0., dataSize);
        fill(dfInDot2->data, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(dfInOutDot0->data[i] - 12.) < eps);
        }

        delete dfInOutDot0;
        delete dfInDot1;
        delete dfInDot2;
        clearExtrapolator(extrapolator);
    }

    /***********************************************************************************************
     * \brief GenMSExtrapolator test case: Test Signal extrapolation with same input and output.
     * \author Michael Andre
     ***********/
    void testSignalGenMSExtrapolateSameIO() {
        int numInput = 3;
        int seqLen = 6;
        bool sumOutput = true;
        double deltaTime = 0.1;
        int dataSize = 4;
        double eps = 1.e-20;
        vector<double> coefDot0(6, 0.);
        coefDot0[0] = 1.;
        vector<double> coefDot1(6, 0.);
        coefDot1[1] = 0.5;
        vector<double> coefDot2(6, 0.);
        coefDot2[2] = 3.;
        vector<double> coefOut(6, 0.);
        coefOut[2] = 1.;
        coefOut[3] = 2.;
        coefOut[4] = 4.;

        Signal *sigInOutDot0 = new Signal("dummyInOut0", dataSize, 1, 1);
        Signal *sigInDot1 = new Signal("dummyIn1", dataSize, 1, 1);
        Signal *sigInDot2 = new Signal("dummyIn2", dataSize, 1, 1);

        extrapolator = new GenMSExtrapolator("dummyExtrapolator", numInput, sumOutput, seqLen,
                deltaTime, &coefDot0, &coefDot1, &coefDot2, &coefOut);

        ConnectionIOSetup::addInputForExtrapolator(extrapolator, sigInOutDot0);
        ConnectionIOSetup::addInputForExtrapolator(extrapolator, sigInDot1);
        ConnectionIOSetup::addInputForExtrapolator(extrapolator, sigInDot2);
        ConnectionIOSetup::addOutputForExtrapolator(extrapolator, sigInOutDot0);

        // step 1
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInOutDot0->array, 1., dataSize);
        fill(sigInDot1->array, (2. / deltaTime), dataSize);
        fill(sigInDot2->array, (2. / (3 * deltaTime * deltaTime)), dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigInOutDot0->array[i] - 1.) < eps);
        }
        // step 2
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInOutDot0->array, 1., dataSize);
        fill(sigInDot1->array, (2. / deltaTime), dataSize);
        fill(sigInDot2->array, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigInOutDot0->array[i] - 2.) < eps);
        }
        // step 3
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInOutDot0->array, 1., dataSize);
        fill(sigInDot1->array, 0., dataSize);
        fill(sigInDot2->array, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigInOutDot0->array[i] - 4.) < eps);
        }
        // step 4
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInOutDot0->array, 0., dataSize);
        fill(sigInDot1->array, 0., dataSize);
        fill(sigInDot2->array, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigInOutDot0->array[i] - 1.) < eps);
        }
        // step 5
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInOutDot0->array, 0., dataSize);
        fill(sigInDot1->array, 0., dataSize);
        fill(sigInDot2->array, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigInOutDot0->array[i] - 4.) < eps);
        }
        // step 6
        CPPUNIT_ASSERT(extrapolator->doExtrapolate == false);
        extrapolator->setDoExtrapolate();
        fill(sigInOutDot0->array, 0., dataSize);
        fill(sigInDot1->array, 0., dataSize);
        fill(sigInDot2->array, 0., dataSize);
        extrapolator->extrapolate();
        for (int i = 0; i < dataSize; i++) {
            CPPUNIT_ASSERT(fabs(sigInOutDot0->array[i] - 12.) < eps);
        }

        delete sigInOutDot0;
        delete sigInDot1;
        delete sigInDot2;
        clearExtrapolator(extrapolator);
    }

CPPUNIT_TEST_SUITE( TestExtrapolator );
        CPPUNIT_TEST( testSimpleExtrapolate);
        CPPUNIT_TEST( testDataFieldGenMSExtrapolateDiffIO);
        CPPUNIT_TEST( testSignalGenMSExtrapolateDiffIO);
        CPPUNIT_TEST( testDataFieldGenMSExtrapolateSameIO);
        CPPUNIT_TEST( testSignalGenMSExtrapolateSameIO);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestExtrapolator);
