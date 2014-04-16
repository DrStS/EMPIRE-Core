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
#include <iostream>
#include <typeinfo>
#include <map>

#include "ClientCode.h"
#include "DataField.h"
#include "FEMesh.h"
#include "MetaDatabase.h"
#include "Emperor.h"
#include "Aitken.h"
#include "ConstantRelaxation.h"
#include "Residual.h"
#include "Connection.h"
#include "CouplingLogicSequence.h"
#include "TimeStepLoop.h"
#include "OptimizationLoop.h"
#include "IterativeCouplingLoop.h"
#include "ConvergenceChecker.h"
#include "DataOutput.h"
#include "ConnectionIO.h"
#include "LinearExtrapolator.h"

namespace EMPIRE {

using namespace std;
/********//**
 * \brief Test the class Emperor
 ***********/
class TestEmperor: public CppUnit::TestFixture {
private:
public:
    void setUp() {

    }
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief Test case: dummy objects of ClientCodes, Connections, Meshes and so on are created.
     *        We write objects of structCouplingLogic to set up the whole global coupling logic (time step
     *        loop, iterative coupling loop, etc.). Finally, we test whether coupling logic settings are
     *        parsed correctly in function initGlobalCouplingLogic().
     ***********/
    void testInitCouplingLogic() {
        const int INTY = 5;
        const double DOUBY = 0.5;
        const string STRING_DUMMY = "dummy";

        // 0. initialize Emperor instance
        Emperor *emperor = new Emperor();
        MetaDatabase::metaDatabase = new MetaDatabase();

        // 1. set up client codes
        const string STRING_MESH = "dummyMesh";
        const string STRING_CLIENT_A = "clientA";
        const string STRING_CLIENT_B = "clientB";
        ClientCode *clientA = new ClientCode(STRING_CLIENT_A);
        ClientCode *clientB = new ClientCode(STRING_CLIENT_B);
        const string STRING_DATAFIELD = "dummyDatafield";

        // dummyMesh and dummyDataField are used only to pass NULL pointer check
        FEMesh *dummyMeshA = new FEMesh(STRING_MESH, INTY, INTY);
        dummyMeshA->addDataField(STRING_DATAFIELD, EMPIRE_DataField_atNode, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        DataField *dummyDataFieldA = dummyMeshA->getDataFieldByName(STRING_DATAFIELD);
        clientA->nameToMeshMap.insert(pair<string, FEMesh*>(dummyMeshA->name, dummyMeshA));
        FEMesh *dummyMeshB = new FEMesh(STRING_MESH, INTY, INTY);
        dummyMeshB->addDataField(STRING_DATAFIELD, EMPIRE_DataField_atNode, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        DataField *dummyDataFieldB = dummyMeshB->getDataFieldByName(STRING_DATAFIELD);
        clientB->nameToMeshMap.insert(pair<string, FEMesh*>(dummyMeshB->name, dummyMeshB));
        emperor->nameToClientCodeMap.insert(pair<string, ClientCode*>(STRING_CLIENT_A, clientA));
        emperor->nameToClientCodeMap.insert(pair<string, ClientCode*>(STRING_CLIENT_B, clientB));

        // set up dummy dataOutput
        structDataOutput settingDataOutput;
        settingDataOutput.name = STRING_DUMMY;
        DataOutput *dataOutput = new DataOutput(settingDataOutput, emperor->nameToClientCodeMap);
        emperor->nameToDataOutputMap.insert(pair<string, DataOutput *>(STRING_DUMMY, dataOutput));

        // 2. set up coupling algorithms
        const string STRING_AITKEN = "aitken";
        Aitken *aitken = new Aitken(STRING_AITKEN, DOUBY);
        Residual *residual = new Residual(INTY);
        aitken->addResidual(residual, INTY);
        //aitken->addOutput(NULL, INTY);
        emperor->nameToCouplingAlgorithmMap.insert(
                pair<string, AbstractCouplingAlgorithm*>(STRING_AITKEN, aitken));

        // 2.5 set up extrapolator
        const string STRING_EXTRAPOLATOR = "linear";
        LinearExtrapolator *linearExtrapolator = new LinearExtrapolator(STRING_EXTRAPOLATOR);
        emperor->nameToExtrapolatorMap.insert(
                pair<string, AbstractExtrapolator*>(STRING_EXTRAPOLATOR, linearExtrapolator));

        // 3. set up connections
        const string STRING_CONNECTION_A = "connectionA";
        const string STRING_CONNECTION_B = "connectionB";
        Connection *connectionA = new Connection(STRING_CONNECTION_A);
        Connection *connectionB = new Connection(STRING_CONNECTION_B);

        emperor->nameToConnetionMap.insert(
                pair<string, Connection*>(STRING_CONNECTION_A, connectionA));
        emperor->nameToConnetionMap.insert(
                pair<string, Connection*>(STRING_CONNECTION_B, connectionB));

        // 4. set up coupling logic
        structCouplingLogic singleConnection1;
        singleConnection1.type = EMPIRE_connection;
        singleConnection1.connectionRef.connectionName = STRING_CONNECTION_A;
        structCouplingLogic singleConnection2;
        singleConnection2.type = EMPIRE_connection;
        singleConnection2.connectionRef.connectionName = STRING_CONNECTION_B;

        structCouplingLogic iterativeCouplingLoop;
        iterativeCouplingLoop.type = EMPIRE_IterativeCouplingLoop;
        iterativeCouplingLoop.sequence.push_back(singleConnection1);
        iterativeCouplingLoop.sequence.push_back(singleConnection2);
        iterativeCouplingLoop.iterativeCouplingLoop.convergenceChecker.maxNumOfIterations = INTY;
        structCouplingLogic::structIterativeCouplingLoop::structConvergenceChecker::structCheckResidual cr;
        cr.residualRef.couplingAlgorithmName = STRING_AITKEN;
        cr.residualRef.index = INTY;
        iterativeCouplingLoop.iterativeCouplingLoop.convergenceChecker.checkResiduals.push_back(cr);
        iterativeCouplingLoop.iterativeCouplingLoop.convergenceObservers.push_back(STRING_CLIENT_A);
        iterativeCouplingLoop.iterativeCouplingLoop.convergenceObservers.push_back(STRING_CLIENT_B);
        iterativeCouplingLoop.iterativeCouplingLoop.couplingAlgorithmRefs.push_back(STRING_AITKEN);

        structCouplingLogic timeStepLoop;
        timeStepLoop.type = EMPIRE_TimeStepLoop;
        timeStepLoop.sequence.push_back(iterativeCouplingLoop);
        timeStepLoop.timeStepLoop.numTimeSteps = INTY;
        timeStepLoop.timeStepLoop.extrapolatorRef.first = true;
        timeStepLoop.timeStepLoop.extrapolatorRef.second = STRING_EXTRAPOLATOR;
        timeStepLoop.timeStepLoop.dataOutputRefs.push_back(STRING_DUMMY);

        structCouplingLogic optimizationLoop;
        optimizationLoop.type = EMPIRE_OptimizationLoop;
        optimizationLoop.optimizationLoop.maxNumOfIterations = INTY;
        optimizationLoop.optimizationLoop.convergenceSignalSender = STRING_CLIENT_A;
        optimizationLoop.optimizationLoop.convergenceSignalReceivers.push_back(STRING_CLIENT_B);
        optimizationLoop.sequence.push_back(timeStepLoop);

        structCouplingLogic coSimulation;
        coSimulation.type = EMPIRE_CouplingLogicSequence;
        coSimulation.sequence.push_back(optimizationLoop);
        // Finally, test initCouplingLogic()
        MetaDatabase::getSingleton()->settingGlobalCouplingLogic = coSimulation;
        emperor->initGlobalCouplingLogic(); // This is actually what is being tested

        CPPUNIT_ASSERT(typeid(*(emperor->globalCouplingLogic)) == typeid(CouplingLogicSequence));
        CPPUNIT_ASSERT(emperor->globalCouplingLogic->couplingLogicSequence.size() == 1);

        AbstractCouplingLogic *opt = emperor->globalCouplingLogic->couplingLogicSequence[0];
        CPPUNIT_ASSERT( typeid(*opt) == typeid(OptimizationLoop));
        OptimizationLoop *opt2 = dynamic_cast<OptimizationLoop*>(opt);
        CPPUNIT_ASSERT( opt2->maxNumOfIterations == INTY);
        CPPUNIT_ASSERT( opt2->convergenceSignalSender == clientA);
        CPPUNIT_ASSERT( opt2->convergenceSignalReceivers.size() == 1);
        CPPUNIT_ASSERT( opt2->convergenceSignalReceivers[0] == clientB);
        CPPUNIT_ASSERT( opt2->dataOutputVec.size() == 0);
        CPPUNIT_ASSERT( opt2->couplingLogicSequence.size() == 1);

        AbstractCouplingLogic *tsl = opt->couplingLogicSequence[0];
        CPPUNIT_ASSERT( typeid(*tsl) == typeid(TimeStepLoop));
        TimeStepLoop *tsl2 = dynamic_cast<TimeStepLoop*>(tsl);
        CPPUNIT_ASSERT( tsl2->numTimeSteps == INTY);
        CPPUNIT_ASSERT( tsl2->dataOutputVec.size() == 1);
        CPPUNIT_ASSERT( tsl2->dataOutputVec[0] == dataOutput);
        CPPUNIT_ASSERT( tsl2->extrapolator == linearExtrapolator);
        CPPUNIT_ASSERT(tsl->couplingLogicSequence.size() == 1);

        AbstractCouplingLogic *icl = tsl->couplingLogicSequence[0];
        CPPUNIT_ASSERT(typeid(*icl) == typeid(IterativeCouplingLoop));
        IterativeCouplingLoop *icl2 = dynamic_cast<IterativeCouplingLoop*>(icl);
        CPPUNIT_ASSERT(icl2->convergenceChecker->MAX_NUM_ITERATIONS == INTY);
        CPPUNIT_ASSERT( icl2->convergenceChecker->checkResiduals.size() == 1);
        CPPUNIT_ASSERT( icl2->convergenceChecker->checkResiduals[0]->couplingAlgorithm == aitken);
        CPPUNIT_ASSERT( icl2->convergenceChecker->checkResiduals[0]->residualIndex == INTY);
        CPPUNIT_ASSERT(icl2->convergenceObserverVec[0] == clientA);
        CPPUNIT_ASSERT(icl2->convergenceObserverVec[1] == clientB);
        CPPUNIT_ASSERT(icl2->couplingAlgorithmVec[0] == aitken);
        CPPUNIT_ASSERT(icl2->couplingLogicSequence.size() == 2);
        CPPUNIT_ASSERT(icl2->couplingLogicSequence[0] == connectionA);
        CPPUNIT_ASSERT(icl2->couplingLogicSequence[1] == connectionB);
        delete emperor; // the destructor is tested (failed if there is segmentation fault)
    }

CPPUNIT_TEST_SUITE( TestEmperor );
        CPPUNIT_TEST( testInitCouplingLogic);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestEmperor);
