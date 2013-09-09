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

#include "MetaDatabase.h"
#include "Message.h"
#include "ticpp.h"
#include "ConnectionIO.h"

using namespace std;

extern string pathToFolderOfFiles;

namespace EMPIRE {
/********//**
 * \brief Test the class TestMetaDatabase
 ***********/
class TestMetaDatabase: public CppUnit::TestFixture {
private:
public:
    void setUp() {
    }
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief Test case: the filling functions inside MetaDatabase are called to set up the
     *        setting structures. The contest inside the setting structures are compared with the
     *        counterparts in the XML input file
     ***********/
    void testParsing() {
        string inputFile(pathToFolderOfFiles);
        inputFile.append("inputFileOfTestMetaDatabase.xml");
        // 0. initialize Emperor instance
        MetaDatabase::metaDatabase = NULL;
        MetaDatabase::init(const_cast<char*>(inputFile.c_str()));

        vector<structClientCode> &settingClientCodesVec =
                MetaDatabase::getSingleton()->settingClientCodeVec;
        vector<structDataOutput> &settingDataOutputVec =
                MetaDatabase::getSingleton()->settingDataOutputVec;
        vector<structMapper> &settingMapperVec = MetaDatabase::getSingleton()->settingMapperVec;
        std::vector<structCouplingAlgorithm> &settingCouplingAlgorithmVec =
                MetaDatabase::getSingleton()->settingCouplingAlgorithmVec;
        std::vector<structExtrapolator> &settingExtrapolatorVec =
                MetaDatabase::getSingleton()->settingExtrapolatorVec;
        std::vector<structConnection> &settingConnectionVec =
                MetaDatabase::getSingleton()->settingConnectionVec;
        structCouplingLogic & settingCouplingLogic =
                MetaDatabase::getSingleton()->settingGlobalCouplingLogic;

        { // check block client codes
            CPPUNIT_ASSERT(settingClientCodesVec.size() == 2);
            structClientCode client0 = settingClientCodesVec[0];
            CPPUNIT_ASSERT(client0.name == "meshClientA");
            CPPUNIT_ASSERT(client0.meshes.size() == 1);
            CPPUNIT_ASSERT(client0.meshes[0].name == "myMesh");
            CPPUNIT_ASSERT(client0.meshes[0].dataFields.size() == 2);
            CPPUNIT_ASSERT(client0.meshes[0].dataFields[0].name == "displacements");
            CPPUNIT_ASSERT(client0.meshes[0].dataFields[0].dimension == EMPIRE_DataField_vector);
            CPPUNIT_ASSERT(client0.meshes[0].dataFields[0].location == EMPIRE_DataField_atNode);
            CPPUNIT_ASSERT(
                    client0.meshes[0].dataFields[0].typeOfQuantity == EMPIRE_DataField_field);
            CPPUNIT_ASSERT(client0.meshes[0].dataFields[1].name == "forces");
            CPPUNIT_ASSERT(
                    client0.meshes[0].dataFields[1].typeOfQuantity == EMPIRE_DataField_fieldIntegral);
            CPPUNIT_ASSERT(client0.signals.size() == 1);
            CPPUNIT_ASSERT(client0.signals[0].name == "signal");
            CPPUNIT_ASSERT(client0.signals[0].size3D[0] == 2);
            CPPUNIT_ASSERT(client0.signals[0].size3D[1] == 2);
            CPPUNIT_ASSERT(client0.signals[0].size3D[2] == 2);

            structClientCode client1 = settingClientCodesVec[1];
            CPPUNIT_ASSERT(client1.name == "meshClientB");
            CPPUNIT_ASSERT(client1.meshes.size() == 1);
            CPPUNIT_ASSERT(client1.meshes[0].name == "myMesh");
            CPPUNIT_ASSERT(client1.meshes[0].dataFields.size() == 3);
            CPPUNIT_ASSERT(client1.meshes[0].dataFields[1].name == "tractionsElem");
            CPPUNIT_ASSERT(client1.meshes[0].dataFields[1].dimension == EMPIRE_DataField_vector);
            CPPUNIT_ASSERT(
                    client1.meshes[0].dataFields[1].location == EMPIRE_DataField_atElemCentroid);
            CPPUNIT_ASSERT(
                    client1.meshes[0].dataFields[1].typeOfQuantity == EMPIRE_DataField_field);
            CPPUNIT_ASSERT(client1.signals.size() == 1);
            CPPUNIT_ASSERT(client1.signals[0].name == "signal");
            CPPUNIT_ASSERT(client1.signals[0].size3D[0] == 1);
            CPPUNIT_ASSERT(client1.signals[0].size3D[1] == 1);
            CPPUNIT_ASSERT(client1.signals[0].size3D[2] == 5);
        }
        { // check block data outputs
            CPPUNIT_ASSERT(settingDataOutputVec.size()==2);
            { // the 1st dataOutput
                structDataOutput dataOutput = settingDataOutputVec[0];
                CPPUNIT_ASSERT(dataOutput.name=="timeStep");
                CPPUNIT_ASSERT(dataOutput.interval==5);
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs.size()==2);
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[0].clientCodeName=="meshClientA");
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[0].meshName=="myMesh");
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[0].dataFieldName=="displacements");
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[1].clientCodeName=="meshClientB");
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[1].meshName=="myMesh");
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[1].dataFieldName=="forces");
                CPPUNIT_ASSERT(dataOutput.signalRefs.size()==2);
                CPPUNIT_ASSERT(dataOutput.signalRefs[0].clientCodeName=="meshClientA");
                CPPUNIT_ASSERT(dataOutput.signalRefs[0].signalName=="signal");
                CPPUNIT_ASSERT(dataOutput.signalRefs[1].clientCodeName=="meshClientB");
                CPPUNIT_ASSERT(dataOutput.signalRefs[1].signalName=="signal");
            }
            { // the 2nd dataOutput
                structDataOutput dataOutput = settingDataOutputVec[1];
                CPPUNIT_ASSERT(dataOutput.name=="iterativeCoupling");
                CPPUNIT_ASSERT(dataOutput.interval==1);
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs.size()==2);
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[0].clientCodeName=="meshClientA");
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[0].meshName=="myMesh");
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[0].dataFieldName=="displacements");
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[1].clientCodeName=="meshClientB");
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[1].meshName=="myMesh");
                CPPUNIT_ASSERT(dataOutput.dataFieldRefs[1].dataFieldName=="forces");
                CPPUNIT_ASSERT(dataOutput.signalRefs.size()==2);
                CPPUNIT_ASSERT(dataOutput.signalRefs[0].clientCodeName=="meshClientA");
                CPPUNIT_ASSERT(dataOutput.signalRefs[0].signalName=="signal");
                CPPUNIT_ASSERT(dataOutput.signalRefs[1].clientCodeName=="meshClientB");
                CPPUNIT_ASSERT(dataOutput.signalRefs[1].signalName=="signal");
            }
        }
        { // check block mappers
            CPPUNIT_ASSERT(settingMapperVec.size() == 3);
            { // 1st mapper
                structMapper settingMapper = settingMapperVec[0];
                CPPUNIT_ASSERT(settingMapper.name == "mortar1");
                CPPUNIT_ASSERT(settingMapper.meshRefA.clientCodeName == "meshClientA");
                CPPUNIT_ASSERT(settingMapper.meshRefB.clientCodeName == "meshClientB");
                CPPUNIT_ASSERT(settingMapper.meshRefA.meshName == "myMesh");
                CPPUNIT_ASSERT(settingMapper.meshRefB.meshName == "myMesh");
                CPPUNIT_ASSERT(settingMapper.type == EMPIRE_MortarMapper);
                CPPUNIT_ASSERT(settingMapper.mortarMapper.oppositeSurfaceNormal == false);
                CPPUNIT_ASSERT(settingMapper.mortarMapper.dual == true);
                CPPUNIT_ASSERT(settingMapper.mortarMapper.enforceConsistency == true);
            }
            { // 2nd mapper
                structMapper settingMapper = settingMapperVec[1];
                CPPUNIT_ASSERT(settingMapper.name == "nn");
                CPPUNIT_ASSERT(settingMapper.meshRefA.clientCodeName == "meshClientA");
                CPPUNIT_ASSERT(settingMapper.meshRefB.clientCodeName == "meshClientB");
                CPPUNIT_ASSERT(settingMapper.meshRefA.meshName == "myMesh");
                CPPUNIT_ASSERT(settingMapper.meshRefB.meshName == "myMesh");
                CPPUNIT_ASSERT(settingMapper.type == EMPIRE_NearestNeighborMapper);
            }
            { // 3rd mapper
                structMapper settingMapper = settingMapperVec[2];
                CPPUNIT_ASSERT(settingMapper.name == "bi");
                CPPUNIT_ASSERT(settingMapper.meshRefA.clientCodeName == "meshClientA");
                CPPUNIT_ASSERT(settingMapper.meshRefB.clientCodeName == "meshClientB");
                CPPUNIT_ASSERT(settingMapper.meshRefA.meshName == "myMesh");
                CPPUNIT_ASSERT(settingMapper.meshRefB.meshName == "myMesh");
                CPPUNIT_ASSERT(settingMapper.type == EMPIRE_BarycentricInterpolationMapper);
            }
        }
        { // check block coupling algorithms
            CPPUNIT_ASSERT(settingCouplingAlgorithmVec.size() == 1);
            structCouplingAlgorithm settingCoupAlg = settingCouplingAlgorithmVec[0];
            CPPUNIT_ASSERT(settingCoupAlg.type == EMPIRE_Aitken);
            CPPUNIT_ASSERT(settingCoupAlg.name == "aitken A disp");
            CPPUNIT_ASSERT(settingCoupAlg.aitken.initialAitkenFactor == 0.3);
        }
        { // check block extrapolators
            CPPUNIT_ASSERT(settingExtrapolatorVec.size() == 1);
            structExtrapolator settingExtrapolator = settingExtrapolatorVec[0];
            CPPUNIT_ASSERT(settingExtrapolator.type == EMPIRE_SimpleExtrapolator);
            CPPUNIT_ASSERT(settingExtrapolator.name == "simple");
        }
        { // check block connections
            CPPUNIT_ASSERT(settingConnectionVec.size() == 3);
            {
                structConnection settingConnection = settingConnectionVec[0];
                CPPUNIT_ASSERT(settingConnection.name == "transfer displacements");
                CPPUNIT_ASSERT(settingConnection.inputs.size() == 1);
                CPPUNIT_ASSERT(settingConnection.inputs[0].type == EMPIRE_ConnectionIO_DataField);
                CPPUNIT_ASSERT(
                        settingConnection.inputs[0].dataFieldRef.clientCodeName == "meshClientA");
                CPPUNIT_ASSERT(settingConnection.inputs[0].dataFieldRef.meshName == "myMesh");
                CPPUNIT_ASSERT(
                        settingConnection.inputs[0].dataFieldRef.dataFieldName == "displacements");
                CPPUNIT_ASSERT(settingConnection.outputs.size() == 1);
                CPPUNIT_ASSERT(settingConnection.outputs[0].type == EMPIRE_ConnectionIO_DataField);
                CPPUNIT_ASSERT(
                        settingConnection.outputs[0].dataFieldRef.clientCodeName == "meshClientB");
                CPPUNIT_ASSERT(settingConnection.outputs[0].dataFieldRef.meshName == "myMesh");
                CPPUNIT_ASSERT(
                        settingConnection.outputs[0].dataFieldRef.dataFieldName == "displacements");
                CPPUNIT_ASSERT(settingConnection.filterSequence.size() == 2);
                {
                    structFilter &settingfilter = settingConnection.filterSequence[0];
                    CPPUNIT_ASSERT(settingfilter.type == EMPIRE_CouplingAlgorithmFilter);
                    CPPUNIT_ASSERT(
                            settingfilter.couplingAlgorithmFilter.couplingAlgorithmName == "aitken A disp");
                    CPPUNIT_ASSERT(settingfilter.inputs.size() == 1);
                    CPPUNIT_ASSERT(settingfilter.inputs[0].type == EMPIRE_ConnectionIO_DataField);
                    CPPUNIT_ASSERT(
                            settingfilter.inputs[0].dataFieldRef.clientCodeName == "meshClientA");
                    CPPUNIT_ASSERT(settingfilter.inputs[0].dataFieldRef.meshName== "myMesh");
                    CPPUNIT_ASSERT(
                            settingfilter.inputs[0].dataFieldRef.dataFieldName == "displacements");
                    CPPUNIT_ASSERT(settingfilter.outputs.size() == 1);
                    CPPUNIT_ASSERT(settingfilter.outputs[0].type == EMPIRE_ConnectionIO_DataField);
                    CPPUNIT_ASSERT(
                            settingfilter.outputs[0].dataFieldRef.clientCodeName == "meshClientA");
                    CPPUNIT_ASSERT(settingfilter.outputs[0].dataFieldRef.meshName== "myMesh");
                    CPPUNIT_ASSERT(
                            settingfilter.outputs[0].dataFieldRef.dataFieldName == "displacements");
                }
                {
                    structFilter &settingfilter = settingConnection.filterSequence[1];
                    CPPUNIT_ASSERT(settingfilter.type == EMPIRE_MappingFilter);
                    CPPUNIT_ASSERT(settingfilter.mappingFilter.mapperName == "mortar1");
                    CPPUNIT_ASSERT(settingfilter.inputs.size() == 1);
                    CPPUNIT_ASSERT(settingfilter.inputs[0].type == EMPIRE_ConnectionIO_DataField);
                    CPPUNIT_ASSERT(
                            settingfilter.inputs[0].dataFieldRef.clientCodeName == "meshClientA");
                    CPPUNIT_ASSERT(settingfilter.inputs[0].dataFieldRef.meshName== "myMesh");
                    CPPUNIT_ASSERT(
                            settingfilter.inputs[0].dataFieldRef.dataFieldName == "displacements");
                    CPPUNIT_ASSERT(settingfilter.outputs.size() == 1);
                    CPPUNIT_ASSERT(settingfilter.outputs[0].type == EMPIRE_ConnectionIO_DataField);
                    CPPUNIT_ASSERT(
                            settingfilter.outputs[0].dataFieldRef.clientCodeName == "meshClientB");
                    CPPUNIT_ASSERT(settingfilter.outputs[0].dataFieldRef.meshName== "myMesh");
                    CPPUNIT_ASSERT(
                            settingfilter.outputs[0].dataFieldRef.dataFieldName == "displacements");
                }
            }
            {
                structConnection settingConnection = settingConnectionVec[1];
                CPPUNIT_ASSERT(settingConnection.name == "transfer forces");
                CPPUNIT_ASSERT(settingConnection.inputs.size() == 1);
                CPPUNIT_ASSERT(settingConnection.inputs[0].type == EMPIRE_ConnectionIO_DataField);
                CPPUNIT_ASSERT(
                        settingConnection.inputs[0].dataFieldRef.clientCodeName == "meshClientB");
                CPPUNIT_ASSERT(settingConnection.inputs[0].dataFieldRef.meshName == "myMesh");
                CPPUNIT_ASSERT(
                        settingConnection.inputs[0].dataFieldRef.dataFieldName == "tractionsElem");
                CPPUNIT_ASSERT(settingConnection.outputs.size() == 1);
                CPPUNIT_ASSERT(settingConnection.outputs[0].type == EMPIRE_ConnectionIO_DataField);
                CPPUNIT_ASSERT(
                        settingConnection.outputs[0].dataFieldRef.clientCodeName == "meshClientA");
                CPPUNIT_ASSERT(settingConnection.outputs[0].dataFieldRef.meshName == "myMesh");
                CPPUNIT_ASSERT(settingConnection.outputs[0].dataFieldRef.dataFieldName == "forces");
                CPPUNIT_ASSERT(settingConnection.filterSequence.size() == 2);

                {
                    structFilter &settingfilter = settingConnection.filterSequence[0];
                    CPPUNIT_ASSERT(settingfilter.type == EMPIRE_LocationFilter);
                    CPPUNIT_ASSERT(settingfilter.inputs.size() == 1);
                    CPPUNIT_ASSERT(settingfilter.inputs[0].type == EMPIRE_ConnectionIO_DataField);
                    CPPUNIT_ASSERT(
                            settingfilter.inputs[0].dataFieldRef.clientCodeName == "meshClientB");
                    CPPUNIT_ASSERT(settingfilter.inputs[0].dataFieldRef.meshName== "myMesh");
                    CPPUNIT_ASSERT(
                            settingfilter.inputs[0].dataFieldRef.dataFieldName == "tractionsElem");
                    CPPUNIT_ASSERT(settingfilter.outputs.size() == 1);
                    CPPUNIT_ASSERT(settingfilter.outputs[0].type == EMPIRE_ConnectionIO_DataField);
                    CPPUNIT_ASSERT(
                            settingfilter.outputs[0].dataFieldRef.clientCodeName == "meshClientB");
                    CPPUNIT_ASSERT(settingfilter.outputs[0].dataFieldRef.meshName== "myMesh");
                    CPPUNIT_ASSERT(
                            settingfilter.outputs[0].dataFieldRef.dataFieldName == "tractionsNode");
                }

                {
                    structFilter &settingfilter = settingConnection.filterSequence[1];
                    CPPUNIT_ASSERT(settingfilter.type == EMPIRE_MappingFilter);
                    CPPUNIT_ASSERT(settingfilter.mappingFilter.mapperName == "mortar1");
                    CPPUNIT_ASSERT(settingfilter.inputs.size() == 1);
                    CPPUNIT_ASSERT(settingfilter.inputs[0].type == EMPIRE_ConnectionIO_DataField);
                    CPPUNIT_ASSERT(
                            settingfilter.inputs[0].dataFieldRef.clientCodeName == "meshClientB");
                    CPPUNIT_ASSERT(settingfilter.inputs[0].dataFieldRef.meshName== "myMesh");
                    CPPUNIT_ASSERT(
                            settingfilter.inputs[0].dataFieldRef.dataFieldName == "tractionsNode");
                    CPPUNIT_ASSERT(settingfilter.outputs.size() == 1);
                    CPPUNIT_ASSERT(settingfilter.outputs[0].type == EMPIRE_ConnectionIO_DataField);
                    CPPUNIT_ASSERT(
                            settingfilter.outputs[0].dataFieldRef.clientCodeName == "meshClientA");
                    CPPUNIT_ASSERT(settingfilter.outputs[0].dataFieldRef.meshName== "myMesh");
                    CPPUNIT_ASSERT(
                            settingfilter.outputs[0].dataFieldRef.dataFieldName == "forces");
                }
            }
            {
                structConnection &settingConnection = settingConnectionVec[2];
                CPPUNIT_ASSERT(settingConnection.name == "transfer signal");
                CPPUNIT_ASSERT(settingConnection.inputs.size() == 1);
                CPPUNIT_ASSERT(settingConnection.inputs[0].type == EMPIRE_ConnectionIO_Signal);
                CPPUNIT_ASSERT(
                        settingConnection.inputs[0].signalRef.clientCodeName == "meshClientA");
                CPPUNIT_ASSERT(settingConnection.inputs[0].signalRef.signalName== "signal");
                CPPUNIT_ASSERT(settingConnection.outputs.size() == 1);
                CPPUNIT_ASSERT(settingConnection.outputs[0].type == EMPIRE_ConnectionIO_Signal);
                CPPUNIT_ASSERT(
                        settingConnection.outputs[0].signalRef.clientCodeName == "meshClientB");
                CPPUNIT_ASSERT(settingConnection.outputs[0].signalRef.signalName== "signal");
                {
                    structFilter &settingfilter = settingConnection.filterSequence[0];
                    CPPUNIT_ASSERT(settingfilter.type == EMPIRE_CopyFilter);
                    CPPUNIT_ASSERT(settingfilter.inputs.size() == 1);
                    CPPUNIT_ASSERT(settingfilter.inputs[0].type == EMPIRE_ConnectionIO_Signal);
                    CPPUNIT_ASSERT(
                            settingfilter.inputs[0].signalRef.clientCodeName == "meshClientA");
                    CPPUNIT_ASSERT(settingfilter.inputs[0].signalRef.signalName== "signal");
                    CPPUNIT_ASSERT(settingfilter.outputs.size() == 1);
                    CPPUNIT_ASSERT(settingfilter.outputs[0].type == EMPIRE_ConnectionIO_Signal);
                    CPPUNIT_ASSERT(
                            settingfilter.outputs[0].signalRef.clientCodeName == "meshClientB");
                    CPPUNIT_ASSERT(settingfilter.outputs[0].signalRef.signalName== "signal");
                }
            }
        }
        // check block coSimulation
        CPPUNIT_ASSERT(settingCouplingLogic.type == EMPIRE_CouplingLogicSequence);

        CPPUNIT_ASSERT(settingCouplingLogic.sequence.size() == 1);
        structCouplingLogic settingTSL = settingCouplingLogic.sequence[0];
        CPPUNIT_ASSERT(settingTSL.type == EMPIRE_TimeStepLoop);
        CPPUNIT_ASSERT(settingTSL.timeStepLoop.numTimeSteps == 5);
        //CPPUNIT_ASSERT(settingTSL.timeStepLoop.extrapolatorRefs.size() == 1);
        //CPPUNIT_ASSERT(settingTSL.timeStepLoop.extrapolatorRefs[0] == "transfer displacements");

        CPPUNIT_ASSERT(settingTSL.timeStepLoop.dataOutputRefs.size() == 1);
        CPPUNIT_ASSERT(settingTSL.timeStepLoop.dataOutputRefs[0] == "timeStep");
        CPPUNIT_ASSERT(settingTSL.timeStepLoop.extrapolatorRefs.size() == 1);
        CPPUNIT_ASSERT(settingTSL.timeStepLoop.extrapolatorRefs[0] == "simple");
        CPPUNIT_ASSERT(settingTSL.sequence.size() == 1);
        structCouplingLogic settingICL = settingTSL.sequence[0];
        {
            CPPUNIT_ASSERT(settingICL.type == EMPIRE_IterativeCouplingLoop);
            string stmp =
                    settingICL.iterativeCouplingLoop.convergenceChecker.dataFieldRef.clientCodeName;
            CPPUNIT_ASSERT(stmp == "meshClientA");
            stmp = settingICL.iterativeCouplingLoop.convergenceChecker.dataFieldRef.meshName;
            CPPUNIT_ASSERT(stmp == "myMesh");
            stmp = settingICL.iterativeCouplingLoop.convergenceChecker.dataFieldRef.dataFieldName;
            CPPUNIT_ASSERT(stmp == "displacements");
            CPPUNIT_ASSERT(
                    settingICL.iterativeCouplingLoop.convergenceChecker.absoluteTolerance==1e-12);
            CPPUNIT_ASSERT(
                    settingICL.iterativeCouplingLoop.convergenceChecker.relativeTolerance==1e-6);
            CPPUNIT_ASSERT(
                    settingICL.iterativeCouplingLoop.convergenceChecker.maxNumOfIterations==100);
            CPPUNIT_ASSERT(settingICL.iterativeCouplingLoop.convergenceObservers.size() == 2);
            stmp = settingICL.iterativeCouplingLoop.convergenceObservers[0];
            CPPUNIT_ASSERT(stmp == "meshClientA");
            stmp = settingICL.iterativeCouplingLoop.convergenceObservers[1];
            CPPUNIT_ASSERT(stmp == "meshClientB");
            CPPUNIT_ASSERT(settingICL.iterativeCouplingLoop.couplingAlgorithmRefs.size() == 1);
            stmp = settingICL.iterativeCouplingLoop.couplingAlgorithmRefs[0];
            CPPUNIT_ASSERT(stmp == "aitken A disp");
            CPPUNIT_ASSERT(settingICL.iterativeCouplingLoop.dataOutputRefs.size() == 1);
            stmp = settingICL.iterativeCouplingLoop.dataOutputRefs[0];
            CPPUNIT_ASSERT(stmp == "iterativeCoupling");
            CPPUNIT_ASSERT(settingICL.sequence.size() == 3);
            CPPUNIT_ASSERT(settingICL.sequence[0].type == EMPIRE_connection);
            stmp = settingICL.sequence[0].connectionRef.connectionName;
            CPPUNIT_ASSERT(stmp == "transfer displacements");
            CPPUNIT_ASSERT(settingICL.sequence[1].type == EMPIRE_connection);
            stmp = settingICL.sequence[1].connectionRef.connectionName;
            CPPUNIT_ASSERT(stmp == "transfer forces");
            CPPUNIT_ASSERT(settingICL.sequence[2].type == EMPIRE_connection);
            stmp = settingICL.sequence[2].connectionRef.connectionName;
            CPPUNIT_ASSERT(stmp == "transfer signal");
        }

        // block general
        CPPUNIT_ASSERT(MetaDatabase::getSingleton()->verbosity == "DeBuG");
        /// Test CompareStringInsensitive
        CPPUNIT_ASSERT(Message::userSetOutputLevel==Message::DEBUG);
        CPPUNIT_ASSERT(MetaDatabase::getSingleton()->serverPortFile == "server.port");
    }

CPPUNIT_TEST_SUITE( TestMetaDatabase );
        CPPUNIT_TEST( testParsing);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestMetaDatabase);
