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
#include <iostream>
#include <string.h>
#include <assert.h>
#include <sstream>

#include "MetaDatabase.h"
#include "ticpp.h"
#include "Message.h"
#include "AuxiliaryFunctions.h"

using namespace std;
using namespace ticpp;

namespace EMPIRE {

MetaDatabase *MetaDatabase::metaDatabase = NULL;

void MetaDatabase::init(char *inputFileName) {
    assert(metaDatabase == NULL);
    metaDatabase = new MetaDatabase(inputFileName);
}

MetaDatabase *MetaDatabase::getSingleton() {
    assert(metaDatabase != NULL);
    return metaDatabase;
}

MetaDatabase::MetaDatabase() {
    // do nothing
}

MetaDatabase::MetaDatabase(char *inputFileName) {
    try {
        inputFile = new Document(inputFileName);
        inputFile->LoadFile();
        /// Fill up data base
        fillServerPortFile();
        fillVerbosity();
        fillSettingClientCodesVec();
        fillSettingDataOutputVec();
        fillSettingMapperVec();
        fillSettingCouplingAlgorithmVec();
        fillSettingExtrapolatorVec();
        fillSettingConnectionVec();
        fillSettingCouplingLogic();
    } catch (ticpp::Exception& ex) {
        cout << "ERROR Parser: " << ex.what() << endl;
        exit(EXIT_FAILURE);
    }
}

MetaDatabase::~MetaDatabase() {
    metaDatabase = NULL;
}

void MetaDatabase::fillServerPortFile() {
    Element *pXMLElement =
            inputFile->FirstChildElement()->FirstChildElement("general")->FirstChildElement(
                    "portFile");
    serverPortFile = pXMLElement->GetText();
}

void MetaDatabase::fillVerbosity() {
    Element *pXMLElement =
            inputFile->FirstChildElement()->FirstChildElement("general")->FirstChildElement(
                    "verbosity");
    verbosity = pXMLElement->GetText();
    ///Set verbosity in Message class for output
    if (AuxiliaryFunctions::CompareStringInsensitive(verbosity, "debug")) {
        Message::userSetOutputLevel = Message::DEBUG;
    } else if (AuxiliaryFunctions::CompareStringInsensitive(verbosity, "info")) {
        Message::userSetOutputLevel = Message::INFO;
    } else if (AuxiliaryFunctions::CompareStringInsensitive(verbosity, "warning")) {
        Message::userSetOutputLevel = Message::WARNING;
    } else if (AuxiliaryFunctions::CompareStringInsensitive(verbosity, "error")) {
        Message::userSetOutputLevel = Message::ERROR;
    } else {
        Message::userSetOutputLevel = Message::INFO;
    }

}

bool MetaDatabase::checkForClientCodeName(std::string clientName) {
    for (int i = 0; i < settingClientCodeVec.size(); i++)
        if (settingClientCodeVec[i].name == clientName)
            return true;
    return false;
}

void MetaDatabase::fillSettingClientCodesVec() {
    assert(settingClientCodeVec.size() == 0);
    ticpp::Element *xmlEMPEROR = inputFile->FirstChildElement("EMPEROR");
    ticpp::Iterator<Element> xmlClientCode("clientCode");
    for (xmlClientCode = xmlClientCode.begin(xmlEMPEROR); xmlClientCode != xmlClientCode.end();
            xmlClientCode++) {
        structClientCode clientCode;
        assert(xmlClientCode->HasAttribute("name"));
        clientCode.name = xmlClientCode->GetAttribute("name");

        ticpp::Iterator<Element> xmlMesh("mesh");

        for (xmlMesh = xmlMesh.begin(xmlClientCode.Get()); xmlMesh != xmlMesh.end(); xmlMesh++) {
            structClientCode::structMesh mesh;
            mesh.name = xmlMesh->GetAttribute("name");
            if (xmlMesh->HasAttribute("type")) {
                string meshType = xmlMesh->GetAttribute("type");
                if (meshType == "FEMesh") {
                    mesh.type = EMPIRE_Mesh_FEMesh;
                } else if (meshType == "IGAMesh") {
                    mesh.type = EMPIRE_Mesh_IGAMesh;
                } else {
                    assert(false);
                }
            } else {
                mesh.type = EMPIRE_Mesh_FEMesh;
            }
            ticpp::Iterator<Element> xmlDataField("dataField");
            for (xmlDataField = xmlDataField.begin(xmlMesh.Get());
                    xmlDataField != xmlDataField.end(); xmlDataField++) {
                structClientCode::structMesh::structDataField dataField;
                assert(xmlDataField->HasAttribute("name"));
                assert(xmlDataField->HasAttribute("dimension"));
                assert(xmlDataField->HasAttribute("location"));
                assert(xmlDataField->HasAttribute("typeOfQuantity"));
                dataField.name = xmlDataField->GetAttribute("name");
                if (xmlDataField->GetAttribute("dimension") == "vector")
                    dataField.dimension = EMPIRE_DataField_vector;
                else if (xmlDataField->GetAttribute("dimension") == "scalar")
                    dataField.dimension = EMPIRE_DataField_scalar;
                else
                    assert(false);
                if (xmlDataField->GetAttribute("location") == "atNode")
                    dataField.location = EMPIRE_DataField_atNode;
                else if (xmlDataField->GetAttribute("location") == "atElemCentroid")
                    dataField.location = EMPIRE_DataField_atElemCentroid;
                else
                    assert(false);
                if (xmlDataField->GetAttribute("typeOfQuantity") == "field")
                    dataField.typeOfQuantity = EMPIRE_DataField_field;
                else if (xmlDataField->GetAttribute("typeOfQuantity") == "fieldIntegral")
                    dataField.typeOfQuantity = EMPIRE_DataField_fieldIntegral;
                else
                    assert(false);

                mesh.dataFields.push_back(dataField);
            }
            clientCode.meshes.push_back(mesh);
        }

        ticpp::Iterator<Element> xmlSignal("signal");
        for (xmlSignal = xmlSignal.begin(xmlClientCode.Get()); xmlSignal != xmlSignal.end();
                xmlSignal++) {
            structClientCode::structSignal signal;
            signal.name = xmlSignal->GetAttribute("name");
            string size = xmlSignal->GetAttribute("size");
            { // get the size
                stringstream ss(size);
                int count = 0;
                int size[3];
                for (int i = 0; i < 3; i++) {
                    if (ss.eof()) {
                        break;
                    } else {
                        ss >> size[i];
                        count++;
                    }
                }
                if (count == 1) {
                    signal.size3D[0] = 1;
                    signal.size3D[1] = 1;
                    signal.size3D[2] = size[0];
                } else if (count == 2) {
                    signal.size3D[0] = 1;
                    signal.size3D[1] = size[0];
                    signal.size3D[2] = size[1];
                } else if (count == 3) {
                    signal.size3D[0] = size[0];
                    signal.size3D[1] = size[1];
                    signal.size3D[2] = size[2];
                } else {
                    assert(false);
                }
            }
            clientCode.signals.push_back(signal);
        }
        settingClientCodeVec.push_back(clientCode);
    }
}

void MetaDatabase::fillSettingDataOutputVec() {
    assert(settingDataOutputVec.size()==0);
    ticpp::Element *xmlEMPEROR = inputFile->FirstChildElement("EMPEROR");
    ticpp::Iterator<Element> xmlDataOutput("dataOutput");
    for (xmlDataOutput = xmlDataOutput.begin(xmlEMPEROR); xmlDataOutput != xmlDataOutput.end();
            xmlDataOutput++) {
        structDataOutput dataOutput;
        dataOutput.name = xmlDataOutput->GetAttribute("name");
        dataOutput.interval = xmlDataOutput->GetAttribute<int>("interval");
        dataOutput.connectionIOs = parseConnectionIORefs(xmlDataOutput.Get());

        settingDataOutputVec.push_back(dataOutput);
    }
}

void MetaDatabase::fillSettingMapperVec() {
    assert(settingMapperVec.size() == 0);
    ticpp::Element *xmlEMPEROR = inputFile->FirstChildElement("EMPEROR");
    ticpp::Iterator<Element> xmlMapper("mapper");
    for (xmlMapper = xmlMapper.begin(xmlEMPEROR); xmlMapper != xmlMapper.end(); xmlMapper++) {
        structMapper mapper;
        mapper.name = xmlMapper->GetAttribute("name");
        ticpp::Element *xmlMeshRefA = xmlMapper->FirstChildElement("meshA")->FirstChildElement(
                "meshRef");
        mapper.meshRefA.clientCodeName = xmlMeshRefA->GetAttribute("clientCodeName");
        mapper.meshRefA.meshName = xmlMeshRefA->GetAttribute("meshName");
        ticpp::Element *xmlMeshRefB = xmlMapper->FirstChildElement("meshB")->FirstChildElement(
                "meshRef");
        mapper.meshRefB.clientCodeName = xmlMeshRefB->GetAttribute("clientCodeName");
        mapper.meshRefB.meshName = xmlMeshRefB->GetAttribute("meshName");
        if (xmlMapper->GetAttribute("type") == "mortarMapper") {
            mapper.type = EMPIRE_MortarMapper;
            ticpp::Element *xmlMortarMapper = xmlMapper->FirstChildElement("mortarMapper");
            if (xmlMortarMapper->GetAttribute("oppositeSurfaceNormal") == "true")
                mapper.mortarMapper.oppositeSurfaceNormal = true;
            else if (xmlMortarMapper->GetAttribute("oppositeSurfaceNormal") == "false")
                mapper.mortarMapper.oppositeSurfaceNormal = false;
            else
                assert(false);
            if (xmlMortarMapper->GetAttribute("dual") == "true")
                mapper.mortarMapper.dual = true;
            else if (xmlMortarMapper->GetAttribute("dual") == "false")
                mapper.mortarMapper.dual = false;
            else
                assert(false);
            if (xmlMortarMapper->GetAttribute("enforceConsistency") == "true")
                mapper.mortarMapper.enforceConsistency = true;
            else if (xmlMortarMapper->GetAttribute("enforceConsistency") == "false")
                mapper.mortarMapper.enforceConsistency = false;
            else
                assert(false);
        } else if (xmlMapper->GetAttribute("type") == "nearestNeighborMapper") {
            mapper.type = EMPIRE_NearestNeighborMapper;
        } else if (xmlMapper->GetAttribute("type") == "barycentricInterpolationMapper") {
            mapper.type = EMPIRE_BarycentricInterpolationMapper;
        } else if (xmlMapper->GetAttribute("type") == "IGAMortarMapper") {
            mapper.type = EMPIRE_IGAMortarMapper;
            ticpp::Element *xmlIGAMortar = xmlMapper->FirstChildElement("IGAMortarMapper");
            mapper.igaMortarMapper.tolProjectionDistance = xmlIGAMortar->GetAttribute<double>(
                    "tolProjectionDistance");
            mapper.igaMortarMapper.numGPsTriangle = xmlIGAMortar->GetAttribute<int>(
                    "numGPsTriangle");
            mapper.igaMortarMapper.numGPsQuad = xmlIGAMortar->GetAttribute<int>("numGPsQuad");
        } else {
            assert(false);
        }
        settingMapperVec.push_back(mapper);
    }
}

void MetaDatabase::fillSettingCouplingAlgorithmVec() {
    assert(settingCouplingAlgorithmVec.size() == 0);
    ticpp::Element *xmlEMPEROR = inputFile->FirstChildElement("EMPEROR");
    ticpp::Iterator<Element> xmlCoupAlg("couplingAlgorithm");
    for (xmlCoupAlg = xmlCoupAlg.begin(xmlEMPEROR); xmlCoupAlg != xmlCoupAlg.end(); xmlCoupAlg++) {
        structCouplingAlgorithm coupAlg;
        coupAlg.name = xmlCoupAlg->GetAttribute("name");

        { // residuals
            ticpp::Iterator<Element> xmlResidual("residual");
            for (xmlResidual = xmlResidual.begin(xmlCoupAlg.Get());
                    xmlResidual != xmlResidual.end(); xmlResidual++) {
                structResidual residual;
                residual.index = xmlResidual->GetAttribute<int>("index");
                ticpp::Iterator<Element> xmlCompenent("component");
                for (xmlCompenent = xmlCompenent.begin(xmlResidual.Get());
                        xmlCompenent != xmlCompenent.end(); xmlCompenent++) {
                    structResidual::structComponent component;
                    component.coefficient = xmlCompenent->GetAttribute<double>("coefficient");
                    component.timeToUpdate = xmlCompenent->GetAttribute<string>("timeToUpdate");

                    structConnectionIO io = parseConnectionIORef(xmlCompenent.Get());
                    component.connectionIO = io;

                    residual.components.push_back(component);
                }
                coupAlg.residuals.push_back(residual);
            }
        }
        { // outputs
            ticpp::Iterator<Element> xmlOutput("output");
            for (xmlOutput = xmlOutput.begin(xmlCoupAlg.Get()); xmlOutput != xmlOutput.end();
                    xmlOutput++) {
                structCouplingAlgorithm::structOutput output;
                output.index = xmlOutput->GetAttribute<int>("index");

                structConnectionIO io = parseConnectionIORef(xmlOutput.Get());
                output.connectionIO = io;

                coupAlg.outputs.push_back(output);
            }
        }
        if (xmlCoupAlg->GetAttribute("type") == "aitken") {
            coupAlg.type = EMPIRE_Aitken;
            ticpp::Element *xmlAitken = xmlCoupAlg->FirstChildElement("aitken");
            double tmpDouble = xmlAitken->GetAttribute<double>("initialRelaxationFactor");
            coupAlg.aitken.initialRelaxationFactor = tmpDouble;
        } else if (xmlCoupAlg->GetAttribute("type") == "constantRelaxation") {
            coupAlg.type = EMPIRE_ConstantRelaxation;
            ticpp::Element *xmlAitken = xmlCoupAlg->FirstChildElement("constantRelaxation");
            double tmpDouble = xmlAitken->GetAttribute<double>("relaxationFactor");
            coupAlg.constantRelaxation.relaxationFactor = tmpDouble;
        } else {
            assert(false);
        }
        settingCouplingAlgorithmVec.push_back(coupAlg);
    }
}

void MetaDatabase::fillSettingExtrapolatorVec() {
    assert(settingExtrapolatorVec.size() == 0);
    ticpp::Element *xmlEMPEROR = inputFile->FirstChildElement("EMPEROR");
    ticpp::Iterator<Element> xmlExtrapolator("extrapolator");
    for (xmlExtrapolator = xmlExtrapolator.begin(xmlEMPEROR);
            xmlExtrapolator != xmlExtrapolator.end(); xmlExtrapolator++) {
        structExtrapolator settingExtrapolator;
        settingExtrapolator.name = xmlExtrapolator->GetAttribute("name");
        if (xmlExtrapolator->GetAttribute("type") == "linearExtrapolator") {
            settingExtrapolator.type = EMPIRE_LinearExtrapolator;
        } else {
            assert(false);
        }
        settingExtrapolator.connectionIOs = parseConnectionIORefs(xmlExtrapolator.Get());
        settingExtrapolatorVec.push_back(settingExtrapolator);
    }
}

void MetaDatabase::fillSettingConnectionVec() {
    assert(settingConnectionVec.size() == 0);
    ticpp::Element *xmlEMPEROR = inputFile->FirstChildElement("EMPEROR");
    ticpp::Iterator<Element> xmlConnection("connection");
    for (xmlConnection = xmlConnection.begin(xmlEMPEROR); xmlConnection != xmlConnection.end();
            xmlConnection++) {
        structConnection connection;
        connection.name = xmlConnection->GetAttribute("name");
        // inputs and outputs
        if (xmlConnection->FirstChildElement("inputAndOutput", false) != NULL) {
            ticpp::Element *xmlIO = xmlConnection->FirstChildElement("inputAndOutput");

            structConnectionIO io = parseConnectionIORef(xmlIO);
            connection.inputs.push_back(io);
            connection.outputs.push_back(io);
        } else {
            ticpp::Iterator<Element> xmlInput("input");
            for (xmlInput = xmlInput.begin(xmlConnection.Get()); xmlInput != xmlInput.end();
                    xmlInput++) {
                structConnectionIO input = parseConnectionIORef(xmlInput.Get());
                connection.inputs.push_back(input);
            }
            ticpp::Iterator<Element> xmlOutput("output");
            for (xmlOutput = xmlOutput.begin(xmlConnection.Get()); xmlOutput != xmlOutput.end();
                    xmlOutput++) {
                structConnectionIO output = parseConnectionIORef(xmlOutput.Get());
                connection.outputs.push_back(output);
            }
        }

        if (xmlConnection->FirstChildElement("sequence", false) != NULL) {
            ticpp::Element *xmlFilters = xmlConnection->FirstChildElement("sequence");
            ticpp::Iterator<Element> xmlFilter("filter");
            for (xmlFilter = xmlFilter.begin(xmlFilters); xmlFilter != xmlFilter.end();
                    xmlFilter++) {
                structFilter filter;
                // inputs and outputs
                if (xmlFilter->FirstChildElement("inputAndOutput", false) != NULL) {
                    ticpp::Element *xmlIO = xmlFilter->FirstChildElement("inputAndOutput");
                    structConnectionIO io = parseConnectionIORef(xmlIO);
                    filter.inputs.push_back(io);
                    filter.outputs.push_back(io);
                } else {
                    ticpp::Iterator<Element> xmlInput("input");
                    for (xmlInput = xmlInput.begin(xmlFilter.Get()); xmlInput != xmlInput.end();
                            xmlInput++) {
                        structConnectionIO input = parseConnectionIORef(xmlInput.Get());
                        filter.inputs.push_back(input);
                    }
                    ticpp::Iterator<Element> xmlOutput("output");
                    for (xmlOutput = xmlOutput.begin(xmlFilter.Get()); xmlOutput != xmlOutput.end();
                            xmlOutput++) {
                        structConnectionIO output = parseConnectionIORef(xmlOutput.Get());
                        filter.outputs.push_back(output);
                    }
                }
                // filter type
                if (xmlFilter->GetAttribute("type") == "mappingFilter") {
                    filter.type = EMPIRE_MappingFilter;
                    filter.mappingFilter.mapperName =
                            xmlFilter->FirstChildElement("mappingFilter")->FirstChildElement(
                                    "mapperRef")->GetAttribute("mapperName");
                } else if (xmlFilter->GetAttribute("type") == "locationFilter") {
                    filter.type = EMPIRE_LocationFilter;
                } else if (xmlFilter->GetAttribute("type") == "scalingFilter") {
                    filter.type = EMPIRE_ScalingFilter;
                    filter.scalingFilter.factor =
                            xmlFilter->FirstChildElement("scalingFilter")->GetAttribute<double>(
                                    "factor");
                } else if (xmlFilter->GetAttribute("type") == "setFilter") {
                    filter.type = EMPIRE_SetFilter;
                    string valueString = xmlFilter->FirstChildElement("setFilter")->GetAttribute<
                            string>("value");
                    std::stringstream ss(valueString);
                    while (!ss.eof()) {
                        double value;
                        ss >> value;
                        filter.setFilter.value.push_back(value);
                    }
                } else if (xmlFilter->GetAttribute("type") == "copyFilter") {
                    filter.type = EMPIRE_CopyFilter;
                } else if (xmlFilter->GetAttribute("type") == "dataFieldIntegrationFilter") {
                    filter.type = EMPIRE_DataFieldIntegrationFilter;
                    ticpp::Element *xmlMeshRef = xmlFilter->FirstChildElement(
                            "dataFieldIntegrationFilter")->FirstChildElement("meshRef");
                    filter.dataFieldIntegrationFilter.meshRef.clientCodeName =
                            xmlMeshRef->GetAttribute("clientCodeName");
                    filter.dataFieldIntegrationFilter.meshRef.meshName = xmlMeshRef->GetAttribute(
                            "meshName");

                } else {
                    assert(false);
                }
                connection.filterSequence.push_back(filter);
            }
        }
        settingConnectionVec.push_back(connection);
    }
}

void parseCouplingLogicBlock(ticpp::Iterator<Element> &xmlCouplingLogicIn,
        structCouplingLogic &couplingLogicIn) { // a global function instead of a member function
    // 1. parse coupling logic setting
    if (xmlCouplingLogicIn->GetAttribute("type") == "timeStepLoop") {
        couplingLogicIn.type = EMPIRE_TimeStepLoop;
        ticpp::Element *xmlTimeStepLoop = xmlCouplingLogicIn->FirstChildElement("timeStepLoop");
        int numTimeSteps = xmlTimeStepLoop->GetAttribute<int>("numTimeSteps");
        couplingLogicIn.timeStepLoop.numTimeSteps = numTimeSteps;

        { // add extrapolator
            ticpp::Element *xmlExtrapolatorRef = xmlTimeStepLoop->FirstChildElement(
                    "extrapolatorRef", false);
            if (xmlExtrapolatorRef != NULL) {
                couplingLogicIn.timeStepLoop.extrapolatorRef.first = true;
                couplingLogicIn.timeStepLoop.extrapolatorRef.second =
                        xmlExtrapolatorRef->GetAttribute("extrapolatorName");
            } else {
                couplingLogicIn.iterativeCouplingLoop.couplingAlgorithmRef.first = false;
            }
        }

        { // add dataOutputs
            ticpp::Iterator<Element> xmlDataOutputRef("dataOutputRef");
            for (xmlDataOutputRef = xmlDataOutputRef.begin(xmlTimeStepLoop);
                    xmlDataOutputRef != xmlDataOutputRef.end(); xmlDataOutputRef++) {
                string dataOutputName = xmlDataOutputRef->GetAttribute("dataOutputName");
                couplingLogicIn.timeStepLoop.dataOutputRefs.push_back(dataOutputName);
            }
        }
    } else if (xmlCouplingLogicIn->GetAttribute("type") == "iterativeCouplingLoop") {
        couplingLogicIn.type = EMPIRE_IterativeCouplingLoop;
        ticpp::Element *xmlIterativeCouplingLoop = xmlCouplingLogicIn->FirstChildElement(
                "iterativeCouplingLoop");
        { // add convergence checker
            ticpp::Element *xmlConvergenceChecker = xmlIterativeCouplingLoop->FirstChildElement(
                    "convergenceChecker");
            couplingLogicIn.iterativeCouplingLoop.convergenceChecker.maxNumOfIterations =
                    xmlConvergenceChecker->GetAttribute<double>("maxNumOfIterations");

            ticpp::Iterator<Element> xmlCheckResidual("checkResidual");
            for (xmlCheckResidual = xmlCheckResidual.begin(xmlConvergenceChecker);
                    xmlCheckResidual != xmlCheckResidual.end(); xmlCheckResidual++) {
                structCouplingLogic::structIterativeCouplingLoop::structConvergenceChecker::structCheckResidual checkResidual;
                checkResidual.relativeTolerance = xmlCheckResidual->GetAttribute<double>(
                        "relativeTolerance");
                checkResidual.absoluteTolerance = xmlCheckResidual->GetAttribute<double>(
                        "absoluteTolerance");
                checkResidual.residualRef.couplingAlgorithmName =
                        xmlCheckResidual->FirstChildElement("residualRef")->GetAttribute<string>(
                                "couplingAlgorithmName");
                checkResidual.residualRef.index =
                        xmlCheckResidual->FirstChildElement("residualRef")->GetAttribute<int>(
                                "index");
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.checkResiduals.push_back(
                        checkResidual);
            }
        }
        { // add Convergence Observer
            ticpp::Iterator<Element> xmlConvergenceObserver("convergenceObserver");
            for (xmlConvergenceObserver = xmlConvergenceObserver.begin(xmlIterativeCouplingLoop);
                    xmlConvergenceObserver != xmlConvergenceObserver.end();
                    xmlConvergenceObserver++) {
                string tmpString =
                        xmlConvergenceObserver->FirstChildElement("clientCodeRef")->GetAttribute(
                                "clientCodeName");
                couplingLogicIn.iterativeCouplingLoop.convergenceObservers.push_back(tmpString);
            }
        }

        { // add coupling algorithm ref
            ticpp::Element *xmlCouplingAlgorithmRef = xmlIterativeCouplingLoop->FirstChildElement(
                    "couplingAlgorithmRef", false);

            if (xmlCouplingAlgorithmRef != NULL) {
                couplingLogicIn.iterativeCouplingLoop.couplingAlgorithmRef.first = true;
                couplingLogicIn.iterativeCouplingLoop.couplingAlgorithmRef.second =
                        xmlCouplingAlgorithmRef->GetAttribute("couplingAlgorithmName");
            } else {
                couplingLogicIn.iterativeCouplingLoop.couplingAlgorithmRef.first = false;
            }
        }
        { // add dataOutputs
            ticpp::Iterator<Element> xmlDataOutputRef("dataOutputRef");
            for (xmlDataOutputRef = xmlDataOutputRef.begin(xmlIterativeCouplingLoop);
                    xmlDataOutputRef != xmlDataOutputRef.end(); xmlDataOutputRef++) {
                string dataOutputName = xmlDataOutputRef->GetAttribute("dataOutputName");
                couplingLogicIn.iterativeCouplingLoop.dataOutputRefs.push_back(dataOutputName);
            }
        }
        // check whether there is coupling algorithm ref when check residuals in the convergence checker
        if (couplingLogicIn.iterativeCouplingLoop.convergenceChecker.checkResiduals.size() > 0) {
            assert(couplingLogicIn.iterativeCouplingLoop.couplingAlgorithmRef.first);
        }
    } else if (xmlCouplingLogicIn->GetAttribute("type") == "connection") {
        couplingLogicIn.type = EMPIRE_connection;
        couplingLogicIn.connectionRef.connectionName = xmlCouplingLogicIn->FirstChildElement(
                "connectionRef")->GetAttribute("connectionName");
    } else {
        assert(false);
    }

    // 2. add child coupling logics
    if (xmlCouplingLogicIn->FirstChildElement("sequence", false) != NULL) {
        ticpp::Element *xmlCouplingLogicSequence = xmlCouplingLogicIn->FirstChildElement(
                "sequence");
        ticpp::Iterator<Element> xmlCouplingLogic("couplingLogic");
        for (xmlCouplingLogic = xmlCouplingLogic.begin(xmlCouplingLogicSequence);
                xmlCouplingLogic != xmlCouplingLogic.end(); xmlCouplingLogic++) {
            structCouplingLogic couplingLogic;
            parseCouplingLogicBlock(xmlCouplingLogic, couplingLogic);
            couplingLogicIn.sequence.push_back(couplingLogic);
        }
    }
}

void MetaDatabase::fillSettingCouplingLogic() {
    assert(settingGlobalCouplingLogic.sequence.size() == 0);
    settingGlobalCouplingLogic.type = EMPIRE_CouplingLogicSequence; // the type of globalCouplingLogic is fixed
    ticpp::Element *xmlGlobalCouplingLogic =
            inputFile->FirstChildElement("EMPEROR")->FirstChildElement("coSimulation");
    ticpp::Element *xmlCouplingLogicSequence = xmlGlobalCouplingLogic->FirstChildElement(
            "sequence");
    ticpp::Iterator<Element> xmlCouplingLogic("couplingLogic");
    for (xmlCouplingLogic = xmlCouplingLogic.begin(xmlCouplingLogicSequence);
            xmlCouplingLogic != xmlCouplingLogic.end(); xmlCouplingLogic++) {
        structCouplingLogic couplingLogic;
        parseCouplingLogicBlock(xmlCouplingLogic, couplingLogic);
        settingGlobalCouplingLogic.sequence.push_back(couplingLogic);
    }
}

structConnectionIO MetaDatabase::parseConnectionIORef(ticpp::Element *xmlElement) {
    ticpp::Element *xmlDataFieldRef = xmlElement->FirstChildElement("dataFieldRef", false);
    ticpp::Element *xmlSignalRef = xmlElement->FirstChildElement("signalRef", false);
    structConnectionIO settingConnectionIO;
    if (xmlDataFieldRef != NULL) {
        settingConnectionIO.type = EMPIRE_ConnectionIO_DataField;
        settingConnectionIO.dataFieldRef.clientCodeName = xmlDataFieldRef->GetAttribute(
                "clientCodeName");
        settingConnectionIO.dataFieldRef.meshName = xmlDataFieldRef->GetAttribute("meshName");
        settingConnectionIO.dataFieldRef.dataFieldName = xmlDataFieldRef->GetAttribute(
                "dataFieldName");
    } else if (xmlSignalRef != NULL) {
        settingConnectionIO.type = EMPIRE_ConnectionIO_Signal;
        settingConnectionIO.signalRef.clientCodeName = xmlSignalRef->GetAttribute("clientCodeName");
        settingConnectionIO.signalRef.signalName = xmlSignalRef->GetAttribute("signalName");
    } else {
        assert(false);
    }
    return settingConnectionIO;
}

std::vector<structConnectionIO> MetaDatabase::parseConnectionIORefs(ticpp::Element *xmlElement) {
    ticpp::Iterator<Element> xmlDataFieldRef("dataFieldRef");
    std::vector<structConnectionIO> settingConnectionIOs;
    for (xmlDataFieldRef = xmlDataFieldRef.begin(xmlElement);
            xmlDataFieldRef != xmlDataFieldRef.end(); xmlDataFieldRef++) {
        structConnectionIO io;
        io.type = EMPIRE_ConnectionIO_DataField;
        io.dataFieldRef.clientCodeName = xmlDataFieldRef->GetAttribute("clientCodeName");
        io.dataFieldRef.meshName = xmlDataFieldRef->GetAttribute("meshName");
        io.dataFieldRef.dataFieldName = xmlDataFieldRef->GetAttribute("dataFieldName");
        settingConnectionIOs.push_back(io);
    }
    ticpp::Iterator<Element> xmlSignalRef("signalRef");
    for (xmlSignalRef = xmlSignalRef.begin(xmlElement); xmlSignalRef != xmlSignalRef.end();
            xmlSignalRef++) {
        structConnectionIO io;
        io.type = EMPIRE_ConnectionIO_Signal;
        io.signalRef.clientCodeName = xmlSignalRef->GetAttribute("clientCodeName");
        io.signalRef.signalName = xmlSignalRef->GetAttribute("signalName");
        settingConnectionIOs.push_back(io);
    }
    return settingConnectionIOs;
}

} /* namespace EMPIRE */
