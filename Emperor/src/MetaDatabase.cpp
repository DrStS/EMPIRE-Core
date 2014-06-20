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
        clientCode.name = xmlClientCode->GetAttribute<string>("name");

        ticpp::Iterator<Element> xmlMesh("mesh");

        for (xmlMesh = xmlMesh.begin(xmlClientCode.Get()); xmlMesh != xmlMesh.end(); xmlMesh++) {
            structClientCode::structMesh mesh;
            mesh.name = xmlMesh->GetAttribute<string>("name");
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
            if (xmlMesh->HasAttribute("triangulateAll")) {
                string tmp = xmlMesh->GetAttribute<string>("triangulateAll");
                bool triangulateAll;
                if (tmp == "true") {
                    triangulateAll = true;
                } else if (tmp == "false") {
                    triangulateAll = false;
                } else {
                    assert(false);
                }
                mesh.triangulateAll = triangulateAll;
            } else {
                mesh.triangulateAll = false;
            }
            ticpp::Iterator<Element> xmlDataField("dataField");
            for (xmlDataField = xmlDataField.begin(xmlMesh.Get());
                    xmlDataField != xmlDataField.end(); xmlDataField++) {
                structClientCode::structMesh::structDataField dataField;
                assert(xmlDataField->HasAttribute("name"));
                assert(xmlDataField->HasAttribute("dimension"));
                assert(xmlDataField->HasAttribute("location"));
                assert(xmlDataField->HasAttribute("typeOfQuantity"));
                dataField.name = xmlDataField->GetAttribute<string>("name");
                if (xmlDataField->GetAttribute<string>("dimension") == "vector")
                    dataField.dimension = EMPIRE_DataField_vector;
                else if (xmlDataField->GetAttribute<string>("dimension") == "scalar")
                    dataField.dimension = EMPIRE_DataField_scalar;
                else
                    assert(false);
                if (xmlDataField->GetAttribute<string>("location") == "atNode")
                    dataField.location = EMPIRE_DataField_atNode;
                else if (xmlDataField->GetAttribute<string>("location") == "atElemCentroid")
                    dataField.location = EMPIRE_DataField_atElemCentroid;
                else
                    assert(false);
                if (xmlDataField->GetAttribute<string>("typeOfQuantity") == "field")
                    dataField.typeOfQuantity = EMPIRE_DataField_field;
                else if (xmlDataField->GetAttribute<string>("typeOfQuantity") == "fieldIntegral")
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
            signal.name = xmlSignal->GetAttribute<string>("name");
            string size = xmlSignal->GetAttribute<string>("size");
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
        dataOutput.name = xmlDataOutput->GetAttribute<string>("name");
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
        mapper.name = xmlMapper->GetAttribute<string>("name");
        ticpp::Element *xmlMeshRefA = xmlMapper->FirstChildElement("meshA")->FirstChildElement(
                "meshRef");
        mapper.meshRefA.clientCodeName = xmlMeshRefA->GetAttribute<string>("clientCodeName");
        mapper.meshRefA.meshName = xmlMeshRefA->GetAttribute<string>("meshName");
        ticpp::Element *xmlMeshRefB = xmlMapper->FirstChildElement("meshB")->FirstChildElement(
                "meshRef");
        mapper.meshRefB.clientCodeName = xmlMeshRefB->GetAttribute<string>("clientCodeName");
        mapper.meshRefB.meshName = xmlMeshRefB->GetAttribute<string>("meshName");
        if (xmlMapper->GetAttribute<string>("type") == "mortarMapper") {
            mapper.type = EMPIRE_MortarMapper;
            ticpp::Element *xmlMortarMapper = xmlMapper->FirstChildElement("mortarMapper");
            if (xmlMortarMapper->GetAttribute<string>("oppositeSurfaceNormal") == "true")
                mapper.mortarMapper.oppositeSurfaceNormal = true;
            else if (xmlMortarMapper->GetAttribute<string>("oppositeSurfaceNormal") == "false")
                mapper.mortarMapper.oppositeSurfaceNormal = false;
            else
                assert(false);
            if (xmlMortarMapper->GetAttribute<string>("dual") == "true")
                mapper.mortarMapper.dual = true;
            else if (xmlMortarMapper->GetAttribute<string>("dual") == "false")
                mapper.mortarMapper.dual = false;
            else
                assert(false);
            if (xmlMortarMapper->GetAttribute<string>("enforceConsistency") == "true")
                mapper.mortarMapper.enforceConsistency = true;
            else if (xmlMortarMapper->GetAttribute<string>("enforceConsistency") == "false")
                mapper.mortarMapper.enforceConsistency = false;
            else
                assert(false);
        } else if (xmlMapper->GetAttribute<string>("type") == "nearestNeighborMapper") {
            mapper.type = EMPIRE_NearestNeighborMapper;
        } else if (xmlMapper->GetAttribute<string>("type") == "barycentricInterpolationMapper") {
            mapper.type = EMPIRE_BarycentricInterpolationMapper;
        } else if (xmlMapper->GetAttribute<string>("type") == "IGAMortarMapper") {
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
        coupAlg.name = xmlCoupAlg->GetAttribute<string>("name");

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
        if (xmlCoupAlg->GetAttribute<string>("type") == "aitken") {
            coupAlg.type = EMPIRE_Aitken;
            ticpp::Element *xmlAitken = xmlCoupAlg->FirstChildElement("aitken");
            double tmpDouble = xmlAitken->GetAttribute<double>("initialRelaxationFactor");
            coupAlg.aitken.initialRelaxationFactor = tmpDouble;
        } else if (xmlCoupAlg->GetAttribute<string>("type") == "constantRelaxation") {
            coupAlg.type = EMPIRE_ConstantRelaxation;
            ticpp::Element *constantRelaxation = xmlCoupAlg->FirstChildElement(
                    "constantRelaxation");
            double tmpDouble = constantRelaxation->GetAttribute<double>("relaxationFactor");
            coupAlg.constantRelaxation.relaxationFactor = tmpDouble;
        } else if (xmlCoupAlg->GetAttribute<string>("type") == "IJCSA") {
            { // interfaceJacobian
                ticpp::Iterator<Element> xmlOutput("interfaceJacobian");
                for (xmlOutput = xmlOutput.begin(xmlCoupAlg.Get()); xmlOutput != xmlOutput.end();
                        xmlOutput++) {

                    structCouplingAlgorithm::structInterfaceJacobian interfaceJacobian;
                    interfaceJacobian.indexRow    = xmlOutput->GetAttribute<unsigned int>("indexRow");
                    interfaceJacobian.indexColumn = xmlOutput->GetAttribute<unsigned int>("indexColumn");

                    if(xmlOutput->FirstChildElement("constantValue", false) != NULL){
                        ticpp::Element *constantValue = xmlOutput->FirstChildElement("constantValue");
                        interfaceJacobian.value       = constantValue->GetAttribute<double>("value");
                        interfaceJacobian.isAutoDiff  = false;
                        interfaceJacobian.isSignal    = false;
                        interfaceJacobian.isConstant  = true;
                    }

                    if(xmlOutput->FirstChildElement("automaticDetermination", false) != NULL){

                    	interfaceJacobian.coefficient   =xmlOutput->FirstChildElement("automaticDetermination")->GetAttribute<double>("coefficient");
                        ticpp::Element *functionInput   =xmlOutput->FirstChildElement("automaticDetermination")->FirstChildElement("functionInput");
                        ticpp::Element *functionOutput  =xmlOutput->FirstChildElement("automaticDetermination")->FirstChildElement("functionOutput");
                        interfaceJacobian.functionInput =parseConnectionIORef(functionInput);
                        interfaceJacobian.functionOutput=parseConnectionIORef(functionOutput);
                        interfaceJacobian.isConstant  = false;
                        interfaceJacobian.isSignal    = false;
                        interfaceJacobian.isAutoDiff  = true;
                    }
                    coupAlg.interfaceJacobians.push_back(interfaceJacobian);
                }
            }
            coupAlg.type = EMPIRE_IJCSA;
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
        settingExtrapolator.name = xmlExtrapolator->GetAttribute<string>("name");
        if (xmlExtrapolator->GetAttribute<string>("type") == "linearExtrapolator") {
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
        connection.name = xmlConnection->GetAttribute<string>("name");
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
                if (xmlFilter->GetAttribute<string>("type") == "mappingFilter") {
                    filter.type = EMPIRE_MappingFilter;
                    filter.mappingFilter.mapperName =
                            xmlFilter->FirstChildElement("mappingFilter")->FirstChildElement(
                                    "mapperRef")->GetAttribute<string>("mapperName");
                } else if (xmlFilter->GetAttribute<string>("type") == "locationFilter") {
                    filter.type = EMPIRE_LocationFilter;
                } else if (xmlFilter->GetAttribute<string>("type") == "scalingFilter") {
                    filter.type = EMPIRE_ScalingFilter;
                    filter.scalingFilter.factor =
                            xmlFilter->FirstChildElement("scalingFilter")->GetAttribute<double>(
                                    "factor");
                } else if (xmlFilter->GetAttribute<string>("type") == "setFilter") {
                    filter.type = EMPIRE_SetFilter;
                    string valueString = xmlFilter->FirstChildElement("setFilter")->GetAttribute<
                            string>("value");
                    std::stringstream ss(valueString);
                    while (!ss.eof()) {
                        double value;
                        ss >> value;
                        filter.setFilter.value.push_back(value);
                    }
                } else if (xmlFilter->GetAttribute<string>("type") == "copyFilter") {
                    filter.type = EMPIRE_CopyFilter;
                    filter.copyFilter.signalOffset = 0;
                    ticpp::Element *pXMLElement = xmlFilter->FirstChildElement("copyFilter", false);
                    if (pXMLElement != NULL) {
                        filter.copyFilter.signalOffset = pXMLElement->GetAttribute<int>(
                                "signalOffset");
                    }

                } else if (xmlFilter->GetAttribute<string>("type")
                        == "dataFieldIntegrationFilter") {
                    filter.type = EMPIRE_DataFieldIntegrationFilter;
                    ticpp::Element *xmlMeshRef = xmlFilter->FirstChildElement(
                            "dataFieldIntegrationFilter")->FirstChildElement("meshRef");
                    filter.dataFieldIntegrationFilter.meshRef.clientCodeName =
                            xmlMeshRef->GetAttribute<string>("clientCodeName");
                    filter.dataFieldIntegrationFilter.meshRef.meshName = xmlMeshRef->GetAttribute<
                            string>("meshName");
                } else if (xmlFilter->GetAttribute<string>("type") == "additionFilter") {
                    filter.type = EMPIRE_AdditionFilter;
                    filter.additionFilter.a =
                            xmlFilter->FirstChildElement("additionFilter")->GetAttribute<double>(
                                    "a");
                    filter.additionFilter.b =
                            xmlFilter->FirstChildElement("additionFilter")->GetAttribute<double>(
                                    "b");
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
    if (xmlCouplingLogicIn->GetAttribute<string>("type") == "timeStepLoop") {
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
                        xmlExtrapolatorRef->GetAttribute<string>("extrapolatorName");
            } else {
                couplingLogicIn.timeStepLoop.extrapolatorRef.first = false;
            }
        }

        { // add dataOutputs
            ticpp::Iterator<Element> xmlDataOutputRef("dataOutputRef");
            for (xmlDataOutputRef = xmlDataOutputRef.begin(xmlTimeStepLoop);
                    xmlDataOutputRef != xmlDataOutputRef.end(); xmlDataOutputRef++) {
                string dataOutputName = xmlDataOutputRef->GetAttribute<string>("dataOutputName");
                couplingLogicIn.timeStepLoop.dataOutputRefs.push_back(dataOutputName);
            }
        }
    } else if (xmlCouplingLogicIn->GetAttribute<string>("type") == "iterativeCouplingLoop") {
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
                        xmlConvergenceObserver->FirstChildElement("clientCodeRef")->GetAttribute<
                                string>("clientCodeName");
                couplingLogicIn.iterativeCouplingLoop.convergenceObservers.push_back(tmpString);
            }
        }

        {   // add coupling algorithm refs
            ticpp::Iterator<Element> xmlCouplingAlgorithmRef("couplingAlgorithmRef");
            for (xmlCouplingAlgorithmRef = xmlCouplingAlgorithmRef.begin(xmlIterativeCouplingLoop);
            		xmlCouplingAlgorithmRef != xmlCouplingAlgorithmRef.end();
            		xmlCouplingAlgorithmRef++) {


                couplingLogicIn.iterativeCouplingLoop.couplingAlgorithmRefs.push_back(
                        xmlCouplingAlgorithmRef->GetAttribute<string>("couplingAlgorithmName"));
            }
        }
        { // add dataOutputs
            ticpp::Iterator<Element> xmlDataOutputRef("dataOutputRef");
            for (xmlDataOutputRef = xmlDataOutputRef.begin(xmlIterativeCouplingLoop);
                    xmlDataOutputRef != xmlDataOutputRef.end(); xmlDataOutputRef++) {
                string dataOutputName = xmlDataOutputRef->GetAttribute<string>("dataOutputName");
                couplingLogicIn.iterativeCouplingLoop.dataOutputRefs.push_back(dataOutputName);
            }
        }
        // check whether there is coupling algorithm ref when check residuals in the convergence checker
        if (couplingLogicIn.iterativeCouplingLoop.convergenceChecker.checkResiduals.size() > 0) {
            assert(couplingLogicIn.iterativeCouplingLoop.couplingAlgorithmRefs.size()>0);
        }
    } else if (xmlCouplingLogicIn->GetAttribute<string>("type") == "connection") {
        couplingLogicIn.type = EMPIRE_connection;
        couplingLogicIn.connectionRef.connectionName = xmlCouplingLogicIn->FirstChildElement(
                "connectionRef")->GetAttribute<string>("connectionName");
    } else if (xmlCouplingLogicIn->GetAttribute<string>("type") == "optimizationLoop") {
        couplingLogicIn.type = EMPIRE_OptimizationLoop;
        ticpp::Element *xmlOptimizationLoop = xmlCouplingLogicIn->FirstChildElement(
                "optimizationLoop");
        couplingLogicIn.optimizationLoop.maxNumOfIterations =
                xmlOptimizationLoop->GetAttribute<int>("maxNumOfIterations");
        { // add convergence signal sender
            ticpp::Element *xmlConvergenceSignalSender = xmlOptimizationLoop->FirstChildElement(
                    "convergenceSignalSender");
            string tmpString =
                    xmlConvergenceSignalSender->FirstChildElement("clientCodeRef")->GetAttribute<
                            string>("clientCodeName");
            couplingLogicIn.optimizationLoop.convergenceSignalSender = tmpString;
        }

        { // add convergence signal receivers
            ticpp::Iterator<Element> xmlConvergenceSignalReceiver("convergenceSignalReceiver");
            for (xmlConvergenceSignalReceiver = xmlConvergenceSignalReceiver.begin(
                    xmlOptimizationLoop);
                    xmlConvergenceSignalReceiver != xmlConvergenceSignalReceiver.end();
                    xmlConvergenceSignalReceiver++) {
                string tmpString =
                        xmlConvergenceSignalReceiver->FirstChildElement("clientCodeRef")->GetAttribute<
                                string>("clientCodeName");
                couplingLogicIn.optimizationLoop.convergenceSignalReceivers.push_back(tmpString);
            }

        }

        { // add dataOutputs
            ticpp::Iterator<Element> xmlDataOutputRef("dataOutputRef");
            for (xmlDataOutputRef = xmlDataOutputRef.begin(xmlOptimizationLoop);
                    xmlDataOutputRef != xmlDataOutputRef.end(); xmlDataOutputRef++) {
                string dataOutputName = xmlDataOutputRef->GetAttribute<string>("dataOutputName");
                couplingLogicIn.optimizationLoop.dataOutputRefs.push_back(dataOutputName);
            }
        }
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
        settingConnectionIO.dataFieldRef.clientCodeName = xmlDataFieldRef->GetAttribute<string>(
                "clientCodeName");
        settingConnectionIO.dataFieldRef.meshName = xmlDataFieldRef->GetAttribute<string>(
                "meshName");
        settingConnectionIO.dataFieldRef.dataFieldName = xmlDataFieldRef->GetAttribute<string>(
                "dataFieldName");
    } else if (xmlSignalRef != NULL) {
        settingConnectionIO.type = EMPIRE_ConnectionIO_Signal;
        settingConnectionIO.signalRef.clientCodeName = xmlSignalRef->GetAttribute<string>(
                "clientCodeName");
        settingConnectionIO.signalRef.signalName = xmlSignalRef->GetAttribute<string>("signalName");
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
        io.dataFieldRef.clientCodeName = xmlDataFieldRef->GetAttribute<string>("clientCodeName");
        io.dataFieldRef.meshName = xmlDataFieldRef->GetAttribute<string>("meshName");
        io.dataFieldRef.dataFieldName = xmlDataFieldRef->GetAttribute<string>("dataFieldName");
        settingConnectionIOs.push_back(io);
    }
    ticpp::Iterator<Element> xmlSignalRef("signalRef");
    for (xmlSignalRef = xmlSignalRef.begin(xmlElement); xmlSignalRef != xmlSignalRef.end();
            xmlSignalRef++) {
        structConnectionIO io;
        io.type = EMPIRE_ConnectionIO_Signal;
        io.signalRef.clientCodeName = xmlSignalRef->GetAttribute<string>("clientCodeName");
        io.signalRef.signalName = xmlSignalRef->GetAttribute<string>("signalName");
        settingConnectionIOs.push_back(io);
    }
    return settingConnectionIOs;
}

} /* namespace EMPIRE */
