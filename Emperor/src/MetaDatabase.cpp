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
        ticpp::Iterator<Element> xmlDataFieldRef("dataFieldRef");
        for (xmlDataFieldRef = xmlDataFieldRef.begin(xmlDataOutput.Get());
                xmlDataFieldRef != xmlDataFieldRef.end(); xmlDataFieldRef++) {
            structDataFieldRef dataFieldRef;
            dataFieldRef.clientCodeName = xmlDataFieldRef->GetAttribute("clientCodeName");
            dataFieldRef.meshName = xmlDataFieldRef->GetAttribute("meshName");
            dataFieldRef.dataFieldName = xmlDataFieldRef->GetAttribute("dataFieldName");
            dataOutput.dataFieldRefs.push_back(dataFieldRef);
        }
        ticpp::Iterator<Element> xmlSignalRef("signalRef");
        for (xmlSignalRef = xmlSignalRef.begin(xmlDataOutput.Get());
                xmlSignalRef != xmlSignalRef.end(); xmlSignalRef++) {
            structSignalRef signalRef;
            signalRef.clientCodeName = xmlSignalRef->GetAttribute("clientCodeName");
            signalRef.signalName = xmlSignalRef->GetAttribute("signalName");
            dataOutput.signalRefs.push_back(signalRef);
        }
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
        } else if (xmlMapper->GetAttribute("type") == "IGAMortarMapper"){
        	mapper.type = EMPIRE_IGAMortarMapper;
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
        if (xmlCoupAlg->GetAttribute("type") == "aitken") {
            coupAlg.type = EMPIRE_Aitken;
            ticpp::Element *xmlAitken = xmlCoupAlg->FirstChildElement("aitken");
            double tmpDouble = xmlAitken->GetAttribute<double>("initialAitkenFactor");
            coupAlg.aitken.initialAitkenFactor = tmpDouble;
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
        structExtrapolator extrapolator;
        extrapolator.name = xmlExtrapolator->GetAttribute("name");
        if (xmlExtrapolator->GetAttribute("type") == "simpleExtrapolator") {
            extrapolator.type = EMPIRE_SimpleExtrapolator;
        } else if (xmlExtrapolator->GetAttribute("type") == "genMSExtrapolator") {
            extrapolator.type = EMPIRE_GenMSExtrapolator;
            ticpp::Element *xmlGenMSExtrapolator = xmlExtrapolator->FirstChildElement(
                    "genMSExtrapolator");

            // number of inputs
            int numInput = xmlGenMSExtrapolator->GetAttribute<int>("numInput");
            assert(numInput >= 1 && numInput <= 3);
            extrapolator.genMSExtrapolator.numInput = numInput;

            // sequence length
            int seqLen = xmlGenMSExtrapolator->GetAttribute<int>("seqLen");
            assert(seqLen > 0);
            extrapolator.genMSExtrapolator.seqLen = seqLen;

            // sum old extrapolated data
            string sumOutput = xmlGenMSExtrapolator->GetAttribute("sumOutput");
            if (sumOutput == "true")
                extrapolator.genMSExtrapolator.sumOutput = true;
            else
                extrapolator.genMSExtrapolator.sumOutput = false;

            // time step
            double deltaTime = xmlGenMSExtrapolator->GetAttribute<double>("deltaTime");
            extrapolator.genMSExtrapolator.deltaTime = deltaTime;

            int index;
            // zeroth derivative coefficients
            vector<double>& coefficientDot0 = extrapolator.genMSExtrapolator.coefficientDot0;
            coefficientDot0.assign(seqLen, 0.);
            ticpp::Iterator<Element> xmlCoefficientDot0("coefficientDot0");
            index = 0;
            for (xmlCoefficientDot0 = xmlCoefficientDot0.begin(xmlGenMSExtrapolator);
                    xmlCoefficientDot0 != xmlCoefficientDot0.end(); xmlCoefficientDot0++) {
                coefficientDot0.at(index) = xmlCoefficientDot0->GetAttribute<double>("value");
                index++;
                if (index == seqLen) // ignore extra coefficients
                    break;
            }

            // first derivative coefficients
            if (numInput >= 2) {
                vector<double>& coefficientDot1 = extrapolator.genMSExtrapolator.coefficientDot1;
                coefficientDot1.assign(seqLen, 0.);
                ticpp::Iterator<Element> xmlCoefficientDot1("coefficientDot1");
                index = 0;
                for (xmlCoefficientDot1 = xmlCoefficientDot1.begin(xmlGenMSExtrapolator);
                        xmlCoefficientDot1 != xmlCoefficientDot1.end(); xmlCoefficientDot1++) {
                    coefficientDot1.at(index) = xmlCoefficientDot1->GetAttribute<double>("value");
                    index++;
                    if (index == seqLen) // ignore extra coefficients
                        break;
                }
                // second derivative coefficients
                if (numInput == 3) {
                    vector<double>& coefficientDot2 = extrapolator.genMSExtrapolator.coefficientDot2;
                    coefficientDot2.assign(seqLen, 0.);
                    ticpp::Iterator<Element> xmlCoefficientDot2("coefficientDot2");
                    index = 0;
                    for (xmlCoefficientDot2 = xmlCoefficientDot2.begin(xmlGenMSExtrapolator);
                            xmlCoefficientDot2 != xmlCoefficientDot2.end(); xmlCoefficientDot2++) {
                        coefficientDot2.at(index) = xmlCoefficientDot2->GetAttribute<double>(
                                "value");
                        index++;
                        if (index == seqLen) // ignore extra coefficients
                            break;
                    }
                }
            }

            // coefficients for summing old extrapolated data.
            if (sumOutput == "true") {
                vector<double>& coefficientOut = extrapolator.genMSExtrapolator.coefficientOut;
                coefficientOut.assign(seqLen, 0.);
                ticpp::Iterator<Element> xmlCoefficientOut("coefficientOut");
                index = 0;
                for (xmlCoefficientOut = xmlCoefficientOut.begin(xmlGenMSExtrapolator);
                        xmlCoefficientOut != xmlCoefficientOut.end(); xmlCoefficientOut++) {
                    coefficientOut.at(index) = xmlCoefficientOut->GetAttribute<double>("value");
                    index++;
                    if (index == seqLen) // ignore extra coefficients
                        break;
                }
            }
        } else {
            assert(false);
        }
        settingExtrapolatorVec.push_back(extrapolator);
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
            ticpp::Element *xmlDataFieldRef = xmlIO->FirstChildElement("dataFieldRef", false);
            ticpp::Element *xmlSignalRef = xmlIO->FirstChildElement("signalRef", false);
            if (xmlDataFieldRef != NULL) {
                structConnectionIO io;
                io.type = EMPIRE_ConnectionIO_DataField;
                io.dataFieldRef.clientCodeName = xmlDataFieldRef->GetAttribute("clientCodeName");
                io.dataFieldRef.meshName = xmlDataFieldRef->GetAttribute("meshName");
                io.dataFieldRef.dataFieldName = xmlDataFieldRef->GetAttribute("dataFieldName");
                connection.inputs.push_back(io);
                connection.outputs.push_back(io);
            } else if (xmlSignalRef != NULL) {
                structConnectionIO io;
                io.type = EMPIRE_ConnectionIO_Signal;
                io.signalRef.clientCodeName = xmlSignalRef->GetAttribute("clientCodeName");
                io.signalRef.signalName = xmlSignalRef->GetAttribute("signalName");
                connection.inputs.push_back(io);
                connection.outputs.push_back(io);
            } else {
                assert(false);
            }
        } else {
            ticpp::Iterator<Element> xmlInput("input");
            for (xmlInput = xmlInput.begin(xmlConnection.Get()); xmlInput != xmlInput.end();
                    xmlInput++) {
                ticpp::Element *xmlDataFieldRef = xmlInput->FirstChildElement("dataFieldRef",
                        false);
                ticpp::Element *xmlSignalRef = xmlInput->FirstChildElement("signalRef", false);
                if (xmlDataFieldRef != NULL) {
                    structConnectionIO input;
                    input.type = EMPIRE_ConnectionIO_DataField;
                    input.dataFieldRef.clientCodeName = xmlDataFieldRef->GetAttribute(
                            "clientCodeName");
                    input.dataFieldRef.meshName = xmlDataFieldRef->GetAttribute("meshName");
                    input.dataFieldRef.dataFieldName = xmlDataFieldRef->GetAttribute(
                            "dataFieldName");
                    connection.inputs.push_back(input);
                } else if (xmlSignalRef != NULL) {
                    structConnectionIO input;
                    input.type = EMPIRE_ConnectionIO_Signal;
                    input.signalRef.clientCodeName = xmlSignalRef->GetAttribute("clientCodeName");
                    input.signalRef.signalName = xmlSignalRef->GetAttribute("signalName");
                    connection.inputs.push_back(input);
                } else {
                    assert(false);
                }
            }
            ticpp::Iterator<Element> xmlOutput("output");
            for (xmlOutput = xmlOutput.begin(xmlConnection.Get()); xmlOutput != xmlOutput.end();
                    xmlOutput++) {
                ticpp::Element *xmlDataFieldRef = xmlOutput->FirstChildElement("dataFieldRef",
                        false);
                ticpp::Element *xmlSignalRef = xmlOutput->FirstChildElement("signalRef", false);
                if (xmlDataFieldRef != NULL) {
                    structConnectionIO output;
                    output.type = EMPIRE_ConnectionIO_DataField;
                    output.dataFieldRef.clientCodeName = xmlDataFieldRef->GetAttribute(
                            "clientCodeName");
                    output.dataFieldRef.meshName = xmlDataFieldRef->GetAttribute("meshName");
                    output.dataFieldRef.dataFieldName = xmlDataFieldRef->GetAttribute(
                            "dataFieldName");
                    connection.outputs.push_back(output);
                } else if (xmlSignalRef != NULL) {
                    structConnectionIO output;
                    output.type = EMPIRE_ConnectionIO_Signal;
                    output.signalRef.clientCodeName = xmlSignalRef->GetAttribute("clientCodeName");
                    output.signalRef.signalName = xmlSignalRef->GetAttribute("signalName");
                    connection.outputs.push_back(output);
                } else {
                    assert(false);
                }
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
                    ticpp::Element *xmlDataFieldRef = xmlIO->FirstChildElement("dataFieldRef",
                            false);
                    ticpp::Element *xmlSignalRef = xmlIO->FirstChildElement("signalRef", false);
                    if (xmlDataFieldRef != NULL) {
                        structConnectionIO io;
                        io.type = EMPIRE_ConnectionIO_DataField;
                        io.dataFieldRef.clientCodeName = xmlDataFieldRef->GetAttribute(
                                "clientCodeName");
                        io.dataFieldRef.meshName = xmlDataFieldRef->GetAttribute("meshName");
                        io.dataFieldRef.dataFieldName = xmlDataFieldRef->GetAttribute(
                                "dataFieldName");
                        filter.inputs.push_back(io);
                        filter.outputs.push_back(io);
                    } else if (xmlSignalRef != NULL) {
                        structConnectionIO io;
                        io.type = EMPIRE_ConnectionIO_Signal;
                        io.signalRef.clientCodeName = xmlSignalRef->GetAttribute("clientCodeName");
                        io.signalRef.signalName = xmlSignalRef->GetAttribute("signalName");
                        filter.inputs.push_back(io);
                        filter.outputs.push_back(io);
                    } else {
                        assert(false);
                    }
                } else {
                    ticpp::Iterator<Element> xmlInput("input");
                    for (xmlInput = xmlInput.begin(xmlFilter.Get()); xmlInput != xmlInput.end();
                            xmlInput++) {
                        ticpp::Element *xmlDataFieldRef = xmlInput->FirstChildElement(
                                "dataFieldRef", false);
                        ticpp::Element *xmlSignalRef = xmlInput->FirstChildElement("signalRef",
                                false);
                        if (xmlDataFieldRef != NULL) {
                            structConnectionIO input;
                            input.type = EMPIRE_ConnectionIO_DataField;
                            input.dataFieldRef.clientCodeName = xmlDataFieldRef->GetAttribute(
                                    "clientCodeName");
                            input.dataFieldRef.meshName = xmlDataFieldRef->GetAttribute("meshName");
                            input.dataFieldRef.dataFieldName = xmlDataFieldRef->GetAttribute(
                                    "dataFieldName");
                            filter.inputs.push_back(input);
                        } else if (xmlSignalRef != NULL) {
                            structConnectionIO input;
                            input.type = EMPIRE_ConnectionIO_Signal;
                            input.signalRef.clientCodeName = xmlSignalRef->GetAttribute(
                                    "clientCodeName");
                            input.signalRef.signalName = xmlSignalRef->GetAttribute("signalName");
                            filter.inputs.push_back(input);
                        } else {
                            assert(false);
                        }
                    }
                    ticpp::Iterator<Element> xmlOutput("output");
                    for (xmlOutput = xmlOutput.begin(xmlFilter.Get()); xmlOutput != xmlOutput.end();
                            xmlOutput++) {
                        ticpp::Element *xmlDataFieldRef = xmlOutput->FirstChildElement(
                                "dataFieldRef", false);
                        ticpp::Element *xmlSignalRef = xmlOutput->FirstChildElement("signalRef",
                                false);
                        if (xmlDataFieldRef != NULL) {
                            structConnectionIO output;
                            output.type = EMPIRE_ConnectionIO_DataField;
                            output.dataFieldRef.clientCodeName = xmlDataFieldRef->GetAttribute(
                                    "clientCodeName");
                            output.dataFieldRef.meshName = xmlDataFieldRef->GetAttribute(
                                    "meshName");
                            output.dataFieldRef.dataFieldName = xmlDataFieldRef->GetAttribute(
                                    "dataFieldName");
                            filter.outputs.push_back(output);
                        } else if (xmlSignalRef != NULL) {
                            structConnectionIO output;
                            output.type = EMPIRE_ConnectionIO_Signal;
                            output.signalRef.clientCodeName = xmlSignalRef->GetAttribute(
                                    "clientCodeName");
                            output.signalRef.signalName = xmlSignalRef->GetAttribute("signalName");
                            filter.outputs.push_back(output);
                        } else {
                            assert(false);
                        }
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
                } else if (xmlFilter->GetAttribute("type") == "couplingAlgorithmFilter") {
                    filter.type = EMPIRE_CouplingAlgorithmFilter;
                    filter.couplingAlgorithmFilter.couplingAlgorithmName =
                            xmlFilter->FirstChildElement("couplingAlgorithmFilter")->FirstChildElement(
                                    "couplingAlgorithmRef")->GetAttribute("couplingAlgorithmName");
                } else if (xmlFilter->GetAttribute("type") == "extrapolatingFilter") {
                    filter.type = EMPIRE_ExtrapolatingFilter;
                    filter.extrapolatingFilter.extrapolatorName =
                            xmlFilter->FirstChildElement("extrapolatingFilter")->FirstChildElement(
                                    "extrapolatorRef")->GetAttribute("extrapolatorName");
                } else if (xmlFilter->GetAttribute("type") == "scalingFilter") {
                    filter.type = EMPIRE_ScalingFilter;
                    filter.scalingFilter.factor =
                            xmlFilter->FirstChildElement("scalingFilter")->GetAttribute<double>(
                                    "factor");
                } else if (xmlFilter->GetAttribute("type") == "copyFilter") {
                    filter.type = EMPIRE_CopyFilter;
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

        { // add extrapolators
            ticpp::Iterator<Element> xmlExtrapolatorRef("extrapolatorRef");
            for (xmlExtrapolatorRef = xmlExtrapolatorRef.begin(xmlTimeStepLoop);
                    xmlExtrapolatorRef != xmlExtrapolatorRef.end(); xmlExtrapolatorRef++) {
                string tmpString = xmlExtrapolatorRef->GetAttribute("extrapolatorName");
                couplingLogicIn.timeStepLoop.extrapolatorRefs.push_back(tmpString);
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
            ticpp::Element *xmlDataFieldRef = xmlConvergenceChecker->FirstChildElement(
                    "dataFieldRef", false);
            ticpp::Element *xmlSignalRef = xmlConvergenceChecker->FirstChildElement("signalRef",
                    false);
            ticpp::Element *xmlCouplingAlgRef = xmlConvergenceChecker->FirstChildElement(
                    "couplingAlgorithmRef", false);
            if (xmlDataFieldRef != NULL && xmlSignalRef == NULL && xmlCouplingAlgRef == NULL) {
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.whichRef =
                        EMPIRE_ConvergenceChecker_dataFieldRef;
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.dataFieldRef.clientCodeName =
                        xmlDataFieldRef->GetAttribute("clientCodeName");
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.dataFieldRef.meshName =
                        xmlDataFieldRef->GetAttribute("meshName");
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.dataFieldRef.dataFieldName =
                        xmlDataFieldRef->GetAttribute("dataFieldName");
            } else if (xmlSignalRef != NULL && xmlDataFieldRef == NULL && xmlCouplingAlgRef == NULL) {
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.whichRef =
                        EMPIRE_ConvergenceChecker_signalRef;
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.signalRef.clientCodeName =
                        xmlSignalRef->GetAttribute("clientCodeName");
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.signalRef.signalName =
                        xmlSignalRef->GetAttribute("signalName");
            } else if (xmlCouplingAlgRef != NULL && xmlSignalRef == NULL && xmlDataFieldRef == NULL) {
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.whichRef =
                        EMPIRE_ConvergenceChecker_couplingAlgorithmRef;
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.couplingAlgorithmRef =
                        xmlCouplingAlgRef->GetAttribute("couplingAlgorithmName");
            } else {
                assert(false);
            }
            bool hasAtleastOne = false;
            if (xmlConvergenceChecker->HasAttribute("absoluteTolerance")) {
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.absoluteTolerance =
                        xmlConvergenceChecker->GetAttribute<double>("absoluteTolerance");
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.hasAbsTol = true;
                hasAtleastOne = hasAtleastOne
                        || couplingLogicIn.iterativeCouplingLoop.convergenceChecker.hasAbsTol;
            }
            if (xmlConvergenceChecker->HasAttribute("relativeTolerance")) {
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.relativeTolerance =
                        xmlConvergenceChecker->GetAttribute<double>("relativeTolerance");
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.hasRelTol = true;
                hasAtleastOne = hasAtleastOne
                        || couplingLogicIn.iterativeCouplingLoop.convergenceChecker.hasRelTol;
            }
            if (xmlConvergenceChecker->HasAttribute("maxNumOfIterations")) {
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.maxNumOfIterations =
                        xmlConvergenceChecker->GetAttribute<double>("maxNumOfIterations");
                couplingLogicIn.iterativeCouplingLoop.convergenceChecker.hasMaxNumOfIters = true;
                hasAtleastOne =
                        hasAtleastOne
                                || couplingLogicIn.iterativeCouplingLoop.convergenceChecker.hasMaxNumOfIters;
            }
            assert(hasAtleastOne);
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

        { // add coupling algorithm refs
            ticpp::Iterator<Element> xmlCouplingAlgorithmRef("couplingAlgorithmRef");
            for (xmlCouplingAlgorithmRef = xmlCouplingAlgorithmRef.begin(xmlIterativeCouplingLoop);
                    xmlCouplingAlgorithmRef != xmlCouplingAlgorithmRef.end();
                    xmlCouplingAlgorithmRef++) {
                string tmpString = xmlCouplingAlgorithmRef->GetAttribute("couplingAlgorithmName");
                couplingLogicIn.iterativeCouplingLoop.couplingAlgorithmRefs.push_back(tmpString);
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

} /* namespace EMPIRE */
