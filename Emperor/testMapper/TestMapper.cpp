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
#include "TestMapper.h"

namespace EMPIRE {

TestMapper::TestMapper() {
}

TestMapper::~TestMapper() {
    for (int i = 0; i < mappers.size(); i++)
        delete mappers[i];
}

void TestMapper::parseInputFile(char *fileName) {
    try {
        // 1. parse the input file
        Document *inputFile = new Document(fileName);
        inputFile->LoadFile();
        /// Fill up data base
        assert(settingMappers.size() == 0);
        ticpp::Element *xmlTestMapper = inputFile->FirstChildElement("testMapper");
        // mappers
        ticpp::Iterator<Element> xmlMapper("mapper");
        for (xmlMapper = xmlMapper.begin(xmlTestMapper); xmlMapper != xmlMapper.end();
                xmlMapper++) {
            // GetAttribute("name"); // throw no exception
            // GetAttribute<string>("name"); // throw exception
            // set up general stuff
            StructMapper mapper;
            mapper.name = xmlMapper->GetAttribute<string>("name");
            ticpp::Element *xmlMeshA = xmlMapper->FirstChildElement("meshA");
            if (xmlMeshA->HasAttribute("GiDMesh")) {
                mapper.meshAFormat = "GiDMesh";
                mapper.meshFileA = xmlMeshA->GetAttribute<string>("GiDMesh");
            } else {
                assert(false);
            }

            ticpp::Element *xmlMeshB = xmlMapper->FirstChildElement("meshB");
            if (xmlMeshB->HasAttribute("GiDMesh")) {
                mapper.meshBFormat = "GiDMesh";
                mapper.meshFileB = xmlMeshB->GetAttribute<string>("GiDMesh");
            } else {
                assert(false);
            }
            // set up mapper
            if (xmlMapper->GetAttribute("type") == "mortarMapper") {
                mapper.type = "mortarMapper";
                ticpp::Element *xmlMortarMapper = xmlMapper->FirstChildElement("mortarMapper");
                if (xmlMortarMapper->GetAttribute<string>("oppositeSurfaceNormal") == "true")
                    mapper.settingMortarMapper.oppositeSurfaceNormal = true;
                else if (xmlMortarMapper->GetAttribute<string>("oppositeSurfaceNormal") == "false")
                    mapper.settingMortarMapper.oppositeSurfaceNormal = false;
                else
                    assert(false);
                if (xmlMortarMapper->GetAttribute<string>("dual") == "true")
                    mapper.settingMortarMapper.dual = true;
                else if (xmlMortarMapper->GetAttribute<string>("dual") == "false")
                    mapper.settingMortarMapper.dual = false;
                else
                    assert(false);
                if (xmlMortarMapper->GetAttribute<string>("enforceConsistency") == "true")
                    mapper.settingMortarMapper.enforceConsistency = true;
                else if (xmlMortarMapper->GetAttribute<string>("enforceConsistency") == "false")
                    mapper.settingMortarMapper.enforceConsistency = false;
                else
                    assert(false);
            } else if (xmlMapper->GetAttribute("type") == "nearestNeighborMapper") {
                mapper.type = "nearestNeighborMapper";
            } else if (xmlMapper->GetAttribute("type") == "barycentricInterpolationMapper") {
                mapper.type = "barycentricInterpolationMapper";
            } else if (xmlMapper->GetAttribute("type") == "nearestElementMapper") {
                mapper.type = "nearestElementMapper";
            } else {
                assert(false);
            }
            // set consistent mapping
            ticpp::Element *xmlConsistentMapping = xmlMapper->FirstChildElement("consistentMapping",
                    false);
            if (xmlConsistentMapping != NULL) {
                mapper.settingConsistentMapping.doConsistentMapping = true;
                mapper.settingConsistentMapping.dataFieldATypeOfQuantity =
                        xmlConsistentMapping->GetAttribute<string>("dataFieldATypeOfQuantity");
                mapper.settingConsistentMapping.dataFieldBTypeOfQuantity =
                        xmlConsistentMapping->GetAttribute<string>("dataFieldBTypeOfQuantity");

                parseXMLDataField(xmlConsistentMapping->FirstChildElement("dataFieldA"),
                        mapper.settingConsistentMapping.dataFieldA);

                ticpp::Element *xmlDataFieldBRef = xmlConsistentMapping->FirstChildElement(
                        "dataFieldBRef", false);
                if (xmlDataFieldBRef != NULL) {
                    mapper.settingConsistentMapping.doErrorCalculation = true;
                    parseXMLDataField(xmlDataFieldBRef,
                            mapper.settingConsistentMapping.dataFieldBRef);
                } else {
                    mapper.settingConsistentMapping.doErrorCalculation = false;
                }
            } else {
                mapper.settingConsistentMapping.doConsistentMapping = false;
            }
            // set conservative mapping
            ticpp::Element *xmlConservativeMapping = xmlMapper->FirstChildElement(
                    "conservativeMapping", false);
            if (xmlConservativeMapping != NULL) {
                mapper.settingConservativeMapping.doConservativeMapping = true;
                mapper.settingConservativeMapping.dataFieldATypeOfQuantity =
                        xmlConservativeMapping->GetAttribute<string>("dataFieldATypeOfQuantity");
                mapper.settingConservativeMapping.dataFieldBTypeOfQuantity =
                        xmlConservativeMapping->GetAttribute<string>("dataFieldBTypeOfQuantity");

                parseXMLDataField(xmlConservativeMapping->FirstChildElement("dataFieldB"),
                        mapper.settingConservativeMapping.dataFieldB);

                ticpp::Element *xmlDataFieldARef = xmlConservativeMapping->FirstChildElement(
                        "dataFieldARef", false);
                if (xmlDataFieldARef != NULL) {
                    mapper.settingConservativeMapping.doErrorCalculation = true;
                    parseXMLDataField(xmlDataFieldARef,
                            mapper.settingConservativeMapping.dataFieldARef);
                } else {
                    mapper.settingConservativeMapping.doErrorCalculation = false;
                }
            } else {
                mapper.settingConservativeMapping.doConservativeMapping = false;
            }
            settingMappers.push_back(mapper);
        }

        // conservation analysises
        ticpp::Iterator<Element> xmlConservationAnalysis("conservationAnalysis");
        for (xmlConservationAnalysis = xmlConservationAnalysis.begin(xmlTestMapper);
                xmlConservationAnalysis != xmlConservationAnalysis.end();
                xmlConservationAnalysis++) {
            StructConservationAnalysis conservationAnalysis;
            conservationAnalysis.name = xmlConservationAnalysis->GetAttribute<string>("name");
            conservationAnalysis.displacementsATypeOfQuantity =
                    xmlConservationAnalysis->GetAttribute<string>("displacementsATypeOfQuantity");
            conservationAnalysis.forcesATypeOfQuantity = xmlConservationAnalysis->GetAttribute<
                    string>("forcesATypeOfQuantity");
            conservationAnalysis.displacementsBTypeOfQuantity =
                    xmlConservationAnalysis->GetAttribute<string>("displacementsBTypeOfQuantity");
            conservationAnalysis.forcesBTypeOfQuantity = xmlConservationAnalysis->GetAttribute<
                    string>("forcesBTypeOfQuantity");

            ticpp::Element *xmlMeshA = xmlConservationAnalysis->FirstChildElement("meshA");
            if (xmlMeshA->HasAttribute("GiDMesh")) {
                conservationAnalysis.meshAFormat = "GiDMesh";
                conservationAnalysis.meshFileA = xmlMeshA->GetAttribute<string>("GiDMesh");
            } else {
                assert(false);
            }

            ticpp::Element *xmlMeshB = xmlConservationAnalysis->FirstChildElement("meshB");
            if (xmlMeshB->HasAttribute("GiDMesh")) {
                conservationAnalysis.meshBFormat = "GiDMesh";
                conservationAnalysis.meshFileB = xmlMeshB->GetAttribute<string>("GiDMesh");
            } else {
                assert(false);
            }

            parseXMLDataField(xmlConservationAnalysis->FirstChildElement("displacementsA"),
                    conservationAnalysis.displacementsA);
            parseXMLDataField(xmlConservationAnalysis->FirstChildElement("forcesA"),
                    conservationAnalysis.forcesA);
            parseXMLDataField(xmlConservationAnalysis->FirstChildElement("displacementsB"),
                    conservationAnalysis.displacementsB);
            parseXMLDataField(xmlConservationAnalysis->FirstChildElement("forcesB"),
                    conservationAnalysis.forcesB);

            settingConservationAnalysises.push_back(conservationAnalysis);
        }

        delete inputFile;
    } catch (ticpp::Exception& ex) {
        cout << "ERROR Parser: " << ex.what() << endl;
        exit(EXIT_FAILURE);
    }
}

void TestMapper::outputSettingToShell() {
    cout << "testMapper {" << endl;
    // output setting for mappers
    for (int i = 0; i < settingMappers.size(); i++) {
        cout << "mapper {" << endl;
        cout << "\t" << "name: " << settingMappers[i].name << endl;
        cout << "\t" << "type: " << settingMappers[i].type << endl;
        cout << "\t" << "meshAFormat: " << settingMappers[i].meshAFormat << endl;
        cout << "\t" << "meshFileA: " << settingMappers[i].meshFileA << endl;
        cout << "\t" << "meshBFormat: " << settingMappers[i].meshBFormat << endl;
        cout << "\t" << "meshFileB: " << settingMappers[i].meshFileB << endl;
        if (settingMappers[i].type == "mortarMapper") {
            cout << "\t" << "mortarMapper {" << endl;
            cout << "\t" << "\t" << "oppositeSurfaceNormal: "
                    << settingMappers[i].settingMortarMapper.oppositeSurfaceNormal << endl;
            cout << "\t" << "\t" << "dual: " << settingMappers[i].settingMortarMapper.dual << endl;
            cout << "\t" << "\t" << "enforceConsistency: "
                    << settingMappers[i].settingMortarMapper.enforceConsistency << endl;
            cout << "\t" << "}" << endl;
        }
        if (settingMappers[i].settingConsistentMapping.doConsistentMapping) {
            cout << "\t" << "consistentMapping {" << endl;
            cout << "\t" << "\t" << "dataFieldATypeOfQuantity: "
                    << settingMappers[i].settingConsistentMapping.dataFieldATypeOfQuantity << endl;
            cout << "\t" << "\t" << "dataFieldBTypeOfQuantity: "
                    << settingMappers[i].settingConsistentMapping.dataFieldBTypeOfQuantity << endl;
            outputXMLDataFieldToShell(settingMappers[i].settingConsistentMapping.dataFieldA,
                    "dataFieldA", "\t\t");

            if (settingMappers[i].settingConsistentMapping.doErrorCalculation) {
                outputXMLDataFieldToShell(settingMappers[i].settingConsistentMapping.dataFieldBRef,
                        "dataFieldBRef", "\t\t");
            }
            cout << "\t" << "}" << endl;
        }
        if (settingMappers[i].settingConservativeMapping.doConservativeMapping) {
            cout << "\t" << "conservativeMapping {" << endl;
            cout << "\t" << "\t" << "dataFieldATypeOfQuantity: "
                    << settingMappers[i].settingConservativeMapping.dataFieldATypeOfQuantity
                    << endl;
            cout << "\t" << "\t" << "dataFieldBTypeOfQuantity: "
                    << settingMappers[i].settingConservativeMapping.dataFieldBTypeOfQuantity
                    << endl;
            outputXMLDataFieldToShell(settingMappers[i].settingConservativeMapping.dataFieldB,
                    "dataFieldB", "\t\t");

            if (settingMappers[i].settingConservativeMapping.doErrorCalculation) {
                outputXMLDataFieldToShell(
                        settingMappers[i].settingConservativeMapping.dataFieldARef, "dataFieldARef",
                        "\t\t");
            }
            cout << "\t" << "}" << endl;
        }
        cout << "}" << endl;
    }

    // output setting for conservation analysises
    for (int i = 0; i < settingConservationAnalysises.size(); i++) {
        cout << "conservationAnalysis {" << endl;
        cout << "\t" << "name: " << settingConservationAnalysises[i].name << endl;
        cout << "\t" << "displacementsATypeOfQuantity: "
                << settingConservationAnalysises[i].displacementsATypeOfQuantity << endl;
        cout << "\t" << "forcesATypeOfQuantity: "
                << settingConservationAnalysises[i].forcesATypeOfQuantity << endl;
        cout << "\t" << "displacementsBTypeOfQuantity: "
                << settingConservationAnalysises[i].displacementsBTypeOfQuantity << endl;
        cout << "\t" << "forcesBTypeOfQuantity: "
                << settingConservationAnalysises[i].forcesBTypeOfQuantity << endl;
        cout << "\t" << "meshAFormat: " << settingConservationAnalysises[i].meshAFormat << endl;
        cout << "\t" << "meshFileA: " << settingConservationAnalysises[i].meshFileA << endl;
        cout << "\t" << "meshBFormat: " << settingConservationAnalysises[i].meshBFormat << endl;
        cout << "\t" << "meshFileB: " << settingConservationAnalysises[i].meshFileB << endl;

        outputXMLDataFieldToShell(settingConservationAnalysises[i].displacementsA, "displacementsA",
                "\t");
        outputXMLDataFieldToShell(settingConservationAnalysises[i].forcesA, "forcesA", "\t");
        outputXMLDataFieldToShell(settingConservationAnalysises[i].displacementsB, "displacementsB",
                "\t");
        outputXMLDataFieldToShell(settingConservationAnalysises[i].forcesB, "forcesB", "\t");
        cout << "}" << endl;
    }

    cout << "}" << endl;
}

void TestMapper::doMapping() {
    for (int i = 0; i < settingMappers.size(); i++) {
        cout << "================================================" << endl;
        cout << "Do mapping: " << settingMappers[i].name << endl;

        // *. read mesh
        int numNodesA;
        int numElemsA;
        int* numNodesPerElemA;
        double *nodeCoorsA;
        int *nodeIDsA;
        int *elemTableA;
        int *elemIDsA;
        if (settingMappers[i].meshAFormat == "GiDMesh") {
            GiDFileIO::readDotMsh(settingMappers[i].meshFileA, numNodesA, numElemsA, nodeCoorsA,
                    nodeIDsA, numNodesPerElemA, elemTableA, elemIDsA);
        } else {
            assert(false);
        }

        DataFieldIntegration *DFI_A = new DataFieldIntegration(numNodesA, numElemsA,
                numNodesPerElemA, nodeCoorsA, nodeIDsA, elemTableA);

        int numNodesB;
        int numElemsB;
        int* numNodesPerElemB;
        double *nodeCoorsB;
        int *nodeIDsB;
        int *elemTableB;
        int *elemIDsB;
        if (settingMappers[i].meshBFormat == "GiDMesh") {
            GiDFileIO::readDotMsh(settingMappers[i].meshFileB, numNodesB, numElemsB, nodeCoorsB,
                    nodeIDsB, numNodesPerElemB, elemTableB, elemIDsB);
        } else {
            assert(false);
        }

        DataFieldIntegration *DFI_B = new DataFieldIntegration(numNodesB, numElemsB,
                numNodesPerElemB, nodeCoorsB, nodeIDsB, elemTableB);

        // *. write mesh
        string meshFileA;
        meshFileA = settingMappers[i].name;
        meshFileA.append("_meshA.msh");
        GiDFileIO::writeDotMsh(meshFileA, numNodesA, numElemsA, nodeCoorsA, nodeIDsA,
                numNodesPerElemA, elemTableA, elemIDsA);
        string meshFileB;
        meshFileB = settingMappers[i].name;
        meshFileB.append("_meshB.msh");
        GiDFileIO::writeDotMsh(meshFileB, numNodesB, numElemsB, nodeCoorsB, nodeIDsB,
                numNodesPerElemB, elemTableB, elemIDsB);

        // *. initialize resultA
        string resultFileA;
        resultFileA = settingMappers[i].name;
        resultFileA.append("_meshA.res");
        GiDFileIO::initDotRes(resultFileA);
        string resultFileB;
        resultFileB = settingMappers[i].name;
        resultFileB.append("_meshB.res");
        GiDFileIO::initDotRes(resultFileB);

        // *. construct mapper
        AbstractMapper *mapper;
        if (settingMappers[i].type == "mortarMapper") {
            bool oppositeSurfaceNormal = settingMappers[i].settingMortarMapper.oppositeSurfaceNormal;
            bool dual = settingMappers[i].settingMortarMapper.dual;
            bool enforceConsistency = settingMappers[i].settingMortarMapper.enforceConsistency;
            mapper = new MortarMapper(numNodesA, numElemsA, numNodesPerElemA, nodeCoorsA, nodeIDsA,
                    elemTableA, numNodesB, numElemsB, numNodesPerElemB, nodeCoorsB, nodeIDsB,
                    elemTableB, oppositeSurfaceNormal, dual, enforceConsistency);
        } else if (settingMappers[i].type == "nearestNeighborMapper") {
            mapper = new NearestNeighborMapper(numNodesA, nodeCoorsA, numNodesB, nodeCoorsB);
        } else if (settingMappers[i].type == "barycentricInterpolationMapper") {
            mapper = new BarycentricInterpolationMapper(numNodesA, nodeCoorsA, numNodesB,
                    nodeCoorsB);
        } else if (settingMappers[i].type == "nearestElementMapper") {
            mapper = new NearestElementMapper(numNodesA, numElemsA, numNodesPerElemA, nodeCoorsA,
                    nodeIDsA, elemTableA, numNodesB, numElemsB, numNodesPerElemB, nodeCoorsB,
                    nodeIDsB, elemTableB);
        } else {
            assert(false);
        }

        // *. consistent mapping
        if (settingMappers[i].settingConsistentMapping.doConsistentMapping) {
            cout << '\t' << "do consistent mapping ... " << endl;
            double *dataFieldA = new double[numNodesA * 3];
            double *dataFieldB = new double[numNodesB * 3];

            initDataField(settingMappers[i].settingConsistentMapping.dataFieldA, numNodesA,
                    numElemsA, nodeCoorsA, nodeIDsA, numNodesPerElemA, elemTableA, elemIDsA,
                    dataFieldA);

            for (int j = 0; j < 3; j++) { // x,y,z
                double *dataFieldA_j = new double[numNodesA];
                double *dataFieldB_j = new double[numNodesB];
                for (int k = 0; k < numNodesA; k++)
                    dataFieldA_j[k] = dataFieldA[k * 3 + j];

                if (settingMappers[i].settingConsistentMapping.dataFieldATypeOfQuantity
                        == "fieldIntegral") { // deIntegrate dataFieldA_j
                    double *tmpDF = new double[numNodesA];
                    DFI_A->deIntegrate(dataFieldA_j, tmpDF);
                    for (int k = 0; k < numNodesA; k++)
                        dataFieldA_j[k] = tmpDF[k];
                    delete[] tmpDF;
                } else {
                    assert(
                            settingMappers[i].settingConsistentMapping.dataFieldATypeOfQuantity == "field");
                }
                mapper->consistentMapping(dataFieldA_j, dataFieldB_j);

                if (settingMappers[i].settingConsistentMapping.dataFieldBTypeOfQuantity
                        == "fieldIntegral") { // integrate dataFieldB_j
                    double *tmpDF = new double[numNodesB];
                    DFI_B->integrate(dataFieldB_j, tmpDF);
                    for (int k = 0; k < numNodesB; k++)
                        dataFieldB_j[k] = tmpDF[k];
                    delete[] tmpDF;
                } else {
                    assert(
                            settingMappers[i].settingConsistentMapping.dataFieldBTypeOfQuantity == "field");
                }

                for (int k = 0; k < numNodesB; k++)
                    dataFieldB[k * 3 + j] = dataFieldB_j[k];
                delete[] dataFieldA_j;
                delete[] dataFieldB_j;
            }

            GiDFileIO::appendNodalDataToDotRes(resultFileA, "\"consistentMappingDataFieldA\"",
                    "\"testMapper\"", 1, "Vector", numNodesA, nodeIDsA, dataFieldA);
            GiDFileIO::appendNodalDataToDotRes(resultFileB, "\"consistentMappingDataFieldB\"",
                    "\"testMapper\"", 1, "Vector", numNodesB, nodeIDsB, dataFieldB);

            if (settingMappers[i].settingConsistentMapping.doErrorCalculation) {
                double *dataFieldBRef = new double[numNodesB * 3];
                initDataField(settingMappers[i].settingConsistentMapping.dataFieldBRef, numNodesB,
                        numElemsB, nodeCoorsB, nodeIDsB, numNodesPerElemB, elemTableB, elemIDsB,
                        dataFieldBRef);
                double *dataFieldBErr = new double[numNodesB * 3];
                for (int k = 0; k < numNodesB * 3; k++) {
                    dataFieldBErr[k] = dataFieldB[k] - dataFieldBRef[k];
                    dataFieldBErr[k] = fabs(dataFieldBErr[k]);
                }
                GiDFileIO::appendNodalDataToDotRes(resultFileB,
                        "\"consistentMappingDataFieldBRef\"", "\"testMapper\"", 1, "Vector",
                        numNodesB, nodeIDsB, dataFieldBRef);
                GiDFileIO::appendNodalDataToDotRes(resultFileB,
                        "\"consistentMappingDataFieldBErr\"", "\"testMapper\"", 1, "Vector",
                        numNodesB, nodeIDsB, dataFieldBErr);

                // comoute L2 norm of error
                double vecL2Norm[3] = { 0.0, 0.0, 0.0 };
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < numNodesB; k++) {
                        vecL2Norm[j] += dataFieldBErr[k * 3 + j] * dataFieldBErr[k * 3 + j];
                    }
                }
                for (int j = 0; j < 3; j++) {
                    vecL2Norm[j] /= (double) numNodesB;
                    vecL2Norm[j] = sqrt(vecL2Norm[j]);
                }
                cout << '\t' << '\t' << "Error L2 norm: " << "(" << vecL2Norm[0] << ", "
                        << vecL2Norm[1] << ", " << vecL2Norm[2] << ")" << endl;

                delete[] dataFieldBRef;
                delete[] dataFieldBErr;
            }

            delete[] dataFieldA;
            delete[] dataFieldB;
        }

        // *. conservative mapping
        if (settingMappers[i].settingConservativeMapping.doConservativeMapping) {
            cout << '\t' << "do conservative mapping ... " << endl;
            double *dataFieldA = new double[numNodesA * 3];
            double *dataFieldB = new double[numNodesB * 3];
            initDataField(settingMappers[i].settingConservativeMapping.dataFieldB, numNodesB,
                    numElemsB, nodeCoorsB, nodeIDsB, numNodesPerElemB, elemTableB, elemIDsB,
                    dataFieldB);

            for (int j = 0; j < 3; j++) { // x,y,z
                double *dataFieldA_j = new double[numNodesA];
                double *dataFieldB_j = new double[numNodesB];
                for (int k = 0; k < numNodesB; k++)
                    dataFieldB_j[k] = dataFieldB[k * 3 + j];

                if (settingMappers[i].settingConservativeMapping.dataFieldBTypeOfQuantity
                        == "field") { // integrate dataFieldB_j
                    double *tmpDF = new double[numNodesB];
                    DFI_B->integrate(dataFieldB_j, tmpDF);
                    for (int k = 0; k < numNodesB; k++)
                        dataFieldB_j[k] = tmpDF[k];
                    delete[] tmpDF;
                } else {
                    assert(
                            settingMappers[i].settingConservativeMapping.dataFieldBTypeOfQuantity == "fieldIntegral");
                }

                mapper->conservativeMapping(dataFieldB_j, dataFieldA_j);

                if (settingMappers[i].settingConservativeMapping.dataFieldATypeOfQuantity
                        == "field") { // deIntegrate dataFieldA_j
                    double *tmpDF = new double[numNodesA];
                    DFI_A->deIntegrate(dataFieldA_j, tmpDF);
                    for (int k = 0; k < numNodesA; k++)
                        dataFieldA_j[k] = tmpDF[k];
                    delete[] tmpDF;
                } else {
                    assert(
                            settingMappers[i].settingConservativeMapping.dataFieldATypeOfQuantity == "fieldIntegral");
                }

                for (int k = 0; k < numNodesA; k++)
                    dataFieldA[k * 3 + j] = dataFieldA_j[k];
                delete[] dataFieldA_j;
                delete[] dataFieldB_j;
            }

            GiDFileIO::appendNodalDataToDotRes(resultFileA, "\"conservativeMappingDataFieldA\"",
                    "\"testMapper\"", 1, "Vector", numNodesA, nodeIDsA, dataFieldA);
            GiDFileIO::appendNodalDataToDotRes(resultFileB, "\"conservativeMappingDataFieldB\"",
                    "\"testMapper\"", 1, "Vector", numNodesB, nodeIDsB, dataFieldB);

            if (settingMappers[i].settingConservativeMapping.doErrorCalculation) {
                double *dataFieldARef = new double[numNodesA * 3];
                initDataField(settingMappers[i].settingConservativeMapping.dataFieldARef, numNodesA,
                        numElemsA, nodeCoorsA, nodeIDsA, numNodesPerElemA, elemTableA, elemIDsA,
                        dataFieldARef);
                double *dataFieldAErr = new double[numNodesA * 3];
                for (int k = 0; k < numNodesA * 3; k++) {
                    dataFieldAErr[k] = dataFieldA[k] - dataFieldARef[k];
                    dataFieldAErr[k] = fabs(dataFieldAErr[k]);
                }
                GiDFileIO::appendNodalDataToDotRes(resultFileA,
                        "\"conservativeMappingDataFieldARef\"", "\"testMapper\"", 1, "Vector",
                        numNodesA, nodeIDsA, dataFieldARef);
                GiDFileIO::appendNodalDataToDotRes(resultFileA,
                        "\"conservativeMappingDataFieldAErr\"", "\"testMapper\"", 1, "Vector",
                        numNodesA, nodeIDsA, dataFieldAErr);

                // comoute L2 norm of error
                double vecL2Norm[3] = { 0.0, 0.0, 0.0 };
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < numNodesA; k++) {
                        vecL2Norm[j] += dataFieldAErr[k * 3 + j] * dataFieldAErr[k * 3 + j];
                    }
                }
                for (int j = 0; j < 3; j++) {
                    vecL2Norm[j] /= (double) numNodesA;
                    vecL2Norm[j] = sqrt(vecL2Norm[j]);
                }
                cout << '\t' << '\t' << "Error L2 norm: " << "(" << vecL2Norm[0] << ", "
                        << vecL2Norm[1] << ", " << vecL2Norm[2] << ")" << endl;

                delete[] dataFieldARef;
                delete[] dataFieldAErr;
            }

            delete[] dataFieldA;
            delete[] dataFieldB;
        }

        delete[] numNodesPerElemA;
        delete[] nodeCoorsA;
        delete[] nodeIDsA;
        delete[] elemTableA;
        delete[] elemIDsA;
        delete[] numNodesPerElemB;
        delete[] nodeCoorsB;
        delete[] nodeIDsB;
        delete[] elemTableB;
        delete[] elemIDsB;

        delete DFI_A;
        delete DFI_B;

        delete mapper;
    }
}

void TestMapper::doConservationAnalysis() {
    for (int i = 0; i < settingConservationAnalysises.size(); i++) {
        cout << "================================================" << endl;
        cout << "Do conservation analysises: " << settingConservationAnalysises[i].name << endl;

        // *. read meshes
        int numNodesA;
        int numElemsA;
        int* numNodesPerElemA;
        double *nodeCoorsA;
        int *nodeIDsA;
        int *elemTableA;
        int *elemIDsA;
        if (settingConservationAnalysises[i].meshAFormat == "GiDMesh") {
            GiDFileIO::readDotMsh(settingConservationAnalysises[i].meshFileA, numNodesA, numElemsA,
                    nodeCoorsA, nodeIDsA, numNodesPerElemA, elemTableA, elemIDsA);
        } else {
            assert(false);
        }

        DataFieldIntegration *DFI_A = new DataFieldIntegration(numNodesA, numElemsA,
                numNodesPerElemA, nodeCoorsA, nodeIDsA, elemTableA);

        int numNodesB;
        int numElemsB;
        int* numNodesPerElemB;
        double *nodeCoorsB;
        int *nodeIDsB;
        int *elemTableB;
        int *elemIDsB;
        if (settingConservationAnalysises[i].meshBFormat == "GiDMesh") {
            GiDFileIO::readDotMsh(settingConservationAnalysises[i].meshFileB, numNodesB, numElemsB,
                    nodeCoorsB, nodeIDsB, numNodesPerElemB, elemTableB, elemIDsB);
        } else {
            assert(false);
        }

        DataFieldIntegration *DFI_B = new DataFieldIntegration(numNodesB, numElemsB,
                numNodesPerElemB, nodeCoorsB, nodeIDsB, elemTableB);

        // *. read dataFields
        double *displacementsA = new double[numNodesA * 3];
        double *forcesA = new double[numNodesA * 3];
        double *displacementsB = new double[numNodesB * 3];
        double *forcesB = new double[numNodesB * 3];
        initDataField(settingConservationAnalysises[i].displacementsA, numNodesA, numElemsA,
                nodeCoorsA, nodeIDsA, numNodesPerElemA, elemTableA, elemIDsA, displacementsA);
        initDataField(settingConservationAnalysises[i].forcesA, numNodesA, numElemsA, nodeCoorsA,
                nodeIDsA, numNodesPerElemA, elemTableA, elemIDsA, forcesA);
        initDataField(settingConservationAnalysises[i].displacementsB, numNodesB, numElemsB,
                nodeCoorsB, nodeIDsB, numNodesPerElemB, elemTableB, elemIDsB, displacementsB);
        initDataField(settingConservationAnalysises[i].forcesB, numNodesB, numElemsB, nodeCoorsB,
                nodeIDsB, numNodesPerElemB, elemTableB, elemIDsB, forcesB);

        // *. convert dataFields
        for (int j = 0; j < 3; j++) { // x,y,z
            double *dataFieldA_j = new double[numNodesA];
            double *dataFieldB_j = new double[numNodesB];

            // convert displacementsA
            for (int k = 0; k < numNodesA; k++)
                dataFieldA_j[k] = displacementsA[k * 3 + j];
            if (settingConservationAnalysises[i].displacementsATypeOfQuantity == "fieldIntegral") { // deIntegrate dataFieldA_j
                double *tmpDF = new double[numNodesA];
                DFI_A->deIntegrate(dataFieldA_j, tmpDF);
                for (int k = 0; k < numNodesA; k++)
                    dataFieldA_j[k] = tmpDF[k];
                delete[] tmpDF;
            } else {
                assert( settingConservationAnalysises[i].displacementsATypeOfQuantity == "field");
            }
            for (int k = 0; k < numNodesA; k++)
                displacementsA[k * 3 + j] = dataFieldA_j[k];
            // convert forcesA
            for (int k = 0; k < numNodesA; k++)
                dataFieldA_j[k] = forcesA[k * 3 + j];
            if (settingConservationAnalysises[i].forcesATypeOfQuantity == "field") { // integrate dataFieldA_j
                double *tmpDF = new double[numNodesA];
                DFI_A->integrate(dataFieldA_j, tmpDF);
                for (int k = 0; k < numNodesA; k++)
                    dataFieldA_j[k] = tmpDF[k];
                delete[] tmpDF;
            } else {
                assert( settingConservationAnalysises[i].forcesATypeOfQuantity == "fieldIntegral");
            }
            for (int k = 0; k < numNodesA; k++)
                forcesA[k * 3 + j] = dataFieldA_j[k];
            // convert displacementsB
            for (int k = 0; k < numNodesB; k++)
                dataFieldB_j[k] = displacementsB[k * 3 + j];
            if (settingConservationAnalysises[i].displacementsBTypeOfQuantity == "fieldIntegral") { // deIntegrate dataFieldB_j
                double *tmpDF = new double[numNodesB];
                DFI_B->deIntegrate(dataFieldB_j, tmpDF);
                for (int k = 0; k < numNodesB; k++)
                    dataFieldB_j[k] = tmpDF[k];
                delete[] tmpDF;
            } else {
                assert( settingConservationAnalysises[i].displacementsBTypeOfQuantity == "field");
            }
            for (int k = 0; k < numNodesB; k++)
                displacementsB[k * 3 + j] = dataFieldB_j[k];
            // convert forcesB
            for (int k = 0; k < numNodesB; k++)
                dataFieldB_j[k] = forcesB[k * 3 + j];
            if (settingConservationAnalysises[i].forcesBTypeOfQuantity == "field") { // integrate dataFieldB_j
                double *tmpDF = new double[numNodesB];
                DFI_B->integrate(dataFieldB_j, tmpDF);
                for (int k = 0; k < numNodesB; k++)
                    dataFieldB_j[k] = tmpDF[k];
                delete[] tmpDF;
            } else {
                assert( settingConservationAnalysises[i].forcesBTypeOfQuantity == "fieldIntegral");
            }
            for (int k = 0; k < numNodesB; k++)
                forcesB[k * 3 + j] = dataFieldB_j[k];

            delete[] dataFieldA_j;
            delete[] dataFieldB_j;
        }

        // *. check conservation of force
        double resultantForceA[3] = { 0.0, 0.0, 0.0 };
        double resultantForceB[3] = { 0.0, 0.0, 0.0 };
        for (int j = 0; j < 3; j++) { // x,y,z
            for (int k = 0; k < numNodesA; k++) {
                resultantForceA[j] += forcesA[k * 3 + j];
            }
            for (int k = 0; k < numNodesB; k++) {
                resultantForceB[j] += forcesB[k * 3 + j];
            }
        }
        double resultantForceDiff[3];
        for (int j = 0; j < 3; j++) { // x,y,z
            resultantForceDiff[j] = resultantForceA[j] - resultantForceB[j];
        }
        double relativeError = MortarMath::computeVectorLength(resultantForceDiff)
                / MortarMath::computeVectorLength(resultantForceB);
        cout << "\t" << "check conservation of force ..." << endl;
        cout << "\t\t" << "resultantForceA: " << "(" << resultantForceA[0] << ", "
                << resultantForceA[1] << ", " << resultantForceA[2] << ")" << endl;
        cout << "\t\t" << "resultantForceB: " << "(" << resultantForceB[0] << ", "
                << resultantForceB[1] << ", " << resultantForceB[2] << ")" << endl;
        cout << "\t\t" << "resultantForceDiff: " << "(" << resultantForceDiff[0] << ", "
                << resultantForceDiff[1] << ", " << resultantForceDiff[2] << ")" << endl;
        cout << "\t\t" << "relativeError: " << relativeError << endl;

        // *. check conservation of energy
        double energyA = 0.0;
        double energyB = 0.0;
        for (int j = 0; j < numNodesA * 3; j++) {
            energyA += displacementsA[j] * forcesA[j];
        }
        for (int j = 0; j < numNodesB * 3; j++) {
            energyB += displacementsB[j] * forcesB[j];
        }
        double energyDiff = energyA - energyB;
        relativeError = energyDiff / (energyA + energyB) * 2.0;
        relativeError = fabs(relativeError);
        cout << "\t" << "check conservation of energy ..." << endl;
        cout << "\t\t" << "energyA: " << energyA << endl;
        cout << "\t\t" << "energyB: " << energyB << endl;
        cout << "\t\t" << "energyDiff: " << energyDiff << endl;
        cout << "\t\t" << "relativeError: " << relativeError << endl;

        delete[] numNodesPerElemA;
        delete[] nodeCoorsA;
        delete[] nodeIDsA;
        delete[] elemTableA;
        delete[] elemIDsA;
        delete[] numNodesPerElemB;
        delete[] nodeCoorsB;
        delete[] nodeIDsB;
        delete[] elemTableB;
        delete[] elemIDsB;

        delete DFI_A;
        delete DFI_B;

        delete[] displacementsA;
        delete[] forcesA;
        delete[] displacementsB;
        delete[] forcesB;
    }
}

void TestMapper::parseXMLDataField(ticpp::Element *xmlDataField, StructDataField &structDataField) {
    structDataField.format = xmlDataField->GetAttribute<string>("format");
    if (structDataField.format == "GiDResult") {
        ticpp::Element *xmlGiDResult = xmlDataField->FirstChildElement("GiDResult");
        structDataField.file = xmlGiDResult->GetAttribute<string>("file");
        structDataField.resultName = xmlGiDResult->GetAttribute<string>("resultName");
        structDataField.analysisName = xmlGiDResult->GetAttribute<string>("analysisName");
        structDataField.stepNum = xmlGiDResult->GetAttribute<int>("stepNum");
    } else if (structDataField.format == "function") {
        ticpp::Element *xmlFunction = xmlDataField->FirstChildElement("function");
        structDataField.function = xmlFunction->GetAttribute<string>("name");
    } else {
        assert(false);
    }
}

void TestMapper::outputXMLDataFieldToShell(StructDataField &structDataField, string dataFieldName,
        string indent) {
    cout << indent << dataFieldName << " { " << endl;
    string format = structDataField.format;
    cout << indent << "\t" << "format: " << format << endl;
    if (format == "function") {
        cout << indent << "\t" << "function: " << structDataField.function << endl;
    } else if (format == "GiDResult") {
        cout << indent << "\t" << "file: " << structDataField.file << endl;
        cout << indent << "\t" << "resultName: " << structDataField.resultName << endl;
        cout << indent << "\t" << "analysisName: " << structDataField.analysisName << endl;
        cout << indent << "\t" << "stepNum: " << structDataField.stepNum << endl;
    } else {
        assert(false);
    }
    cout << indent << "} " << endl;
}

void TestMapper::initDataField(StructDataField &structDataField, int numNodes, int numElems,
        double *nodeCoors, int *nodeIDs, int* numNodesPerElem, int *elemTable, int *elemIDs,
        double *dataField) {
    if (structDataField.format == "function") {
        string function = structDataField.function;
        const int DUMMY = 0;
        FieldCreator fc(numNodes, DUMMY, DUMMY, nodeCoors, NULL, NULL, function);
        fc.create(dataField);
    } else if (structDataField.format == "GiDResult") {
        string file = structDataField.file;
        string resultName = structDataField.resultName;
        string analysisName = structDataField.analysisName;
        int stepNum = structDataField.stepNum;
        GiDFileIO::readNodalDataFromDotRes(file, resultName, analysisName, stepNum, "vector",
                numNodes, nodeIDs, dataField);
    } else {
        assert(false);
    }

}

} /* namespace EMPIRE */
