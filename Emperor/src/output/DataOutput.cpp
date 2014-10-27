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
#include "DataOutput.h"
#include "Emperor.h"
#include "MetaDataStructures.h"
#include "GiDFileIO.h"
#include "MatlabIGAFileIO.h"
#include "AbstractMesh.h"
#include "FEMesh.h"
#include "IGAMesh.h"
#include "DataField.h"
#include "ClientCode.h"
#include "Signal.h"

#include <assert.h>
#include <iostream>
#include <map>
#include <fstream>

using namespace std;

namespace EMPIRE {

DataOutput::DataOutput(const structDataOutput &_settingDataOutput,
        std::map<std::string, ClientCode*> &_nameToClientCodeMap) :
        settingDataOutput(_settingDataOutput), nameToClientCodeMap(_nameToClientCodeMap) {
}

DataOutput::~DataOutput() {
}

void DataOutput::init(std::string rearPart) {
    dataOutputName = settingDataOutput.name;
    dataOutputName.append(rearPart);
    writeMeshes();
    initDataFieldFiles();
    initSignalFiles();
}

void DataOutput::writeCurrentStep(int step) {
    writeDataFields(step);
    writeSignals(step);
}

void DataOutput::writeMeshes() {
    vector<structDataFieldRef> dataFieldRefs;
    for (int i = 0; i < settingDataOutput.connectionIOs.size(); i++) {
        if (settingDataOutput.connectionIOs[i].type == EMPIRE_ConnectionIO_DataField)
            dataFieldRefs.push_back(settingDataOutput.connectionIOs[i].dataFieldRef);
    }
    // get rid of repeated meshes
    map<string, AbstractMesh*> meshFileNameToMeshMap;
    for (int i = 0; i < dataFieldRefs.size(); i++) {
        const string UNDERSCORE = "_";
        string clientCodeName = dataFieldRefs[i].clientCodeName;
        string meshName = dataFieldRefs[i].meshName;
        string meshFileName = dataOutputName + UNDERSCORE + clientCodeName + UNDERSCORE + meshName
                + ".msh";

        assert(nameToClientCodeMap.find(clientCodeName) != nameToClientCodeMap.end());
        AbstractMesh *mesh = nameToClientCodeMap[clientCodeName]->getMeshByName(meshName);
        meshFileNameToMeshMap.insert(pair<string, AbstractMesh*>(meshFileName, mesh));
    }
    // write meshes
    for (map<string, AbstractMesh*>::iterator it = meshFileNameToMeshMap.begin();
            it != meshFileNameToMeshMap.end(); it++) {
        string meshFileName = it->first;
        if (it->second->type == EMPIRE_Mesh_FEMesh) {
            FEMesh *mesh = dynamic_cast<FEMesh*>(it->second);
            if (mesh->triangulate() != NULL)
                mesh = mesh->triangulate();
            GiDFileIO::writeDotMsh(meshFileName, mesh->numNodes, mesh->numElems, mesh->nodes,
                    mesh->nodeIDs, mesh->numNodesPerElem, mesh->elems, mesh->elemIDs);
        } else if (it->second->type == EMPIRE_Mesh_IGAMesh) {
            IGAMesh* igaMesh = dynamic_cast<IGAMesh*>(it->second);
            MatlabIGAFileIO::writeIGAMesh(igaMesh);
        } else
            assert(0);
    }
}

void DataOutput::initDataFieldFiles() {
    vector<structDataFieldRef> dataFieldRefs;
    for (int i = 0; i < settingDataOutput.connectionIOs.size(); i++) {
        if (settingDataOutput.connectionIOs[i].type == EMPIRE_ConnectionIO_DataField)
            dataFieldRefs.push_back(settingDataOutput.connectionIOs[i].dataFieldRef);
    }
    // get rid of repeated data field files (one data field file corresponds to multiple data fields of a single mesh)
    map<string, FEMesh*> dataFieldFileNameToMeshMap;
    for (int i = 0; i < dataFieldRefs.size(); i++) {
        const string UNDERSCORE = "_";
        string clientCodeName = dataFieldRefs[i].clientCodeName;
        string meshName = dataFieldRefs[i].meshName;
        string dataFieldFileName = dataOutputName + UNDERSCORE + clientCodeName + UNDERSCORE
                + meshName + ".res";

        assert(nameToClientCodeMap.find(clientCodeName) != nameToClientCodeMap.end());
        AbstractMesh *mesh = nameToClientCodeMap[clientCodeName]->getMeshByName(meshName);
        if (mesh->type == EMPIRE_Mesh_FEMesh) {
            FEMesh *feMesh = dynamic_cast<FEMesh*>(mesh);
            dataFieldFileNameToMeshMap.insert(pair<string, FEMesh*>(dataFieldFileName, feMesh));
        } else if (mesh->type == EMPIRE_Mesh_IGAMesh) {

        } else
            assert(0);
    }
    for (map<string, FEMesh*>::iterator it = dataFieldFileNameToMeshMap.begin();
            it != dataFieldFileNameToMeshMap.end(); it++) {
        string dataFieldFileName = it->first;
        FEMesh *mesh = it->second;

        GiDFileIO::initDotRes(dataFieldFileName); // TODO enable output write data on element center of hybrid mesh
    }
}

void DataOutput::writeDataFields(int step) {
    int interval = settingDataOutput.interval;
    if (step % interval == 0) {
        vector<structDataFieldRef> dataFieldRefs;
        for (int i = 0; i < settingDataOutput.connectionIOs.size(); i++) {
            if (settingDataOutput.connectionIOs[i].type == EMPIRE_ConnectionIO_DataField)
                dataFieldRefs.push_back(settingDataOutput.connectionIOs[i].dataFieldRef);
        }

        for (int i = 0; i < dataFieldRefs.size(); i++) {
            const structDataFieldRef &dataFieldRef = dataFieldRefs[i];
            string clientCodeName = dataFieldRef.clientCodeName;
            string meshName = dataFieldRef.meshName;
            string dataFieldName = dataFieldRef.dataFieldName;
            AbstractMesh *mesh = nameToClientCodeMap[clientCodeName]->getMeshByName(meshName);
            if (mesh->type == EMPIRE_Mesh_FEMesh) {
                FEMesh *feMesh = dynamic_cast<FEMesh*>(mesh);
                const string UNDERSCORE = "_";
                string dataFieldFileName = dataOutputName + UNDERSCORE + clientCodeName + UNDERSCORE
                        + meshName + ".res";
                DataField *dataField = feMesh->getDataFieldByName(dataFieldName);
                bool atNode = (dataField->location == EMPIRE_DataField_atNode ? true : false);
                int *locationIDs = (atNode ? feMesh->nodeIDs : feMesh->elemIDs);
                string type;
                if (dataField->dimension == EMPIRE_DataField_vector)
                    type = "Vector";
                else if (dataField->dimension == EMPIRE_DataField_scalar)
                    type = "Scalar";
                else
                    assert(false);
                string tmpdataFieldName = "\"" + dataFieldName + "\"";
                if (atNode) {
                    GiDFileIO::appendNodalDataToDotRes(dataFieldFileName, tmpdataFieldName,
                            "\"EMPIRE_CoSimulation\"", step, type, dataField->numLocations,
                            locationIDs, dataField->data);
                } else {
                    if (feMesh->triangulate() == NULL) {
                        GiDFileIO::appendElementalDataToDotRes(dataFieldFileName, tmpdataFieldName,
                                "\"EMPIRE_CoSimulation\"", step, type, dataField->numLocations,
                                locationIDs, feMesh->numNodesPerElem, dataField->data);
                    } else {
                        assert(false); // writing out data field on element certeroid of triangulated mesh is not implemented yet
                    }
                }
            } else if (mesh->type == EMPIRE_Mesh_IGAMesh) {
                DataField *dataField = mesh->getDataFieldByName(dataFieldName);
                MatlabIGAFileIO::writeVectorFieldOnCPs(dataFieldName, step, dataField);
            } else {
                assert(0);
            }
        }
    }
}

void DataOutput::initSignalFiles() {
    vector<structSignalRef> signalRefs;
    for (int i = 0; i < settingDataOutput.connectionIOs.size(); i++) {
        if (settingDataOutput.connectionIOs[i].type == EMPIRE_ConnectionIO_Signal)
            signalRefs.push_back(settingDataOutput.connectionIOs[i].signalRef);
    }
    for (int i = 0; i < signalRefs.size(); i++) {
        const string UNDERSCORE = "_";
        string clientCodeName = signalRefs[i].clientCodeName;
        string signalName = signalRefs[i].signalName;
        string signalFileName = dataOutputName + UNDERSCORE + clientCodeName + UNDERSCORE
                + signalName + ".csv";
        fstream signalFile;
        signalFile.open(signalFileName.c_str(), ios_base::out);
        assert(!signalFile.fail());
        assert(nameToClientCodeMap.find(clientCodeName) != nameToClientCodeMap.end());
        const Signal *signal = nameToClientCodeMap[clientCodeName]->getSignalByName(signalName);

        signalFile << "Time";
        for (int j = 0; j < signal->size; j++) {
            signalFile << '\t' << "signal[" << j << "]";
        }
        signalFile << endl;

        signalFile.close();
    }
}

void DataOutput::writeSignals(int step) {
    vector<structSignalRef> signalRefs;
    for (int i = 0; i < settingDataOutput.connectionIOs.size(); i++) {
        if (settingDataOutput.connectionIOs[i].type == EMPIRE_ConnectionIO_Signal)
            signalRefs.push_back(settingDataOutput.connectionIOs[i].signalRef);
    }
    for (int i = 0; i < signalRefs.size(); i++) {
        const string UNDERSCORE = "_";
        string clientCodeName = signalRefs[i].clientCodeName;
        string signalName = signalRefs[i].signalName;
        string signalFileName = dataOutputName + UNDERSCORE + clientCodeName + UNDERSCORE
                + signalName + ".csv";
        fstream signalFile;
        signalFile.open(signalFileName.c_str(), ios_base::out | ios_base::app);
        assert(!signalFile.fail());
        assert(nameToClientCodeMap.find(clientCodeName) != nameToClientCodeMap.end());
        const Signal *signal = nameToClientCodeMap[clientCodeName]->getSignalByName(signalName);

        signalFile << step;
        for (int j = 0; j < signal->size; j++) {
            signalFile << '\t' << signal->array[j];
        }
        signalFile << endl;

        signalFile.close();
    }
}

} /* namespace EMPIRE */
