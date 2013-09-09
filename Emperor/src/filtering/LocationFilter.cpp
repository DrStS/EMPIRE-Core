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
#include <assert.h>
#include <map>

#include "LocationFilter.h"
#include "AbstractMesh.h"
#include "FEMesh.h"
#include "ClientCode.h"
#include "DataField.h"
#include "ConnectionIO.h"

namespace EMPIRE {

using namespace std;

LocationFilter::LocationFilter() :
        AbstractFilter(), mesh(NULL), feMesh(NULL) {
}

LocationFilter::~LocationFilter() {
    if (caseNum == 1) {
        const int numNodes = feMesh->numNodes;
        for (int i = 0; i < numNodes; i++)
            delete nodePosToElemTable[i];
        delete[] nodePosToElemTable;
    } else {
        assert(false);
    }
}

void LocationFilter::filtering() {
    switch (caseNum) {
    case 1:
        filterDataFieldCase1();
        break;
    default:
        assert(false);
        break;
    }

}

void LocationFilter::init() {
    assert(inputVec.size() == 1);
    assert(outputVec.size() == 1);
    assert(inputVec[0]->type == outputVec[0]->type);
    assert(inputVec[0]->type == EMPIRE_ConnectionIO_DataField);
    const DataField *inDataField = inputVec[0]->dataField;
    DataField *outDataField = outputVec[0]->dataField;
    assert(inDataField->dimension == outDataField->dimension);
    assert(inDataField->typeOfQuantity == outDataField->typeOfQuantity);

    assert(mesh==NULL);
    assert(inputVec[0]->mesh == outputVec[0]->mesh);
    mesh = inputVec[0]->mesh;

    feMesh = dynamic_cast<FEMesh*>(mesh);

    if (inDataField->typeOfQuantity == EMPIRE_DataField_field) {
        if ((inDataField->location == EMPIRE_DataField_atElemCentroid)
                && (outDataField->location == EMPIRE_DataField_atNode)) {
            caseNum = 1;
            computeNodePosToElemTable(); // mesh should be set before here!!!
        }
    } else {
        assert(false);
    }
}

void LocationFilter::filterDataFieldCase1() {
    const DataField *inDataField = inputVec[0]->dataField;
    DataField *outDataField = outputVec[0]->dataField;
    const int dimension = inDataField->dimension;
    const double *inData = inDataField->data;
    double *outData = outDataField->data;

    const int numNodes = feMesh->numNodes;
    for (int i = 0; i < numNodes; i++) {
        vector<int> *elemsContainingMe = nodePosToElemTable[i];
        double weight = 1.0 / double(elemsContainingMe->size());
        for (int j = 0; j < dimension; j++) {
            double sum = 0.0;
            for (int k = 0; k < elemsContainingMe->size(); k++)
                sum += inData[elemsContainingMe->at(k) * dimension + j];
            outData[i * dimension + j] = weight * sum;
        }
    }
}

void LocationFilter::computeNodePosToElemTable() {
    const int numNodes = feMesh->numNodes;
    const int numElems = feMesh->numElems;
    const int *nodeIDs = feMesh->nodeIDs;
    const int *elems = feMesh->elems;

    nodePosToElemTable = new vector<int>*[numNodes];
    for (int i = 0; i < numNodes; i++)
        nodePosToElemTable[i] = new vector<int>;

    // a map linking a nodeID to its position in nodes is required
    map<int, int> *nodeIDToNodePosMap = new map<int, int>;
    for (int i = 0; i < numNodes; i++)
        nodeIDToNodePosMap->insert(nodeIDToNodePosMap->end(), pair<int, int>(nodeIDs[i], i));

    int count = 0;
    for (int i = 0; i < numElems; i++) {
        int numNodesThisElem = feMesh->numNodesPerElem[i];
        for (int j = 0; j < numNodesThisElem; j++) {
            int nodeID = elems[count + j];
            int nodePos = nodeIDToNodePosMap->at(nodeID);
            nodePosToElemTable[nodePos]->push_back(i);
        }
        count += numNodesThisElem;
    }assert(count == feMesh->elemsArraySize);
    delete nodeIDToNodePosMap;
}

} /* namespace EMPIRE */
