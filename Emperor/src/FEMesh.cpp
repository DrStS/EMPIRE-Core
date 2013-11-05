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
#include <iostream>
#include <map>
#include "FEMesh.h"
#include "Message.h"
#include "DataField.h"
#include "TriangulatorAdaptor.h"

namespace EMPIRE {

using namespace std;

FEMesh::FEMesh(std::string _name, int _numNodes, int _numElems) :
        AbstractMesh(_name, EMPIRE_Mesh_FEMesh), numNodes(_numNodes), numElems(_numElems) {
    boundingBox.isComputed = false;
    nodes = new double[numNodes * 3];
    nodeIDs = new int[numNodes];
    numNodesPerElem = new int[numElems];
    elems = NULL;

    elemIDs = new int[numElems];
    for (int i = 0; i < numElems; i++)
        elemIDs[i] = i + 1; // set element ID by hand instead of receiving from the client, element ID starts from 1.

    tobeTriangulated = false;
    triangulatedMesh = NULL;
}

FEMesh::~FEMesh() {
    delete[] nodes;
    delete[] nodeIDs;
    delete[] numNodesPerElem;
    delete[] elems;
    delete[] elemIDs;
    if (triangulatedMesh != NULL)
        delete triangulatedMesh;
}

void FEMesh::initElems() {
    assert(elems==NULL);
    elemsArraySize = 0;

    for (int i = 0; i < numElems; i++) {
        int numNodesThisElem = numNodesPerElem[i];
        elemsArraySize += numNodesThisElem;
    }
    elems = new int[elemsArraySize];
    for (int i = 0; i < numElems; i++) {
        int numNodesThisElem = numNodesPerElem[i];
        if (numNodesThisElem > 3) {
            tobeTriangulated = true;
            break;
        }
    }
}

void FEMesh::addDataField(string dataFieldName, EMPIRE_DataField_location location,
        EMPIRE_DataField_dimension dimension, EMPIRE_DataField_typeOfQuantity typeOfQuantity) {
    int numLocatiions = -1;
    if (location == EMPIRE_DataField_atNode)
        numLocatiions = numNodes;
    if (location == EMPIRE_DataField_atElemCentroid)
        numLocatiions = numElems;
    assert(nameToDataFieldMap.find(dataFieldName) == nameToDataFieldMap.end());
    DataField *dataField = new DataField(dataFieldName, location, numLocatiions, dimension,
            typeOfQuantity);
    nameToDataFieldMap.insert(pair<string, DataField*>(dataFieldName, dataField));
}

FEMesh *FEMesh::triangulate() {
    if (!tobeTriangulated)
        return NULL;
    if (triangulatedMesh != NULL)
        return triangulatedMesh;

    map<int, int> *nodeIDToNodePosMap = new map<int, int>;
    for (int i = 0; i < numNodes; i++)
        nodeIDToNodePosMap->insert(nodeIDToNodePosMap->end(), pair<int, int>(nodeIDs[i], i));

    vector<int> *elemsTri = new vector<int>;
    vector<int> *numNodesPerElemTri = new vector<int>;
    int count = 0;
    for (int i = 0; i < numElems; i++) {
        int numNodesThisElem = numNodesPerElem[i];
        if (numNodesThisElem == 3) {
            numNodesPerElemTri->push_back(numNodesThisElem);
            for (int j = 0; j < numNodesThisElem; j++)
                elemsTri->push_back(elems[count + j]);
        } else {
            // if more than 4 nodes, use third party triangulation algorithm
            TriangulatorAdaptor *triangulator = new TriangulatorAdaptor();
            for (int j = 0; j < numNodesThisElem; j++) {
                int nodeID = elems[count + j];
                int nodePos = nodeIDToNodePosMap->at(nodeID);
                double x = nodes[nodePos * 3 + 0];
                double y = nodes[nodePos * 3 + 1];
                double z = nodes[nodePos * 3 + 2];
                triangulator->addPoint(x, y, z);
            }
            int numTriangles = numNodesThisElem - 2;
            int triangleIndexes[numTriangles * 3];
            triangulator->triangulate(triangleIndexes);
            for (int j = 0; j < numTriangles; j++) {
                numNodesPerElemTri->push_back(3);
                for (int k = 0; k < 3; k++) {
                    int nodeID = elems[count + triangleIndexes[j * 3 + k]];
                    elemsTri->push_back(nodeID);
                }
            }
            delete triangulator;
        }
        count += numNodesThisElem;
    }
    delete nodeIDToNodePosMap;
    assert(count == elemsArraySize);

    triangulatedMesh = new FEMesh(name + "_triangulated", numNodes, numNodesPerElemTri->size());
    for (int i = 0; i < numNodes; i++) {
        triangulatedMesh->nodeIDs[i] = nodeIDs[i];
    }
    for (int i = 0; i < numNodes * 3; i++) {
        triangulatedMesh->nodes[i] = nodes[i];
    }
    for (int i = 0; i < numNodesPerElemTri->size(); i++) {
        triangulatedMesh->numNodesPerElem[i] = numNodesPerElemTri->at(i);
    }
    triangulatedMesh->initElems();
    for (int i = 0; i < triangulatedMesh->elemsArraySize; i++) {
        triangulatedMesh->elems[i] = elemsTri->at(i);
    }

    delete numNodesPerElemTri;
    delete elemsTri;
    assert(triangulatedMesh->tobeTriangulated == false);
    return triangulatedMesh;
}

void FEMesh::computeBoundingBox() {
    if (boundingBox.isComputed)
        return;
    boundingBox.xmin = nodes[0 * 3 + 0];
    boundingBox.xmax = nodes[0 * 3 + 0];
    boundingBox.ymin = nodes[0 * 3 + 1];
    boundingBox.ymax = nodes[0 * 3 + 1];
    boundingBox.zmin = nodes[0 * 3 + 2];
    boundingBox.zmax = nodes[0 * 3 + 2];
    for (int i = 1; i < numNodes; i++) {
        double x = nodes[i * 3 + 0];
        double y = nodes[i * 3 + 1];
        double z = nodes[i * 3 + 2];
        if (x < boundingBox.xmin)
            boundingBox.xmin = x;
        else if (x > boundingBox.xmax)
            boundingBox.xmax = x;
        if (y < boundingBox.ymin)
            boundingBox.ymin = y;
        else if (y > boundingBox.ymax)
            boundingBox.ymax = y;
        if (z < boundingBox.zmin)
            boundingBox.zmin = z;
        else if (z > boundingBox.zmax)
            boundingBox.zmax = z;
    }
    boundingBox.isComputed = true;
}

void revertSurfaceNormalOfFEMesh(FEMesh *mesh) {
    int count = 0;
    for (int i = 0; i < mesh->numElems; i++) {
        int numNodesThisElem = mesh->numNodesPerElem[i];
        int tmpElem[numNodesThisElem];
        int *elems = &(mesh->elems[count]);
        for (int j = 0; j < numNodesThisElem; j++)
            tmpElem[j] = elems[numNodesThisElem - j - 1];
        for (int j = 0; j < numNodesThisElem; j++)
            elems[j] = tmpElem[j];
        count += numNodesThisElem;
    }
}

Message &operator<<(Message &message, FEMesh &mesh) {
    message << "\t+" << "FEMesh name: " << mesh.name << endl;
    message << "\t\t+" << "no. of nodes: " << mesh.numNodes << endl;
    message << "\t\t+" << "no. of elements: " << mesh.numElems << endl;
    message << "\t\t+" << "Nodes:" << endl;
    for (int i = 0; i < mesh.numNodes; i++) {
        message << "\t\t+" << "\t" << mesh.nodeIDs[i];
        for (int j = 0; j < 3; j++) {
            message << "\t" << mesh.nodes[i * 3 + j];
        }
        message << endl;
    }
    message << "\t\t+" << "Elements:" << endl;
    int count = 0;
    for (int i = 0; i < mesh.numElems; i++) {
        message << "\t\t+" << '\t' << i + 1;
        for (int j = 0; j < mesh.numNodesPerElem[i]; j++) {
            message << "\t" << mesh.elems[count + j];
        }
        count += mesh.numNodesPerElem[i];
        message << endl;
    }
    message() << "\t+" << "---------------------------------" << endl;
    return message;
}

} /* namespace EMPIRE */
