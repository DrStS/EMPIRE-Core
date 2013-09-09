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
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <list>
#include <map>
#include <assert.h>
#include <vector>
#include "stdlib.h"

#include "GiDFileIO.h"

using namespace std;

namespace GiDFileIO {
/********//**
 * \brief Class Node stores node data
 ***********/
class Node {
public:
    Node(int _id, double *_coordinate) :
            id(_id), coordinate(_coordinate) {
    }
    virtual ~Node() {
        delete[] coordinate;
    }
    const int id;
    const double *coordinate;
};
/********//**
 * \brief Class Element stores element data
 ***********/
class Element {
public:
    Element(int _id, int *_nodeIds, int _numberOfNodes) :
            id(_id), nodeIds(_nodeIds), numberOfNodes(_numberOfNodes) {
    }
    virtual ~Element() {
        delete[] nodeIds;
    }
    const int id;
    const int * nodeIds;
    const int numberOfNodes;
};
/***********************************************************************************************
 * \brief Compare two elements by their number of nodes
 * \param[in] first first element
 * \param[in] second second element
 * \return true if the first element has fewer nodes. Otherwise false.
 * \author Michael Andre
 ***********/
bool compare_elem_by_num_node(Element *first, Element *second) {
    return (first->numberOfNodes < second->numberOfNodes) ? true : false;
}
/***********************************************************************************************
 * \brief Compare two elements by their ids
 * \param[in] first first element
 * \param[in] second second element
 * \return true if the first element has a smaller id
 * \author Michael Andre
 ***********/
bool compare_elem_by_id(Element *first, Element *second) {
    return (first->id < second->id) ? true : false;
}
/***********************************************************************************************
 * \brief Compare two nodes by their ids
 * \param[in] first first node
 * \param[in] second second second
 * \return true if the first node has a smaller id
 * \author Michael Andre
 ***********/
bool compare_node_by_id(Node *first, Node *second) {
    return (first->id < second->id) ? true : false;
}
/***********************************************************************************************
 * \brief Convert the string to a lower case string
 * \param[in] str the string
 * \author Tianyang Wang
 ***********/
void convertToLowerCase(string &str) {
    for (unsigned i = 0; i < str.size(); i++) {
        str[i] = tolower(str[i]);
    }
}
/***********************************************************************************************
 * \brief Output the header line of a gid result file
 * \param[in] outputStream output stream
 * \author Tianyang Wang
 ***********/
void outputHeaderLine(ostream &outputStream) {
    outputStream << "GiD Post Results File 1.0" << endl << endl;
}
/***********************************************************************************************
 * \brief Output the gauss points block of a gid result file
 * \param[in] outputStream output stream
 * \author Tianyang Wang
 ***********/
void outputGaussPoints(ostream &outputStream) {
    // write GaussPoints for triangles in the mesh
    outputStream << "GaussPoints " << "\"GP_triangles\"" << " ElemType " << "Triangle" << " "
            << "\"mesh_triangles\"" << endl; // use space instead of tab, otherwise GiD would fail to parse
    outputStream << "  " << "Number Of Gauss Points: 1" << endl;
    outputStream << "  " << "Natural Coordinates: internal" << endl;
    outputStream << "end gausspoints" << endl << endl;
    // write GaussPoints for quads in the mesh
    outputStream << "GaussPoints " << "\"GP_quads\"" << " ElemType " << "Quadrilateral" << " "
            << "\"mesh_quads\"" << endl; // use space instead of tab, otherwise GiD would fail to parse
    outputStream << "  " << "Number Of Gauss Points: 1" << endl;
    outputStream << "  " << "Natural Coordinates: internal" << endl;
    outputStream << "end gausspoints" << endl << endl;
}

void readDotMsh(string fileName, int &numberOfMeshNodes, int &numberOfElements,
        double *&meshNodeCoordinates, int *&meshNodeIds, int *&numberOfNodesPerElement,
        int *&elementNodeTables, int *&elementIds) {
    ifstream inputStream(fileName.c_str());
    if (!inputStream) {
        cerr << "readDotMsh: GiD .msh file \"" << fileName << "\" cannot be found" << '\n';
        exit(EXIT_FAILURE);
    }
    const int XYZ = 3;
    bool inNodes = false, inElements = false, inMesh = false;
    int numberOfNodesThisElement = 0;
    map<int, Node *> *nodeAtId = new map<int, Node *>;
    list<Node *> *meshNodes = new list<Node *>;
    list<Element *> *elements = new list<Element *>;
    string elementType, textLine, meshName, lineToken;

    while (true) { /* parse line */
        getline(inputStream, textLine, '\n');

        istringstream lineStream(textLine);
        lineStream >> skipws >> lineToken;
        convertToLowerCase(lineToken);

        if (lineToken.size() == 0 || lineToken[0] == '#') {
            // ignore comments and empty lines
        } else if (inNodes) {
            if (lineToken.compare("end") == 0) {
                inNodes = false;
            } else { // add node
                int id;
                istringstream lineStream(textLine);
                lineStream >> id;
                double *coordinate = new double[XYZ];
                for (int i = 0; i < XYZ; i++) {
                    lineStream >> coordinate[i];
                }assert(nodeAtId->find(id) == nodeAtId->end());
                // enforce unique node ids
                (*nodeAtId)[id] = new Node(id, coordinate);
            }
        } else if (inElements) {
            if (lineToken.compare("end") == 0) {
                inElements = false;
                inMesh = false;
            } else {
                { // 1. check the number of entries in this line
                    istringstream lineStream(textLine);
                    int count = 0;
                    while (lineStream >> lineToken)
                        count++;
                    // there may be material number
                    assert(
                            (count == (numberOfNodesThisElement + 1)) || (count == (numberOfNodesThisElement + 2)));
                }
                { // 2. parse the elements
                    istringstream lineStream(textLine);
                    int id;
                    lineStream >> id;
                    int *nodeIds = new int[numberOfNodesThisElement];
                    for (int i = 0; i < numberOfNodesThisElement; i++) {
                        lineStream >> nodeIds[i];
                    }
                    elements->push_back(new Element(id, nodeIds, numberOfNodesThisElement));
                }
            }
        } else if (lineToken.compare("mesh") == 0) { // new mesh
            int dimension;
            inMesh = true;
            int count = 0;
            { // get number of words in the line
                istringstream lineStream(textLine);
                while (lineStream >> lineToken)
                    count++;
            }
            { // extract header
                const int headerLineWordsWithMeshName = 8;
                istringstream lineStream(textLine);
                lineStream >> lineToken; // skip the word "mesh"
                meshName.clear(); // ignore previous mesh names
                if (count == headerLineWordsWithMeshName) { // there is a mesh name
                    do {
                        lineStream >> lineToken;
                        meshName.append(" ").append(lineToken); // here a space is added to the beginning
                    } while (meshName[meshName.size() - 1] != '\"'); // parse the name with spaces in
                    meshName.erase(meshName.begin()); // remove the space at the beginning
                } else {
                    assert(count == 7);
                    meshName = "\"\"";
                }

                lineStream >> lineToken; // skip the keyword "dimension"
                convertToLowerCase(lineToken);
                assert(lineToken.compare("dimension") == 0);
                lineStream >> dimension;
                assert(dimension == XYZ);

                lineStream >> lineToken; // skip the keyword "ElemType"
                convertToLowerCase(lineToken);
                assert(lineToken.compare("elemtype") == 0);
                lineStream >> elementType;

                lineStream >> lineToken; // skip the keyword "Nnode"
                convertToLowerCase(lineToken);
                assert(lineToken.compare("nnode") == 0);
                lineStream >> numberOfNodesThisElement;

                if (elementType == "Triangle")
                    assert(numberOfNodesThisElement == 3);
                else if (elementType == "Quadrilateral")
                    assert(numberOfNodesThisElement == 4);
                else
                    assert(false);
            }
        } else if (lineToken.compare("coordinates") == 0) {
            inNodes = true;
        } else if (lineToken.compare("elements") == 0) {
            inElements = true;
        }
        // terminate if failbit, eofbit or badbit is set
        if (!inputStream.good())
            break;
    } /* parse line */

    assert(inMesh == false);
    assert(inNodes == false);
    assert(inElements == false);

    numberOfElements = elements->size();

    { // number of nodes per element
        int i = 0;
        numberOfNodesPerElement = new int[numberOfElements];
        list<Element *>::iterator it;
        for (it = elements->begin(); it != elements->end(); ++it, ++i)
            numberOfNodesPerElement[i] = (*it)->numberOfNodes;
    }

    { // element Ids
        int i = 0;
        elementIds = new int[numberOfElements];
        list<Element *>::iterator it;
        for (it = elements->begin(); it != elements->end(); ++it, ++i)
            elementIds[i] = (*it)->id;
    }

    { // element tables
        int size = 0;
        list<Element *>::iterator it;
        for (it = elements->begin(); it != elements->end(); ++it) {
            size += (*it)->numberOfNodes;
        }
        elementNodeTables = new int[size];

        int i_g = 0;
        for (it = elements->begin(); it != elements->end(); ++it) {
            for (int i = 0; i < (*it)->numberOfNodes; i++)
                elementNodeTables[i_g++] = (*it)->nodeIds[i];
        }
    }

    { // get nodes with element association
        int id;
        list<Element *>::iterator it;
        for (it = elements->begin(); it != elements->end(); ++it) {
            for (int i = 0; i < (*it)->numberOfNodes; i++) {
                id = (*it)->nodeIds[i];
                assert(nodeAtId->find(id) != nodeAtId->end());
                // node must exist
                meshNodes->push_back((*nodeAtId)[id]);
            }
        }
        meshNodes->sort(); // remove duplicates
        meshNodes->unique();
        meshNodes->sort(compare_node_by_id);
    }

    numberOfMeshNodes = meshNodes->size();

    { // node coordinates
        int i = 0;
        meshNodeCoordinates = new double[numberOfMeshNodes * XYZ];
        list<Node *>::iterator it;
        for (it = meshNodes->begin(); it != meshNodes->end(); ++it, ++i)
            for (int j = 0; j < XYZ; j++)
                meshNodeCoordinates[i * XYZ + j] = (*it)->coordinate[j];
    }

    { // node Ids
        int i = 0;
        meshNodeIds = new int[numberOfMeshNodes];
        list<Node *>::iterator it;
        for (it = meshNodes->begin(); it != meshNodes->end(); ++it, ++i)
            meshNodeIds[i] = (*it)->id;
    }

    { // delete nodes
        for (map<int, Node *>::iterator it = nodeAtId->begin(); it != nodeAtId->end(); ++it)
            delete (it->second);
        delete nodeAtId;
        delete meshNodes;
    }

    { // delete elements
        for (list<Element *>::iterator it = elements->begin(); it != elements->end(); ++it)
            delete (*it);
        delete elements;
    }
}

void writeDotMsh(std::string fileName, int numberOfMeshNodes,
        int numberOfElements, const double *meshNodeCoordinates, const int *meshNodeIds,
        const int *numberOfNodesPerElement, const int *elementNodeTables, const int *elementIds) {
    const int XYZ = 3;

    map<int, Node *> *nodeAtId = new map<int, Node *>;
    map<int, bool> *writeNodeAtId = new map<int, bool>; // if not written yet, true; if written already, false
    list<Element *> *elements = new list<Element *>;

    { // extract elements
        int i_g = 0;
        for (int i = 0; i < numberOfElements; i++) {
            int *elementNodeTable = new int[numberOfNodesPerElement[i]];
            for (int j = 0; j < numberOfNodesPerElement[i]; j++)
                elementNodeTable[j] = elementNodeTables[i_g++];

            elements->push_back(
                    new Element(elementIds[i], elementNodeTable, numberOfNodesPerElement[i]));
        }
    }

    { // extract nodes
        for (int i = 0; i < numberOfMeshNodes; i++) {
            double *coordinate = new double[XYZ];
            for (int j = 0; j < XYZ; j++)
                coordinate[j] = meshNodeCoordinates[i * XYZ + j];

            (*nodeAtId)[meshNodeIds[i]] = new Node(meshNodeIds[i], coordinate);
            (*writeNodeAtId)[meshNodeIds[i]] = true;
        }
    }

    // order elements by number of nodes
    elements->sort(compare_elem_by_num_node);

    ofstream outputStream(fileName.c_str());

    while (elements->size() > 0) {
        string elementType;
        int numberOfNodesThisElement = (*elements->begin())->numberOfNodes;
        list<Node *> *nodeSublist = new list<Node *>;
        list<Element *> *elementSublist = new list<Element *>;

        if (numberOfNodesThisElement == 3)
            elementType = "Triangle";
        else if (numberOfNodesThisElement == 4)
            elementType = "Quadrilateral";
        else
            assert(false);

        { // extract element sublist
            list<Element *>::iterator it;
            for (it = elements->begin(); it != elements->end(); ++it) {
                if ((*it)->numberOfNodes != numberOfNodesThisElement)
                    break;
            }
            elementSublist->splice(elementSublist->end(), *elements, elements->begin(), it);
            elementSublist->sort(compare_elem_by_id);
        }

        { // extract node sublist
            list<Element *>::iterator it;
            for (it = elementSublist->begin(); it != elementSublist->end(); ++it) {
                for (int i = 0; i < numberOfNodesThisElement; i++) {
                    int id = (*it)->nodeIds[i];
                    assert(writeNodeAtId->find(id) != writeNodeAtId->end());
                    // node must exist

                    if ((*writeNodeAtId)[id]) {
                        (*writeNodeAtId)[id] = false;
                        nodeSublist->push_back((*nodeAtId)[id]);
                    }
                }
            }
        }
        nodeSublist->sort(compare_node_by_id);

        { // write header
            string meshName;
            if (numberOfNodesThisElement == 3)
                meshName = "\"mesh_triangles\"";
            else if (numberOfNodesThisElement == 4)
                meshName = "\"mesh_quads\"";
            else
                assert(false);
            outputStream << "MESH " << meshName << " dimension " << XYZ << " Elemtype ";
            outputStream << elementType << " Nnode " << numberOfNodesThisElement << endl;
        }

        { // write nodes
            Node *node;
            list<Node *>::iterator it;
            outputStream << "Coordinates" << endl;
            for (it = nodeSublist->begin(); it != nodeSublist->end(); it++) {
                node = *it;
                outputStream << node->id;
                for (int i = 0; i < XYZ; i++)
                    outputStream << '\t' << node->coordinate[i];
                outputStream << endl;
            }
            outputStream << "End Coordinates" << endl;
        }

        { // write elements
            Element *element;
            list<Element *>::iterator it;
            outputStream << "Elements" << endl;
            for (it = elementSublist->begin(); it != elementSublist->end(); it++) {
                element = *it;
                outputStream << element->id;
                for (int i = 0; i < numberOfNodesThisElement; i++)
                    outputStream << '\t' << element->nodeIds[i];
                outputStream << endl;
            }
            outputStream << "End Elements" << endl;
        }
        { // delete elements in sublist
            list<Element *>::iterator it;
            for (it = elementSublist->begin(); it != elementSublist->end(); it++)
                delete (*it);
            delete elementSublist;
        }
        delete nodeSublist;
    }

    { // delete rest of elements
        for (list<Element *>::iterator it = elements->begin(); it != elements->end(); it++)
            delete (*it);
        delete elements;
    }
    { // delete nodes
        for (map<int, Node *>::iterator it = nodeAtId->begin(); it != nodeAtId->end(); it++)
            delete (it->second);
        delete nodeAtId;
    }
    delete writeNodeAtId;
}

void initDotRes(const string fileName) {
    fstream dotResFile;
    dotResFile.open(fileName.c_str(), ios_base::out);
    assert(!dotResFile.fail());
    outputHeaderLine(dotResFile);
    outputGaussPoints(dotResFile);
    dotResFile.close();
}

void appendNodalDataToDotRes(std::string fileName, std::string resultName,
        std::string analysisName, int stepNum, std::string type,
        int numberOfNodes, const int *nodeIds, const double *data) {
    fstream dotResFile;
    dotResFile.open(fileName.c_str(), ios_base::out | ios_base::app);
    assert(!dotResFile.fail());

    assert(resultName.size()!=0);
    if (resultName[0] != '\"')
        resultName.insert(0, "\"");
    if (resultName[resultName.size()-1] != '\"')
        resultName.append("\"");
    assert(analysisName.size()!=0);
    if (analysisName[0] != '\"')
        analysisName.insert(0, "\"");
    if (analysisName[analysisName.size()-1] != '\"')
        analysisName.append("\"");

    convertToLowerCase(type);
    int dimension = 0;
    if (type == "vector")
        dimension = 3;
    else if (type == "scalar")
        dimension = 1;
    else
        assert(false);

    dotResFile << "Result\t" << resultName << '\t' << analysisName << '\t' << stepNum << '\t'
            << type << '\t' << "OnNodes" << endl;
    dotResFile << "Values" << endl;

    for (int i = 0; i < numberOfNodes; i++) {
        dotResFile << '\t' << nodeIds[i];
        for (int j = 0; j < dimension; j++)
            dotResFile << '\t' << data[i * dimension + j];
        dotResFile << endl;
    }
    dotResFile << "End Values" << endl << endl;

    dotResFile.close();
}

void appendElementalDataToDotRes(std::string fileName, std::string resultName,
        std::string analysisName, int stepNum, std::string type,
        int numberOfElements, const int *elemIds, const int *numberOfNodesPerElement,
        const double *data) {
    fstream dotResFile;
    dotResFile.open(fileName.c_str(), ios_base::out | ios_base::app);
    assert(!dotResFile.fail());

    assert(resultName.size()!=0);
    if (resultName[0] != '\"')
        resultName.insert(0, "\"");
    if (resultName[resultName.size()-1] != '\"')
        resultName.append("\"");
    assert(analysisName.size()!=0);
    if (analysisName[0] != '\"')
        analysisName.insert(0, "\"");
    if (analysisName[analysisName.size()-1] != '\"')
        analysisName.append("\"");

    convertToLowerCase(type);
    int dimension = 0;
    if (type == "vector")
        dimension = 3;
    else if (type == "scalar")
        dimension = 1;
    else
        assert(false);

    vector<int> *triangles = new vector<int>;
    vector<int> *quads = new vector<int>;
    for (int i = 0; i < numberOfElements; i++) {
        if (numberOfNodesPerElement[i] == 3)
            triangles->push_back(i);
        else if (numberOfNodesPerElement[i] == 4)
            quads->push_back(i);
        else
            assert(false);
    }

    if (triangles->size() != 0) { // output data on triangular element
        dotResFile << "Result\t" << resultName << '\t' << analysisName << '\t' << stepNum << '\t'
                << type << '\t' << "OnGaussPoints" << '\t' << "\"GP_triangles\"" << endl;
        dotResFile << "Values" << endl;

        for (int i = 0; i < triangles->size(); i++) {
            int pos = triangles->at(i);
            dotResFile << '\t' << elemIds[pos];
            for (int j = 0; j < dimension; j++)
                dotResFile << '\t' << data[pos * dimension + j];
            dotResFile << endl;
        }
        dotResFile << "End Values" << endl << endl;
    }

    if (quads->size() != 0) { // output data on quad element
        dotResFile << "Result\t" << resultName << '\t' << analysisName << '\t' << stepNum << '\t'
                << type << '\t' << "OnGaussPoints" << '\t' << "\"GP_quads\"" << endl;
        dotResFile << "Values" << endl;

        for (int i = 0; i < quads->size(); i++) {
            int pos = quads->at(i);
            dotResFile << '\t' << elemIds[pos];
            for (int j = 0; j < dimension; j++)
                dotResFile << '\t' << data[pos * dimension + j];
            dotResFile << endl;
        }
        dotResFile << "End Values" << endl << endl;
    }
    dotResFile.close();
    delete triangles;
    delete quads;
}

void readNodalDataFromDotRes(std::string fileName, std::string resultName,
        std::string analysisName, int stepNum, std::string type,
        int numberOfNodes, const int *nodeIds, double *data) {
    ifstream dotResFile(fileName.c_str());
    assert(dotResFile);

    assert(resultName.size()!=0);
    if (resultName[0] != '\"')
        resultName.insert(0, "\"");
    if (resultName[resultName.size()-1] != '\"')
        resultName.append("\"");
    assert(analysisName.size()!=0);
    if (analysisName[0] != '\"')
        analysisName.insert(0, "\"");
    if (analysisName[analysisName.size()-1] != '\"')
        analysisName.append("\"");


    int dimension = 0;
    convertToLowerCase(type);
    if (type == "vector")
        dimension = 3;
    else if (type == "scalar")
        dimension = 1;
    else
        assert(false);

    map<int, int> *nodeIDToPosMap = new map<int, int>;
    for (int i = 0; i < numberOfNodes; i++) {
        nodeIDToPosMap->insert(nodeIDToPosMap->end(), pair<int, int>(nodeIds[i], i));
    }

    while (true) { /* parse line */
        if (!dotResFile.good()) { // terminate if failbit, eofbit or badbit is set
            cerr << "result with the following parameters is not found in " << fileName << endl;
            cerr << "result name: " << resultName << endl;
            cerr << "analysis name: " << analysisName << endl;
            cerr << "step number: " << stepNum << endl;
            cerr << "type: " << type << endl;
            cerr << "Exit!" << endl;
            exit(EXIT_FAILURE);
        }
        string textLine;
        string lineToken;
        getline(dotResFile, textLine, '\n');

        istringstream lineStream(textLine);
        lineStream >> skipws >> lineToken;
        convertToLowerCase(lineToken);
        if (lineToken == "result") {
            string resultNameInFile;
            do {
                lineStream >> lineToken;
                resultNameInFile.append(" ").append(lineToken); // here a space is added to the beginning
            } while (resultNameInFile[resultNameInFile.size() - 1] != '\"'); // parse the name with spaces in
            resultNameInFile.erase(resultNameInFile.begin()); // remove the space at the beginning

            string analysisNameInFile;
            do {
                lineStream >> lineToken;
                analysisNameInFile.append(" ").append(lineToken); // here a space is added to the beginning
            } while (analysisNameInFile[analysisNameInFile.size() - 1] != '\"'); // parse the name with spaces in
            analysisNameInFile.erase(analysisNameInFile.begin()); // remove the space at the beginning

            int stepNumInFile;
            lineStream >> stepNumInFile;

            string typeInFile;
            lineStream >> typeInFile;
            convertToLowerCase(typeInFile);

            string locationInFile;
            lineStream >> locationInFile;
            convertToLowerCase(locationInFile);

            if (stepNumInFile > stepNum) { // save time if it is beyond the interesting time step
                cerr << "result with the following parameters is not found in " << fileName << endl;
                cerr << "result name: " << resultName << endl;
                cerr << "analysis name: " << analysisName << endl;
                cerr << "step number: " << stepNum << endl;
                cerr << "type: " << type << endl;
                cerr << "Exit!" << endl;
                exit(EXIT_FAILURE);
            }

            if (resultName == resultNameInFile && analysisName == analysisNameInFile
                    && stepNum == stepNumInFile && type == typeInFile
                    && locationInFile == "onnodes") { // read nodal data, do not allow empty lines or comments
                getline(dotResFile, textLine, '\n');
                istringstream lineStream1(textLine);
                lineStream1 >> skipws >> lineToken;
                convertToLowerCase(lineToken);
                assert(lineToken == "values");
                for (int i = 0; i < numberOfNodes; i++) {
                    getline(dotResFile, textLine, '\n');
                    istringstream lineStream2(textLine);
                    int nodeID;
                    lineStream2 >> skipws >> nodeID;
                    assert(nodeIDToPosMap->find(nodeID) != nodeIDToPosMap->end());
                    for (int j = 0; j < dimension; j++) {
                        int nodePos = nodeIDToPosMap->at(nodeID);
                        lineStream2 >> data[nodePos * dimension + j];
                    }
                }
                getline(dotResFile, textLine, '\n');
                istringstream lineStream3(textLine);
                lineStream3 >> skipws >> lineToken;
                convertToLowerCase(lineToken);
                assert(lineToken == "end");

                return;
            } else {
                continue;
            }

        }
    }
    delete nodeIDToPosMap;
}

void readElementalDataFromDotRes(std::string fileName, std::string resultName,
        std::string analysisName, int stepNum, std::string type,
        int numberOfElements, const int *elemIds, const int *numberOfNodesPerElement,
        double *data) {
    ifstream dotResFile(fileName.c_str());
    assert(dotResFile);

    assert(resultName.size()!=0);
    if (resultName[0] != '\"')
        resultName.insert(0, "\"");
    if (resultName[resultName.size()-1] != '\"')
        resultName.append("\"");
    assert(analysisName.size()!=0);
    if (analysisName[0] != '\"')
        analysisName.insert(0, "\"");
    if (analysisName[analysisName.size()-1] != '\"')
        analysisName.append("\"");

    int dimension = 0;
    convertToLowerCase(type);
    if (type == "vector")
        dimension = 3;
    else if (type == "scalar")
        dimension = 1;
    else
        assert(false);

    map<int, int> *elemIDToPosMapTriangle = new map<int, int>;
    for (int i = 0; i < numberOfElements; i++) {
        if (numberOfNodesPerElement[i] == 3)
            elemIDToPosMapTriangle->insert(elemIDToPosMapTriangle->end(),
                    pair<int, int>(elemIds[i], i));
    }

    map<int, int> *elemIDToPosMapQuad = new map<int, int>;
    for (int i = 0; i < numberOfElements; i++) {
        if (numberOfNodesPerElement[i] == 4)
            elemIDToPosMapQuad->insert(elemIDToPosMapQuad->end(), pair<int, int>(elemIds[i], i));
    }

    bool hasTriangles = false;
    bool hasQuads = false;

    while (true) { /* parse line */
        if (!dotResFile.good()) { // terminate if failbit, eofbit or badbit is set
            if (!hasTriangles && !hasQuads) {
                cerr << "result with the following parameters is not found in " << fileName
                        << endl;
                cerr << "result name: " << resultName << endl;
                cerr << "analysis name: " << analysisName << endl;
                cerr << "step number: " << stepNum << endl;
                cerr << "type: " << type << endl;
                cerr << "Exit!" << endl;
                exit(EXIT_FAILURE);
            } else {
                return;
            }
        }
        string textLine;
        string lineToken;
        getline(dotResFile, textLine, '\n');

        istringstream lineStream(textLine);
        lineStream >> skipws >> lineToken;
        convertToLowerCase(lineToken);
        if (lineToken == "result") {
            string resultNameInFile;
            do {
                lineStream >> lineToken;
                resultNameInFile.append(" ").append(lineToken); // here a space is added to the beginning
            } while (resultNameInFile[resultNameInFile.size() - 1] != '\"'); // parse the name with spaces in
            resultNameInFile.erase(resultNameInFile.begin()); // remove the space at the beginning

            string analysisNameInFile;
            do {
                lineStream >> lineToken;
                analysisNameInFile.append(" ").append(lineToken); // here a space is added to the beginning
            } while (analysisNameInFile[analysisNameInFile.size() - 1] != '\"'); // parse the name with spaces in
            analysisNameInFile.erase(analysisNameInFile.begin()); // remove the space at the beginning

            int stepNumInFile;
            lineStream >> stepNumInFile;

            string typeInFile;
            lineStream >> typeInFile;
            convertToLowerCase(typeInFile);

            string locationInFile;
            lineStream >> locationInFile;
            convertToLowerCase(locationInFile);

            if (stepNumInFile > stepNum) { // save time if it is beyond the interesting time step
                if (!hasTriangles && !hasQuads) {
                    cerr << "result with the following parameters is not found in " << fileName
                            << endl;
                    cerr << "result name: " << resultName << endl;
                    cerr << "analysis name: " << analysisName << endl;
                    cerr << "step number: " << stepNum << endl;
                    cerr << "type: " << type << endl;
                    cerr << "Exit!" << endl;
                    exit(EXIT_FAILURE);
                } else {
                    return;
                }
            }

            if (resultName == resultNameInFile && analysisName == analysisNameInFile
                    && stepNum == stepNumInFile && type == typeInFile
                    && locationInFile == "ongausspoints") { // read nodal data, do not allow empty lines or comments
                string triangleOrQuad;
                lineStream >> triangleOrQuad;
                if (triangleOrQuad == "\"GP_triangles\"") {
                    getline(dotResFile, textLine, '\n');
                    istringstream lineStream1(textLine);
                    lineStream1 >> skipws >> lineToken;
                    convertToLowerCase(lineToken);
                    assert(lineToken == "values");
                    for (int i = 0; i < elemIDToPosMapTriangle->size(); i++) {
                        getline(dotResFile, textLine, '\n');
                        istringstream lineStream2(textLine);
                        int elemID;
                        lineStream2 >> skipws >> elemID;
                        assert(
                                elemIDToPosMapTriangle->find(elemID) != elemIDToPosMapTriangle->end());
                        for (int j = 0; j < dimension; j++) {
                            int elemPos = elemIDToPosMapTriangle->at(elemID);
                            lineStream2 >> data[elemPos * dimension + j];
                        }
                    }
                    getline(dotResFile, textLine, '\n');
                    istringstream lineStream3(textLine);
                    lineStream3 >> skipws >> lineToken;
                    convertToLowerCase(lineToken);
                    assert(lineToken == "end");
                    hasTriangles = true;
                } else if (triangleOrQuad == "\"GP_quads\"") {
                    getline(dotResFile, textLine, '\n');
                    istringstream lineStream1(textLine);
                    lineStream1 >> skipws >> lineToken;
                    convertToLowerCase(lineToken);
                    assert(lineToken == "values");
                    for (int i = 0; i < elemIDToPosMapQuad->size(); i++) {
                        getline(dotResFile, textLine, '\n');
                        istringstream lineStream2(textLine);
                        int elemID;
                        lineStream2 >> skipws >> elemID;
                        assert(elemIDToPosMapQuad->find(elemID) != elemIDToPosMapQuad->end());
                        for (int j = 0; j < dimension; j++) {
                            int elemPos = elemIDToPosMapQuad->at(elemID);
                            lineStream2 >> data[elemPos * dimension + j];
                        }
                    }
                    getline(dotResFile, textLine, '\n');
                    istringstream lineStream3(textLine);
                    lineStream3 >> skipws >> lineToken;
                    convertToLowerCase(lineToken);
                    assert(lineToken == "end");
                    hasQuads = true;
                } else {
                    assert(false);
                }
            } else {
                continue;
            }

        }
    }
    delete elemIDToPosMapTriangle;
    delete elemIDToPosMapQuad;
}

} /* namespace GiDFileIO */
