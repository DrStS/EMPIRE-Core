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
#include "stdlib.h"

#include "GmshFileIO.h"
#include "meshIO.h"

namespace EMPIRE {

using namespace std;

void GmshFileIO::readDotMsh(string fileName, int &numberOfMeshNodes, int &numberOfElements, double *&meshNodeCoordinates,
			    int *&meshNodeIds, int *&numberOfNodesPerElement, int *&elementNodeTables, int *&elementIds) {
  ifstream inputStream(fileName.c_str());
  if (!inputStream) {
    cerr << "GmshFileIO::readDotMsh: Gmsh .msh file \"" << fileName << "\" cannot be found"
	 << '\n';
    exit(EXIT_FAILURE);
  }
  const int XYZ = 3;
  bool inNodes = false, inElements = false;
  int numberOfNodes;
  int numberOfNodesThisElement = 0;
  map<int, MeshIO::Node *> nodeAtId;
  list<MeshIO::Node *> meshNodes;
  list<MeshIO::Element *> elements;
  string textLine, lineToken;
    
  while (true) { /* parse line */
    getline(inputStream, textLine, '\n');
    // terminate if failbit, eofbit or badbit is set
    if (!inputStream.good())
      break;

    istringstream lineStream(textLine);
    lineStream >> skipws >> lineToken;

    if (lineToken.size() == 0) { 
      // ignore empty lines
    }
    else if (inNodes) { 
      if (lineToken.compare("$EndNodes") == 0) {
	inNodes = false;
      } 
      else { // add one node
	int id;
	istringstream lineStream(textLine);
	lineStream >> id;
	double *coordinate = new double[XYZ];
	for (int i = 0; i < XYZ; i++) {
	  lineStream >> coordinate[i];
	}
	assert(nodeAtId.find(id) == nodeAtId.end()); // enforce unique node ids
	nodeAtId[id] = new MeshIO::Node(id, coordinate);
      }
    }
    else if (inElements) { // add element
      if (lineToken.compare("$EndElements") == 0) {
	inElements = false;
      } 
      else { // add one element
	istringstream lineStream(textLine);
	int id, type, numTags, tag;
	lineStream >> id >> type >> numTags;

	for (int i = 0; i < numTags; i++) // ignore tags
	  lineStream >> tag;

	switch (type) {
	case 2:
	  numberOfNodesThisElement = 3; // triangle
	  break;
	case 3:
	  numberOfNodesThisElement = 4; // quad
	  break;
	default:
	  assert(false);
	}

	int *nodeIdsThisElement = new int[numberOfNodesThisElement];
	for (int i = 0; i < numberOfNodesThisElement; i++) { 
	  lineStream >> nodeIdsThisElement[i];
	}
	elements.push_back(new MeshIO::Element(id, nodeIdsThisElement, numberOfNodesThisElement));
      }
    }
    else if (lineToken.compare("$MeshFormat") == 0) {
      getline(inputStream, textLine, '\n');
      istringstream lineStream(textLine);
      double versionNumber;
      int fileType, dataSize;

      lineStream >> versionNumber >> fileType >> dataSize;
      assert(versionNumber == 2.2);
      assert(fileType == 0); // ascii file format
      assert(dataSize == 8); // dataSize = sizeof(double)
    } 
    else if (lineToken.compare("$Nodes") == 0) {
      inNodes = true;
      getline(inputStream, textLine, '\n');
      istringstream lineStream(textLine);
      lineStream >> numberOfNodes;
    }
    else if (lineToken.compare("$Elements") == 0) {
      inElements = true;
      getline(inputStream, textLine, '\n');
      istringstream lineStream(textLine);
      lineStream >> numberOfElements;
    }
  } /* parse line */

  assert(nodeAtId.size() == numberOfNodes);
  assert(elements.size() == numberOfElements);
  assert(inNodes == false);
  assert(inElements == false);

  { // number of nodes per element
    int i = 0;
    numberOfNodesPerElement = new int[numberOfElements];
    list<MeshIO::Element *>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it, ++i)
      numberOfNodesPerElement[i] = (*it)->numberOfNodes;
  }

  { // element Ids
    int i = 0;
    elementIds = new int[numberOfElements];
    list<MeshIO::Element *>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it, ++i)
      elementIds[i] = (*it)->id;
  }

  { // element tables
    int size = 0;
    list<MeshIO::Element *>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it) {
      size += (*it)->numberOfNodes;
    }
    elementNodeTables = new int[size];

    int i_g = 0;
    for (it = elements.begin(); it != elements.end(); ++it) {
      for (int i=0; i < (*it)->numberOfNodes; i++)
	elementNodeTables[i_g++] = (*it)->nodeIds[i];
    }
  }

  { // get mesh nodes
    int id;
    list<MeshIO::Element *>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it) {
      for (int i=0; i < (*it)->numberOfNodes; i++) {
	id = (*it)->nodeIds[i];
	assert(nodeAtId.find(id) != nodeAtId.end()); // node must exist
	meshNodes.push_back(nodeAtId[id]);
      }
    }
    meshNodes.sort();  // remove duplicates
    meshNodes.unique();
    meshNodes.sort(MeshIO::compare_node_by_id);
  }

  numberOfMeshNodes = meshNodes.size();

  { // node coordinates
    int i_g = 0;
    meshNodeCoordinates = new double[numberOfMeshNodes * XYZ];
    list<MeshIO::Node *>::iterator it;
    for (it = meshNodes.begin(); it != meshNodes.end(); ++it) 
      for (int i=0; i < XYZ; i++) 
	meshNodeCoordinates[i_g++] = (*it)->coordinate[i];
  }

  { // node Ids
    int i_g = 0;
    meshNodeIds = new int[numberOfMeshNodes];
    list<MeshIO::Node *>::iterator it;
    for (it = meshNodes.begin(); it != meshNodes.end(); ++it) 
      meshNodeIds[i_g++] = (*it)->id;
  }

  { // delete nodes
    map<int, MeshIO::Node *>::iterator it;
    MeshIO::Node *node;
    for (it = nodeAtId.begin(); it != nodeAtId.end(); ++it)
      delete (it->second);
    nodeAtId.clear();
    meshNodes.clear();
  }

  { // delete elements
    list<MeshIO::Element *>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it) 
      delete (*it);
    elements.clear();
  }
}

void GmshFileIO::writeDotMsh(string fileName, int numberOfMeshNodes, int numberOfElements, double *&meshNodeCoordinates,
			     int *meshNodeIds, int *numberOfNodesPerElement, int *elementNodeTables, int *elementIds) {
  const int XYZ = 3;
  list<MeshIO::Node *> meshNodes;
  list<MeshIO::Element *> elements;
  
  { // extract elements
    int i_g = 0;
    for (int i=0; i<numberOfElements; i++) {
      int *elementNodeTable = new int[numberOfNodesPerElement[i]];
      for (int j = 0; j < numberOfNodesPerElement[i]; j++)
	elementNodeTable[j] = elementNodeTables[i_g++];
      
      elements.push_back(new MeshIO::Element(elementIds[i], elementNodeTable, numberOfNodesPerElement[i]));
    }
  }

  { // extract nodes
    for (int i=0; i<numberOfMeshNodes; i++) {
      double *coordinate = new double[XYZ];
      for (int j=0; j<XYZ; j++)
	coordinate[j] = meshNodeCoordinates[i*XYZ + j];

      meshNodes.push_back(new MeshIO::Node(meshNodeIds[i], coordinate));
    }
  }

  // order nodes by id
  meshNodes.sort(MeshIO::compare_node_by_id);

  // order elements by number of nodes
  elements.sort(MeshIO::compare_elem_by_num_node);

  ofstream outputStream(fileName.c_str());
  if (!outputStream) {
    cerr << "GmshFileIO::writeDotMsh: Gmsh .msh file \"" << fileName << "\" cannot be created" << endl;
    exit(EXIT_FAILURE);
  }
    
  // write header
  outputStream << "$MeshFormat"    << endl;
  outputStream << "2.2 0 8"        << endl;
  outputStream << "$EndMeshFormat" << endl;

  { // write nodes
    list<MeshIO::Node *>::iterator it;
    outputStream << "$Nodes" << endl;
    outputStream << numberOfMeshNodes << endl;
    
    for (it = meshNodes.begin(); it != meshNodes.end(); it++) {
      outputStream << (*it)->id;
      for (int i=0; i < XYZ; i++) 
	outputStream << '\t' << (*it)->coordinate[i];
      outputStream << endl;
    }
    
    outputStream << "$EndNodes" << endl;
  }
  
  { // write elements
    list<MeshIO::Element *>::iterator it;
    double type;
    int numberOfNodesThisElement;
    int numberOfTags = 0;
    outputStream << "$Elements" << endl;
    outputStream << numberOfElements << endl;

    for (it = elements.begin(); it != elements.end(); it++) {
      if ((*it)->numberOfNodes == 3)
	type = 2;
      else if ((*it)->numberOfNodes == 4)
	type = 3;
      else
	assert(false);

      outputStream << (*it)->id << '\t' << type << '\t' << numberOfTags;
      numberOfNodesThisElement = (*it)->numberOfNodes;

      for (int i=0; i < numberOfNodesThisElement; i++)
	outputStream << '\t' << (*it)->nodeIds[i];
      outputStream << endl;
    }

    outputStream << "$EndElements" << endl;
  }

  { // delete elements
    list<MeshIO::Element *>::iterator it;
    for (it = elements.begin(); it != elements.end(); it++) 
      delete (*it);
  }
  { // delete nodes
    list<MeshIO::Node *>::iterator it;
    for (it = meshNodes.begin(); it != meshNodes.end(); it++)
      delete (*it);
  }
}

} /* namespace EMPIRE */
