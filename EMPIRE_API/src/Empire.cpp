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
#include "Empire.h"
#include "ClientMetaDatabase.h"
#include "ClientCommunication.h"
#include <assert.h>
#include <string.h>
#include <iostream>
#include "EMPIRE_API.h"
using namespace std;

namespace EMPIRE {

Empire::Empire() {

}

Empire::~Empire() {
}

void Empire::initEnvironment(char *inputFileName) {
    ClientMetaDatabase::init(inputFileName);
}

string Empire::getUserDefinedText(string elementName) {
    return ClientMetaDatabase::getSingleton()->getUserDefinedText(elementName);
}

void Empire::connect() {
    ClientCommunication::getSingleton()->connect();
}

void Empire::disconnect() {
    ClientCommunication::getSingleton()->disconnect();
    //ClientMetaDatabase::getSingleton()->free();
    //ClientCommunication::getSingleton()->free();
}

void Empire::sendMesh(int numNodes, int numElems, double *nodes, int *nodeIDs, int *numNodesPerElem,
        int *elems) {
    const int BUFFER_SIZE = 2;
    const int DIMENSION = 3;
    int meshInfo[BUFFER_SIZE] = { numNodes, numElems };
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(BUFFER_SIZE, meshInfo);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(numNodes * DIMENSION, nodes);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(numNodes, nodeIDs);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(numElems, numNodesPerElem);
    int count = 0;
    for (int i = 0; i < numElems; i++)
        count += numNodesPerElem[i];
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(count, elems);
}

void Empire::sendIGAMesh(int _numPatches, int _numControlPoints, double* _globalControlPoints, int* _controlPointID){
	cout << "Sending IGA Mesh......................" << endl;
	const int BUFFER_SIZE = 2;
	int meshInfo[BUFFER_SIZE] = { _numPatches, _numControlPoints};
	ClientCommunication::getSingleton()->sendToServerBlocking<int>(BUFFER_SIZE,meshInfo);
	ClientCommunication::getSingleton()->sendToServerBlocking<double>(_numControlPoints * 4, _globalControlPoints);
	ClientCommunication::getSingleton()->sendToServerBlocking<int>(_numControlPoints, _controlPointID);
}

void Empire::sendIGAPatch(int _pDegree,	int _uNoKnots, double* _uKnotVector, int _qDegree, int _vNoKnots,
		double* _vKnotVector, int _uNoControlPoints, int _vNoControlPoints, int* _controlPointNetID) {
	cout << "sendIGAPatch....................." << endl;
    const int BUFFER_SIZE = 6;
    int meshInfo[BUFFER_SIZE] = { _pDegree, _uNoKnots, _qDegree, _vNoKnots,_uNoControlPoints, _vNoControlPoints };
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(BUFFER_SIZE,meshInfo);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_uNoKnots,_uKnotVector);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(_vNoKnots,_vKnotVector);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(_uNoControlPoints * _vNoControlPoints, _controlPointNetID);
}

void Empire::sendDataField(int sizeOfArray, double *dataField) {
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(1, &sizeOfArray);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(sizeOfArray, dataField);
}

void Empire::recvDataField(int sizeOfArray, double *dataField) {
    int sizeOfArrayRecv = 0;
    ClientCommunication::getSingleton()->receiveFromServerBlocking<int>(1, &sizeOfArrayRecv);
    assert(sizeOfArray == sizeOfArrayRecv);
    ClientCommunication::getSingleton()->receiveFromServerBlocking<double>(sizeOfArray, dataField);
}

void Empire::sendSignal_double(char *name, int sizeOfArray, double *signal) {
    char nameBuffer[EMPIRE_API_NAME_STRING_LENGTH];
    strcpy(nameBuffer, name);
    ClientCommunication::getSingleton()->sendToServerBlocking<char>(EMPIRE_API_NAME_STRING_LENGTH,
            nameBuffer);
    ClientCommunication::getSingleton()->sendToServerBlocking<int>(1, &sizeOfArray);
    ClientCommunication::getSingleton()->sendToServerBlocking<double>(sizeOfArray, signal);
}

void Empire::recvSignal_double(char *name, int sizeOfArray, double *signal) {
    char nameRecv[EMPIRE_API_NAME_STRING_LENGTH];
    ClientCommunication::getSingleton()->receiveFromServerBlocking<char>(
            EMPIRE_API_NAME_STRING_LENGTH, nameRecv);
    if (strcmp(name, nameRecv) != 0)
        assert(false);
    int sizeOfArrayRecv = 0;
    ClientCommunication::getSingleton()->receiveFromServerBlocking<int>(1, &sizeOfArrayRecv);
    assert(sizeOfArray == sizeOfArrayRecv);
    ClientCommunication::getSingleton()->receiveFromServerBlocking<double>(sizeOfArray, signal);
}

int Empire::recvConvergenceSignal() {
    int signal;
    ClientCommunication::getSingleton()->receiveFromServerBlocking<int>(1, &signal);
    return signal;
}

void Empire::printDataField(char *name, int sizeOfArray, double *dataField) {
    std::cout << name << std::endl;
    for (int i = 0; i < sizeOfArray; i++)
        std::cout << dataField[i] << '\t';
    std::cout << std::endl;
}
}/* namespace EMPIRE */

