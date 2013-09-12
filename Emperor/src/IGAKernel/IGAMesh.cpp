// Inclusion of standard libraries
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// Inclusion of user defined libraries
#include "IGAPatchSurface.h"
#include "IGAControlPoint.h"
#include "IGAMesh.h"
#include "DataField.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

IGAMesh::IGAMesh(std::string _name, int _numControlPoints,
		double* _globalControlPoints, int* _controlPointID) :
		AbstractMesh(_name, EMPIRE_Mesh_IGAMesh), controlPointID(
				_controlPointID), numControlPoints(_numControlPoints) {

	cout << "numControlPoints: " << numControlPoints << endl;
	for (int i = 0; i < numControlPoints; i++) {
		globalControlPoints.push_back(new IGAControlPoint(controlPointID[i], &_globalControlPoints[i * 4]));
		cpMap.insert(pair<int, int>(controlPointID[i], i));
	}

}

IGAMesh::~IGAMesh() {

	for (int i = 0; i < globalControlPoints.size(); i++)
		delete globalControlPoints[i];
	for (int i = 0; i < surfacePatches.size(); i++)
		delete surfacePatches[i];

}

void IGAMesh::addPatch(int _pDegree, int _uNoKnots, double* _uKnotVector,
		int _qDegree, int _vNoKnots, double* _vKnotVector,
		int _uNoControlPoints, int _vNoControlPoints, int* _controlPointNetID) {

	std::string patchName = name + " Patch";
	int IDBasis = 0;

	int numCPs = _uNoControlPoints * _vNoControlPoints;
	IGAControlPoint **cpNet;
	cpNet = new IGAControlPoint*[numCPs];

	for (int i = 0; i < numCPs; i++){
		if (cpMap.find(_controlPointNetID[i]) == cpMap.end()){
			cout << "i = " << i << "   Control Point ID:" << _controlPointNetID[i] << endl;
			cout << "Control Point " << _controlPointNetID[i] << "doesn't exist! " << endl;
			assert(cpMap.find(_controlPointNetID[i]) != cpMap.end());
		}
		cpNet[i] = globalControlPoints[cpMap[_controlPointNetID[i]]];
	}

	surfacePatches.push_back(
			new IGAPatchSurface(patchName, IDBasis, _pDegree, _uNoKnots,
					_uKnotVector, _qDegree, _vNoKnots, _vKnotVector,
					_uNoControlPoints, _vNoControlPoints, cpNet));
//	surfacePatches[surfacePatches.size()-1]->printSelf();
}

void IGAMesh::computeBoundingBox() {
	if (boundingBox.isComputed)
		return;
	boundingBox.xmin = globalControlPoints[0]->getX();
	boundingBox.xmax = globalControlPoints[0]->getX();
	boundingBox.ymin = globalControlPoints[0]->getY();
	boundingBox.ymax = globalControlPoints[0]->getY();
	boundingBox.zmin = globalControlPoints[0]->getZ();
	boundingBox.zmax = globalControlPoints[0]->getZ();
	for (int i = 1; i < globalControlPoints.size(); i++) {
		double x = globalControlPoints[i]->getX();
		double y = globalControlPoints[i]->getY();
		double z = globalControlPoints[i]->getZ();
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

void IGAMesh::addDataField(string dataFieldName,
		EMPIRE_DataField_location location,
		EMPIRE_DataField_dimension dimension,
		EMPIRE_DataField_typeOfQuantity typeOfQuantity) {

	int numLocatiions = -1;
	if (location == EMPIRE_DataField_atNode)
		numLocatiions = globalControlPoints.size();
	else
		assert(false);

	assert(nameToDataFieldMap.find(dataFieldName) == nameToDataFieldMap.end());
	DataField *dataField = new DataField(dataFieldName, location, numLocatiions,
			dimension, typeOfQuantity);
	nameToDataFieldMap.insert(
			pair<string, DataField*>(dataFieldName, dataField));

}


Message &operator<<(Message &message, IGAMesh &mesh) {
	message << "\t" << "IGA Mesh name: " << mesh.name << endl;

	message << "\t\tNumber of Patches: " << mesh.getSurfacePatches().size() << endl;

	for (int k = 0; k < mesh.getSurfacePatches().size(); k++){

		message << "\t\t\tPatch" << k << ":" << endl;
		message << "\t\t\tpDegree:  "	<< mesh.getSurfacePatches()[k]->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree() << endl;
		message << "\t\t\tqDegree:  "	<< mesh.getSurfacePatches()[k]->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree() << endl;
		message << "\t\t\tKnots Vector U: \t";

		for (int i = 0; i < mesh.getSurfacePatches()[k]->getIGABasis()->getUBSplineBasis1D()->getNoKnots(); i++)
			message << mesh.getSurfacePatches()[k]->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i]	<< "  ";
		message << endl;

		message << "\t\t\tKnots Vector V: \t";
		for (int i = 0; i < mesh.getSurfacePatches()[k]->getIGABasis()->getVBSplineBasis1D()->getNoKnots();	i++)
			message << mesh.getSurfacePatches()[k]->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i]	<< "  ";
		message << endl;

		message << "\t\t\t" << "number of control points U: "	<< mesh.getSurfacePatches()[k]->getUNoControlPoints() << endl;
		message << "\t\t\t" << "number of control points V: "	<< mesh.getSurfacePatches()[k]->getVNoControlPoints() << endl;

		message << "\t\t\tControl Points Net: " << endl;
		int count = 0;
		for (int i = 0; i < mesh.getSurfacePatches()[k]->getUNoControlPoints(); i++) {
			cout << "\t\t\t";
			for (int j = 0; j < mesh.getSurfacePatches()[k]->getVNoControlPoints(); j++) {
				message << mesh.getSurfacePatches()[k]->getControlPointNet()[count]->getX() << ", "
						<< mesh.getSurfacePatches()[k]->getControlPointNet()[count]->getY() << ", "
						<< mesh.getSurfacePatches()[k]->getControlPointNet()[count]->getZ() << "\t";
				count++;
			}
			message << endl;
		}
		message() << "\t" << "---------------------------------" << endl;
	}
	return message;
}

}/* namespace EMPIRE */

