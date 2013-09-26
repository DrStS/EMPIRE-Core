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

IGAMesh::IGAMesh(std::string _name, int _numControlPoints, double* _globalControlPoints,
        int* _controlPointID) :
        AbstractMesh(_name, EMPIRE_Mesh_IGAMesh), controlPointID(_controlPointID), numControlPoints(
                _numControlPoints) {

    for (int i = 0; i < numControlPoints; i++) {
        globalControlPoints.push_back(
                new IGAControlPoint(controlPointID[i], &_globalControlPoints[i * 4]));
        mapControlPointIDToIndex.insert(pair<int, int>(controlPointID[i], i));
    }

}

IGAMesh::~IGAMesh() {

    for (int i = 0; i < globalControlPoints.size(); i++)
        delete globalControlPoints[i];
    for (int i = 0; i < surfacePatches.size(); i++)
        delete surfacePatches[i];

}

void IGAMesh::addPatch(int _pDegree, int _uNoKnots, double* _uKnotVector, int _qDegree,
        int _vNoKnots, double* _vKnotVector, int _uNoControlPoints, int _vNoControlPoints,
        int* _controlPointNetID) {

    std::string patchName = name + " Patch";
    int IDBasis = 0;

    int numCPs = _uNoControlPoints * _vNoControlPoints;
    IGAControlPoint **cpNet;
    cpNet = new IGAControlPoint*[numCPs];

    for (int i = 0; i < numCPs; i++) {

        assert(
                mapControlPointIDToIndex.find(_controlPointNetID[i]) != mapControlPointIDToIndex.end());
        cpNet[i] = globalControlPoints[mapControlPointIDToIndex[_controlPointNetID[i]]];
    }

    surfacePatches.push_back(
            new IGAPatchSurface(IDBasis, _pDegree, _uNoKnots, _uKnotVector, _qDegree, _vNoKnots,
                    _vKnotVector, _uNoControlPoints, _vNoControlPoints, cpNet));
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

void IGAMesh::addDataField(string _dataFieldName, EMPIRE_DataField_location _location,
        EMPIRE_DataField_dimension _dimension, EMPIRE_DataField_typeOfQuantity _typeOfQuantity) {

    int numLocatiions = -1;
    if (_location == EMPIRE_DataField_atNode)
        numLocatiions = globalControlPoints.size();
    else
        assert(false);

    assert(nameToDataFieldMap.find(_dataFieldName) == nameToDataFieldMap.end());
    DataField *dataField = new DataField(_dataFieldName, _location, numLocatiions, _dimension,
            _typeOfQuantity);
    nameToDataFieldMap.insert(pair<string, DataField*>(_dataFieldName, dataField));

}

Message &operator<<(Message & _message, IGAMesh & _mesh) {
    _message << "\t" << "IGA Mesh name: " << _mesh.name << endl;

    _message << "\t\tNumber of Patches: " << _mesh.getSurfacePatches().size() << endl;

    for (int k = 0; k < _mesh.getSurfacePatches().size(); k++) {

        _message << "\t\t\tPatch" << k << ":" << endl;
        _message << "\t\t\tpDegree:  "
                << _mesh.getSurfacePatches()[k]->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
                << endl;
        _message << "\t\t\tqDegree:  "
                << _mesh.getSurfacePatches()[k]->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
                << endl;
        _message << "\t\t\tKnots Vector U: \t";

        for (int i = 0;
                i < _mesh.getSurfacePatches()[k]->getIGABasis()->getUBSplineBasis1D()->getNoKnots();
                i++)
            _message
                    << _mesh.getSurfacePatches()[k]->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i]
                    << "  ";
        _message << endl;

        _message << "\t\t\tKnots Vector V: \t";
        for (int i = 0;
                i < _mesh.getSurfacePatches()[k]->getIGABasis()->getVBSplineBasis1D()->getNoKnots();
                i++)
            _message
                    << _mesh.getSurfacePatches()[k]->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i]
                    << "  ";
        _message << endl;

        _message << "\t\t\t" << "number of control points U: "
                << _mesh.getSurfacePatches()[k]->getUNoControlPoints() << endl;
        _message << "\t\t\t" << "number of control points V: "
                << _mesh.getSurfacePatches()[k]->getVNoControlPoints() << endl;

        _message << "\t\t\tControl Points Net: " << endl;
        int count = 0;
        for (int i = 0; i < _mesh.getSurfacePatches()[k]->getUNoControlPoints(); i++) {
            cout << "\t\t\t";
            for (int j = 0; j < _mesh.getSurfacePatches()[k]->getVNoControlPoints(); j++) {
                _message << _mesh.getSurfacePatches()[k]->getControlPointNet()[count]->getX()
                        << ", " << _mesh.getSurfacePatches()[k]->getControlPointNet()[count]->getY()
                        << ", " << _mesh.getSurfacePatches()[k]->getControlPointNet()[count]->getZ()
                        << "\t";
                count++;
            }
            _message << endl;
        }
        _message() << "\t" << "---------------------------------" << endl;
    }
    return _message;
}

}/* namespace EMPIRE */

