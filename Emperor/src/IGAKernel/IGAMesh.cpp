/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Munich
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

IGAMesh::IGAMesh(std::string _name, int _numNodes) :
        AbstractMesh(_name, EMPIRE_Mesh_IGAMesh), numNodes(_numNodes), untrimmedNumNodes(-1) {
}

IGAMesh::~IGAMesh() {
    for (int i = 0; i < surfacePatches.size(); i++)
        delete surfacePatches[i];
}

IGAPatchSurface* IGAMesh::addPatch(int _pDegree, int _uNoKnots, double* _uKnotVector, int _qDegree,
        int _vNoKnots, double* _vKnotVector, int _uNoControlPoints, int _vNoControlPoints,
        double* _controlPointNet, int* _dofIndexNet) {

    std::string patchName = name + " Patch";
    int IDBasis = 0;

    int numCPs = _uNoControlPoints * _vNoControlPoints;
    IGAControlPoint **cpNet;
    cpNet = new IGAControlPoint*[numCPs];

    for (int i = 0; i < numCPs; i++) {
        if (_dofIndexNet[i] < numNodes && _dofIndexNet[i] >= 0)
            cpNet[i] = new IGAControlPoint(_dofIndexNet[i], &_controlPointNet[i * 4]);
        else {
            ERROR_OUT() << "DOF " << _dofIndexNet[i] << " has not been defined" << endl;
            exit(-1);
        }
    }

    surfacePatches.push_back(
            new IGAPatchSurface(IDBasis, _pDegree, _uNoKnots, _uKnotVector, _qDegree, _vNoKnots,
                    _vKnotVector, _uNoControlPoints, _vNoControlPoints, cpNet));
    return surfacePatches.back();
    
}

void IGAMesh::computeBoundingBox() {
    if (boundingBox.isComputed)
        return;
    boundingBox.xmin = surfacePatches[0]->getControlPointNet()[0]->getX();
    boundingBox.xmax = surfacePatches[0]->getControlPointNet()[0]->getX();
    boundingBox.ymin = surfacePatches[0]->getControlPointNet()[0]->getY();
    boundingBox.ymax = surfacePatches[0]->getControlPointNet()[0]->getY();
    boundingBox.zmin = surfacePatches[0]->getControlPointNet()[0]->getZ();
    boundingBox.zmax = surfacePatches[0]->getControlPointNet()[0]->getZ();

    for (int patchCount = 0; patchCount < surfacePatches.size(); patchCount++) {
        IGAPatchSurface* patch = surfacePatches[patchCount];
        for (int cpCount = 0; cpCount < patch->getNoControlPoints(); cpCount++) {
            double x = patch->getControlPointNet()[cpCount]->getX();
            double y = patch->getControlPointNet()[cpCount]->getY();
            double z = patch->getControlPointNet()[cpCount]->getZ();
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
    }

    boundingBox.isComputed = true;
}

void IGAMesh::addDataField(string _dataFieldName, EMPIRE_DataField_location _location,
        EMPIRE_DataField_dimension _dimension, EMPIRE_DataField_typeOfQuantity _typeOfQuantity) {

    int numLocations = -1;
    if (_location == EMPIRE_DataField_atNode)
        numLocations = numNodes;
    else
        assert(false);

    assert(nameToDataFieldMap.find(_dataFieldName) == nameToDataFieldMap.end());
    DataField *dataField = new DataField(_dataFieldName, _location, numLocations, _dimension,
            _typeOfQuantity);
    nameToDataFieldMap.insert(pair<string, DataField*>(_dataFieldName, dataField));

}

int IGAMesh::getUntrimmedNumNodes() {
	if(untrimmedNumNodes==-1){
		std::set<int> CPids;
		for(int i=0;i<surfacePatches.size();i++){
			surfacePatches[i]->getUntrimmedCPindexes(CPids);
		}
		untrimmedNumNodes=CPids.size();
	}
	return untrimmedNumNodes;
}

Message &operator<<(Message & _message, const IGAMesh & _mesh) {
    _message <<endl;
    _message << "\t" << "---------------------------------Start Mesh" << endl;
    _message << "\t" << "IGA Mesh name: " << _mesh.name << endl;
    _message << "\t\tNumber of Patches: " << _mesh.getSurfacePatches().size() << endl;
    //_message << "\t" << "---------------------------------" << endl;
    for (int k = 0; k < _mesh.getSurfacePatches().size(); k++) {
        _message << "\t" << "---------------------------------Start Patch" << endl;
        _message << "\tPatch[" << k << "]:" << endl;
        _message << *(_mesh.getSurfacePatches()[k]);
//        _message << "\t\t\tpDegree:  "
//                << _mesh.getSurfacePatches()[k]->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
//                << endl;
//        _message << "\t\t\tqDegree:  "
//                << _mesh.getSurfacePatches()[k]->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
//                << endl;
//        _message << "\t\t\tKnots Vector U: \t";
//
//        for (int i = 0;
//                i < _mesh.getSurfacePatches()[k]->getIGABasis()->getUBSplineBasis1D()->getNoKnots();
//                i++)
//            _message
//                    << _mesh.getSurfacePatches()[k]->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i]
//                    << "  ";
//        _message << endl;
//
//        _message << "\t\t\tKnots Vector V: \t";
//        for (int i = 0;
//                i < _mesh.getSurfacePatches()[k]->getIGABasis()->getVBSplineBasis1D()->getNoKnots();
//                i++)
//            _message
//                    << _mesh.getSurfacePatches()[k]->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i]
//                    << "  ";
//        _message << endl;
//
//        _message << "\t\t\t" << "number of control points U: "
//                << _mesh.getSurfacePatches()[k]->getUNoControlPoints() << endl;
//        _message << "\t\t\t" << "number of control points V: "
//                << _mesh.getSurfacePatches()[k]->getVNoControlPoints() << endl;
//
//        _message << "\t\t\tControl Points Net: " << endl;
//        int count = 0;
//        for (int i = 0; i < _mesh.getSurfacePatches()[k]->getUNoControlPoints(); i++) {
//            _message << "\t\t\t";
//            for (int j = 0; j < _mesh.getSurfacePatches()[k]->getVNoControlPoints(); j++) {
//                _message << _mesh.getSurfacePatches()[k]->getControlPointNet()[count]->getX()
//                        << ", " << _mesh.getSurfacePatches()[k]->getControlPointNet()[count]->getY()
//                        << ", " << _mesh.getSurfacePatches()[k]->getControlPointNet()[count]->getZ()
//                        << "\t \t";
//                count++;
//            }
//            _message << endl;
//        }
    }
    _message << "\t" << "---------------------------------End Mesh" << endl;
    return _message;
}

}/* namespace EMPIRE */

