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
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <list>
#include <map>
#include <assert.h>
#include <vector>
#include "stdlib.h"

#include "MatlabIGAFileIO.h"
#include "IGAMesh.h"
#include "IGAPatchSurface.h"
#include "BSplineBasis2D.h"
#include "BSplineBasis1D.h"
#include "IGAControlPoint.h"
#include "DataField.h"

using namespace std;
namespace EMPIRE {
namespace MatlabIGAFileIO {
/********//**
 * \brief Class Node stores node data
 ***********/

void writeIGAMesh(IGAMesh* igaMesh) {
    // Matlab visualization

    ofstream myfile;

    myfile.open("SurfacePatches.m");
    myfile.precision(14);
    myfile << std::dec;

    std::vector<IGAPatchSurface*> surfacePatches = igaMesh->getSurfacePatches();
    myfile << "numNodesIGA =  " << igaMesh->getNumNodes() << ";" << endl;
    int numPatches = surfacePatches.size();
    IGAPatchSurface* patch;
    for (int patchCount = 0; patchCount < numPatches; patchCount++) {
        patch = surfacePatches[patchCount];
        myfile << "%% Patch: " << patchCount + 1 << endl;
        myfile << "surfacePatch(" << patchCount + 1 << ").p = "
                << patch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree() << ";" << endl;
        myfile << "surfacePatch(" << patchCount + 1 << ").q = "
                << patch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree() << ";" << endl;
        myfile << "surfacePatch(" << patchCount + 1 << ").Xi = [";

        for (int i = 0; i < patch->getIGABasis()->getUBSplineBasis1D()->getNoKnots(); i++)
            myfile << patch->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] << "  ";
        myfile << "];" << endl;
        myfile << "surfacePatch(" << patchCount + 1 << ").Eta = [";
        for (int i = 0; i < patch->getIGABasis()->getVBSplineBasis1D()->getNoKnots(); i++)
            myfile << patch->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] << "  ";
        myfile << "];" << endl;

        int cpCount = 0;
        for (int vCount = 0; vCount < patch->getVNoControlPoints(); vCount++) {
            for (int uCount = 0; uCount < patch->getUNoControlPoints(); uCount++) {
                IGAControlPoint* cp = patch->getControlPointNet()[cpCount];
                myfile << "surfacePatch(" << patchCount + 1 << ").CP(" << uCount + 1 << ", "
                        << vCount + 1 << ", 1:4) = [" << cp->getX() << ", " << cp->getY() << ", "
                        << cp->getZ() << ", " << cp->getW() << "]; " << endl;
                cpCount++;
            }
        }

        myfile << "surfacePatch(" << patchCount + 1 << ").dof = [";
        for (int i = 0; i < patch->getNoControlPoints(); i++)
            myfile << patch->getControlPointNet()[i]->getDofIndex() + 1 << " ";
        myfile << "];" << endl;
    }

    myfile.close();



}

void writeVectorFieldOnCPs(string _dataFieldName, int _step, DataField* _dataField) {
    ofstream myfile;

    string fileName = _dataFieldName+"IGA.m";
    myfile.open(fileName.c_str() , ios::app);
    myfile.precision(14);
    myfile << std::dec;
    myfile << "timestep:" << _step << "\n";
    for (int i = 0; i < _dataField->numLocations * 3; i++)
        myfile << _dataField->data[i] << "\n";
    myfile.close();

}
} /* namespace GiDFileIO */
} /* namespace EMPIRE */
