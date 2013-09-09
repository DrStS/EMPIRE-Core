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
#include "TriangulatorAdaptor.h"
#include <math.h>
#include <assert.h>
#include <map>
#include "polypartition.h"
#include <iostream>

using namespace std;

namespace EMPIRE {

TriangulatorAdaptor::TriangulatorAdaptor() {
}

TriangulatorAdaptor::~TriangulatorAdaptor() {
    for (int i = 0; i < polygon.size(); i++)
        delete[] polygon[i];
}

void TriangulatorAdaptor::addPoint(double x, double y, double z) {
    double *point = new double[3];
    point[0] = x;
    point[1] = y;
    point[2] = z;
    polygon.push_back(point);
}

void TriangulatorAdaptor::triangulate(int *triangleIndexes) {
    int XX, YY;
    from3DTo2D(XX, YY);

    map<pair<double, double> , int> coorToPosMap; // coordinate to position map
    for (int i = 0; i < polygon.size(); i++) {
        coorToPosMap.insert(
                pair<pair<double, double> , int>(
                        pair<double, double>(polygon[i][XX], polygon[i][YY]), i));
    }

    TPPLPoly poly;
    poly.Init(polygon.size());
    for (int i = 0; i < polygon.size(); i++) {
        poly[i].x = polygon[i][XX];
        poly[i].y = polygon[i][YY];
    }

    if (poly.GetOrientation() == TPPL_CW) {
        isClockwise = true;
        poly.Invert(); // polygon has to be counter-clockwise
    } else {
        isClockwise = false;
    }

    list<TPPLPoly> triangles;

    TPPLPartition triangulator; // see http://code.google.com/p/polypartition/
    int tmp = triangulator.Triangulate_OPT(&poly, &triangles);
    assert (tmp == 1);
    assert(triangles.size() +2 == polygon.size());
    int count = 0;
    for (list<TPPLPoly>::iterator it = triangles.begin(); it != triangles.end(); it++) {
        if (!isClockwise) {
            for (int j = 0; j < 3; j++) {
                triangleIndexes[count * 3 + j] = coorToPosMap[pair<double, double>((*it)[j].x,
                        (*it)[j].y)];
            }
        } else {
            for (int j = 0; j < 3; j++) {
                triangleIndexes[count * 3 + 2 - j] = coorToPosMap[pair<double, double>((*it)[j].x,
                        (*it)[j].y)];
            }
        }
        count++;
    }
    assert(count == triangles.size());
}

void TriangulatorAdaptor::from3DTo2D(int &XX, int &YY) {
    assert(polygon.size() >= 4);
    // determine whether to project to XY or XZ or YZ plane
    // by checking the size of projected triangle
    // the plane where the size of the projected triangle is the largest, then use that plane
    int flag = 0;
    double areaMax = 0.0;
    { // the 1st triangle
        double dx1 = polygon[0][0] - polygon[1][0];
        double dy1 = polygon[0][1] - polygon[1][1];
        double dz1 = polygon[0][2] - polygon[1][2];
        double dx2 = polygon[0][0] - polygon[2][0];
        double dy2 = polygon[0][1] - polygon[2][1];
        double dz2 = polygon[0][2] - polygon[2][2];
        double areaXY = fabs(dx1*dy2 - dx2*dy1);
        double areaYZ = fabs(dy1*dz2 - dy2*dz1);
        double areaXZ = fabs(dx1*dz2 - dx2*dz1);
        if (areaXY >= areaYZ && areaXY >= areaXZ) {
            flag = 2;
            areaMax = areaXY;
        } else if (areaYZ >= areaXY && areaYZ >= areaXZ) {
            flag = 0;
            areaMax = areaYZ;
        } else {
            assert(areaXZ >= areaXY && areaXZ >= areaYZ);
            flag = 1;
            areaMax = areaXZ;
        }
    }
    { // the 2nd triangle
        double dx1 = polygon[0][0] - polygon[1][0];
        double dy1 = polygon[0][1] - polygon[1][1];
        double dz1 = polygon[0][2] - polygon[1][2];
        double dx2 = polygon[0][0] - polygon[3][0];
        double dy2 = polygon[0][1] - polygon[3][1];
        double dz2 = polygon[0][2] - polygon[3][2];
        double areaXY = fabs(dx1*dy2 - dx2*dy1);
        double areaYZ = fabs(dy1*dz2 - dy2*dz1);
        double areaXZ = fabs(dx1*dz2 - dx2*dz1);
        if (areaXY > areaMax) {
            flag = 2;
        } else if (areaYZ > areaMax) {
            flag = 0;
        } else if (areaXZ > areaMax) {
            flag = 1;
        }
    }
    XX = (flag +1) % 3;
    YY = (flag +2) % 3;
}

} /* namespace EMPIRE */
