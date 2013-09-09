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
#include <iostream>
#include <stdlib.h>

#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

#include "MappingFilter.h"
#include "DataField.h"
#include "FEMesh.h"
#include "LocationFilter.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"

#ifdef USE_INTEL_MKL
#define LAPACK
#endif

#ifdef LAPACK_GCC
#define LAPACK
#endif

using namespace std;

namespace EMPIRE {
/********//**
 ***********************************************************************************************
 * \brief Test the class LocationFilter
 ***********/
class TestLocationFilter: public CppUnit::TestFixture {
private:
    FEMesh *mesh;
    DataField *atElemCentroidScalar;
    DataField *atNodeScalar;
    DataField *atElemCentroidVector;
    DataField *atNodeVector;

public:
    /********//**
     ************************************************************************************************
     * \brief Set up the mesh, create one scalar field (1), one vector field (1, 0, 2).
     *        Both fields locate on element centroid.
     ***********/
    void setUp() {
        /*
         * 6---5---4
         * |   | \ |  mesh
         * 1---2---3
         */
        int numNodes = 6;
        int numElems = 3;
        mesh = new FEMesh("", numNodes, numElems);
        mesh->numNodesPerElem[0] = 4;
        mesh->numNodesPerElem[1] = 3;
        mesh->numNodesPerElem[2] = 3;
        mesh->initElems();
        for (int i = 0; i < numNodes; i++)
            mesh->nodeIDs[i] = i + 1;
        mesh->nodes[0 * 3 + 0] = 0;
        mesh->nodes[0 * 3 + 1] = 0;
        mesh->nodes[0 * 3 + 2] = 0;
        mesh->nodes[1 * 3 + 0] = 1;
        mesh->nodes[1 * 3 + 1] = 0;
        mesh->nodes[1 * 3 + 2] = 0;
        mesh->nodes[2 * 3 + 0] = 2;
        mesh->nodes[2 * 3 + 1] = 0;
        mesh->nodes[2 * 3 + 2] = 0;
        mesh->nodes[3 * 3 + 0] = 2;
        mesh->nodes[3 * 3 + 1] = 1;
        mesh->nodes[3 * 3 + 2] = 0;
        mesh->nodes[4 * 3 + 0] = 1;
        mesh->nodes[4 * 3 + 1] = 1;
        mesh->nodes[4 * 3 + 2] = 0;
        mesh->nodes[5 * 3 + 0] = 0;
        mesh->nodes[5 * 3 + 1] = 1;
        mesh->nodes[5 * 3 + 2] = 0;

        mesh->elems[0] = 1;
        mesh->elems[1] = 2;
        mesh->elems[2] = 5;
        mesh->elems[3] = 6;

        mesh->elems[4] = 2;
        mesh->elems[5] = 3;
        mesh->elems[6] = 5;

        mesh->elems[7] = 3;
        mesh->elems[8] = 4;
        mesh->elems[9] = 5;

        atElemCentroidScalar = new DataField("b1", EMPIRE_DataField_atElemCentroid, mesh->numElems,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        atNodeScalar = new DataField("", EMPIRE_DataField_atNode, mesh->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);

        atElemCentroidVector = new DataField("", EMPIRE_DataField_atElemCentroid, mesh->numElems,
                EMPIRE_DataField_vector, EMPIRE_DataField_field);
        atNodeVector = new DataField("", EMPIRE_DataField_atNode, mesh->numNodes,
                EMPIRE_DataField_vector, EMPIRE_DataField_field);
        for (int i = 0; i < mesh->numElems; i++)
            atElemCentroidScalar->data[i] = 1.0;
        for (int i = 0; i < mesh->numElems; i++) {
            atElemCentroidVector->data[i * EMPIRE_DataField_vector + 0] = 1.0;
            atElemCentroidVector->data[i * EMPIRE_DataField_vector + 1] = 0.0;
            atElemCentroidVector->data[i * EMPIRE_DataField_vector + 2] = 2.0;
        }

    }
    void tearDown() {
        delete mesh;
        delete atElemCentroidScalar, atElemCentroidVector, atNodeScalar, atNodeVector;
    }
    /********//**
     ***********************************************************************************************
     * \brief Test case: Call ElemCentroidToNodeFilter to filter both fields to node.
     ***********/
    void testFiltering() {
        LocationFilter *filterScalar = new LocationFilter();
        ConnectionIOSetup::setupIOForFilter(filterScalar, mesh, atElemCentroidScalar, mesh,
                atNodeScalar);
        LocationFilter *filterVector = new LocationFilter();
        ConnectionIOSetup::setupIOForFilter(filterVector, mesh, atElemCentroidVector, mesh,
                atNodeVector);
        filterScalar->filtering();
        filterVector->filtering();
        for (int i = 0; i < mesh->numNodes; i++)
            CPPUNIT_ASSERT(atNodeScalar->data[i] == 1.0);
        for (int i = 0; i < mesh->numElems; i++) {
            CPPUNIT_ASSERT(atNodeVector->data[i* EMPIRE_DataField_vector + 0] == 1.0);
            CPPUNIT_ASSERT(atNodeVector->data[i* EMPIRE_DataField_vector + 1] == 0.0);
            CPPUNIT_ASSERT(atNodeVector->data[i* EMPIRE_DataField_vector + 2] == 2.0);
        }
        delete filterScalar;
        delete filterVector;
    }

CPPUNIT_TEST_SUITE( TestLocationFilter );
        CPPUNIT_TEST( testFiltering);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestLocationFilter);
