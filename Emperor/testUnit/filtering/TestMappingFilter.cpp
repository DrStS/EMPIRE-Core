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
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

#include "MappingFilter.h"
#include "DataField.h"
#include "FEMesh.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"
#include "AbstractMapper.h"
#include "MapperAdapter.h"
#include "Message.h"
#include <math.h>
#include <iostream>

#ifdef USE_INTEL_MKL
#define LAPACK
#endif

#ifdef LAPACK_GCC
#define LAPACK
#endif

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test the class MappingFilter and the mappers (consistent mapping and conservative mapping)
 ***********/
class TestMappingFilter: public CppUnit::TestFixture {
private:
    FEMesh *meshA;
    FEMesh *meshB;
    DataField *a1;
    DataField *b1;
    DataField *a2;
    DataField *b2;
    static const double EPS = 1E-8;

public:
    /***********************************************************************************************
     * \brief Set up two meshes and data fields on both meshes. In fact, we do consistent mapping
     *        from a1 to b1, do conservative mapping (forces to forces) from b2 to a2, and do conservative
     *        mapping (pressures to forces) from b3 to a3.
     ***********/
    void setUp() {
        int numNodesPerElem = 4;
        /*
         * 4-------3
         * |       | mesh A
         * 1-------2
         */
        int numNodesA = 4;
        int numElemsA = 1;
        meshA = new FEMesh("", numNodesA, numElemsA);
        for (int i = 0; i < numElemsA; i++)
            meshA->numNodesPerElem[i] = numNodesPerElem;
        meshA->initElems();

        for (int i = 0; i < numNodesA; i++)
            meshA->nodeIDs[i] = i + 1;
        meshA->nodes[0 * 3 + 0] = 0;
        meshA->nodes[0 * 3 + 1] = 0;
        meshA->nodes[0 * 3 + 2] = 0;
        meshA->nodes[1 * 3 + 0] = 2;
        meshA->nodes[1 * 3 + 1] = 0;
        meshA->nodes[1 * 3 + 2] = 0;
        meshA->nodes[2 * 3 + 0] = 2;
        meshA->nodes[2 * 3 + 1] = 1;
        meshA->nodes[2 * 3 + 2] = 0;
        meshA->nodes[3 * 3 + 0] = 0;
        meshA->nodes[3 * 3 + 1] = 1;
        meshA->nodes[3 * 3 + 2] = 0;
        meshA->elems[0] = 1;
        meshA->elems[1] = 2;
        meshA->elems[2] = 3;
        meshA->elems[3] = 4;
        /*
         * 6---5---4
         * |   |   | mesh B
         * 1---2---3
         */
        int numNodesB = 6;
        int numElemsB = 2;
        meshB = new FEMesh("", numNodesB, numElemsB);
        for (int i = 0; i < numElemsB; i++)
            meshB->numNodesPerElem[i] = numNodesPerElem;
        meshB->initElems();

        for (int i = 0; i < numNodesB; i++)
            meshB->nodeIDs[i] = i + 1;
        meshB->nodes[0 * 3 + 0] = 0;
        meshB->nodes[0 * 3 + 1] = 0;
        meshB->nodes[0 * 3 + 2] = 0;
        meshB->nodes[1 * 3 + 0] = 1;
        meshB->nodes[1 * 3 + 1] = 0;
        meshB->nodes[1 * 3 + 2] = 0;
        meshB->nodes[2 * 3 + 0] = 2;
        meshB->nodes[2 * 3 + 1] = 0;
        meshB->nodes[2 * 3 + 2] = 0;
        meshB->nodes[3 * 3 + 0] = 2;
        meshB->nodes[3 * 3 + 1] = 1;
        meshB->nodes[3 * 3 + 2] = 0;
        meshB->nodes[4 * 3 + 0] = 1;
        meshB->nodes[4 * 3 + 1] = 1;
        meshB->nodes[4 * 3 + 2] = 0;
        meshB->nodes[5 * 3 + 0] = 0;
        meshB->nodes[5 * 3 + 1] = 1;
        meshB->nodes[5 * 3 + 2] = 0;
        meshB->elems[0 * 4 + 0] = 1;
        meshB->elems[0 * 4 + 1] = 2;
        meshB->elems[0 * 4 + 2] = 5;
        meshB->elems[0 * 4 + 3] = 6;
        meshB->elems[1 * 4 + 0] = 2;
        meshB->elems[1 * 4 + 1] = 3;
        meshB->elems[1 * 4 + 2] = 4;
        meshB->elems[1 * 4 + 3] = 5;

        a1 = new DataField("a1", EMPIRE_DataField_atNode, meshA->numNodes, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        b1 = new DataField("b1", EMPIRE_DataField_atNode, meshB->numNodes, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        a2 = new DataField("a2", EMPIRE_DataField_atNode, meshA->numNodes, EMPIRE_DataField_scalar,
                EMPIRE_DataField_fieldIntegral);
        b2 = new DataField("b2", EMPIRE_DataField_atNode, meshB->numNodes, EMPIRE_DataField_scalar,
                EMPIRE_DataField_fieldIntegral);
        for (int i = 0; i < meshA->numNodes; i++)
            a1->data[i] = 1.0;
        for (int i = 0; i < meshB->numNodes; i++)
            b2->data[i] = 1.0;
        b2->data[1] = 2.0;
        b2->data[4] = 2.0; // on nodes 2 and 5, the forces are equal to 2.0
    }
    void tearDown() {
        delete meshA;
        delete meshB;
        delete a1, b1, a2, b2;
    }
    /***********************************************************************************************
     * \brief Test case: Test consistent and conservative property of all mappers
     ***********/
    void testMappers() {
        // <math.h> abs/fabs is for floating point numbers, <stdlib.h> abs is for integers!!!
        for (int i = 0; i < 2; i++) {
            bool dual;
            if (i == 0)
                dual = false;
            else if (i == 1)
                dual = true;
#ifndef USE_INTEL_MKL
            if (!dual)
                continue;
#endif
#ifndef LAPACK
            continue;
#endif
            bool oppositeSurfaceNormal = false;
            bool enforceConsistency = false;

            MapperAdapter *mapper = new MapperAdapter("testMortar", meshA, meshB);
            mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);

            {
                AbstractFilter *filterConsistant = new MappingFilter(mapper);
                ConnectionIOSetup::setupIOForFilter(filterConsistant, meshA, a1, meshB, b1);
                filterConsistant->filtering();
                for (int i = 0; i < meshB->numNodes; i++)
                    CPPUNIT_ASSERT(fabs(b1->data[i] - 1.0) < EPS);
                delete filterConsistant;
            }
            {
                AbstractFilter *filterConservative = new MappingFilter(mapper);
                ConnectionIOSetup::setupIOForFilter(filterConservative, meshB, b2, meshA, a2);
                filterConservative->filtering();
                for (int i = 0; i < meshA->numNodes; i++)
                    CPPUNIT_ASSERT(fabs(a2->data[i] - 2.0) < EPS);
                delete filterConservative;
            }
            delete mapper;
        }
    }
CPPUNIT_TEST_SUITE( TestMappingFilter );
        CPPUNIT_TEST( testMappers);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestMappingFilter);
