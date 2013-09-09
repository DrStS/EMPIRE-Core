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
#include "MapperAdapter.h"
#include "MapperAdapter.h"
#include "Message.h"
#include <math.h>
#include <iostream>

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test mortar mapper. Also some test cases used to cause problems can be put here.
 ***********/
class TestMortarMapper: public CppUnit::TestFixture {
private:

public:
    void setUp() {

    }
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief compareC_BBandC_BA
     ***********/
    void compareC_BBandC_BA() {

    }
    /***********************************************************************************************
     * \brief rowsum of C_BA_dual is negative
     ***********/
    void problem1() {
        FEMesh *meshA;
        FEMesh *meshB;
        DataField *a1;
        DataField *b1;
        static const double EPS = 1E-8;
        int numNodesPerElem = 4;
        int numNodesA = 4;
        int numElemsA = 1;
        meshA = new FEMesh("", numNodesA, numElemsA);
        for (int i = 0; i < numElemsA; i++)
            meshA->numNodesPerElem[i] = numNodesPerElem;
        meshA->initElems();

        for (int i = 0; i < numNodesA; i++)
            meshA->nodeIDs[i] = i + 1;

        meshA->nodes[0 * 3 + 0] = 145.5446;
        meshA->nodes[0 * 3 + 1] = 3840.0;
        meshA->nodes[0 * 3 + 2] = 1.811688;

        meshA->nodes[1 * 3 + 0] = 170.3338;
        meshA->nodes[1 * 3 + 1] = 3840.0;
        meshA->nodes[1 * 3 + 2] = 2.97447;

        meshA->nodes[2 * 3 + 0] = 169.7244;
        meshA->nodes[2 * 3 + 1] = 3840.0;
        meshA->nodes[2 * 3 + 2] = 33.02235;

        meshA->nodes[3 * 3 + 0] = 145.5336;
        meshA->nodes[3 * 3 + 1] = 3840.0;
        meshA->nodes[3 * 3 + 2] = 30.04398;

        meshA->elems[0] = 1;
        meshA->elems[1] = 2;
        meshA->elems[2] = 3;
        meshA->elems[3] = 4;

        int numNodesB = 4;
        int numElemsB = 1;
        meshB = new FEMesh("", numNodesB, numElemsB);
        for (int i = 0; i < numElemsB; i++)
            meshB->numNodesPerElem[i] = numNodesPerElem;
        meshB->initElems();

        for (int i = 0; i < numNodesB; i++)
            meshB->nodeIDs[i] = i + 1;

        meshB->nodes[0 * 3 + 0] = 178.1;
        meshB->nodes[0 * 3 + 1] = 3840.0;
        meshB->nodes[0 * 3 + 2] = 0.4;

        meshB->nodes[1 * 3 + 0] = 111.4;
        meshB->nodes[1 * 3 + 1] = 3840.0;
        meshB->nodes[1 * 3 + 2] = -1.8;

        meshB->nodes[2 * 3 + 0] = 109.4;
        meshB->nodes[2 * 3 + 1] = 3840.0;
        meshB->nodes[2 * 3 + 2] = 26.6;

        meshB->nodes[3 * 3 + 0] = 176.8;
        meshB->nodes[3 * 3 + 1] = 3840.0;
        meshB->nodes[3 * 3 + 2] = 34.6;

        meshB->elems[0] = 1;
        meshB->elems[1] = 2;
        meshB->elems[2] = 3;
        meshB->elems[3] = 4;

        a1 = new DataField("a1", EMPIRE_DataField_atNode, meshA->numNodes, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        b1 = new DataField("b1", EMPIRE_DataField_atNode, meshB->numNodes, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        for (int i = 0; i < meshA->numNodes; i++)
            a1->data[i] = 1.0;

        bool oppositeSurfaceNormal = true;
        bool dual = true;
        bool enforceConsistency = true;
        MapperAdapter *mapper = new MapperAdapter("testMortar", meshA, meshB);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);

        AbstractFilter *filterConsistant = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistant, meshA, a1, meshB, b1);
        filterConsistant->filtering();
        for (int i = 0; i < meshB->numNodes; i++)
            CPPUNIT_ASSERT(fabs(b1->data[i] - 1.0) < EPS);

        delete filterConsistant;
        delete mapper;

        delete meshA;
        delete meshB;
        delete a1, b1;
    }
    /***********************************************************************************************
     * \brief Test the memory leak of the constructor by calling it 1,000,000 times
     *        This function should not be put into the test suite except when you really want to check
     *        the memory leak.
     ***********/
    void testMemoryLeakOfConstructor() {
        FEMesh *meshA;
        FEMesh *meshB;
        int numNodesPerElem = 4;
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

        meshA->nodes[1 * 3 + 0] = 1;
        meshA->nodes[1 * 3 + 1] = 0;
        meshA->nodes[1 * 3 + 2] = 0;

        meshA->nodes[2 * 3 + 0] = 2;
        meshA->nodes[2 * 3 + 1] = 2;
        meshA->nodes[2 * 3 + 2] = 0;

        meshA->nodes[3 * 3 + 0] = 0;
        meshA->nodes[3 * 3 + 1] = 1;
        meshA->nodes[3 * 3 + 2] = 0;

        meshA->elems[0] = 1;
        meshA->elems[1] = 2;
        meshA->elems[2] = 3;
        meshA->elems[3] = 4;

        int numNodesB = 4;
        int numElemsB = 1;
        meshB = new FEMesh("", numNodesB, numElemsB);
        for (int i = 0; i < numElemsB; i++)
            meshB->numNodesPerElem[i] = numNodesPerElem;
        meshB->initElems();

        for (int i = 0; i < numNodesB; i++)
            meshB->nodeIDs[i] = i + 1;

        meshB->nodes[0 * 3 + 0] = 0;
        meshB->nodes[0 * 3 + 1] = 0;
        meshB->nodes[0 * 3 + 2] = 0;

        meshB->nodes[1 * 3 + 0] = 1.5;
        meshB->nodes[1 * 3 + 1] = 0;
        meshB->nodes[1 * 3 + 2] = 0;

        meshB->nodes[2 * 3 + 0] = 1.5;
        meshB->nodes[2 * 3 + 1] = 1.5;
        meshB->nodes[2 * 3 + 2] = 0;

        meshB->nodes[3 * 3 + 0] = 0;
        meshB->nodes[3 * 3 + 1] = 1.5;
        meshB->nodes[3 * 3 + 2] = 0;

        meshB->elems[0] = 1;
        meshB->elems[1] = 2;
        meshB->elems[2] = 3;
        meshB->elems[3] = 4;

        bool oppositeSurfaceNormal = false;
        bool dual = true;
        bool enforceConsistency = true;

        for (int i = 0; i < 10000000; i++) {
            MapperAdapter *mapper = new MapperAdapter("testMortar", meshA, meshB);
            mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
            delete mapper;
        }

        delete meshA;
        delete meshB;
    }
CPPUNIT_TEST_SUITE( TestMortarMapper );
        CPPUNIT_TEST( compareC_BBandC_BA);
        CPPUNIT_TEST( problem1);
//CPPUNIT_TEST( testMemoryLeakOfConstructor); // test memory leak, comment it except when checking memory leak
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestMortarMapper);
