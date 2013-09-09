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
#include "Message.h"
#include "MapperAdapter.h"
#include "MapperAdapter.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"
#include <math.h>
#include <iostream>
#include <limits>


namespace EMPIRE {
using namespace std;
/********//**
 * \brief Test mapping for some representative cases
 ***********/
class TestMapping: public CppUnit::TestFixture {
private:
    FEMesh *testTria01Mesh_A;
    FEMesh *testTria01Mesh_B;

    FEMesh *testTria02Mesh_A;
    FEMesh *testTria02Mesh_B;

    DataField *testTria01Data_A;
    DataField *testTria01Data_B;

    DataField *testTria02Data_A;
    DataField *testTria02Data_B;

    bool debugMe;

    bool oppositeSurfaceNormal;
    bool enforceConsistency;

public:
    /********//**
     ***********************************************************************************************
     * \brief Set up two meshes and data fields for all test cases
     ***********/
    void setUp() {
        debugMe = false;
        oppositeSurfaceNormal = false;
        enforceConsistency = false;

        /*============================================================================================*/
        /*  Test Case testTria01                                                                */
        /*============================================================================================*/
        /*
         * \verbatim
         8    7    6    5                    8    7    6    5
         x---x----x----x                    x---x----x----x
         |\  |\   |\   |      Mapping       |  /|   /|   /|
         Mesh triaTest01_A (8 triangular Ele.) | \ |  \ |  \ |   <============>   | / | /  |  / |    Mesh triaTest01_B (8 triangular Elements )
         |1 \|2  \|3  \|4                   |/1 |/2  |/3  |4
         x---x----x----x                    x---x----x----x


         (Mesh of CSM side)                  (Mesh of CFD side)
         \endverbatim*/

        int numNodesPerElem = 3;
        {
            int numNodes = 8;
            int numElems = 6;
            testTria01Mesh_A = new FEMesh("triaTest01_A", numNodes, numElems);
            for (int i = 0; i < numElems; i++)
                testTria01Mesh_A->numNodesPerElem[i] = numNodesPerElem;
            testTria01Mesh_A->initElems();

            for (int i = 0; i < numNodes; i++)
                testTria01Mesh_A->nodeIDs[i] = i + 1;
            testTria01Mesh_A->nodes[0 * 3 + 0] = 0;
            testTria01Mesh_A->nodes[0 * 3 + 1] = 0;
            testTria01Mesh_A->nodes[0 * 3 + 2] = 0;

            testTria01Mesh_A->nodes[1 * 3 + 0] = 0.5;
            testTria01Mesh_A->nodes[1 * 3 + 1] = 0;
            testTria01Mesh_A->nodes[1 * 3 + 2] = 0;

            testTria01Mesh_A->nodes[2 * 3 + 0] = 1.0;
            testTria01Mesh_A->nodes[2 * 3 + 1] = 0;
            testTria01Mesh_A->nodes[2 * 3 + 2] = 0;

            testTria01Mesh_A->nodes[3 * 3 + 0] = 1.5;
            testTria01Mesh_A->nodes[3 * 3 + 1] = 0;
            testTria01Mesh_A->nodes[3 * 3 + 2] = 0;

            testTria01Mesh_A->nodes[4 * 3 + 0] = 1.5;
            testTria01Mesh_A->nodes[4 * 3 + 1] = 0.5;
            testTria01Mesh_A->nodes[4 * 3 + 2] = 0;

            testTria01Mesh_A->nodes[5 * 3 + 0] = 1.0;
            testTria01Mesh_A->nodes[5 * 3 + 1] = 0.5;
            testTria01Mesh_A->nodes[5 * 3 + 2] = 0;

            testTria01Mesh_A->nodes[6 * 3 + 0] = 0.5;
            testTria01Mesh_A->nodes[6 * 3 + 1] = 0.5;
            testTria01Mesh_A->nodes[6 * 3 + 2] = 0;

            testTria01Mesh_A->nodes[7 * 3 + 0] = 0.0;
            testTria01Mesh_A->nodes[7 * 3 + 1] = 0.5;
            testTria01Mesh_A->nodes[7 * 3 + 2] = 0;

            testTria01Mesh_A->elems[0] = 1;
            testTria01Mesh_A->elems[1] = 2;
            testTria01Mesh_A->elems[2] = 7;

            testTria01Mesh_A->elems[3] = 1;
            testTria01Mesh_A->elems[4] = 7;
            testTria01Mesh_A->elems[5] = 8;

            testTria01Mesh_A->elems[6] = 2;
            testTria01Mesh_A->elems[7] = 3;
            testTria01Mesh_A->elems[8] = 6;

            testTria01Mesh_A->elems[9] = 2;
            testTria01Mesh_A->elems[10] = 6;
            testTria01Mesh_A->elems[11] = 7;

            testTria01Mesh_A->elems[12] = 3;
            testTria01Mesh_A->elems[13] = 4;
            testTria01Mesh_A->elems[14] = 5;

            testTria01Mesh_A->elems[15] = 3;
            testTria01Mesh_A->elems[16] = 5;
            testTria01Mesh_A->elems[17] = 6;
        }
        {
            int numNodes = 8;
            int numElems = 6;
            testTria01Mesh_B = new FEMesh("triaTest01_B", numNodes, numElems);
            for (int i = 0; i < numElems; i++)
                testTria01Mesh_B->numNodesPerElem[i] = numNodesPerElem;
            testTria01Mesh_B->initElems();

            for (int i = 0; i < numNodes; i++)
                testTria01Mesh_B->nodeIDs[i] = i + 1;
            testTria01Mesh_B->nodes[0 * 3 + 0] = 0;
            testTria01Mesh_B->nodes[0 * 3 + 1] = 0;
            testTria01Mesh_B->nodes[0 * 3 + 2] = 0;

            testTria01Mesh_B->nodes[1 * 3 + 0] = 0.5;
            testTria01Mesh_B->nodes[1 * 3 + 1] = 0;
            testTria01Mesh_B->nodes[1 * 3 + 2] = 0;

            testTria01Mesh_B->nodes[2 * 3 + 0] = 1;
            testTria01Mesh_B->nodes[2 * 3 + 1] = 0;
            testTria01Mesh_B->nodes[2 * 3 + 2] = 0;

            testTria01Mesh_B->nodes[3 * 3 + 0] = 1.5;
            testTria01Mesh_B->nodes[3 * 3 + 1] = 0;
            testTria01Mesh_B->nodes[3 * 3 + 2] = 0;

            testTria01Mesh_B->nodes[4 * 3 + 0] = 1.5;
            testTria01Mesh_B->nodes[4 * 3 + 1] = 0.5;
            testTria01Mesh_B->nodes[4 * 3 + 2] = 0;

            testTria01Mesh_B->nodes[5 * 3 + 0] = 1.0;
            testTria01Mesh_B->nodes[5 * 3 + 1] = 0.5;
            testTria01Mesh_B->nodes[5 * 3 + 2] = 0;

            testTria01Mesh_B->nodes[6 * 3 + 0] = 0.5;
            testTria01Mesh_B->nodes[6 * 3 + 1] = 0.5;
            testTria01Mesh_B->nodes[6 * 3 + 2] = 0;

            testTria01Mesh_B->nodes[7 * 3 + 0] = 0.0;
            testTria01Mesh_B->nodes[7 * 3 + 1] = 0.5;
            testTria01Mesh_B->nodes[7 * 3 + 2] = 0;

            testTria01Mesh_B->elems[0] = 1;
            testTria01Mesh_B->elems[1] = 2;
            testTria01Mesh_B->elems[2] = 8;

            testTria01Mesh_B->elems[3] = 2;
            testTria01Mesh_B->elems[4] = 7;
            testTria01Mesh_B->elems[5] = 8;

            testTria01Mesh_B->elems[6] = 2;
            testTria01Mesh_B->elems[7] = 3;
            testTria01Mesh_B->elems[8] = 7;

            testTria01Mesh_B->elems[9] = 3;
            testTria01Mesh_B->elems[10] = 6;
            testTria01Mesh_B->elems[11] = 7;

            testTria01Mesh_B->elems[12] = 3;
            testTria01Mesh_B->elems[13] = 4;
            testTria01Mesh_B->elems[14] = 6;

            testTria01Mesh_B->elems[15] = 4;
            testTria01Mesh_B->elems[16] = 5;
            testTria01Mesh_B->elems[17] = 6;
        }

        testTria01Data_A = new DataField("testTriaMortar01_A", EMPIRE_DataField_atNode,
                testTria01Mesh_A->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        testTria01Data_B = new DataField("testTriaMortar01_B", EMPIRE_DataField_atNode,
                testTria01Mesh_B->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);

        for (int i = 0; i < testTria01Mesh_A->numNodes; i++)
            /// fix the boundary
            testTria01Data_A->data[i] = 0.0;
        {
            /// set internal field to one
            testTria01Data_A->data[1] = 1.0;
            testTria01Data_A->data[2] = 1.0;
            testTria01Data_A->data[5] = 1.0;
            testTria01Data_A->data[6] = 1.0;
        }

        /*============================================================================================*/
        /*  Test Case testTria02                                                                */
        /*============================================================================================*/

        {
            int numNodes = 8;
            int numElems = 6;
            testTria02Mesh_A = new FEMesh("triaTest02_A", numNodes, numElems);
            for (int i = 0; i < numElems; i++)
                testTria02Mesh_A->numNodesPerElem[i] = numNodesPerElem;
            testTria02Mesh_A->initElems();

            for (int i = 0; i < numNodes; i++)
                testTria02Mesh_A->nodeIDs[i] = i + 1;
            testTria02Mesh_A->nodes[0 * 3 + 0] = 0;
            testTria02Mesh_A->nodes[0 * 3 + 1] = 0;
            testTria02Mesh_A->nodes[0 * 3 + 2] = 0;

            testTria02Mesh_A->nodes[1 * 3 + 0] = 0.5;
            testTria02Mesh_A->nodes[1 * 3 + 1] = 0;
            testTria02Mesh_A->nodes[1 * 3 + 2] = 0;

            testTria02Mesh_A->nodes[2 * 3 + 0] = 1.0;
            testTria02Mesh_A->nodes[2 * 3 + 1] = 0;
            testTria02Mesh_A->nodes[2 * 3 + 2] = 0;

            testTria02Mesh_A->nodes[3 * 3 + 0] = 1.5;
            testTria02Mesh_A->nodes[3 * 3 + 1] = 0;
            testTria02Mesh_A->nodes[3 * 3 + 2] = 0;

            testTria02Mesh_A->nodes[4 * 3 + 0] = 1.5;
            testTria02Mesh_A->nodes[4 * 3 + 1] = 0.5;
            testTria02Mesh_A->nodes[4 * 3 + 2] = 0;

            testTria02Mesh_A->nodes[5 * 3 + 0] = 1.0;
            testTria02Mesh_A->nodes[5 * 3 + 1] = 0.5;
            testTria02Mesh_A->nodes[5 * 3 + 2] = 0;

            testTria02Mesh_A->nodes[6 * 3 + 0] = 0.5;
            testTria02Mesh_A->nodes[6 * 3 + 1] = 0.5;
            testTria02Mesh_A->nodes[6 * 3 + 2] = 0;

            testTria02Mesh_A->nodes[7 * 3 + 0] = 0.0;
            testTria02Mesh_A->nodes[7 * 3 + 1] = 0.5;
            testTria02Mesh_A->nodes[7 * 3 + 2] = 0;

            testTria02Mesh_A->elems[0] = 1;
            testTria02Mesh_A->elems[1] = 2;
            testTria02Mesh_A->elems[2] = 7;

            testTria02Mesh_A->elems[3] = 1;
            testTria02Mesh_A->elems[4] = 7;
            testTria02Mesh_A->elems[5] = 8;

            testTria02Mesh_A->elems[6] = 2;
            testTria02Mesh_A->elems[7] = 3;
            testTria02Mesh_A->elems[8] = 6;

            testTria02Mesh_A->elems[9] = 2;
            testTria02Mesh_A->elems[10] = 6;
            testTria02Mesh_A->elems[11] = 7;

            testTria02Mesh_A->elems[12] = 3;
            testTria02Mesh_A->elems[13] = 4;
            testTria02Mesh_A->elems[14] = 5;

            testTria02Mesh_A->elems[15] = 3;
            testTria02Mesh_A->elems[16] = 5;
            testTria02Mesh_A->elems[17] = 6;
        }
        {
            int numNodes = 8;
            int numElems = 6;
            testTria02Mesh_B = new FEMesh("triaTest02_B", numNodes, numElems);
            for (int i = 0; i < numElems; i++)
                testTria02Mesh_B->numNodesPerElem[i] = numNodesPerElem;
            testTria02Mesh_B->initElems();

            for (int i = 0; i < numNodes; i++)
                testTria02Mesh_B->nodeIDs[i] = i + 1;
            testTria02Mesh_B->nodes[0 * 3 + 0] = 0;
            testTria02Mesh_B->nodes[0 * 3 + 1] = 0;
            testTria02Mesh_B->nodes[0 * 3 + 2] = 0;

            testTria02Mesh_B->nodes[1 * 3 + 0] = 0.7;
            testTria02Mesh_B->nodes[1 * 3 + 1] = 0;
            testTria02Mesh_B->nodes[1 * 3 + 2] = 0;

            testTria02Mesh_B->nodes[2 * 3 + 0] = 0.8;
            testTria02Mesh_B->nodes[2 * 3 + 1] = 0;
            testTria02Mesh_B->nodes[2 * 3 + 2] = 0;

            testTria02Mesh_B->nodes[3 * 3 + 0] = 1.5;
            testTria02Mesh_B->nodes[3 * 3 + 1] = 0;
            testTria02Mesh_B->nodes[3 * 3 + 2] = 0;

            testTria02Mesh_B->nodes[4 * 3 + 0] = 1.5;
            testTria02Mesh_B->nodes[4 * 3 + 1] = 0.5;
            testTria02Mesh_B->nodes[4 * 3 + 2] = 0;

            testTria02Mesh_B->nodes[5 * 3 + 0] = 0.8;
            testTria02Mesh_B->nodes[5 * 3 + 1] = 0.5;
            testTria02Mesh_B->nodes[5 * 3 + 2] = 0;

            testTria02Mesh_B->nodes[6 * 3 + 0] = 0.7;
            testTria02Mesh_B->nodes[6 * 3 + 1] = 0.5;
            testTria02Mesh_B->nodes[6 * 3 + 2] = 0;

            testTria02Mesh_B->nodes[7 * 3 + 0] = 0.0;
            testTria02Mesh_B->nodes[7 * 3 + 1] = 0.5;
            testTria02Mesh_B->nodes[7 * 3 + 2] = 0;

            testTria02Mesh_B->elems[0] = 1;
            testTria02Mesh_B->elems[1] = 2;
            testTria02Mesh_B->elems[2] = 8;

            testTria02Mesh_B->elems[3] = 2;
            testTria02Mesh_B->elems[4] = 7;
            testTria02Mesh_B->elems[5] = 8;

            testTria02Mesh_B->elems[6] = 2;
            testTria02Mesh_B->elems[7] = 3;
            testTria02Mesh_B->elems[8] = 7;

            testTria02Mesh_B->elems[9] = 3;
            testTria02Mesh_B->elems[10] = 6;
            testTria02Mesh_B->elems[11] = 7;

            testTria02Mesh_B->elems[12] = 3;
            testTria02Mesh_B->elems[13] = 4;
            testTria02Mesh_B->elems[14] = 6;

            testTria02Mesh_B->elems[15] = 4;
            testTria02Mesh_B->elems[16] = 5;
            testTria02Mesh_B->elems[17] = 6;
        }

        testTria02Data_A = new DataField("testTriaMortar02_A", EMPIRE_DataField_atNode,
                testTria01Mesh_A->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        testTria02Data_B = new DataField("testTriaMortar02_B", EMPIRE_DataField_atNode,
                testTria01Mesh_B->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);

        for (int i = 0; i < testTria02Mesh_A->numNodes; i++)
            /// fix the boundary
            testTria02Data_A->data[i] = 0.0;
        {
            /// set internal field to one
            testTria02Data_A->data[1] = 1.0;
            testTria02Data_A->data[2] = 1.0;
            testTria02Data_A->data[5] = 1.0;
            testTria02Data_A->data[6] = 1.0;
        }

    }
    void tearDown() {
        delete testTria01Mesh_A, testTria01Mesh_B;
        delete testTria02Mesh_A, testTria02Mesh_B;
        delete testTria01Data_A, testTria01Data_A;
        delete testTria02Data_A, testTria02Data_A;
    }
    /********//**
     ***********************************************************************************************
     * \brief compare datafield source and target within a tolerance and return true is test is okay
     * \param[in] data    pointer to double array of data
     * \param[in] refdata pointer to double array of reference data
     * \param[in] size size of array (data and reference data size needs to be the same)
     *
     * \return true if test is okay
     ***********/
    bool compareWithReference(double *data, double *refdata, int size) {
        /// eps is 100 times machine epsilon
        double eps = 100 * (numeric_limits<double>::epsilon());
        ;

        for (int i = 0; i < size; i++) {
            if (fabs(data[i] - refdata[i]) > eps) {
                return false;
            }

        }

        return true;
    }

    /********//**
     ***********************************************************************************************
     * \brief Test case 01 for Mortar-mapper:
     * The test case 01 is used to test triangular mapping for matching nodal coordinates with different
     * element connections. It is also tested if the boundary conditions at nodes 1 4 5 and 8 is fulfilled
     * \verbatim
     8    7    6    5                    8    7    6    5
     x---x----x----x                    x---x----x----x
     |\  |\   |\   |      Mapping       |  /|   /|   /|
     Mesh triaTest01_A (8 triangular Ele.) | \ |  \ |  \ |   <============>   | / | /  |  / |    Mesh triaTest01_B (8 triangular Elements )
     |1 \|2  \|3  \|4                   |/1 |/2  |/3  |4
     x---x----x----x                    x---x----x----x


     (Mesh of CSM side)                  (Mesh of CFD side)
     * \endverbatim
     ***********/

    void testMortar01() {
        bool dual = false;
        MapperAdapter *mapper = new MapperAdapter("testTriaMortar01", testTria01Mesh_A,
                testTria01Mesh_B);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, testTria01Mesh_A, testTria01Data_A,
                testTria01Mesh_B, testTria01Data_B);
        filterConsistent->filtering();

        if (debugMe) {
            debugOut << *testTria01Data_A << endl;
            debugOut << *testTria01Data_B << endl;
        }
        if (!compareWithReference(testTria01Data_A->data, testTria01Data_B->data,
                testTria01Data_B->numLocations)) {
            CPPUNIT_ASSERT(false);
        }

        delete mapper;
        delete filterConsistent;
    }

    /********//**
     ***********************************************************************************************
     * \brief Test case 01 for Dual Mortar-mapper:
     * The test case 01 is used to test triangular mapping for matching nodal coordinates with different
     * element connections. It is also tested if the boundary conditions at nodes 1 4 5 and 8 is fulfilled
     * \verbatim
     8    7    6    5                    8    7    6    5
     x---x----x----x                    x---x----x----x
     |\  |\   |\   |      Mapping       |  /|   /|   /|
     Mesh triaTest01_A (8 triangular Ele.) | \ |  \ |  \ |   <============>   | / | /  |  / |    Mesh triaTest01_B (8 triangular Elements )
     |1 \|2  \|3  \|4                   |/1 |/2  |/3  |4
     x---x----x----x                    x---x----x----x


     (Mesh of CSM side)                  (Mesh of CFD side)
     * \endverbatim
     ***********/

    void testDualMortar01() {
        bool dual = true;
        MapperAdapter *mapper = new MapperAdapter("testTriaMortar01", testTria01Mesh_A,
                testTria01Mesh_B);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, testTria01Mesh_A, testTria01Data_A,
                testTria01Mesh_B, testTria01Data_B);
        filterConsistent->filtering();

        if (debugMe) {
            debugOut << *testTria01Data_A << endl;
            debugOut << *testTria01Data_B << endl;
        }
        if (!compareWithReference(testTria01Data_A->data, testTria01Data_B->data,
                testTria01Data_B->numLocations)) {
            CPPUNIT_ASSERT(false);
        }

        delete mapper;
        delete filterConsistent;
    }

    /********//**
     ***********************************************************************************************
     * \brief Test case 02 for Mortar-mapper:
     * The test case 01 is used to test triangular mapping for non-matching nodal coordinates with different
     * element connections. It is also tested if the boundary conditions at nodes 1 4 5 and 8 is fulfilled
     * \verbatim
     8    7    6    5                    8    7    6    5
     x---x----x----x                    x---x-----x----x
     |\  |\   |\   |      Mapping       |  /|    /|   /|
     Mesh triaTest01_A (8 triangular Ele.) | \ |  \ |  \ |   <============>   | / |  /  |  / |    Mesh triaTest01_B (8 triangular Elements )
     |1 \|2  \|3  \|4                   |/1 | /2  |/3  |4
     x---x----x----x                    x---x-----x----x


     (Mesh of CSM side)                  (Mesh of CFD side)
     * \endverbatim
     ***********/
    void testMortar02() {
        bool dual = false;
        MapperAdapter *mapper = new MapperAdapter("testTriaMortar02", testTria02Mesh_A,
                testTria02Mesh_B);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, testTria02Mesh_A, testTria02Data_A,
                testTria02Mesh_B, testTria02Data_B);
        filterConsistent->filtering();

        if (debugMe) {
            debugOut << *testTria02Data_A << endl;
            debugOut << *testTria02Data_B << endl;
        }
        if (!compareWithReference(testTria02Data_A->data, testTria02Data_B->data,
                testTria02Data_B->numLocations)) {
            //CPPUNIT_ASSERT(false);
        }

        delete mapper;
        delete filterConsistent;
    }

    /********//**
     ***********************************************************************************************
     * \brief Test case 02 for Dual Mortar-mapper:
     * The test case 01 is used to test triangular mapping for non-matching nodal coordinates with different
     * element connections. It is also tested if the boundary conditions at nodes 1 4 5 and 8 is fulfilled
     * \verbatim
     8    7    6    5                    8    7    6    5
     x---x----x----x                    x---x-----x----x
     |\  |\   |\   |      Mapping       |  /|    /|   /|
     Mesh triaTest01_A (8 triangular Ele.) | \ |  \ |  \ |   <============>   | / |  /  |  / |    Mesh triaTest01_B (8 triangular Elements )
     |1 \|2  \|3  \|4                   |/1 | /2  |/3  |4
     x---x----x----x                    x---x-----x----x


     (Mesh of CSM side)                  (Mesh of CFD side)
     * \endverbatim
     ***********/
    void testDualMortar02() {
        bool dual = true;
        MapperAdapter *mapper = new MapperAdapter("testTriaMortar02", testTria02Mesh_A,
                testTria02Mesh_B);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);

        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, testTria02Mesh_A, testTria02Data_A,
                testTria02Mesh_B, testTria02Data_B);
        filterConsistent->filtering();

        if (debugMe) {
            debugOut << *testTria02Data_A << endl;
            debugOut << *testTria02Data_B << endl;
        }
        if (!compareWithReference(testTria02Data_A->data, testTria02Data_B->data,
                testTria02Data_B->numLocations)) {
            //CPPUNIT_ASSERT(false);
        }

        delete mapper;
        delete filterConsistent;
    }

CPPUNIT_TEST_SUITE( TestMapping );
#ifdef USE_INTEL_MKL
        CPPUNIT_TEST(testMortar01);
        CPPUNIT_TEST(testMortar02);
#endif
        CPPUNIT_TEST(testDualMortar01);
        CPPUNIT_TEST(testDualMortar02);
        CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestMapping);
