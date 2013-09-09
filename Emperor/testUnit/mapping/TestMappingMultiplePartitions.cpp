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
#include <math.h>
#include <iostream>

#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

#include "MappingFilter.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"
#include "DataField.h"
#include "FEMesh.h"
#include "MapperAdapter.h"
#include "MapperAdapter.h"
#include "Message.h"
#include "GiDFileIO.h"

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test mapping on multiple partitions (output to shell to check the results instead of making assertions)
 *        (There is no problem when the data field is constant, otherwise, we can observe problems)
 *        There are four test cases.
 *        The 1st test case map on single-partition meshes, which gives the reference result.
 *        In the 2nd test case, mesh A has 2 partitions. The displacements on B are equal to the reference result.
 *        In the 3rd test case, mesh B has 3 partitions. The displacements on B among partitions are non-consistent (which is very wrong!).
 *        In the 4th test case, both meshes are partitioned. The displacements on B are wrong as in the 3rd test case.
 *        In all test cases, the sum of forces on A is equal to the reference result.
 *        The conclusion is that there is non-consistency among partitions when mapping displacements.
 ***********/
class TestMappingMultiplePartitions: public CppUnit::TestFixture {
private:
    FEMesh *meshQuadA;
    FEMesh *meshQuadB;
    FEMesh *meshQuadAMult;
    FEMesh *meshQuadBMult;
    bool debugMe;
    bool constantDataField;
    bool oppositeSurfaceNormal;
    bool dual;
    bool enforceConsistency;

public:
    /***********************************************************************************************
     * \brief Set up two meshes and data fields on both meshes. In fact, we do consistent mapping
     *        from df_A to df_B.
     ***********/
    void setUp() {
        debugMe = false;
        oppositeSurfaceNormal = false;
        dual = true;
        enforceConsistency = false;
        constantDataField = false;
        int numNodesPerElemQuad = 4;
        {
            /*
             * 4-----5-----6
             * |     |     | mesh quad A
             * 1-----2-----3
             */
            int numNodesQuad1 = 6;
            int numElemsQuad1 = 2;
            meshQuadA = new FEMesh("", numNodesQuad1, numElemsQuad1);
            for (int i = 0; i < numElemsQuad1; i++)
                meshQuadA->numNodesPerElem[i] = numNodesPerElemQuad;
            meshQuadA->initElems();

            for (int i = 0; i < numNodesQuad1; i++)
                meshQuadA->nodeIDs[i] = i + 1;

            meshQuadA->nodes[0 * 3 + 0] = 0;
            meshQuadA->nodes[0 * 3 + 1] = 0;
            meshQuadA->nodes[0 * 3 + 2] = 0;

            meshQuadA->nodes[1 * 3 + 0] = 1.5;
            meshQuadA->nodes[1 * 3 + 1] = 0;
            meshQuadA->nodes[1 * 3 + 2] = 0;

            meshQuadA->nodes[2 * 3 + 0] = 3;
            meshQuadA->nodes[2 * 3 + 1] = 0;
            meshQuadA->nodes[2 * 3 + 2] = 0;

            meshQuadA->nodes[3 * 3 + 0] = 0;
            meshQuadA->nodes[3 * 3 + 1] = 1;
            meshQuadA->nodes[3 * 3 + 2] = 0;

            meshQuadA->nodes[4 * 3 + 0] = 1.5;
            meshQuadA->nodes[4 * 3 + 1] = 1;
            meshQuadA->nodes[4 * 3 + 2] = 0;

            meshQuadA->nodes[5 * 3 + 0] = 3;
            meshQuadA->nodes[5 * 3 + 1] = 1;
            meshQuadA->nodes[5 * 3 + 2] = 0;

            meshQuadA->elems[0 * 4 + 0] = 1;
            meshQuadA->elems[0 * 4 + 1] = 2;
            meshQuadA->elems[0 * 4 + 2] = 5;
            meshQuadA->elems[0 * 4 + 3] = 4;

            meshQuadA->elems[1 * 4 + 0] = 2;
            meshQuadA->elems[1 * 4 + 1] = 3;
            meshQuadA->elems[1 * 4 + 2] = 6;
            meshQuadA->elems[1 * 4 + 3] = 5;
        }
        {
            /*       (8)
             * 4-----5-----6
             * |     |     | mesh quad A multiple partitions
             * 1-----2-----3
             *       (7)
             */
            int numNodesQuad1 = 8;
            int numElemsQuad1 = 2;
            meshQuadAMult = new FEMesh("", numNodesQuad1, numElemsQuad1);
            for (int i = 0; i < numElemsQuad1; i++)
                meshQuadAMult->numNodesPerElem[i] = numNodesPerElemQuad;
            meshQuadAMult->initElems();

            for (int i = 0; i < numNodesQuad1; i++)
                meshQuadAMult->nodeIDs[i] = i + 1;

            meshQuadAMult->nodes[0 * 3 + 0] = 0;
            meshQuadAMult->nodes[0 * 3 + 1] = 0;
            meshQuadAMult->nodes[0 * 3 + 2] = 0;

            meshQuadAMult->nodes[1 * 3 + 0] = 1.5;
            meshQuadAMult->nodes[1 * 3 + 1] = 0;
            meshQuadAMult->nodes[1 * 3 + 2] = 0;

            meshQuadAMult->nodes[2 * 3 + 0] = 3;
            meshQuadAMult->nodes[2 * 3 + 1] = 0;
            meshQuadAMult->nodes[2 * 3 + 2] = 0;

            meshQuadAMult->nodes[3 * 3 + 0] = 0;
            meshQuadAMult->nodes[3 * 3 + 1] = 1;
            meshQuadAMult->nodes[3 * 3 + 2] = 0;

            meshQuadAMult->nodes[4 * 3 + 0] = 1.5;
            meshQuadAMult->nodes[4 * 3 + 1] = 1;
            meshQuadAMult->nodes[4 * 3 + 2] = 0;

            meshQuadAMult->nodes[5 * 3 + 0] = 3;
            meshQuadAMult->nodes[5 * 3 + 1] = 1;
            meshQuadAMult->nodes[5 * 3 + 2] = 0;

            meshQuadAMult->nodes[6 * 3 + 0] = meshQuadAMult->nodes[1 * 3 + 0];
            meshQuadAMult->nodes[6 * 3 + 1] = meshQuadAMult->nodes[1 * 3 + 1];
            meshQuadAMult->nodes[6 * 3 + 2] = meshQuadAMult->nodes[1 * 3 + 2];

            meshQuadAMult->nodes[7 * 3 + 0] = meshQuadAMult->nodes[4 * 3 + 0];
            meshQuadAMult->nodes[7 * 3 + 1] = meshQuadAMult->nodes[4 * 3 + 1];
            meshQuadAMult->nodes[7 * 3 + 2] = meshQuadAMult->nodes[4 * 3 + 2];

            meshQuadAMult->elems[0 * 4 + 0] = 1;
            meshQuadAMult->elems[0 * 4 + 1] = 2;
            meshQuadAMult->elems[0 * 4 + 2] = 5;
            meshQuadAMult->elems[0 * 4 + 3] = 4;

            meshQuadAMult->elems[1 * 4 + 0] = 7;
            meshQuadAMult->elems[1 * 4 + 1] = 3;
            meshQuadAMult->elems[1 * 4 + 2] = 6;
            meshQuadAMult->elems[1 * 4 + 3] = 8;
        }
        {
            /*
             * 5---6---7---8
             * |   |   |   | mesh quad B
             * 1---2---3---4
             */
            int numNodesQuad2 = 8;
            int numElemsQuad2 = 3;
            meshQuadB = new FEMesh("", numNodesQuad2, numElemsQuad2);
            for (int i = 0; i < numElemsQuad2; i++)
                meshQuadB->numNodesPerElem[i] = numNodesPerElemQuad;
            meshQuadB->initElems();

            for (int i = 0; i < numNodesQuad2; i++)
                meshQuadB->nodeIDs[i] = i + 1;

            meshQuadB->nodes[0 * 3 + 0] = 0;
            meshQuadB->nodes[0 * 3 + 1] = 0;
            meshQuadB->nodes[0 * 3 + 2] = 0;

            meshQuadB->nodes[1 * 3 + 0] = 1;
            meshQuadB->nodes[1 * 3 + 1] = 0;
            meshQuadB->nodes[1 * 3 + 2] = 0;

            meshQuadB->nodes[2 * 3 + 0] = 2;
            meshQuadB->nodes[2 * 3 + 1] = 0;
            meshQuadB->nodes[2 * 3 + 2] = 0;

            meshQuadB->nodes[3 * 3 + 0] = 3;
            meshQuadB->nodes[3 * 3 + 1] = 0;
            meshQuadB->nodes[3 * 3 + 2] = 0;

            meshQuadB->nodes[4 * 3 + 0] = 0;
            meshQuadB->nodes[4 * 3 + 1] = 1;
            meshQuadB->nodes[4 * 3 + 2] = 0;

            meshQuadB->nodes[5 * 3 + 0] = 1;
            meshQuadB->nodes[5 * 3 + 1] = 1;
            meshQuadB->nodes[5 * 3 + 2] = 0;

            meshQuadB->nodes[6 * 3 + 0] = 2;
            meshQuadB->nodes[6 * 3 + 1] = 1;
            meshQuadB->nodes[6 * 3 + 2] = 0;

            meshQuadB->nodes[7 * 3 + 0] = 3;
            meshQuadB->nodes[7 * 3 + 1] = 1;
            meshQuadB->nodes[7 * 3 + 2] = 0;

            meshQuadB->elems[0 * 4 + 0] = 1;
            meshQuadB->elems[0 * 4 + 1] = 2;
            meshQuadB->elems[0 * 4 + 2] = 6;
            meshQuadB->elems[0 * 4 + 3] = 5;

            meshQuadB->elems[1 * 4 + 0] = 2;
            meshQuadB->elems[1 * 4 + 1] = 3;
            meshQuadB->elems[1 * 4 + 2] = 7;
            meshQuadB->elems[1 * 4 + 3] = 6;

            meshQuadB->elems[2 * 4 + 0] = 3;
            meshQuadB->elems[2 * 4 + 1] = 4;
            meshQuadB->elems[2 * 4 + 2] = 8;
            meshQuadB->elems[2 * 4 + 3] = 7;
        }
        {
            /*
             *   (11)  (12)
             * 5---6---7---8
             * |   |   |   | mesh quad B multiple partitions
             * 1---2---3---4
             *    (9) (10)
             */
            int numNodesQuad2 = 12;
            int numElemsQuad2 = 3;
            meshQuadBMult = new FEMesh("", numNodesQuad2, numElemsQuad2);
            for (int i = 0; i < numElemsQuad2; i++)
                meshQuadBMult->numNodesPerElem[i] = numNodesPerElemQuad;
            meshQuadBMult->initElems();

            for (int i = 0; i < numNodesQuad2; i++)
                meshQuadBMult->nodeIDs[i] = i + 1;

            meshQuadBMult->nodes[0 * 3 + 0] = 0;
            meshQuadBMult->nodes[0 * 3 + 1] = 0;
            meshQuadBMult->nodes[0 * 3 + 2] = 0;

            meshQuadBMult->nodes[1 * 3 + 0] = 1;
            meshQuadBMult->nodes[1 * 3 + 1] = 0;
            meshQuadBMult->nodes[1 * 3 + 2] = 0;

            meshQuadBMult->nodes[2 * 3 + 0] = 2;
            meshQuadBMult->nodes[2 * 3 + 1] = 0;
            meshQuadBMult->nodes[2 * 3 + 2] = 0;

            meshQuadBMult->nodes[3 * 3 + 0] = 3;
            meshQuadBMult->nodes[3 * 3 + 1] = 0;
            meshQuadBMult->nodes[3 * 3 + 2] = 0;

            meshQuadBMult->nodes[4 * 3 + 0] = 0;
            meshQuadBMult->nodes[4 * 3 + 1] = 1;
            meshQuadBMult->nodes[4 * 3 + 2] = 0;

            meshQuadBMult->nodes[5 * 3 + 0] = 1;
            meshQuadBMult->nodes[5 * 3 + 1] = 1;
            meshQuadBMult->nodes[5 * 3 + 2] = 0;

            meshQuadBMult->nodes[6 * 3 + 0] = 2;
            meshQuadBMult->nodes[6 * 3 + 1] = 1;
            meshQuadBMult->nodes[6 * 3 + 2] = 0;

            meshQuadBMult->nodes[7 * 3 + 0] = 3;
            meshQuadBMult->nodes[7 * 3 + 1] = 1;
            meshQuadBMult->nodes[7 * 3 + 2] = 0;

            meshQuadBMult->nodes[8 * 3 + 0] = meshQuadBMult->nodes[1 * 3 + 0];
            meshQuadBMult->nodes[8 * 3 + 1] = meshQuadBMult->nodes[1 * 3 + 1];
            meshQuadBMult->nodes[8 * 3 + 2] = meshQuadBMult->nodes[1 * 3 + 2];

            meshQuadBMult->nodes[9 * 3 + 0] = meshQuadBMult->nodes[2 * 3 + 0];
            meshQuadBMult->nodes[9 * 3 + 1] = meshQuadBMult->nodes[2 * 3 + 1];
            meshQuadBMult->nodes[9 * 3 + 2] = meshQuadBMult->nodes[2 * 3 + 2];

            meshQuadBMult->nodes[10 * 3 + 0] = meshQuadBMult->nodes[5 * 3 + 0];
            meshQuadBMult->nodes[10 * 3 + 1] = meshQuadBMult->nodes[5 * 3 + 1];
            meshQuadBMult->nodes[10 * 3 + 2] = meshQuadBMult->nodes[5 * 3 + 2];

            meshQuadBMult->nodes[11 * 3 + 0] = meshQuadBMult->nodes[6 * 3 + 0];
            meshQuadBMult->nodes[11 * 3 + 1] = meshQuadBMult->nodes[6 * 3 + 1];
            meshQuadBMult->nodes[11 * 3 + 2] = meshQuadBMult->nodes[6 * 3 + 2];

            meshQuadBMult->elems[0 * 4 + 0] = 1;
            meshQuadBMult->elems[0 * 4 + 1] = 2;
            meshQuadBMult->elems[0 * 4 + 2] = 6;
            meshQuadBMult->elems[0 * 4 + 3] = 5;

            meshQuadBMult->elems[1 * 4 + 0] = 9;
            meshQuadBMult->elems[1 * 4 + 1] = 3;
            meshQuadBMult->elems[1 * 4 + 2] = 7;
            meshQuadBMult->elems[1 * 4 + 3] = 11;

            meshQuadBMult->elems[2 * 4 + 0] = 10;
            meshQuadBMult->elems[2 * 4 + 1] = 4;
            meshQuadBMult->elems[2 * 4 + 2] = 8;
            meshQuadBMult->elems[2 * 4 + 3] = 12;
        }
    }
    void tearDown() {
        delete meshQuadA;
        delete meshQuadB;
        delete meshQuadAMult;
        delete meshQuadBMult;
    }
    /***********************************************************************************************
     * \brief Test case 1: mapping from meshQuadA to meshQuadB
     ***********/
    void testQuadAToQuadB() {
        DataField *df_A1;
        DataField *df_B1;
        df_A1 = new DataField("testQuadAToQuadB-df_A1(ref)", EMPIRE_DataField_atNode,
                meshQuadA->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        df_B1 = new DataField("testQuadAToQuadB-df_B1(ref)", EMPIRE_DataField_atNode,
                meshQuadB->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        if (constantDataField) { // df_A1 is constant
            df_A1->data[0] = 1;
            df_A1->data[1] = 1;
            df_A1->data[2] = 1;
            df_A1->data[3] = 1;
            df_A1->data[4] = 1;
            df_A1->data[5] = 1;
        } else { // df_A1 is symmetric along x-direction
            df_A1->data[0] = 0;
            df_A1->data[1] = 1;
            df_A1->data[2] = 0;
            df_A1->data[3] = 0;
            df_A1->data[4] = 1;
            df_A1->data[5] = 0;
        }
        DataField *df_A2;
        DataField *df_B2;
        df_A2 = new DataField("testQuadAToQuadB-df_A2(ref)", EMPIRE_DataField_atNode,
                meshQuadA->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        df_B2 = new DataField("testQuadAToQuadB-df_B2(ref)", EMPIRE_DataField_atNode,
                meshQuadB->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        if (constantDataField) { // df_B2 is integrated from constant
            df_B2->data[0] = 1;
            df_B2->data[1] = 2;
            df_B2->data[2] = 2;
            df_B2->data[3] = 1;
            df_B2->data[4] = 1;
            df_B2->data[5] = 2;
            df_B2->data[6] = 2;
            df_B2->data[7] = 1;
        } else { // df_B2 is symmetric along x-direction
            df_B2->data[0] = 0;
            df_B2->data[1] = 1;
            df_B2->data[2] = 4;
            df_B2->data[3] = 9;
            df_B2->data[4] = 0;
            df_B2->data[5] = 1;
            df_B2->data[6] = 4;
            df_B2->data[7] = 9;
        }
        MapperAdapter *mapper = new MapperAdapter("testQuadAToQuadB", meshQuadA, meshQuadB);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        { // consistent mapping
            AbstractFilter *filterConsistent = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConsistent, meshQuadA, df_A1, meshQuadB,
                    df_B1);
            filterConsistent->filtering();
            if (debugMe)
                infoOut << *df_B1;
            delete filterConsistent;
        }
        { // conservative mapping
            AbstractFilter *filterConservative = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConservative, meshQuadB, df_B2, meshQuadA,
                    df_A2);
            filterConservative->filtering();
            if (debugMe)
                infoOut << *df_A2;
            delete filterConservative;
        }
        delete df_A1, df_B1, df_A2, df_B2;
        delete mapper;
    }
    /***********************************************************************************************
     * \brief Test case 2: mapping from meshQuadAMult to meshQuadB
     ***********/
    void testQuadAMultToQuadB() {
        DataField *df_A1;
        DataField *df_B1;
        df_A1 = new DataField("testQuadAMultToQuadB-df_A1", EMPIRE_DataField_atNode,
                meshQuadAMult->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        df_B1 = new DataField("testQuadAMultToQuadB-df_B1", EMPIRE_DataField_atNode,
                meshQuadB->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        if (constantDataField) { // df_A1 is constant
            df_A1->data[0] = 1;
            df_A1->data[1] = 1;
            df_A1->data[2] = 1;
            df_A1->data[3] = 1;
            df_A1->data[4] = 1;
            df_A1->data[5] = 1;
            df_A1->data[6] = df_A1->data[1];
            df_A1->data[7] = df_A1->data[4];
        } else { // df_A1 is symmetric along x-direction
            df_A1->data[0] = 0;
            df_A1->data[1] = 1;
            df_A1->data[2] = 0;
            df_A1->data[3] = 0;
            df_A1->data[4] = 1;
            df_A1->data[5] = 0;
            df_A1->data[6] = df_A1->data[1];
            df_A1->data[7] = df_A1->data[4];
        }
        DataField *df_A2;
        DataField *df_B2;
        df_A2 = new DataField("testQuadAMultToQuadB-df_A2", EMPIRE_DataField_atNode,
                meshQuadAMult->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        df_B2 = new DataField("testQuadAMultToQuadB-df_B2", EMPIRE_DataField_atNode,
                meshQuadB->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        if (constantDataField) { // df_B2 is integrated from constant
            df_B2->data[0] = 1;
            df_B2->data[1] = 2;
            df_B2->data[2] = 2;
            df_B2->data[3] = 1;
            df_B2->data[4] = 1;
            df_B2->data[5] = 2;
            df_B2->data[6] = 2;
            df_B2->data[7] = 1;
        } else { // df_B2 is symmetric along x-direction
            df_B2->data[0] = 0;
            df_B2->data[1] = 1;
            df_B2->data[2] = 4;
            df_B2->data[3] = 9;
            df_B2->data[4] = 0;
            df_B2->data[5] = 1;
            df_B2->data[6] = 4;
            df_B2->data[7] = 9;
        }
        MapperAdapter *mapper = new MapperAdapter("testQuadAMultToQuadB", meshQuadAMult, meshQuadB);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        { // consistent mapping
            AbstractFilter *filterConsistent = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConsistent, meshQuadAMult, df_A1, meshQuadB,
                    df_B1);
            filterConsistent->filtering();
            if (debugMe)
                infoOut << *df_B1;
            delete filterConsistent;
        }
        { // conservative mapping
            AbstractFilter *filterConservative = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConservative, meshQuadB, df_B2, meshQuadAMult,
                    df_A2);
            filterConservative->filtering();
            if (debugMe)
                infoOut << *df_A2;
            delete filterConservative;
        }
        delete df_A1, df_B1, df_A2, df_B2;
        delete mapper;
    }
    /***********************************************************************************************
     * \brief Test case 3: mapping from meshQuadA to meshQuadBMult
     ***********/
    void testQuadAToQuadBMult() {
        DataField *df_A1;
        DataField *df_B1;
        df_A1 = new DataField("testQuadAToQuadBMult-df_A1", EMPIRE_DataField_atNode,
                meshQuadA->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        df_B1 = new DataField("testQuadAToQuadBMult-df_B1", EMPIRE_DataField_atNode,
                meshQuadBMult->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        if (constantDataField) { // df_A1 is constant
            df_A1->data[0] = 1;
            df_A1->data[1] = 1;
            df_A1->data[2] = 1;
            df_A1->data[3] = 1;
            df_A1->data[4] = 1;
            df_A1->data[5] = 1;
        } else { // df_A1 is symmetric along x-direction
            df_A1->data[0] = 0;
            df_A1->data[1] = 1;
            df_A1->data[2] = 0;
            df_A1->data[3] = 0;
            df_A1->data[4] = 1;
            df_A1->data[5] = 0;
        }
        DataField *df_A2;
        DataField *df_B2;
        df_A2 = new DataField("testQuadAToQuadBMult-df_A2", EMPIRE_DataField_atNode,
                meshQuadA->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        df_B2 = new DataField("testQuadAToQuadBMult-df_B2", EMPIRE_DataField_atNode,
                meshQuadBMult->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        if (constantDataField) { // df_B2 is integrated from constant
            df_B2->data[0] = 1;
            df_B2->data[1] = 2.0 / 2.0;
            df_B2->data[2] = 2.0 / 2.0;
            df_B2->data[3] = 1;
            df_B2->data[4] = 1;
            df_B2->data[5] = 2.0 / 2.0;
            df_B2->data[6] = 2.0 / 2.0;
            df_B2->data[7] = 1;
            df_B2->data[8] = df_B2->data[1];
            df_B2->data[9] = df_B2->data[2];
            df_B2->data[10] = df_B2->data[5];
            df_B2->data[11] = df_B2->data[6];
        } else { // df_B2 is symmetric along x-direction
            df_B2->data[0] = 0;
            df_B2->data[1] = 1.0 / 2.0;
            df_B2->data[2] = 4.0 / 2.0;
            df_B2->data[3] = 9;
            df_B2->data[4] = 0;
            df_B2->data[5] = 1.0 / 2.0;
            df_B2->data[6] = 4.0 / 2.0;
            df_B2->data[7] = 9;
            df_B2->data[8] = df_B2->data[1];
            df_B2->data[9] = df_B2->data[2];
            df_B2->data[10] = df_B2->data[5];
            df_B2->data[11] = df_B2->data[6];
        }
        MapperAdapter *mapper = new MapperAdapter("testQuadAToQuadBMult", meshQuadA, meshQuadBMult);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        { // consistent mapping
            AbstractFilter *filterConsistent = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConsistent, meshQuadA, df_A1, meshQuadBMult,
                    df_B1);

            filterConsistent->filtering();
            if (debugMe)
                infoOut << *df_B1;
            delete filterConsistent;
        }
        { // conservative mapping
            AbstractFilter *filterConservative = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConservative, meshQuadBMult, df_B2, meshQuadA,
                    df_A2);
            filterConservative->filtering();
            if (debugMe)
                infoOut << *df_A2;
            delete filterConservative;
        }
        delete df_A1, df_B1, df_A2, df_B2;
        delete mapper;
    }
    /***********************************************************************************************
     * \brief Test case 4: mapping from meshQuadAMult to meshQuadBMult
     ***********/
    void testQuadAMultToQuadBMult() {
        DataField *df_A1;
        DataField *df_B1;
        df_A1 = new DataField("testQuadAMultToQuadBMult-df_A1", EMPIRE_DataField_atNode,
                meshQuadAMult->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        df_B1 = new DataField("testQuadAMultToQuadBMult-df_B1", EMPIRE_DataField_atNode,
                meshQuadBMult->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        if (constantDataField) { // df_A1 is constant
            df_A1->data[0] = 1;
            df_A1->data[1] = 1;
            df_A1->data[2] = 1;
            df_A1->data[3] = 1;
            df_A1->data[4] = 1;
            df_A1->data[5] = 1;
            df_A1->data[6] = df_A1->data[1];
            df_A1->data[7] = df_A1->data[4];
        } else { // df_A1 is symmetric along x-direction
            df_A1->data[0] = 0;
            df_A1->data[1] = 1;
            df_A1->data[2] = 0;
            df_A1->data[3] = 0;
            df_A1->data[4] = 1;
            df_A1->data[5] = 0;
            df_A1->data[6] = df_A1->data[1];
            df_A1->data[7] = df_A1->data[4];
        }
        DataField *df_A2;
        DataField *df_B2;
        df_A2 = new DataField("testQuadAMultToQuadBMult-df_A2", EMPIRE_DataField_atNode,
                meshQuadAMult->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        df_B2 = new DataField("testQuadAMultToQuadBMult-df_B2", EMPIRE_DataField_atNode,
                meshQuadBMult->numNodes, EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        if (constantDataField) { // df_B2 is integrated from constant
            df_B2->data[0] = 1;
            df_B2->data[1] = 2.0 / 2.0;
            df_B2->data[2] = 2.0 / 2.0;
            df_B2->data[3] = 1;
            df_B2->data[4] = 1;
            df_B2->data[5] = 2.0 / 2.0;
            df_B2->data[6] = 2.0 / 2.0;
            df_B2->data[7] = 1;
            df_B2->data[8] = df_B2->data[1];
            df_B2->data[9] = df_B2->data[2];
            df_B2->data[10] = df_B2->data[5];
            df_B2->data[11] = df_B2->data[6];
        } else { // df_B2 is symmetric along x-direction
            df_B2->data[0] = 0;
            df_B2->data[1] = 1.0 / 2.0;
            df_B2->data[2] = 4.0 / 2.0;
            df_B2->data[3] = 9;
            df_B2->data[4] = 0;
            df_B2->data[5] = 1.0 / 2.0;
            df_B2->data[6] = 4.0 / 2.0;
            df_B2->data[7] = 9;
            df_B2->data[8] = df_B2->data[1];
            df_B2->data[9] = df_B2->data[2];
            df_B2->data[10] = df_B2->data[5];
            df_B2->data[11] = df_B2->data[6];
        }
        MapperAdapter *mapper = new MapperAdapter("testQuadAMultToQuadBMult", meshQuadAMult,
                meshQuadBMult);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        { // consistent mapping
            AbstractFilter *filterConsistent = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConsistent, meshQuadAMult, df_A1,
                    meshQuadBMult, df_B1);

            filterConsistent->filtering();
            if (debugMe)
                infoOut << *df_B1;
            delete filterConsistent;
        }
        { // conservative mapping
            AbstractFilter *filterConservative = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConservative, meshQuadBMult, df_B2,
                    meshQuadAMult, df_A2);
            filterConservative->filtering();
            if (debugMe)
                infoOut << *df_A2;
            delete filterConservative;
        }
        delete df_A1, df_B1, df_A2, df_B2;
        delete mapper;
    }
CPPUNIT_TEST_SUITE( TestMappingMultiplePartitions );
        CPPUNIT_TEST( testQuadAToQuadB);
        CPPUNIT_TEST( testQuadAMultToQuadB);
        CPPUNIT_TEST( testQuadAToQuadBMult);
        CPPUNIT_TEST( testQuadAMultToQuadBMult);

    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestMappingMultiplePartitions);
