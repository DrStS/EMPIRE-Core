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
#include "DataField.h"
#include "FEMesh.h"
#include "MapperAdapter.h"
#include "MapperAdapter.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"
#include "Message.h"
#include "GiDFileIO.h"

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test the symmetry of mortar mapper when mapping on a stripe (quasi 2D mesh). The data is symmetric
 *        along the stripe's length direction. The first four test cases show that only when mapping to triangle
 *        gives non-symmetric result. Interpretation: the result is independent of the element type of meshA, because
 *        both triangle and quad elements can represent symmetrical data. But when we integrate a symmetric data field
 *        on triangles (C_BA*u_A, B is triangular mesh), we cannot get a symmetrical result. The last test case shows
 *        that when map from quad to triangle, when the points are overlapped, we get symmetric result.
 ***********/
class TestMappingUnSymmetric: public CppUnit::TestFixture {
private:
    FEMesh *meshTri1;
    FEMesh *meshTri2;
    FEMesh *meshQuad1;
    FEMesh *meshQuad2;
    bool debugMe;

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
        int numNodesPerElemTri = 3;
        {
            /*
             * 4-----5-----6
             * |  \  |  \  | mesh triangle 1 (four triangles)
             * 1-----2-----3
             */
            int numNodesTri1 = 6;
            int numElemsTri1 = 4;
            meshTri1 = new FEMesh("", numNodesTri1, numElemsTri1);
            for (int i = 0; i < numElemsTri1; i++)
                meshTri1->numNodesPerElem[i] = numNodesPerElemTri;
            meshTri1->initElems();

            for (int i = 0; i < numNodesTri1; i++)
                meshTri1->nodeIDs[i] = i + 1;

            meshTri1->nodes[0 * 3 + 0] = 0;
            meshTri1->nodes[0 * 3 + 1] = 0;
            meshTri1->nodes[0 * 3 + 2] = 0;

            meshTri1->nodes[1 * 3 + 0] = 1.5;
            meshTri1->nodes[1 * 3 + 1] = 0;
            meshTri1->nodes[1 * 3 + 2] = 0;

            meshTri1->nodes[2 * 3 + 0] = 3;
            meshTri1->nodes[2 * 3 + 1] = 0;
            meshTri1->nodes[2 * 3 + 2] = 0;

            meshTri1->nodes[3 * 3 + 0] = 0;
            meshTri1->nodes[3 * 3 + 1] = 1;
            meshTri1->nodes[3 * 3 + 2] = 0;

            meshTri1->nodes[4 * 3 + 0] = 1.5;
            meshTri1->nodes[4 * 3 + 1] = 1;
            meshTri1->nodes[4 * 3 + 2] = 0;

            meshTri1->nodes[5 * 3 + 0] = 3;
            meshTri1->nodes[5 * 3 + 1] = 1;
            meshTri1->nodes[5 * 3 + 2] = 0;

            meshTri1->elems[0 * 3 + 0] = 1;
            meshTri1->elems[0 * 3 + 1] = 2;
            meshTri1->elems[0 * 3 + 2] = 4;

            meshTri1->elems[1 * 3 + 0] = 2;
            meshTri1->elems[1 * 3 + 1] = 5;
            meshTri1->elems[1 * 3 + 2] = 4;

            meshTri1->elems[2 * 3 + 0] = 2;
            meshTri1->elems[2 * 3 + 1] = 3;
            meshTri1->elems[2 * 3 + 2] = 5;

            meshTri1->elems[3 * 3 + 0] = 3;
            meshTri1->elems[3 * 3 + 1] = 6;
            meshTri1->elems[3 * 3 + 2] = 5;

            /*GiDFileIO::writeDotMsh("meshTri1.msh", meshTri1->numNodes, meshTri1->numElems,
             meshTri1->nodes, meshTri1->nodeIDs, meshTri1->numNodesPerElem, meshTri1->elems,
             meshTri1->elemIDs);*/
        }
        {
            /*
             * 5---6---7---8
             * | \ | \ | \ | mesh triangle 2 (six triangles)
             * 1---2---3---4
             */
            int numNodesTri2 = 8;
            int numElemsTri2 = 6;
            meshTri2 = new FEMesh("", numNodesTri2, numElemsTri2);
            for (int i = 0; i < numElemsTri2; i++)
                meshTri2->numNodesPerElem[i] = numNodesPerElemTri;
            meshTri2->initElems();

            for (int i = 0; i < numNodesTri2; i++)
                meshTri2->nodeIDs[i] = i + 1;

            meshTri2->nodes[0 * 3 + 0] = 0;
            meshTri2->nodes[0 * 3 + 1] = 0;
            meshTri2->nodes[0 * 3 + 2] = 0;

            meshTri2->nodes[1 * 3 + 0] = 1;
            meshTri2->nodes[1 * 3 + 1] = 0;
            meshTri2->nodes[1 * 3 + 2] = 0;

            meshTri2->nodes[2 * 3 + 0] = 2;
            meshTri2->nodes[2 * 3 + 1] = 0;
            meshTri2->nodes[2 * 3 + 2] = 0;

            meshTri2->nodes[3 * 3 + 0] = 3;
            meshTri2->nodes[3 * 3 + 1] = 0;
            meshTri2->nodes[3 * 3 + 2] = 0;

            meshTri2->nodes[4 * 3 + 0] = 0;
            meshTri2->nodes[4 * 3 + 1] = 1;
            meshTri2->nodes[4 * 3 + 2] = 0;

            meshTri2->nodes[5 * 3 + 0] = 1;
            meshTri2->nodes[5 * 3 + 1] = 1;
            meshTri2->nodes[5 * 3 + 2] = 0;

            meshTri2->nodes[6 * 3 + 0] = 2;
            meshTri2->nodes[6 * 3 + 1] = 1;
            meshTri2->nodes[6 * 3 + 2] = 0;

            meshTri2->nodes[7 * 3 + 0] = 3;
            meshTri2->nodes[7 * 3 + 1] = 1;
            meshTri2->nodes[7 * 3 + 2] = 0;

            meshTri2->elems[0 * 3 + 0] = 1;
            meshTri2->elems[0 * 3 + 1] = 2;
            meshTri2->elems[0 * 3 + 2] = 5;

            meshTri2->elems[1 * 3 + 0] = 2;
            meshTri2->elems[1 * 3 + 1] = 6;
            meshTri2->elems[1 * 3 + 2] = 5;

            meshTri2->elems[2 * 3 + 0] = 2;
            meshTri2->elems[2 * 3 + 1] = 3;
            meshTri2->elems[2 * 3 + 2] = 6;

            meshTri2->elems[3 * 3 + 0] = 3;
            meshTri2->elems[3 * 3 + 1] = 7;
            meshTri2->elems[3 * 3 + 2] = 6;

            meshTri2->elems[4 * 3 + 0] = 3;
            meshTri2->elems[4 * 3 + 1] = 4;
            meshTri2->elems[4 * 3 + 2] = 7;

            meshTri2->elems[5 * 3 + 0] = 4;
            meshTri2->elems[5 * 3 + 1] = 8;
            meshTri2->elems[5 * 3 + 2] = 7;

            /*GiDFileIO::writeDotMsh("meshTri2.msh", meshTri2->numNodes, meshTri2->numElems,
             meshTri2->nodes, meshTri2->nodeIDs, meshTri2->numNodesPerElem, meshTri2->elems,
             meshTri2->elemIDs);*/
        }
        int numNodesPerElemQuad = 4;
        {
            /*
             * 4-----5-----6
             * |     |     | mesh quad 1
             * 1-----2-----3
             */
            int numNodesQuad1 = 6;
            int numElemsQuad1 = 2;
            meshQuad1 = new FEMesh("", numNodesQuad1, numElemsQuad1);
            for (int i = 0; i < numElemsQuad1; i++)
                meshQuad1->numNodesPerElem[i] = numNodesPerElemQuad;
            meshQuad1->initElems();

            for (int i = 0; i < numNodesQuad1; i++)
                meshQuad1->nodeIDs[i] = i + 1;

            meshQuad1->nodes[0 * 3 + 0] = 0;
            meshQuad1->nodes[0 * 3 + 1] = 0;
            meshQuad1->nodes[0 * 3 + 2] = 0;

            meshQuad1->nodes[1 * 3 + 0] = 1.5;
            meshQuad1->nodes[1 * 3 + 1] = 0;
            meshQuad1->nodes[1 * 3 + 2] = 0;

            meshQuad1->nodes[2 * 3 + 0] = 3;
            meshQuad1->nodes[2 * 3 + 1] = 0;
            meshQuad1->nodes[2 * 3 + 2] = 0;

            meshQuad1->nodes[3 * 3 + 0] = 0;
            meshQuad1->nodes[3 * 3 + 1] = 1;
            meshQuad1->nodes[3 * 3 + 2] = 0;

            meshQuad1->nodes[4 * 3 + 0] = 1.5;
            meshQuad1->nodes[4 * 3 + 1] = 1;
            meshQuad1->nodes[4 * 3 + 2] = 0;

            meshQuad1->nodes[5 * 3 + 0] = 3;
            meshQuad1->nodes[5 * 3 + 1] = 1;
            meshQuad1->nodes[5 * 3 + 2] = 0;

            meshQuad1->elems[0 * 4 + 0] = 1;
            meshQuad1->elems[0 * 4 + 1] = 2;
            meshQuad1->elems[0 * 4 + 2] = 5;
            meshQuad1->elems[0 * 4 + 3] = 4;

            meshQuad1->elems[1 * 4 + 0] = 2;
            meshQuad1->elems[1 * 4 + 1] = 3;
            meshQuad1->elems[1 * 4 + 2] = 6;
            meshQuad1->elems[1 * 4 + 3] = 5;

            /*GiDFileIO::writeDotMsh("meshQuad1.msh", meshQuad1->numNodes, meshQuad1->numElems,
             meshQuad1->nodes, meshQuad1->nodeIDs, meshQuad1->numNodesPerElem,
             meshQuad1->elems, meshQuad1->elemIDs);*/
        }
        {
            /*
             * 5---6---7---8
             * |   |   |   | mesh quad 2
             * 1---2---3---4
             */
            int numNodesQuad2 = 8;
            int numElemsQuad2 = 3;
            meshQuad2 = new FEMesh("", numNodesQuad2, numElemsQuad2);
            for (int i = 0; i < numElemsQuad2; i++)
                meshQuad2->numNodesPerElem[i] = numNodesPerElemQuad;
            meshQuad2->initElems();

            for (int i = 0; i < numNodesQuad2; i++)
                meshQuad2->nodeIDs[i] = i + 1;

            meshQuad2->nodes[0 * 3 + 0] = 0;
            meshQuad2->nodes[0 * 3 + 1] = 0;
            meshQuad2->nodes[0 * 3 + 2] = 0;

            meshQuad2->nodes[1 * 3 + 0] = 1;
            meshQuad2->nodes[1 * 3 + 1] = 0;
            meshQuad2->nodes[1 * 3 + 2] = 0;

            meshQuad2->nodes[2 * 3 + 0] = 2;
            meshQuad2->nodes[2 * 3 + 1] = 0;
            meshQuad2->nodes[2 * 3 + 2] = 0;

            meshQuad2->nodes[3 * 3 + 0] = 3;
            meshQuad2->nodes[3 * 3 + 1] = 0;
            meshQuad2->nodes[3 * 3 + 2] = 0;

            meshQuad2->nodes[4 * 3 + 0] = 0;
            meshQuad2->nodes[4 * 3 + 1] = 1;
            meshQuad2->nodes[4 * 3 + 2] = 0;

            meshQuad2->nodes[5 * 3 + 0] = 1;
            meshQuad2->nodes[5 * 3 + 1] = 1;
            meshQuad2->nodes[5 * 3 + 2] = 0;

            meshQuad2->nodes[6 * 3 + 0] = 2;
            meshQuad2->nodes[6 * 3 + 1] = 1;
            meshQuad2->nodes[6 * 3 + 2] = 0;

            meshQuad2->nodes[7 * 3 + 0] = 3;
            meshQuad2->nodes[7 * 3 + 1] = 1;
            meshQuad2->nodes[7 * 3 + 2] = 0;

            meshQuad2->elems[0 * 4 + 0] = 1;
            meshQuad2->elems[0 * 4 + 1] = 2;
            meshQuad2->elems[0 * 4 + 2] = 6;
            meshQuad2->elems[0 * 4 + 3] = 5;

            meshQuad2->elems[1 * 4 + 0] = 2;
            meshQuad2->elems[1 * 4 + 1] = 3;
            meshQuad2->elems[1 * 4 + 2] = 7;
            meshQuad2->elems[1 * 4 + 3] = 6;

            meshQuad2->elems[2 * 4 + 0] = 3;
            meshQuad2->elems[2 * 4 + 1] = 4;
            meshQuad2->elems[2 * 4 + 2] = 8;
            meshQuad2->elems[2 * 4 + 3] = 7;

            /*GiDFileIO::writeDotMsh("meshQuad2.msh", meshQuad2->numNodes, meshQuad2->numElems,
             meshQuad2->nodes, meshQuad2->nodeIDs, meshQuad2->numNodesPerElem,
             meshQuad2->elems, meshQuad2->elemIDs);*/
        }
    }
    void tearDown() {
        delete meshTri1;
        delete meshTri2;
        delete meshQuad1;
        delete meshQuad2;
    }
    /***********************************************************************************************
     * \brief Test case: mapping from meshQuad1 to meshQuad2
     ***********/
    void testQuad1ToQuad2() {
        DataField *df_A;
        DataField *df_B;
        df_A = new DataField("df_A", EMPIRE_DataField_atNode, meshQuad1->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        df_B = new DataField("df_B", EMPIRE_DataField_atNode, meshQuad2->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        { // df_A is symmetric along x-direction
            df_A->data[0] = 0;
            df_A->data[1] = 1;
            df_A->data[2] = 4;
            df_A->data[3] = 0;
            df_A->data[4] = 1;
            df_A->data[5] = 4;
        }
        MapperAdapter *mapper = new MapperAdapter("", meshQuad1, meshQuad2);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, meshQuad1, df_A, meshQuad2, df_B);
        filterConsistent->filtering();
        if (debugMe)
            infoOut << *df_B;
        {
            // fractional numbers are computed by:
            // http://www.mindspring.com/~alanh/fracs.html
            double EPS = 1e-10;
            CPPUNIT_ASSERT(fabs(df_B->data[0] - 0.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[1] - 7.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[2] - 23.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[3] - 4.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[4] - 0.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[5] - 7.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[6] - 23.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[7] - 4.0) < EPS);
        }
        delete mapper;
        delete filterConsistent;
        delete df_A, df_B;
    }
    /***********************************************************************************************
     * \brief Test case: mapping from meshTri1 to meshQuad2
     ***********/
    void testTri1ToQuad2() {
        DataField *df_A;
        DataField *df_B;
        df_A = new DataField("df_A", EMPIRE_DataField_atNode, meshTri1->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        df_B = new DataField("df_B", EMPIRE_DataField_atNode, meshQuad2->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        { // df_A is symmetric along x-direction
            df_A->data[0] = 0;
            df_A->data[1] = 1;
            df_A->data[2] = 4;
            df_A->data[3] = 0;
            df_A->data[4] = 1;
            df_A->data[5] = 4;
        }
        MapperAdapter *mapper = new MapperAdapter("", meshTri1, meshQuad2);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, meshTri1, df_A, meshQuad2, df_B);

        filterConsistent->filtering();
        if (debugMe)
            infoOut << *df_B;
        {
            // fractional numbers are computed by:
            // http://www.mindspring.com/~alanh/fracs.html
            double EPS = 1e-10;
            CPPUNIT_ASSERT(fabs(df_B->data[0] - 0.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[1] - 7.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[2] - 23.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[3] - 4.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[4] - 0.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[5] - 7.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[6] - 23.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[7] - 4.0) < EPS);
        }
        delete mapper;
        delete filterConsistent;
        delete df_A, df_B;
    }
    /***********************************************************************************************
     * \brief Test case: mapping from meshQuad1 to meshTri2
     ***********/
    void testQuad1ToTri2() {
        DataField *df_A;
        DataField *df_B;
        df_A = new DataField("df_A", EMPIRE_DataField_atNode, meshQuad1->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        df_B = new DataField("df_B", EMPIRE_DataField_atNode, meshTri2->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        { // df_A is symmetric along x-direction
            df_A->data[0] = 0;
            df_A->data[1] = 1;
            df_A->data[2] = 4;
            df_A->data[3] = 0;
            df_A->data[4] = 1;
            df_A->data[5] = 4;
        }
        MapperAdapter *mapper = new MapperAdapter("", meshQuad1, meshTri2);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, meshQuad1, df_A, meshTri2, df_B);
        filterConsistent->filtering();
        if (debugMe)
            infoOut << *df_B;
        {
            // fractional numbers are computed by:
            // http://www.mindspring.com/~alanh/fracs.html
            double EPS = 1e-10;
            CPPUNIT_ASSERT(fabs(df_B->data[0] - 0.0) < EPS);
            //CPPUNIT_ASSERT(fabs(df_B->data[1] - 7.0/12.0) < EPS);
            //CPPUNIT_ASSERT(fabs(df_B->data[2] - 23.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[3] - 4.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[4] - 0.0) < EPS);
            //CPPUNIT_ASSERT(fabs(df_B->data[5] - 7.0/12.0) < EPS);
            //CPPUNIT_ASSERT(fabs(df_B->data[6] - 23.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[7] - 4.0) < EPS);
        }

        delete mapper;
        delete filterConsistent;
        delete df_A, df_B;
    }
    /***********************************************************************************************
     * \brief Test case: mapping from meshTri1 to meshTri2
     ***********/
    void testTri1ToTri2() {
        DataField *df_A;
        DataField *df_B;
        df_A = new DataField("df_A", EMPIRE_DataField_atNode, meshTri1->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        df_B = new DataField("df_B", EMPIRE_DataField_atNode, meshTri2->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        { // df_A is symmetric along x-direction
            df_A->data[0] = 0;
            df_A->data[1] = 1;
            df_A->data[2] = 4;
            df_A->data[3] = 0;
            df_A->data[4] = 1;
            df_A->data[5] = 4;
        }
        MapperAdapter *mapper = new MapperAdapter("", meshTri1, meshTri2);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, meshTri1, df_A, meshTri2, df_B);
        filterConsistent->filtering();
        if (debugMe)
            infoOut << *df_B;
        {
            // fractional numbers are computed by:
            // http://www.mindspring.com/~alanh/fracs.html
            double EPS = 1e-10;
            CPPUNIT_ASSERT(fabs(df_B->data[0] - 0.0) < EPS);
            //CPPUNIT_ASSERT(fabs(df_B->data[1] - 7.0/12.0) < EPS);
            //CPPUNIT_ASSERT(fabs(df_B->data[2] - 23.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[3] - 4.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[4] - 0.0) < EPS);
            //CPPUNIT_ASSERT(fabs(df_B->data[5] - 7.0/12.0) < EPS);
            //CPPUNIT_ASSERT(fabs(df_B->data[6] - 23.0/12.0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[7] - 4.0) < EPS);
        }

        delete mapper;
        delete filterConsistent;
        delete df_A, df_B;
    }
    /***********************************************************************************************
     * \brief Test case: mapping from meshQuad2 to meshTri2
     ***********/
    void testQuad2ToTri2() {
        DataField *df_A;
        DataField *df_B;
        df_A = new DataField("df_A", EMPIRE_DataField_atNode, meshQuad2->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        df_B = new DataField("df_B", EMPIRE_DataField_atNode, meshTri2->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        { // df_A is symmetric along x-direction
            df_A->data[0] = 0;
            df_A->data[1] = 1;
            df_A->data[2] = 4;
            df_A->data[3] = 9;
            df_A->data[4] = 0;
            df_A->data[5] = 1;
            df_A->data[6] = 4;
            df_A->data[7] = 9;
        }
        MapperAdapter *mapper = new MapperAdapter("", meshQuad2, meshTri2);
        mapper->initMortarMapper(oppositeSurfaceNormal, dual, enforceConsistency);
        AbstractFilter *filterConsistent = new MappingFilter(mapper);
        ConnectionIOSetup::setupIOForFilter(filterConsistent, meshQuad2, df_A, meshTri2, df_B);
        filterConsistent->filtering();
        if (debugMe)
            infoOut << *df_B;
        {
            // fractional numbers are computed by:
            // http://www.mindspring.com/~alanh/fracs.html
            double EPS = 1e-10;
            CPPUNIT_ASSERT(fabs(df_B->data[0] - 0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[1] - 1) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[2] - 4) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[3] - 9) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[4] - 0) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[5] - 1) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[6] - 4) < EPS);
            CPPUNIT_ASSERT(fabs(df_B->data[7] - 9) < EPS);
        }

        delete mapper;
        delete filterConsistent;
        delete df_A, df_B;
    }

CPPUNIT_TEST_SUITE( TestMappingUnSymmetric );
        CPPUNIT_TEST( testQuad1ToQuad2);
        CPPUNIT_TEST( testTri1ToQuad2);
        CPPUNIT_TEST( testQuad1ToTri2);
        CPPUNIT_TEST( testTri1ToTri2);
        CPPUNIT_TEST( testQuad2ToTri2);

    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestMappingUnSymmetric);
