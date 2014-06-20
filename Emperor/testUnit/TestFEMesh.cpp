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
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"
#include "FEMesh.h"
#include "DataField.h"
#include "Message.h"
#include <iostream>

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test the class Mesh
 ***********/
class TestFEMesh: public CppUnit::TestFixture {
public:
    void setUp() {
    }
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief Test case: test creation of mesh
     ***********/
    void testMeshCreation() {
        { // 1. mesh with one triangle
            int numNodes = 3;
            int numElems = 1;
            FEMesh *mesh = new FEMesh("triangle", numNodes, numElems);
            CPPUNIT_ASSERT(mesh->name == "triangle");
            CPPUNIT_ASSERT(mesh->numNodes == 3);
            CPPUNIT_ASSERT(mesh->nodes != NULL);
            CPPUNIT_ASSERT(mesh->nodeIDs != NULL);
            CPPUNIT_ASSERT(mesh->numElems == 1);
            CPPUNIT_ASSERT(mesh->numNodesPerElem != NULL);
            CPPUNIT_ASSERT(mesh->elems == NULL);
            mesh->numNodesPerElem[0] = 3;
            mesh->initElems();
            CPPUNIT_ASSERT(mesh->tobeTriangulated == false);
            CPPUNIT_ASSERT(mesh->elems != NULL);
            CPPUNIT_ASSERT(mesh->elemsArraySize == 3);
            delete mesh;
        }
        { // 2. mesh with one quad
            int numNodes = 4;
            int numElems = 1;
            FEMesh *mesh = new FEMesh("quad", numNodes, numElems);
            mesh->numNodesPerElem[0] = 4;
            mesh->initElems();
            CPPUNIT_ASSERT(mesh->elemsArraySize == 4);
            CPPUNIT_ASSERT(mesh->tobeTriangulated == false);
            delete mesh;
        }
        { // 3. mesh with one triangle and one quad and one hexagon
            /*
             *      /\
             *     /__\
             *     |__|
             *     /  \
             *     \__/
             */
            int numNodes = 9;
            int numElems = 3;
            FEMesh *mesh = new FEMesh("dummy", numNodes, numElems);
            mesh->numNodesPerElem[0] = 3;
            mesh->numNodesPerElem[1] = 4;
            mesh->numNodesPerElem[2] = 6;
            mesh->initElems();
            CPPUNIT_ASSERT(mesh->elemsArraySize == 13);
            CPPUNIT_ASSERT(mesh->tobeTriangulated == true);
            delete mesh;
        }
    }
    /***********************************************************************************************
     * \brief Test case: Test the data fields
     ***********/
    void testDataField() {
        int numNodes = 3;
        int numElems = 1;
        FEMesh *mesh = new FEMesh("dummy", numNodes, numElems);
        mesh->addDataField("pressure", EMPIRE_DataField_atElemCentroid, EMPIRE_DataField_vector,
                EMPIRE_DataField_field);
        mesh->addDataField("temperature", EMPIRE_DataField_atNode, EMPIRE_DataField_scalar,
                EMPIRE_DataField_field);
        CPPUNIT_ASSERT(mesh->getDataFieldByName("pressure")->name == "pressure");
        CPPUNIT_ASSERT(
                mesh->getDataFieldByName("pressure")->location == EMPIRE_DataField_atElemCentroid);
        CPPUNIT_ASSERT(mesh->getDataFieldByName("pressure")->dimension == EMPIRE_DataField_vector);
        CPPUNIT_ASSERT(mesh->getDataFieldByName("temperature")->name == "temperature");
        CPPUNIT_ASSERT(
                mesh->getDataFieldByName("temperature")->location == EMPIRE_DataField_atNode);
        CPPUNIT_ASSERT(
                mesh->getDataFieldByName("temperature")->dimension == EMPIRE_DataField_scalar);
        delete mesh;
    }
    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/
    void testRevertSurfaceNormal() {
        int numNodes = 3;
        int numElems = 2;
        FEMesh *mesh = new FEMesh("dummy", numNodes, numElems);
        mesh->numNodesPerElem[0] = 3;
        mesh->numNodesPerElem[1] = 3;
        mesh->initElems();
        mesh->elems[0] = 1;
        mesh->elems[1] = 2;
        mesh->elems[2] = 3;
        mesh->elems[3] = 3;
        mesh->elems[4] = 2;
        mesh->elems[5] = 1;
        revertSurfaceNormalOfFEMesh(mesh);
        CPPUNIT_ASSERT(mesh->elems[0] == 3);
        CPPUNIT_ASSERT(mesh->elems[1] == 2);
        CPPUNIT_ASSERT(mesh->elems[2] == 1);
        CPPUNIT_ASSERT(mesh->elems[3] == 1);
        CPPUNIT_ASSERT(mesh->elems[4] == 2);
        CPPUNIT_ASSERT(mesh->elems[5] == 3);
        delete mesh;
    }
    /***********************************************************************************************
     * \brief Test case: Test mesh triangulation
     ***********/
    void testTriangulation() {
        /*
         *       0
         *      /\
         *    1/__\2
         *    3|__|4
         *    5/  \6
         *    7\__/8
         */
        { // do NOT triangulate all
            int numNodes = 9;
            int numElems = 3;
            FEMesh *mesh = new FEMesh("dummy", numNodes, numElems);
            for (int i = 0; i < numNodes; i++)
                mesh->nodeIDs[i] = i;
            { // set up node coordinates
                mesh->nodes[0 * 3 + 0] = 0;
                mesh->nodes[0 * 3 + 1] = 4;
                mesh->nodes[0 * 3 + 2] = 0;
                mesh->nodes[1 * 3 + 0] = -1;
                mesh->nodes[1 * 3 + 1] = 3;
                mesh->nodes[1 * 3 + 2] = 0;
                mesh->nodes[2 * 3 + 0] = 1;
                mesh->nodes[2 * 3 + 1] = 3;
                mesh->nodes[2 * 3 + 2] = 0;
                mesh->nodes[3 * 3 + 0] = -1;
                mesh->nodes[3 * 3 + 1] = 2;
                mesh->nodes[3 * 3 + 2] = 0;
                mesh->nodes[4 * 3 + 0] = 1;
                mesh->nodes[4 * 3 + 1] = 2;
                mesh->nodes[4 * 3 + 2] = 0;
                mesh->nodes[5 * 3 + 0] = -2;
                mesh->nodes[5 * 3 + 1] = 1;
                mesh->nodes[5 * 3 + 2] = 0;
                mesh->nodes[6 * 3 + 0] = 2;
                mesh->nodes[6 * 3 + 1] = 1;
                mesh->nodes[6 * 3 + 2] = 0;
                mesh->nodes[7 * 3 + 0] = -1;
                mesh->nodes[7 * 3 + 1] = 0;
                mesh->nodes[7 * 3 + 2] = 0;
                mesh->nodes[8 * 3 + 0] = 1;
                mesh->nodes[8 * 3 + 1] = 0;
                mesh->nodes[8 * 3 + 2] = 0;
            }

            mesh->numNodesPerElem[0] = 3;
            mesh->numNodesPerElem[1] = 4;
            mesh->numNodesPerElem[2] = 6;

            mesh->initElems();
            CPPUNIT_ASSERT(mesh->tobeTriangulated == true);

            { // triangle
                mesh->elems[0] = 0;
                mesh->elems[1] = 1;
                mesh->elems[2] = 2;
            }
            { // quad
                mesh->elems[3] = 1;
                mesh->elems[4] = 3;
                mesh->elems[5] = 4;
                mesh->elems[6] = 2;
            }
            { // hexagon
                mesh->elems[7] = 3;
                mesh->elems[8] = 5;
                mesh->elems[9] = 7;
                mesh->elems[10] = 8;
                mesh->elems[11] = 6;
                mesh->elems[12] = 4;
            }

            FEMesh *triangulated = mesh->triangulate();
            { // elemsTri should be defined
                CPPUNIT_ASSERT(triangulated != NULL);
                CPPUNIT_ASSERT(triangulated->numNodes == 9);
                CPPUNIT_ASSERT(triangulated->numElems == 6);
                CPPUNIT_ASSERT(triangulated->numNodesPerElem[0] == 3);
                CPPUNIT_ASSERT(triangulated->numNodesPerElem[1] == 4);
                CPPUNIT_ASSERT(triangulated->numNodesPerElem[2] == 3);
                CPPUNIT_ASSERT(triangulated->numNodesPerElem[3] == 3);
                CPPUNIT_ASSERT(triangulated->numNodesPerElem[4] == 3);
                CPPUNIT_ASSERT(triangulated->numNodesPerElem[5] == 3);

                // output to shell to check whether it is correct or not
                // infoOut << triangulated;
            }

            delete mesh;
        }
        { // triangulate all
            int numNodes = 9;
            int numElems = 3;
            FEMesh *mesh = new FEMesh("dummy", numNodes, numElems);
            mesh->triangulateAll = true;
            for (int i = 0; i < numNodes; i++)
                mesh->nodeIDs[i] = i;
            { // set up node coordinates
                mesh->nodes[0 * 3 + 0] = 0;
                mesh->nodes[0 * 3 + 1] = 4;
                mesh->nodes[0 * 3 + 2] = 0;
                mesh->nodes[1 * 3 + 0] = -1;
                mesh->nodes[1 * 3 + 1] = 3;
                mesh->nodes[1 * 3 + 2] = 0;
                mesh->nodes[2 * 3 + 0] = 1;
                mesh->nodes[2 * 3 + 1] = 3;
                mesh->nodes[2 * 3 + 2] = 0;
                mesh->nodes[3 * 3 + 0] = -1;
                mesh->nodes[3 * 3 + 1] = 2;
                mesh->nodes[3 * 3 + 2] = 0;
                mesh->nodes[4 * 3 + 0] = 1;
                mesh->nodes[4 * 3 + 1] = 2;
                mesh->nodes[4 * 3 + 2] = 0;
                mesh->nodes[5 * 3 + 0] = -2;
                mesh->nodes[5 * 3 + 1] = 1;
                mesh->nodes[5 * 3 + 2] = 0;
                mesh->nodes[6 * 3 + 0] = 2;
                mesh->nodes[6 * 3 + 1] = 1;
                mesh->nodes[6 * 3 + 2] = 0;
                mesh->nodes[7 * 3 + 0] = -1;
                mesh->nodes[7 * 3 + 1] = 0;
                mesh->nodes[7 * 3 + 2] = 0;
                mesh->nodes[8 * 3 + 0] = 1;
                mesh->nodes[8 * 3 + 1] = 0;
                mesh->nodes[8 * 3 + 2] = 0;
            }

            mesh->numNodesPerElem[0] = 3;
            mesh->numNodesPerElem[1] = 4;
            mesh->numNodesPerElem[2] = 6;

            mesh->initElems();
            CPPUNIT_ASSERT(mesh->tobeTriangulated == true);

            { // triangle
                mesh->elems[0] = 0;
                mesh->elems[1] = 1;
                mesh->elems[2] = 2;
            }
            { // quad
                mesh->elems[3] = 1;
                mesh->elems[4] = 3;
                mesh->elems[5] = 4;
                mesh->elems[6] = 2;
            }
            { // hexagon
                mesh->elems[7] = 3;
                mesh->elems[8] = 5;
                mesh->elems[9] = 7;
                mesh->elems[10] = 8;
                mesh->elems[11] = 6;
                mesh->elems[12] = 4;
            }

            FEMesh *triangulated = mesh->triangulate();
            { // elemsTri should be defined
                CPPUNIT_ASSERT(triangulated != NULL);
                CPPUNIT_ASSERT(triangulated->numNodes == 9);
                CPPUNIT_ASSERT(triangulated->numElems == 7);
                for (int i=0; i<7; i++) {
                    CPPUNIT_ASSERT(triangulated->numNodesPerElem[i] == 3);
                }
            }

            delete mesh;
        }

    }

    /***********************************************************************************************
     * \brief Test case: Test triangle orientation after mesh triangulation
     ***********/
    void testTriangulation2() {
        /*
         *       0
         *      /\
         *    1/  \2
         *    3|__|4
         *    5/  \6
         *    7\__/8
         */
        int numNodes = 9;
        int numElems = 2;
        FEMesh *mesh = new FEMesh("dummy", numNodes, numElems);
        for (int i = 0; i < numNodes; i++)
            mesh->nodeIDs[i] = i;
        { // set up node coordinates
            mesh->nodes[0 * 3 + 0] = 0;
            mesh->nodes[0 * 3 + 1] = 4;
            mesh->nodes[0 * 3 + 2] = 0;
            mesh->nodes[1 * 3 + 0] = -1;
            mesh->nodes[1 * 3 + 1] = 3;
            mesh->nodes[1 * 3 + 2] = 0;
            mesh->nodes[2 * 3 + 0] = 1;
            mesh->nodes[2 * 3 + 1] = 3;
            mesh->nodes[2 * 3 + 2] = 0;
            mesh->nodes[3 * 3 + 0] = -1;
            mesh->nodes[3 * 3 + 1] = 2;
            mesh->nodes[3 * 3 + 2] = 0;
            mesh->nodes[4 * 3 + 0] = 1;
            mesh->nodes[4 * 3 + 1] = 2;
            mesh->nodes[4 * 3 + 2] = 0;
            mesh->nodes[5 * 3 + 0] = -2;
            mesh->nodes[5 * 3 + 1] = 1;
            mesh->nodes[5 * 3 + 2] = 0;
            mesh->nodes[6 * 3 + 0] = 2;
            mesh->nodes[6 * 3 + 1] = 1;
            mesh->nodes[6 * 3 + 2] = 0;
            mesh->nodes[7 * 3 + 0] = -1;
            mesh->nodes[7 * 3 + 1] = 0;
            mesh->nodes[7 * 3 + 2] = 0;
            mesh->nodes[8 * 3 + 0] = 1;
            mesh->nodes[8 * 3 + 1] = 0;
            mesh->nodes[8 * 3 + 2] = 0;
        }

        mesh->numNodesPerElem[0] = 5;
        mesh->numNodesPerElem[1] = 6;

        mesh->initElems();
        CPPUNIT_ASSERT(mesh->tobeTriangulated == true);

        { // pentagon
            mesh->elems[0] = 0;
            mesh->elems[1] = 1;
            mesh->elems[2] = 3;
            mesh->elems[3] = 4;
            mesh->elems[4] = 2;
        }
        { // hexagon
            mesh->elems[5] = 3;
            mesh->elems[6] = 4;
            mesh->elems[7] = 6;
            mesh->elems[8] = 8;
            mesh->elems[9] = 7;
            mesh->elems[10] = 5;
        }

        FEMesh *triangulated = mesh->triangulate();
        {
            CPPUNIT_ASSERT(triangulated != NULL);
            // output to shell to check whether it is correct or not
            // infoOut << triangulated;
        }

        delete mesh;
    }

    /***********************************************************************************************
     * \brief Test case: test bounding box computation
     ***********/
    void testBoundingBox() {
        int numNodes = 8;
        int numElems = 6;
        FEMesh *mesh = new FEMesh("dummy", numNodes, numElems);
        for (int i = 0; i < numElems; i++)
            mesh->numNodesPerElem[i] = 4;
        mesh->initElems();

        mesh->nodes[0 * 3 + 0] = 0;
        mesh->nodes[0 * 3 + 1] = 0;
        mesh->nodes[0 * 3 + 2] = 0;

        mesh->nodes[1 * 3 + 0] = 1;
        mesh->nodes[1 * 3 + 1] = 0;
        mesh->nodes[1 * 3 + 2] = 0;

        mesh->nodes[2 * 3 + 0] = 1;
        mesh->nodes[2 * 3 + 1] = 1;
        mesh->nodes[2 * 3 + 2] = 0;

        mesh->nodes[3 * 3 + 0] = 0;
        mesh->nodes[3 * 3 + 1] = 1;
        mesh->nodes[3 * 3 + 2] = 0;

        mesh->nodes[4 * 3 + 0] = 0;
        mesh->nodes[4 * 3 + 1] = 0;
        mesh->nodes[4 * 3 + 2] = 1;

        mesh->nodes[5 * 3 + 0] = 1;
        mesh->nodes[5 * 3 + 1] = 0;
        mesh->nodes[5 * 3 + 2] = 1;

        mesh->nodes[6 * 3 + 0] = 1;
        mesh->nodes[6 * 3 + 1] = 1;
        mesh->nodes[6 * 3 + 2] = 1;

        mesh->nodes[7 * 3 + 0] = 0;
        mesh->nodes[7 * 3 + 1] = 1;
        mesh->nodes[7 * 3 + 2] = 1;

        mesh->computeBoundingBox();
        CPPUNIT_ASSERT(mesh->boundingBox.xmin == 0.0);
        CPPUNIT_ASSERT(mesh->boundingBox.xmax == 1.0);
        CPPUNIT_ASSERT(mesh->boundingBox.ymin == 0.0);
        CPPUNIT_ASSERT(mesh->boundingBox.ymax == 1.0);
        CPPUNIT_ASSERT(mesh->boundingBox.zmin == 0.0);
        CPPUNIT_ASSERT(mesh->boundingBox.zmax == 1.0);

        //infoOut << *mesh;
        delete mesh;
    }

    CPPUNIT_TEST_SUITE(TestFEMesh);
    CPPUNIT_TEST(testMeshCreation);
    CPPUNIT_TEST(testDataField);
    CPPUNIT_TEST(testRevertSurfaceNormal);
    CPPUNIT_TEST(testTriangulation);
    CPPUNIT_TEST(testTriangulation2);
    CPPUNIT_TEST(testBoundingBox);CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestFEMesh);
