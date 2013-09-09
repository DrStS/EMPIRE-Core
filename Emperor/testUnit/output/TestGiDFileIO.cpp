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

#include "GiDFileIO.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <string>

using namespace std;

extern string pathToFolderOfFiles;

namespace EMPIRE {
/********//**
 * \brief Test the class GiDFileIO
 ***********/
class TestGiDFileIO: public CppUnit::TestFixture {
public:
    void setUp() {
    }
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief Test case: Test read/write in GiD .msh format.
     * \author Michael Andre, Tianyang Wang
     ***********/
    void readWriteGiDMsh() {
        string inputMeshFile(pathToFolderOfFiles);
        inputMeshFile.append("hybrid.msh");
        { // check read
            int numNodes, numElems;
            double *nodeCoors;
            int *nodeIDs, *numNodesPerElem, *elemTable, *elemIDs;

            GiDFileIO::readDotMsh(inputMeshFile, numNodes, numElems, nodeCoors, nodeIDs,
                    numNodesPerElem, elemTable, elemIDs);

            const int nodeIDsRef[] = { 1, 2, 3, 5, 6 };
            const double nodeCoorsRef[] = { 0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0.5,
                    1.5, 0. };
            const int numNodesPerElemRef[] = { 4, 3 };
            const int elemTableRef[] = { 1, 2, 3, 5, 5, 3, 6 };
            const int elemIDsRef[] = { 1, 2 };

            CPPUNIT_ASSERT((numNodes == 5));
            CPPUNIT_ASSERT((numElems == 2));
            for (int i = 0; i < numNodes; i++)
                CPPUNIT_ASSERT((nodeIDs[i] == nodeIDsRef[i]));
            for (int i = 0; i < numNodes * 3; i++)
                CPPUNIT_ASSERT((nodeCoors[i] == nodeCoorsRef[i]));
            for (int i = 0; i < numElems; i++)
                CPPUNIT_ASSERT((elemIDs[i] == elemIDsRef[i]));
            for (int i = 0; i < numElems; i++)
                CPPUNIT_ASSERT((numNodesPerElem[i] == numNodesPerElemRef[i]));
            int elemTableSize = 0;
            for (int i = 0; i < numElems; i++)
                elemTableSize += numNodesPerElem[i];
            for (int i = 0; i < elemTableSize; i++)
                CPPUNIT_ASSERT((elemTable[i] == elemTableRef[i]));

            string outputMeshFile("GiDFileIO_unittest_output.msh");
            GiDFileIO::writeDotMsh(outputMeshFile, numNodes, numElems, nodeCoors, nodeIDs,
                    numNodesPerElem, elemTable, elemIDs);
            delete[] nodeCoors, nodeIDs, numNodesPerElem, elemTable, elemIDs;
        }
        { // check write
            string inputMeshFile("GiDFileIO_unittest_output.msh");
            int numNodes, numElems;
            double *nodeCoors;
            int *nodeIDs, *numNodesPerElem, *elemTable, *elemIDs;

            GiDFileIO::readDotMsh(inputMeshFile, numNodes, numElems, nodeCoors, nodeIDs,
                    numNodesPerElem, elemTable, elemIDs);

            const int nodeIDsRef[] = { 1, 2, 3, 5, 6 };
            const double nodeCoorsRef[] = { 0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0.5,
                    1.5, 0. };
            const int numNodesPerElemRef[] = { 3, 4 };
            const int elemTableRef[] = { 5, 3, 6, 1, 2, 3, 5 };
            const int elemIDsRef[] = { 2, 1 };

            CPPUNIT_ASSERT((numNodes == 5));
            CPPUNIT_ASSERT((numElems == 2));
            for (int i = 0; i < numNodes; i++)
                CPPUNIT_ASSERT((nodeIDs[i] == nodeIDsRef[i]));
            for (int i = 0; i < numNodes * 3; i++)
                CPPUNIT_ASSERT((nodeCoors[i] == nodeCoorsRef[i]));
            for (int i = 0; i < numElems; i++)
                CPPUNIT_ASSERT((elemIDs[i] == elemIDsRef[i]));
            for (int i = 0; i < numElems; i++)
                CPPUNIT_ASSERT((numNodesPerElem[i] == numNodesPerElemRef[i]));
            int elemTableSize = 0;
            for (int i = 0; i < numElems; i++)
                elemTableSize += numNodesPerElem[i];
            for (int i = 0; i < elemTableSize; i++)
                CPPUNIT_ASSERT((elemTable[i] == elemTableRef[i]));
            delete[] nodeCoors, nodeIDs, numNodesPerElem, elemTable, elemIDs;
        }
    }
    /***********************************************************************************************
     * \brief Test case: Test read/write in GiD .res format.
     * \author Tianyang Wang
     ***********/
    void readWriteGiDRes() {
        string inputMeshFile(pathToFolderOfFiles);
        inputMeshFile.append("hybrid.msh");
        int numNodes, numElems;
        double *nodeCoors;
        int *nodeIDs, *numNodesPerElem, *elemTable, *elemIDs;
        GiDFileIO::readDotMsh(inputMeshFile, numNodes, numElems, nodeCoors, nodeIDs,
                numNodesPerElem, elemTable, elemIDs);

        { // test read
            string inputResFile(pathToFolderOfFiles);
            inputResFile.append("hybrid.res");
            string resultName1 = "\"nodal data\"";
            string analysisName1 = "\"EMPIRE_CoSimulation\"";
            int stepNum1 = 5;
            string type1 = "vector";
            double *nodalData = new double[numNodes * 3];

            string resultName2 = "\"elemental data\"";
            string analysisName2 = "\"EMPIRE_CoSimulation\"";
            int stepNum2 = 6;
            string type2 = "vector";
            double *elementalData = new double[numElems * 3];

            GiDFileIO::readNodalDataFromDotRes(inputResFile, resultName1, analysisName1, stepNum1,
                    type1, numNodes, nodeIDs, nodalData);

            GiDFileIO::readElementalDataFromDotRes(inputResFile, resultName2, analysisName2,
                    stepNum2, type2, numElems, elemIDs, numNodesPerElem, elementalData);

            for (int i = 0; i < numNodes; i++) {
                CPPUNIT_ASSERT(nodalData[i*3] == double(nodeIDs[i]));
                CPPUNIT_ASSERT(nodalData[i*3+1] == 0.0);
                CPPUNIT_ASSERT(nodalData[i*3+2] == 0.0);
            }
            for (int i = 0; i < numElems; i++) {
                CPPUNIT_ASSERT(elementalData[i*3] == double(elemIDs[i]) * 10.0);
                CPPUNIT_ASSERT(elementalData[i*3+1] == 0.0);
                CPPUNIT_ASSERT(elementalData[i*3+2] == 0.0);
            }

            string outputResFile("GiDFileIO_unittest_output.res");
            GiDFileIO::initDotRes(outputResFile);
            GiDFileIO::appendNodalDataToDotRes(outputResFile, resultName1, analysisName1, stepNum1,
                    type1, numNodes, nodeIDs, nodalData);
            GiDFileIO::appendElementalDataToDotRes(outputResFile, resultName2, analysisName2,
                    stepNum2, type2, numElems, elemIDs, numNodesPerElem, elementalData);

            delete[] nodalData, elementalData;
        }
        { // test write
            string inputResFile("GiDFileIO_unittest_output.res");
            string resultName1 = "\"nodal data\"";
            string analysisName1 = "\"EMPIRE_CoSimulation\"";
            int stepNum1 = 5;
            string type1 = "vector";
            double *nodalData = new double[numNodes * 3];

            string resultName2 = "\"elemental data\"";
            string analysisName2 = "\"EMPIRE_CoSimulation\"";
            int stepNum2 = 6;
            string type2 = "vector";
            double *elementalData = new double[numElems * 3];

            GiDFileIO::readNodalDataFromDotRes(inputResFile, resultName1, analysisName1, stepNum1,
                    type1, numNodes, nodeIDs, nodalData);

            GiDFileIO::readElementalDataFromDotRes(inputResFile, resultName2, analysisName2,
                    stepNum2, type2, numElems, elemIDs, numNodesPerElem, elementalData);

            for (int i = 0; i < numNodes; i++) {
                CPPUNIT_ASSERT(nodalData[i*3] == double(nodeIDs[i]));
                CPPUNIT_ASSERT(nodalData[i*3+1] == 0.0);
                CPPUNIT_ASSERT(nodalData[i*3+2] == 0.0);
            }
            for (int i = 0; i < numElems; i++) {
                CPPUNIT_ASSERT(elementalData[i*3] == double(elemIDs[i]) * 10.0);
                CPPUNIT_ASSERT(elementalData[i*3+1] == 0.0);
                CPPUNIT_ASSERT(elementalData[i*3+2] == 0.0);
            }
            delete[] nodalData, elementalData;
        }
        delete[] nodeCoors, nodeIDs, numNodesPerElem, elemTable, elemIDs;
    }
CPPUNIT_TEST_SUITE( TestGiDFileIO );
        CPPUNIT_TEST( readWriteGiDMsh);
        CPPUNIT_TEST( readWriteGiDRes);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestGiDFileIO);
