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

#include "GmshFileIO.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

namespace EMPIRE {
using namespace std;
/********//**
 * \brief Test the class GmshFileIO
 ***********/
class TestGmshFileIO: public CppUnit::TestFixture {
public:
  void setUp() {
  }
  void tearDown() {
  }
  /***********************************************************************************************
   * \brief Test case: Test read/write for Gmsh mesh in .msh format.
   * \author Michael Andre
   ***********/
  void readWriteTestMsh() {
    const int XYZ = 3;
    { // generate input test mesh in .msh format
      ofstream fout("GmshFileIO_unit_input_mesh.msh");
      CPPUNIT_ASSERT(fout);
      fout  << "$Comments"      << endl;
      fout  << "A mesh generated in Gmsh's .msh format by the unit test suite." << endl;
      fout  << "$EndComments"   << endl;
      fout  << "$MeshFormat"    << endl;
      fout  << "2.2 \t 0 \t 8"  << endl;
      fout  << "$EndMeshFormat" << endl; 
      fout  << "$Nodes"         << endl;
      fout  <<      6           << endl;
      fout  << "    1      "    << '\t' << 0.        << '\t' << 0.        << '\t' << 0.        << endl;
      fout  << "    2      "    << '\t' << 1.        << '\t' << 0.        << '\t' << 0.        << endl;
      fout  << "    3      "    << '\t' << 1.        << '\t' << 1.        << '\t' << 0.        << endl;
      fout  << "    4      "    << '\t' << 1000.     << '\t' << 1000.     << '\t' << 1000.     << endl;
      fout  << "    5      "    << '\t' << 0.        << '\t' << 1.        << '\t' << 0.        << endl;
      fout  << "    6      "    << '\t' << 0.5       << '\t' << 1.5       << '\t' << 0.        << endl;
      fout  << "$EndNodes"      << endl;
      fout  << "$Elements"      << endl;
      fout  <<      2           << endl;
      fout  << "    1      "    << '\t' <<     3     << '\t' <<     0     << '\t' 
	    <<      2           << '\t' <<     1     << '\t' <<     3     << '\t' <<     5     << endl;
      fout  << "    2      "    << '\t' <<     2     << '\t' <<     0     << '\t' 
            <<      3           << '\t' <<     5     << '\t' <<     6     << endl;
      fout  << "$EndElements"   << endl;
      fout.close();
    }
    { // check read
      string inputMeshFile("GmshFileIO_unit_input_mesh.msh");
      int numNodes, numElems;
      double *nodeCoors;
      int *nodeIDs, *numNodesPerElem, *elemTable, *elemIDs;

      GmshFileIO::readDotMsh(inputMeshFile, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable, elemIDs);

      const int nodeIDsRef[] = {1, 2, 3, 5, 6};
      const double nodeCoorsRef[] = {0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0.5, 1.5, 0.};
      const int numNodesPerElemRef[] = {4, 3};
      const int elemTableRef[] = { 2, 1, 3, 5, 3, 5, 6 };
      const int elemIDsRef[] = {1, 2};

      CPPUNIT_ASSERT((numNodes == 5));
      CPPUNIT_ASSERT((numElems == 2));
      for (int i=0; i<numNodes; i++)
	CPPUNIT_ASSERT((nodeIDs[i] == nodeIDsRef[i]));
      for (int i=0; i<numNodes*XYZ; i++)
	CPPUNIT_ASSERT((nodeCoors[i] == nodeCoorsRef[i]));
      for (int i=0; i<numElems; i++)
	CPPUNIT_ASSERT((elemIDs[i] == elemIDsRef[i]));
      for (int i=0; i<numElems; i++)
	CPPUNIT_ASSERT((numNodesPerElem[i] == numNodesPerElemRef[i]));
      int elemTableSize = 0;
      for (int i=0; i<numElems; i++)
	elemTableSize += numNodesPerElem[i];
      for (int i=0; i<elemTableSize; i++)
	CPPUNIT_ASSERT((elemTable[i] == elemTableRef[i]));
      
      string outputMeshFile("GmshFileIO_unit_output_mesh.msh");
      GmshFileIO::writeDotMsh(outputMeshFile, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable, elemIDs);
    }
    { // check write
      string inputMeshFile("GmshFileIO_unit_output_mesh.msh");
      int numNodes, numElems;
      double *nodeCoors;
      int *nodeIDs, *numNodesPerElem, *elemTable, *elemIDs;
      
      GmshFileIO::readDotMsh(inputMeshFile, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable, elemIDs);

      const int nodeIDsRef[] = {1, 2, 3, 5, 6};
      const double nodeCoorsRef[] = {0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0.5, 1.5, 0.};
      const int numNodesPerElemRef[] = {3, 4};
      const int elemTableRef[] = { 3, 5, 6, 2, 1, 3, 5 };
      const int elemIDsRef[] = {2, 1};

      CPPUNIT_ASSERT((numNodes == 5));
      CPPUNIT_ASSERT((numElems == 2));
      for (int i=0; i<numNodes; i++)
	CPPUNIT_ASSERT((nodeIDs[i] == nodeIDsRef[i]));
      for (int i=0; i<numNodes*XYZ; i++)
	CPPUNIT_ASSERT((nodeCoors[i] == nodeCoorsRef[i]));
      for (int i=0; i<numElems; i++)
	CPPUNIT_ASSERT((elemIDs[i] == elemIDsRef[i]));
      for (int i=0; i<numElems; i++)
	CPPUNIT_ASSERT((numNodesPerElem[i] == numNodesPerElemRef[i]));
      int elemTableSize = 0;
      for (int i=0; i<numElems; i++)
	elemTableSize += numNodesPerElem[i];
      for (int i=0; i<elemTableSize; i++)
	CPPUNIT_ASSERT((elemTable[i] == elemTableRef[i]));    
    }

    system("rm GmshFileIO_unit_input_mesh.msh");
    system("rm GmshFileIO_unit_output_mesh.msh");
  }

  CPPUNIT_TEST_SUITE( TestGmshFileIO );
  CPPUNIT_TEST( readWriteTestMsh );
  CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestGmshFileIO);
