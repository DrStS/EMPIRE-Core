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

#include <string>
#include <map>
#include <iostream>

#include "ClientCode.h"
#include "DataField.h"
#include "FEMesh.h"
#include "Signal.h"

namespace EMPIRE {
using namespace std;
/********//**
 * \brief Test the class ClientCode
 ***********/
class TestClientCode: public CppUnit::TestFixture {
private:
    ClientCode *client;
    FEMesh *mesh1;
    FEMesh *mesh2;
    Signal *signal1;
    Signal *signal2;
public:
    void setUp() {
        std::string name("testClient");
        client = new ClientCode(name);
        {
            int numNodes1 = 3;
            int numElems1 = 1;
            int numNodesPerElem1 = 3;
            mesh1 = new FEMesh("dummy1", numNodes1, numElems1);
            client->nameToMeshMap.insert(pair<string, AbstractMesh*>(mesh1->name, mesh1));
        }
        {
            int numNodes2 = 4;
            int numElems2 = 1;
            int numNodesPerElem2 = 4;
            mesh2 = new FEMesh("dummy2", numNodes2, numElems2);
            client->nameToMeshMap.insert(pair<string, AbstractMesh*>(mesh2->name, mesh2));
        }
        {
            signal1 = new Signal("signal1", 1, 1, 5);
            client->nameToSignalMap.insert(pair<string, Signal*>(signal1->name, signal1));
            signal2 = new Signal("signal2", 1, 1, 5);
            client->nameToSignalMap.insert(pair<string, Signal*>(signal2->name, signal2));
        }
    }
    void tearDown() {
        delete client;
        //delete mesh; //cause segmentation fault
    }
    /***********************************************************************************************
     * \brief Test case: Test the name of the mesh
     ***********/
    void testName() {
        //CPPUNIT_ASSERT(client->getName() == "testClient");
    }
    /***********************************************************************************************
     * \brief Test case: Test the mesh pointer
     ***********/
    void testMesh() {
        /*CPPUNIT_ASSERT(client->nameToMeshMap.size()==2);
        CPPUNIT_ASSERT(client->getMeshByName(mesh1->name)== mesh1);
        CPPUNIT_ASSERT(client->getMeshByName(mesh2->name)== mesh2);
        CPPUNIT_ASSERT(client->getSignalByName(signal1->name)== signal1);
        CPPUNIT_ASSERT(client->getSignalByName(signal2->name)== signal2);*/
    }

CPPUNIT_TEST_SUITE( TestClientCode );
        CPPUNIT_TEST( testName);
        CPPUNIT_TEST( testMesh);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestClientCode);
