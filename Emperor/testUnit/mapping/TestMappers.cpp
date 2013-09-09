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
#include <vector>

#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

#include "MappingFilter.h"
#include "DataField.h"
#include "FEMesh.h"
#include "MapperAdapter.h"
#include "NearestNeighborMapper.h"
#include "NearestElementMapper.h"
#include "BarycentricInterpolationMapper.h"
#include "MapperAdapter.h"
#include "ConnectionIO.h"
#include "ConnectionIOSetup.h"
#include "Message.h"
#include "GiDFileIO.h"

using namespace std;

namespace EMPIRE {

class TestMappers: public CppUnit::TestFixture {
private:
    FEMesh *meshQuadA;
    FEMesh *meshQuadACopy;
    FEMesh *meshQuadB;
    bool debugMe;

public:
    void setUp() {
        debugMe = false;
        int numNodesPerElemQuad = 4;
        {
            /*
             * 4-----5-----6
             * |     |     | mesh quad 1
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
            /*
             * 4-----5-----6
             * |     |     | copy of mesh quad 1
             * 1-----2-----3
             */
            int numNodesQuad1 = 6;
            int numElemsQuad1 = 2;
            meshQuadACopy = new FEMesh("", numNodesQuad1, numElemsQuad1);
            for (int i = 0; i < numElemsQuad1; i++)
                meshQuadACopy->numNodesPerElem[i] = numNodesPerElemQuad;
            meshQuadACopy->initElems();

            for (int i = 0; i < numNodesQuad1; i++)
                meshQuadACopy->nodeIDs[i] = meshQuadA->nodeIDs[i];

            for (int i = 0; i < 6 * 3; i++)
                meshQuadACopy->nodes[i] = meshQuadA->nodes[i];
            for (int i = 0; i < 2 * 4; i++)
                meshQuadACopy->elems[i] = meshQuadA->elems[i];
        }
        {
            /*
             * 5---6---7---8
             * |   |   |   | mesh quad 2
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
    }
    void tearDown() {
        delete meshQuadA;
        delete meshQuadB;
        delete meshQuadACopy;
    }
    /***********************************************************************************************
     * \brief Test case: mapping between matching meshes
     ***********/
    void testMappingOnMatchingMeshes() {
        DataField *d_A, *d_B, *f_A, *f_B;
        d_A = new DataField("d_A", EMPIRE_DataField_atNode, meshQuadA->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        d_B = new DataField("d_B", EMPIRE_DataField_atNode, meshQuadACopy->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        f_A = new DataField("d_A", EMPIRE_DataField_atNode, meshQuadA->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        f_B = new DataField("d_B", EMPIRE_DataField_atNode, meshQuadACopy->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        d_A->data[0] = 0;
        d_A->data[1] = 1;
        d_A->data[2] = 2;
        d_A->data[3] = 3;
        d_A->data[4] = 4;
        d_A->data[5] = 5;

        f_B->data[0] = 5;
        f_B->data[1] = 4;
        f_B->data[2] = 3;
        f_B->data[3] = 2;
        f_B->data[4] = 1;
        f_B->data[5] = 0;

        vector<MapperAdapter *> mappers;
        { // NearestNeighborMapper
            MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadACopy);
            mapper->initNearestNeighborMapper();
            mappers.push_back(mapper);
        }
        { // BarycentricInterpolationMapper
            MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadACopy);
            mapper->initBarycentricInterpolationMapper();
            mappers.push_back(mapper);
        }
        { // NearestElementMapper
            MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadACopy);
            mapper->initNearestElementMapper();
            mappers.push_back(mapper);
        }
        { // MortarMapper normal
#ifdef USE_INTEL_MKL
        MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadACopy);
        mapper->initMortarMapper(false, false, false);
        mappers.push_back(mapper);
#endif
        }
        { // MortarMapper dual
#ifdef LAPACK
        MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadACopy);
        mapper->initMortarMapper(false, true, false);
        mappers.push_back(mapper);
#endif
        }

        for (int i = 0; i < mappers.size(); i++) {
            { // consistent mapping
                AbstractFilter *filterConsistent = new MappingFilter(mappers[i]);
                ConnectionIOSetup::setupIOForFilter(filterConsistent, meshQuadA, d_A, meshQuadACopy,
                        d_B);
                filterConsistent->filtering();
                if (debugMe)
                    infoOut << *d_B;

                const double EPS = 1e-10;
                CPPUNIT_ASSERT(fabs(d_B->data[0] - 0.0) < EPS);
                CPPUNIT_ASSERT(fabs(d_B->data[1] - 1.0) < EPS);
                CPPUNIT_ASSERT(fabs(d_B->data[2] - 2.0) < EPS);
                CPPUNIT_ASSERT(fabs(d_B->data[3] - 3.0) < EPS);
                CPPUNIT_ASSERT(fabs(d_B->data[4] - 4.0) < EPS);
                CPPUNIT_ASSERT(fabs(d_B->data[5] - 5.0) < EPS);

                delete filterConsistent;
            }

            { // conservative mapping
                AbstractFilter *filterConservative = new MappingFilter(mappers[i]);
                ConnectionIOSetup::setupIOForFilter(filterConservative, meshQuadACopy, f_B,
                        meshQuadA, f_A);
                filterConservative->filtering();
                if (debugMe)
                    infoOut << *f_A;

                const double EPS = 1e-10;
                CPPUNIT_ASSERT(fabs(f_A->data[0] - 5.0) < EPS);
                CPPUNIT_ASSERT(fabs(f_A->data[1] - 4.0) < EPS);
                CPPUNIT_ASSERT(fabs(f_A->data[2] - 3.0) < EPS);
                CPPUNIT_ASSERT(fabs(f_A->data[3] - 2.0) < EPS);
                CPPUNIT_ASSERT(fabs(f_A->data[4] - 1.0) < EPS);
                CPPUNIT_ASSERT(fabs(f_A->data[5] - 0.0) < EPS);

                delete filterConservative;
            }
        }

        for (int i = 0; i < mappers.size(); i++) {
            delete mappers[i];
        }

        delete d_A, d_B, f_A, f_B;
    }

    /***********************************************************************************************
     * \brief Test case: mapping constant field to check consistency
     ***********/
    void testConsistency() {
        DataField *d_A;
        DataField *d_B;
        d_A = new DataField("d_A", EMPIRE_DataField_atNode, meshQuadA->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        d_B = new DataField("d_B", EMPIRE_DataField_atNode, meshQuadB->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);

        d_A->data[0] = 1;
        d_A->data[1] = 1;
        d_A->data[2] = 1;
        d_A->data[3] = 1;
        d_A->data[4] = 1;
        d_A->data[5] = 1;

        vector<MapperAdapter *> mappers;
        { // NearestNeighborMapper
            MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadB);
            mapper->initNearestNeighborMapper();
            mappers.push_back(mapper);
        }
        { // BarycentricInterpolationMapper
            MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadB);
            mapper->initBarycentricInterpolationMapper();
            mappers.push_back(mapper);
        }
        { // NearestElementMapper
            MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadB);
            mapper->initNearestElementMapper();
            mappers.push_back(mapper);
        }
        { // MortarMapper normal
#ifdef USE_INTEL_MKL
        MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadB);
        mapper->initMortarMapper(false, false, false);
        mappers.push_back(mapper);
#endif
        }
        { // MortarMapper dual
#ifdef LAPACK
        MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadB);
        mapper->initMortarMapper(false, true, false);
        mappers.push_back(mapper);
#endif
        }

        for (int i = 0; i < mappers.size(); i++) {
            { // consistent mapping
                AbstractFilter *filterConsistent = new MappingFilter(mappers[i]);
                ConnectionIOSetup::setupIOForFilter(filterConsistent, meshQuadA, d_A, meshQuadB,
                        d_B);
                filterConsistent->filtering();
                if (debugMe)
                    infoOut << *d_B;

                const double EPS = 1e-10;
                for (int i = 0; i < d_B->numLocations; i++)
                    CPPUNIT_ASSERT(fabs(d_B->data[i] - 1.0) < EPS);

                delete filterConsistent;
            }
        }
        for (int i = 0; i < mappers.size(); i++) {
            delete mappers[i];
        }

        delete d_A, d_B;
    }

    /***********************************************************************************************
     * \brief Test case: Check the conservation of force and energy
     ***********/
    void testConservation() {
        DataField *d_A;
        DataField *d_B;
        DataField *f_A;
        DataField *f_B;
        d_A = new DataField("d_A", EMPIRE_DataField_atNode, meshQuadA->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        d_B = new DataField("d_B", EMPIRE_DataField_atNode, meshQuadB->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_field);
        f_A = new DataField("d_A", EMPIRE_DataField_atNode, meshQuadA->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        f_B = new DataField("d_B", EMPIRE_DataField_atNode, meshQuadB->numNodes,
                EMPIRE_DataField_scalar, EMPIRE_DataField_fieldIntegral);
        // random d_A
        d_A->data[0] = 1;
        d_A->data[1] = 9;
        d_A->data[2] = 5;
        d_A->data[3] = 4;
        d_A->data[4] = 6;
        d_A->data[5] = 3;

        // random f_B
        f_B->data[0] = 3;
        f_B->data[1] = 5;
        f_B->data[2] = 8;
        f_B->data[3] = 6;
        f_B->data[4] = 7;
        f_B->data[5] = 2;
        f_B->data[6] = 6;
        f_B->data[7] = 3;

        vector<MapperAdapter *> mappers;
        { // NearestNeighborMapper
            MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadB);
            mapper->initNearestNeighborMapper();
            mappers.push_back(mapper);
        }
        { // BarycentricInterpolationMapper
            MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadB);
            mapper->initBarycentricInterpolationMapper();
            mappers.push_back(mapper);
        }
        { // NearestElementMapper
            MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadB);
            mapper->initNearestElementMapper();
            mappers.push_back(mapper);
        }
        { // MortarMapper normal
#ifdef USE_INTEL_MKL
        MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadB);
        mapper->initMortarMapper(false, false, false);
        mappers.push_back(mapper);
#endif
        }
        { // MortarMapper dual
#ifdef LAPACK
        MapperAdapter *mapper = new MapperAdapter("", meshQuadA, meshQuadB);
        mapper->initMortarMapper(false, true, false);
        mappers.push_back(mapper);
#endif
        }

        for (int i = 0; i < mappers.size(); i++) {
            { // consistent mapping
                AbstractFilter *filterConsistent = new MappingFilter(mappers[i]);
                ConnectionIOSetup::setupIOForFilter(filterConsistent, meshQuadA, d_A, meshQuadB,
                        d_B);
                filterConsistent->filtering();
                if (debugMe)
                    infoOut << *d_B;

                delete filterConsistent;
            }
            { // conservative mapping
                AbstractFilter *filterConservative = new MappingFilter(mappers[i]);
                ConnectionIOSetup::setupIOForFilter(filterConservative, meshQuadB, f_B, meshQuadA,
                        f_A);
                filterConservative->filtering();
                if (debugMe)
                    infoOut << *f_A;

                delete filterConservative;
            }
            const double EPS = 1e-10;
            { // check force conservation
                double sumForceA = 0.0;
                for (int i = 0; i < f_A->numLocations; i++)
                    sumForceA += f_A->data[i];
                double sumForceB = 0.0;
                for (int i = 0; i < f_B->numLocations; i++)
                    sumForceB += f_B->data[i];
                CPPUNIT_ASSERT(fabs(sumForceA - sumForceB) < EPS);
            }
            { // check energy conservation
                double sumEnergyA = 0.0;
                for (int i = 0; i < f_A->numLocations; i++)
                    sumEnergyA += f_A->data[i] * d_A->data[i];
                double sumEnergyB = 0.0;
                for (int i = 0; i < f_B->numLocations; i++)
                    sumEnergyB += f_B->data[i] * d_B->data[i];
                CPPUNIT_ASSERT(fabs(sumEnergyA - sumEnergyB) < EPS);
            }
        }
        for (int i = 0; i < mappers.size(); i++) {
            delete mappers[i];
        }

        delete d_A, d_B, f_A, f_B;
    }

CPPUNIT_TEST_SUITE( TestMappers );
        CPPUNIT_TEST( testMappingOnMatchingMeshes);
        CPPUNIT_TEST( testConsistency);
        CPPUNIT_TEST( testConservation);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestMappers);
