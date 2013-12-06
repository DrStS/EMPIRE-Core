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

#include <math.h>
#include <iostream>
#include <vector>
#include <assert.h>
#include "IGAMortarMath.h"

namespace EMPIRE {
using namespace std;
using namespace IGAMortarMath;



/********//**
 * \brief Test functions in MortarMath.h
 ***********/
class TestIGAMortarMath: public CppUnit::TestFixture {
public:
    void setUp() {
    }

    void tearDown() {
    }

    void testsolve3x3LinearSystem(){
    	double A[9] = {1,     1,     1,
    					1,     0,     0,
    					0,     1,     0};
    	double b[3] = {1,2,3};
    	solve3x3LinearSystem(A, b, 1e-14);
    	CPPUNIT_ASSERT(fabs(b[0]-2)<1e-12);
    	CPPUNIT_ASSERT(fabs(b[1]-3)<1e-12);
    	CPPUNIT_ASSERT(fabs(b[2]+4)<1e-12);
    }


    void testComputeIntersectionBetweenLineAndTriangle(){
    	double X0[3] = {0.4, 0.5, 0.0};
    	double n[3] = {0.0, -1.0, 1.0};
    	double X[9] = {0.0, 0.0, 0.0,
    				   0.0, 1.0, 1.0,
    				   1.0, 0.0, 1.0};
    	double b[2];
    	computeIntersectionBetweenLineAndTriangle(X, X0, n, b);
    	CPPUNIT_ASSERT(fabs(b[0]-(0.55))<1e-12);
    	CPPUNIT_ASSERT(fabs(b[1]-(0.05))<1e-12);
    }


    void testComputeIntersectionBetweenLineAndQuad(){
    	double X0[3] = {0.4,0.6,0.0};
    	double n[3] = {0.0,1.0,1.0};
    	double X[12] = {0.0, 0.0, 1.0,
    			    1.0, 0.0, 0.0,
    			    1.0, 1.0, 0.0,
    			    0.0, 1.0, 0.0};
    	double b[2];
    	computeIntersectionBetweenLineAndQuad(X, X0, n, b);
    	CPPUNIT_ASSERT(fabs(b[0]-(-0.2))<1e-12);
    	CPPUNIT_ASSERT(fabs(b[1]-(0.5))<1e-12);
    }

CPPUNIT_TEST_SUITE( TestIGAMortarMath );
		CPPUNIT_TEST(testsolve3x3LinearSystem);
        CPPUNIT_TEST(testComputeIntersectionBetweenLineAndTriangle);
        CPPUNIT_TEST(testComputeIntersectionBetweenLineAndQuad);
        //CPPUNIT_TEST(testComputeLocalCoorInQuad);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestIGAMortarMath);
