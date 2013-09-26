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
#include <math.h>
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"


#include "AuxiliaryParameters.h"
#include "MathLibrary.h"

namespace EMPIRE {
using namespace std;
using namespace MathLibrary;

/********//**
 * \brief This class manages tests the MathLibrary of EMPIRE
 **************************************************************************************************/
class TestMathLibrary: public CppUnit::TestFixture {
private:
    double *vecA;
    double *vecB;
    vector<double> *vecC;
    vector<double> *vecD;
    SparseMatrix<double> *sparseMatSymm;
    SparseMatrix<double> *sparseMat;
public:
    /***********************************************************************************************
     * \brief Set up some test vectors and matrices
     * \author Stefan Sicklinger
     ***********/
    void setUp() {
        vecA = new double[4];
        vecB = new double[4];
        // Setup vectors
        vecA[0] = 1.1;
        vecA[1] = 2.2;
        vecA[2] = 3.3;
        vecA[3] = 4.4;
        vecB[0] = 5.5;
        vecB[1] = 6.6;
        vecB[2] = 7.7;
        vecB[3] = 8.8;
        vecC = new vector<double> (4);
        vecD = new vector<double> (4);
        (*vecC)[0] = 1.1;
        (*vecC)[1] = 2.2;
        (*vecC)[2] = 3.3;
        (*vecC)[3] = 4.4;
        (*vecD)[0] = 5.5;
        (*vecD)[1] = 6.6;
        (*vecD)[2] = 7.7;
        (*vecD)[3] = 8.8;

        // Setup sparse matrices
        // This matrix is positive definite
        sparseMatSymm = new SparseMatrix<double> (4,true);
        (*sparseMatSymm)(0,0)=1.38E2;
        (*sparseMatSymm)(1,1)=1.3E1;
        (*sparseMatSymm)(2,2)=9.1E5;
        (*sparseMatSymm)(3,3)=1.5E3;
        (*sparseMatSymm)(0,3)=1.8E2;
        (*sparseMatSymm)(1,2)=8.36E1;
        sparseMat = new SparseMatrix<double> (4,false);
        (*sparseMat)(0,0)=1.38E2;
        (*sparseMat)(1,1)=1.3E1;
        (*sparseMat)(2,2)=9.1E5;
        (*sparseMat)(3,3)=1.5E3;
        (*sparseMat)(0,3)=1.8E2;
        (*sparseMat)(1,2)=8.36E1;
    }
    /***********************************************************************************************
     * \brief Delete test vectors and matrices
     * \author Stefan Sicklinger
     ***********/
    void tearDown() {
        delete vecA;
        delete vecB;
        delete vecC;
        delete vecD;
        delete sparseMatSymm;
        delete sparseMat;
    }
    /***********************************************************************************************
     * \brief Test dense dot product function
     * \author Stefan Sicklinger
     ***********/
    void testDenseDotProduct() {

        double result =  fabs( computeDenseDotProduct(vecA, vecB, 4) - 84.70);
        CPPUNIT_ASSERT(result <= AuxiliaryParameters::machineEpsilon*100);
        result =  fabs( computeDenseDotProduct(*vecC, *vecD) - 84.70);
        CPPUNIT_ASSERT(result <= AuxiliaryParameters::machineEpsilon*100);
    }
    /***********************************************************************************************
     * \brief Test sparse matrix multiplication with vector
     * \author Stefan Sicklinger
     ***********/
    void testSparseMatrixVectorProduct() {

        vector<double> result =(*sparseMat)*(*vecC);
        CPPUNIT_ASSERT(fabs(result[0] - 9.43800000e+02) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result[1] - 3.04480000e+02) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result[2] - 3.00300000e+06) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result[3] - 6.60000000e+03) < 1000000*AuxiliaryParameters::machineEpsilon);


        double result2[4];
        (*sparseMat).mulitplyVec(vecA,result2,4);
        CPPUNIT_ASSERT(fabs(result2[0] - 9.43800000e+02) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result2[1] - 3.04480000e+02) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result2[2] - 3.00300000e+06) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result2[3] - 6.60000000e+03) < 1000000*AuxiliaryParameters::machineEpsilon);


        result =(*sparseMatSymm)*(*vecC);
        CPPUNIT_ASSERT(fabs(result[0] - 9.43800000e+02) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result[1] - 3.04480000e+02) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result[2] - 3.00318392e+06) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result[3] - 6.79800000e+03) < 1000000*AuxiliaryParameters::machineEpsilon);


        (*sparseMatSymm).mulitplyVec(vecA,result2,4);
        CPPUNIT_ASSERT(fabs(result2[0] - 9.43800000e+02) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result2[1] - 3.04480000e+02) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result2[2] - 3.00318392e+06) < 100000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(result2[3] - 6.79800000e+03) < 1000000*AuxiliaryParameters::machineEpsilon);

    }
    /***********************************************************************************************
     * \brief Test sparse matrix direct solver
     * \author Stefan Sicklinger
     ***********/
    void testSparseDirectSolver() {

        // Test symmetric positive definite call
        double* solution;
        solution = new double[4]; // holds the solution
        (*sparseMatSymm).factorize();
        (*sparseMatSymm).solve(solution,vecA);
        CPPUNIT_ASSERT(fabs(solution[0] - 4.9140893470790373e-03) < 1000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(solution[1] - 1.6930747279417244e-01) < 1000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(solution[2] + 1.1927587610541555e-05) < 1000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(solution[3] - 2.3436426116838494e-03) < 1000*AuxiliaryParameters::machineEpsilon);
        (*sparseMatSymm).cleanPardiso();

        // Test asymmetriccall
        (*sparseMat).factorize();
        (*sparseMat).solve(solution,vecA);
        CPPUNIT_ASSERT(fabs(solution[0] - 4.1449275362318849e-03) < 1000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(solution[1] - 1.6920744885883349e-01) < 1000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(solution[2] - 3.6263736263736258e-06) < 1000*AuxiliaryParameters::machineEpsilon);
        CPPUNIT_ASSERT(fabs(solution[3] - 2.9333333333333334e-03) < 1000*AuxiliaryParameters::machineEpsilon);
        (*sparseMat).cleanPardiso();

    }

    void testSparseDirectSolver4Leakage() {

    	(*sparseMatSymm).factorize();
    	(*sparseMat).factorize();
    	for (int i = 0; i < 1e10; i++){
    		// Test symmetric positive definite call
    		double* solution;
    		solution = new double[4]; // holds the solution
    		(*sparseMatSymm).solve(solution,vecA);
    		CPPUNIT_ASSERT(fabs(solution[0] - 4.9140893470790373e-03) < 1000*AuxiliaryParameters::machineEpsilon);
    		CPPUNIT_ASSERT(fabs(solution[1] - 1.6930747279417244e-01) < 1000*AuxiliaryParameters::machineEpsilon);
    		CPPUNIT_ASSERT(fabs(solution[2] + 1.1927587610541555e-05) < 1000*AuxiliaryParameters::machineEpsilon);
    		CPPUNIT_ASSERT(fabs(solution[3] - 2.3436426116838494e-03) < 1000*AuxiliaryParameters::machineEpsilon);

    		// Test asymmetriccall
    		(*sparseMat).solve(solution,vecA);
    		CPPUNIT_ASSERT(fabs(solution[0] - 4.1449275362318849e-03) < 1000*AuxiliaryParameters::machineEpsilon);
    		CPPUNIT_ASSERT(fabs(solution[1] - 1.6920744885883349e-01) < 1000*AuxiliaryParameters::machineEpsilon);
    		CPPUNIT_ASSERT(fabs(solution[2] - 3.6263736263736258e-06) < 1000*AuxiliaryParameters::machineEpsilon);
    		CPPUNIT_ASSERT(fabs(solution[3] - 2.9333333333333334e-03) < 1000*AuxiliaryParameters::machineEpsilon);

    		delete[] solution;
    	}
        (*sparseMatSymm).cleanPardiso();
        (*sparseMat).cleanPardiso();

    }


    CPPUNIT_TEST_SUITE( TestMathLibrary );
    CPPUNIT_TEST(testDenseDotProduct);
    CPPUNIT_TEST(testSparseMatrixVectorProduct);
    CPPUNIT_TEST(testSparseDirectSolver);
//    CPPUNIT_TEST(testSparseDirectSolver4Leakage);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestMathLibrary);
