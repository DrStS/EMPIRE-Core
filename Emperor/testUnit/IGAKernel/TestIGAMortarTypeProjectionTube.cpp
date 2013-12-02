/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Munich
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
// inclusion of standard libraries   (only if really necessary here in *.h)
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"
#include <iostream>
#include <string>
#include <math.h>

// Inclusion of user-defined libraries
#include "IGAPatchSurface.h"

using namespace std;

namespace EMPIRE {

/********//**
 * \brief Test the IGA mortar-type projection for the case of a curved shell structure (tube-like set up)
 ***********/

class TestIGAMortarTypeProjectionTube: public CppUnit::TestFixture {

private:
    IGAPatchSurface* igaPatch2D;
    double Tol;
    double relTol;
    double TolDeriv;

public:
    void setUp() {
        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the functional values
        Tol = 1e-15;

        // Assign a relaxed tolerance value (corresponding to maximum accuracy provided by MATLAB) for the Newton-Rapson iteration error (accumulative error appears here)
        relTol = 1e-14;

        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the derivative functional values
        TolDeriv = 1e-13;

        // Provide an id for the basis
        int id_basis = 3;

        // Provide an id for the patch itself
        int id_patch = 1;

        // The polynomial degrees
        int p = 2;
        int q = 4;

        // Number of knots in both directions
        int uNoKnots = 17;
        int vNoKnots = 14;

        // The knot vectors in each direction
        double* uKnotVector = new double[uNoKnots];
        for (int i = 0; i < 3; i++)
            uKnotVector[i] = 0.0;
        for (int i = 3; i < 4; i++)
            uKnotVector[i] = 0.083333333333333;
        for (int i = 4; i < 5; i++)
            uKnotVector[i] = 0.166666666666667;
        for (int i = 5; i < 6; i++)
            uKnotVector[i] = 0.250000000000000;
        for (int i = 6; i < 7; i++)
            uKnotVector[i] = 0.333333333333333;
        for (int i = 7; i < 8; i++)
            uKnotVector[i] = 0.416666666666667;
        for (int i = 8; i < 9; i++)
            uKnotVector[i] = 0.500000000000000;
        for (int i = 9; i < 10; i++)
            uKnotVector[i] = 0.583333333333333;
        for (int i = 10; i < 11; i++)
            uKnotVector[i] = 0.666666666666667;
        for (int i = 11; i < 12; i++)
            uKnotVector[i] = 0.750000000000000;
        for (int i = 12; i < 13; i++)
            uKnotVector[i] = 0.833333333333333;
        for (int i = 13; i < 14; i++)
            uKnotVector[i] = 0.916666666666667;
        for (int i = 14; i < 17; i++)
            uKnotVector[i] = 1.0;

        double* vKnotVector = new double[vNoKnots];
        for (int i = 0; i < 5; i++)
            vKnotVector[i] = 0.0;
        for (int i = 5; i < 6; i++)
            vKnotVector[i] = 0.200000000000000;
        for (int i = 6; i < 7; i++)
            vKnotVector[i] = 0.400000000000000;
        for (int i = 7; i < 8; i++)
            vKnotVector[i] = 0.600000000000000;
        for (int i = 8; i < 9; i++)
            vKnotVector[i] = 0.800000000000000;
        for (int i = 9; i < 14; i++)
            vKnotVector[i] = 1.0;

        // The Control Point net
        int uNoControlPoints = uNoKnots - p - 1;
        int vNoControlPoints = vNoKnots - q - 1;
        IGAControlPoint** controlPointNet = new IGAControlPoint*[uNoControlPoints * vNoControlPoints];
        double* controlPointWeights = new double[uNoControlPoints * vNoControlPoints];

        // Control Points for a NURBS patch

        // First row
        controlPointNet[0] = new IGAControlPoint(0, 0.0, -1.250000000000000, 10.000000000000000,
                1.000000000000000);
        controlPointNet[1] = new IGAControlPoint(1, 0.0, -1.500000000000000, 9.600000000000001,
                1.000000000000000);
        controlPointNet[2] = new IGAControlPoint(2, 0.0, -1.700000000000000, 8.800000000000001,
                1.000000000000000);
        controlPointNet[3] = new IGAControlPoint(3, 0.0, -1.280000000000000, 7.792000000000000,
                1.000000000000000);
        controlPointNet[4] = new IGAControlPoint(4, 0.0, 0.0, 7.292800000000001, 1.000000000000000);
        controlPointNet[5] = new IGAControlPoint(5, 0.0, 1.280000000000000, 7.792000000000001,
                1.000000000000000);
        controlPointNet[6] = new IGAControlPoint(6, 0.0, 1.700000000000000, 8.800000000000001,
                1.000000000000000);
        controlPointNet[7] = new IGAControlPoint(7, 0.0, 1.500000000000000, 9.600000000000000,
                1.000000000000000);
        controlPointNet[8] = new IGAControlPoint(8, 0.0, 1.250000000000000, 10.000000000000000,
                1.000000000000000);

        // Second row
        controlPointNet[9] = new IGAControlPoint(9, 0.301998945769793, -1.250000000000000,
                10.000000000000002, 0.975592231765546);
        controlPointNet[10] = new IGAControlPoint(10, 0.301998945769793, -1.500000000000000,
                9.600000000000001, 0.975592231765546);
        controlPointNet[11] = new IGAControlPoint(11, 0.301998945769793, -1.700000000000000,
                8.799999999999999, 0.975592231765546);
        controlPointNet[12] = new IGAControlPoint(12, 0.301998945769793, -1.280000000000000,
                7.792000000000000, 0.975592231765546);
        controlPointNet[13] = new IGAControlPoint(13, 0.301998945769793, 0.0, 7.292800000000000,
                0.975592231765546);
        controlPointNet[14] = new IGAControlPoint(14, 0.301998945769793, 1.280000000000000,
                7.79200000000000, 0.975592231765546);
        controlPointNet[15] = new IGAControlPoint(15, 0.301998945769793, 1.700000000000000,
                8.800000000000001, 0.975592231765546);
        controlPointNet[16] = new IGAControlPoint(16, 0.301998945769793, 1.500000000000000,
                9.600000000000001, 0.975592231765546);
        controlPointNet[17] = new IGAControlPoint(17, 0.301998945769793, 1.250000000000000,
                10.000000000000002, 0.975592231765546);

        // Third row
        controlPointNet[18] = new IGAControlPoint(18, 0.988929951704084, -1.250000000000000,
                10.000000000000002, 0.934912618041455);
        controlPointNet[19] = new IGAControlPoint(19, 0.988929951704084, -1.500000000000000,
                9.600000000000001, 0.934912618041455);
        controlPointNet[20] = new IGAControlPoint(20, 0.988929951704084, -1.500000000000000,
                8.800000000000002, 0.934912618041455);
        controlPointNet[21] = new IGAControlPoint(21, 0.988929951704084, -1.280000000000000,
                7.792000000000001, 0.934912618041455);
        controlPointNet[22] = new IGAControlPoint(22, 0.988929951704084, 0.0, 7.292800000000001,
                0.934912618041455);
        controlPointNet[23] = new IGAControlPoint(23, 0.988929951704084, 1.280000000000000,
                7.792000000000000, 0.934912618041455);
        controlPointNet[24] = new IGAControlPoint(24, 0.988929951704084, 1.700000000000000,
                8.799999999999999, 0.934912618041455);
        controlPointNet[25] = new IGAControlPoint(25, 0.988929951704084, 1.500000000000000,
                9.600000000000001, 0.934912618041455);
        controlPointNet[26] = new IGAControlPoint(26, 0.988929951704084, 1.250000000000000,
                10.000000000000002, 0.934912618041455);

        // Fourth row
        controlPointNet[27] = new IGAControlPoint(27, 1.767766952966369, -1.250000000000000,
                10.000000000000002, 0.902368927062182);
        controlPointNet[28] = new IGAControlPoint(28, 1.767766952966369, -1.500000000000000,
                9.600000000000003, 0.902368927062182);
        controlPointNet[29] = new IGAControlPoint(29, 1.767766952966369, -1.700000000000000,
                8.800000000000002, 0.902368927062182);
        controlPointNet[30] = new IGAControlPoint(30, 1.767766952966369, -1.280000000000000,
                7.792000000000001, 0.902368927062182);
        controlPointNet[31] = new IGAControlPoint(31, 1.767766952966369, 0.0, 7.292799999999999,
                0.902368927062182);
        controlPointNet[32] = new IGAControlPoint(32, 1.767766952966369, 1.280000000000000,
                7.792000000000000, 0.902368927062182);
        controlPointNet[33] = new IGAControlPoint(33, 1.767766952966369, 1.700000000000000,
                8.800000000000001, 0.902368927062182);
        controlPointNet[34] = new IGAControlPoint(34, 1.767766952966369, 1.500000000000000,
                9.600000000000003, 0.902368927062182);
        controlPointNet[35] = new IGAControlPoint(35, 1.767766952966369, 1.250000000000000,
                10.000000000000002, 0.902368927062182);

        // Fifth row
        controlPointNet[36] = new IGAControlPoint(36, 2.627078017762149, -1.250000000000000,
                10.000000000000002, 0.877961158827728);
        controlPointNet[37] = new IGAControlPoint(37, 2.627078017762149, -1.500000000000000,
                9.600000000000000, 0.877961158827728);
        controlPointNet[38] = new IGAControlPoint(38, 2.627078017762149, -1.700000000000000,
                8.800000000000001, 0.877961158827728);
        controlPointNet[39] = new IGAControlPoint(39, 2.627078017762149, -1.280000000000000,
                7.792000000000001, 0.877961158827728);
        controlPointNet[40] = new IGAControlPoint(40, 2.627078017762149, 0.0, 7.292800000000002,
                0.877961158827728);
        controlPointNet[41] = new IGAControlPoint(41, 2.627078017762149, 1.280000000000000,
                7.792000000000000, 0.877961158827728);
        controlPointNet[42] = new IGAControlPoint(42, 2.627078017762149, 1.700000000000000,
                8.800000000000001, 0.877961158827728);
        controlPointNet[43] = new IGAControlPoint(43, 2.627078017762149, 1.500000000000000,
                9.600000000000001, 0.877961158827728);
        controlPointNet[44] = new IGAControlPoint(44, 2.627078017762149, 1.250000000000000,
                10.000000000000002, 0.877961158827728);

        // Sixth row
        controlPointNet[45] = new IGAControlPoint(45, 3.549361143684567, -1.250000000000000,
                10.000000000000002, 0.861689313338092);
        controlPointNet[46] = new IGAControlPoint(46, 3.549361143684567, -1.500000000000000,
                9.600000000000003, 0.861689313338092);
        controlPointNet[47] = new IGAControlPoint(47, 3.549361143684567, -1.700000000000000,
                8.800000000000001, 0.861689313338092);
        controlPointNet[48] = new IGAControlPoint(48, 3.549361143684567, -1.280000000000000,
                7.792000000000000, 0.861689313338092);
        controlPointNet[49] = new IGAControlPoint(49, 3.549361143684567, 0.0, 7.292800000000001,
                0.861689313338092);
        controlPointNet[50] = new IGAControlPoint(50, 3.549361143684567, 1.280000000000000,
                7.792000000000002, 0.861689313338092);
        controlPointNet[51] = new IGAControlPoint(51, 3.549361143684567, 1.700000000000000,
                8.800000000000002, 0.861689313338092);
        controlPointNet[52] = new IGAControlPoint(52, 3.549361143684567, 1.500000000000000,
                9.600000000000003, 0.861689313338092);
        controlPointNet[53] = new IGAControlPoint(53, 3.549361143684567, 1.250000000000000,
                10.000000000000002, 0.861689313338092);

        // Seventh row
        controlPointNet[54] = new IGAControlPoint(54, 4.511844635310912, -1.250000000000000,
                10.000000000000000, 0.853553390593274);
        controlPointNet[55] = new IGAControlPoint(55, 4.511844635310912, -1.500000000000000,
                9.600000000000000, 0.853553390593274);
        controlPointNet[56] = new IGAControlPoint(56, 4.511844635310912, -1.700000000000000,
                8.799999999999999, 0.853553390593274);
        controlPointNet[57] = new IGAControlPoint(57, 4.511844635310912, -1.280000000000000,
                7.792000000000000, 0.853553390593274);
        controlPointNet[58] = new IGAControlPoint(58, 4.511844635310912, 0.0, 7.292800000000000,
                0.853553390593274);
        controlPointNet[59] = new IGAControlPoint(59, 4.511844635310912, 1.280000000000000,
                7.792000000000001, 0.853553390593274);
        controlPointNet[60] = new IGAControlPoint(60, 4.511844635310912, 1.700000000000000,
                8.800000000000002, 0.853553390593274);
        controlPointNet[61] = new IGAControlPoint(61, 4.511844635310912, 1.500000000000000,
                9.600000000000001, 0.853553390593274);
        controlPointNet[62] = new IGAControlPoint(62, 4.511844635310912, 1.250000000000000,
                10.000000000000000, 0.853553390593274);

        // Eighth row
        controlPointNet[63] = new IGAControlPoint(63, 5.488155364689087, -1.250000000000000,
                10.000000000000000, 0.853553390593274);
        controlPointNet[64] = new IGAControlPoint(64, 5.488155364689087, -1.500000000000000,
                9.600000000000000, 0.853553390593274);
        controlPointNet[65] = new IGAControlPoint(65, 5.488155364689087, -1.700000000000000,
                8.799999999999999, 0.853553390593274);
        controlPointNet[66] = new IGAControlPoint(66, 5.488155364689087, -1.280000000000000,
                7.792000000000000, 0.853553390593274);
        controlPointNet[67] = new IGAControlPoint(67, 5.488155364689087, 0.0, 7.292800000000000,
                0.853553390593274);
        controlPointNet[68] = new IGAControlPoint(68, 5.488155364689087, 1.280000000000000,
                7.792000000000001, 0.853553390593274);
        controlPointNet[69] = new IGAControlPoint(69, 5.488155364689087, 1.700000000000000,
                8.800000000000002, 0.853553390593274);
        controlPointNet[70] = new IGAControlPoint(70, 5.488155364689087, 1.500000000000000,
                9.600000000000001, 0.853553390593274);
        controlPointNet[71] = new IGAControlPoint(71, 5.488155364689087, 1.250000000000000,
                10.000000000000000, 0.853553390593274);

        // Nineth row
        controlPointNet[72] = new IGAControlPoint(72, 6.450638856315432, -1.250000000000000,
                10.000000000000002, 0.861689313338092);
        controlPointNet[73] = new IGAControlPoint(73, 6.450638856315432, -1.500000000000000,
                9.600000000000003, 0.861689313338092);
        controlPointNet[74] = new IGAControlPoint(74, 6.450638856315432, -1.700000000000000,
                8.800000000000001, 0.861689313338092);
        controlPointNet[75] = new IGAControlPoint(75, 6.450638856315432, -1.280000000000000,
                7.792000000000000, 0.861689313338092);
        controlPointNet[76] = new IGAControlPoint(76, 6.450638856315432, 0.0, 7.292800000000001,
                0.861689313338092);
        controlPointNet[77] = new IGAControlPoint(77, 6.450638856315432, 1.280000000000000,
                7.792000000000002, 0.861689313338092);
        controlPointNet[78] = new IGAControlPoint(78, 6.450638856315432, 1.700000000000000,
                8.800000000000002, 0.861689313338092);
        controlPointNet[79] = new IGAControlPoint(79, 6.450638856315432, 1.500000000000000,
                9.600000000000003, 0.861689313338092);
        controlPointNet[80] = new IGAControlPoint(80, 6.450638856315432, 1.250000000000000,
                10.000000000000002, 0.861689313338092);

        // 10th row
        controlPointNet[81] = new IGAControlPoint(81, 7.372921982237848, -1.250000000000000,
                10.000000000000000, 0.877961158827728);
        controlPointNet[82] = new IGAControlPoint(82, 7.372921982237848, -1.500000000000000,
                9.600000000000000, 0.877961158827728);
        controlPointNet[83] = new IGAControlPoint(83, 7.372921982237848, -1.700000000000000,
                8.800000000000001, 0.877961158827728);
        controlPointNet[84] = new IGAControlPoint(84, 7.372921982237848, -1.280000000000000,
                7.792000000000001, 0.877961158827728);
        controlPointNet[85] = new IGAControlPoint(85, 7.372921982237848, 0.0, 7.292800000000002,
                0.877961158827728);
        controlPointNet[86] = new IGAControlPoint(86, 7.372921982237848, 1.280000000000000,
                7.792000000000001, 0.877961158827728);
        controlPointNet[87] = new IGAControlPoint(87, 7.372921982237848, 1.700000000000000,
                8.800000000000001, 0.877961158827728);
        controlPointNet[88] = new IGAControlPoint(88, 7.372921982237848, 1.500000000000000,
                9.600000000000001, 0.877961158827728);
        controlPointNet[89] = new IGAControlPoint(89, 7.372921982237848, 1.250000000000000,
                10.000000000000000, 0.877961158827728);

        // 11th row
        controlPointNet[90] = new IGAControlPoint(90, 8.232233047033629, -1.250000000000000,
                10.000000000000000, 0.902368927062183);
        controlPointNet[91] = new IGAControlPoint(91, 8.232233047033629, -1.500000000000000,
                9.600000000000000, 0.902368927062183);
        controlPointNet[92] = new IGAControlPoint(92, 8.232233047033629, -1.700000000000000,
                8.800000000000001, 0.902368927062183);
        controlPointNet[93] = new IGAControlPoint(93, 8.232233047033629, -1.280000000000000,
                7.792000000000000, 0.902368927062183);
        controlPointNet[94] = new IGAControlPoint(94, 8.232233047033629, 0.0, 7.292800000000001,
                0.902368927062183);
        controlPointNet[95] = new IGAControlPoint(95, 8.232233047033629, 1.280000000000000,
                7.792000000000000, 0.902368927062183);
        controlPointNet[96] = new IGAControlPoint(96, 8.232233047033629, 1.700000000000000,
                8.800000000000001, 0.902368927062183);
        controlPointNet[97] = new IGAControlPoint(97, 8.232233047033629, 1.500000000000000,
                9.600000000000000, 0.902368927062183);
        controlPointNet[98] = new IGAControlPoint(98, 8.232233047033629, 1.250000000000000,
                10.000000000000000, 0.902368927062183);

        // 12th row
        controlPointNet[99] = new IGAControlPoint(99, 9.011070048295913, -1.250000000000000,
                10.000000000000000, 0.934912618041455);
        controlPointNet[100] = new IGAControlPoint(100, 9.011070048295913, -1.500000000000000,
                9.600000000000001, 0.934912618041455);
        controlPointNet[101] = new IGAControlPoint(101, 9.011070048295913, -1.700000000000000,
                8.800000000000001, 0.934912618041455);
        controlPointNet[102] = new IGAControlPoint(102, 9.011070048295913, -1.280000000000000,
                7.792000000000002, 0.934912618041455);
        controlPointNet[103] = new IGAControlPoint(103, 9.011070048295913, 0.0, 7.292800000000001,
                0.934912618041455);
        controlPointNet[104] = new IGAControlPoint(104, 9.011070048295913, 1.280000000000000,
                7.792000000000001, 0.934912618041455);
        controlPointNet[105] = new IGAControlPoint(105, 9.011070048295913, 1.700000000000000,
                8.799999999999999, 0.934912618041455);
        controlPointNet[106] = new IGAControlPoint(106, 9.011070048295913, 1.500000000000000,
                9.600000000000000, 0.934912618041455);
        controlPointNet[107] = new IGAControlPoint(107, 9.011070048295913, 1.250000000000000,
                10.000000000000000, 0.934912618041455);

        // 13th row
        controlPointNet[108] = new IGAControlPoint(108, 9.698001054230206, -1.250000000000000,
                9.999999999999998, 0.975592231765546);
        controlPointNet[109] = new IGAControlPoint(109, 9.698001054230206, -1.500000000000000,
                9.600000000000000, 0.975592231765546);
        controlPointNet[110] = new IGAControlPoint(110, 9.698001054230206, -1.700000000000000,
                8.799999999999999, 0.975592231765546);
        controlPointNet[111] = new IGAControlPoint(111, 9.698001054230206, -1.280000000000000,
                7.792000000000001, 0.975592231765546);
        controlPointNet[112] = new IGAControlPoint(112, 9.698001054230206, 0.0, 7.292800000000001,
                0.975592231765546);
        controlPointNet[113] = new IGAControlPoint(113, 9.698001054230206, 1.280000000000000,
                7.792000000000001, 0.975592231765546);
        controlPointNet[114] = new IGAControlPoint(114, 9.698001054230206, 1.700000000000000,
                8.799999999999999, 0.975592231765546);
        controlPointNet[115] = new IGAControlPoint(115, 9.698001054230206, 1.500000000000000,
                9.600000000000000, 0.975592231765546);
        controlPointNet[116] = new IGAControlPoint(116, 9.698001054230206, 1.250000000000000,
                9.999999999999998, 0.975592231765546);

        // 14th row
        controlPointNet[117] = new IGAControlPoint(117, 10.000000000000000, -1.250000000000000,
                10.000000000000000, 1.000000000000000);
        controlPointNet[118] = new IGAControlPoint(118, 10.000000000000000, -1.500000000000000,
                9.600000000000001, 1.000000000000000);
        controlPointNet[119] = new IGAControlPoint(119, 10.000000000000000, -1.700000000000000,
                8.800000000000001, 1.000000000000000);
        controlPointNet[120] = new IGAControlPoint(120, 10.000000000000000, -1.280000000000000,
                7.792000000000000, 1.000000000000000);
        controlPointNet[121] = new IGAControlPoint(121, 10.000000000000000, 0.0, 7.292800000000001,
                1.000000000000000);
        controlPointNet[122] = new IGAControlPoint(122, 10.000000000000000, 1.280000000000000,
                7.792000000000001, 1.000000000000000);
        controlPointNet[123] = new IGAControlPoint(123, 10.000000000000000, 1.700000000000000,
                8.800000000000001, 1.000000000000000);
        controlPointNet[124] = new IGAControlPoint(124, 10.000000000000000, 1.500000000000000,
                9.600000000000000, 1.000000000000000);
        controlPointNet[125] = new IGAControlPoint(125, 10.000000000000000, 1.250000000000000,
                10.000000000000000, 1.000000000000000);

        for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
            controlPointWeights[i] = controlPointNet[i]->getW();

        // Test just one object of the class (that works pretty also)
        igaPatch2D = new IGAPatchSurface(id_basis, p, uNoKnots, uKnotVector, q, vNoKnots,
                vKnotVector, uNoControlPoints, vNoControlPoints, controlPointNet);

        // nurbsBasis1D->printPolynomialDegree();
        // nurbsBasis1D->printNoKnots();
        // nurbsBasis1D->printKnotVector();
        // nurbsBasis1D->printNoBasisFunctions();
        // nurbsBasis1D->printControlPointNet();
        // NurbsBasis1D nurbsBasis1DCopy = *nurbsBasis1D;
        // delete nurbsBasis1D;

        // Make test on memory leakage (that works pretty fine)
        /* for (int j = 1; j < 1e9; j++) {

         // Copy the knot vectors
         double* uKnotVectorCopy = new double[uNoKnots];
         for (int i = 0; i < uNoKnots; i++)
         uKnotVectorCopy[i] = uKnotVector[i];
         double* vKnotVectorCopy = new double[vNoKnots];
         for (int i = 0; i < vNoKnots; i++)
         vKnotVectorCopy[i] = vKnotVector[i];

         // Copy the Control Point net and the Control Point weights
         IGAControlPoint* controlPointNetCopy = new IGAControlPoint[uNoControlPoints * vNoControlPoints];
         for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
         controlPointNetCopy[i] = controlPointNet[i];

         // Create an object of the class IGAPatch2D
         IGAPatch2D* igaPatch2DCopy = new IGAPatch2D(id_patch, id_basis, p, uNoKnots, uKnotVectorCopy, q,
         vNoKnots, vKnotVectorCopy, uNoControlPoints, vNoControlPoints,
         controlPointNetCopy);

         // Free the memory on the heap from the pointer
         delete igaPatch2DCopy;
         }*/
    }

    void tearDown() {
        delete igaPatch2D;
    }

    /***********************************************************************************************
     * \brief Test case: Test the projection of an arbitrary point on the 2D IGA patch
     ***********/
    void testProjectionOnIGAPatch() {
        // Initial guesses for the Newton-Rapson iteration
        double u = .5;
        double v = .2;

        // The vertex to be projected onto the NURBS patch
        double vertex[] = { 4.609773274121440, -1.500088888888889, 9.093680000000003 };

        // Flag on the convergence of the Newton-Rapson iterations
        bool flag = 1;

        // Compute the orthogonal projection of the point on the NURBS patch
        flag = igaPatch2D->computePointProjectionOnPatch(u, v, vertex);

        /*cout << endl;
         cout << "The projected on the NURBS surface point:" << endl;
         for (int i = 0; i < 3; i++)
         cout << vertex[i] << " ";
         cout << endl;*/

        // Compare the values with the ones from MATLAB
        // On the return flag
        CPPUNIT_ASSERT(flag==1);

        // On the Cartesian components of the orthogonal projection
        double orthogonalProjection[] = { 4.609773274121439, -1.550611048956046, 9.096336045960706 };

        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(fabs(orthogonalProjection[i] - vertex[i])<=relTol);

        // On the surface parametric values of the orthogonal projection
        CPPUNIT_ASSERT(fabs(u - 0.466666666666667)<=relTol);
        CPPUNIT_ASSERT(fabs(v - 0.115890578620562)<=relTol);
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the projection of a point on the IGA patch for leakage
     ***********/
    void testProjectionOnIGAPatch4Leakage() {
        // Initial guesses for the Newton-Rapson iteration
        double u = .5;
        double v = .2;

        // The vertex to be projected onto the NURBS patch
        double vertex[] = { 4.609773274121440, -1.500088888888889, 9.093680000000003 };

        // Flag on the convergence of the Newton-Rapson iterations
        bool flag = 1;

        // Compute the orthogonal projection of the point on the NURBS patch iteratively
        int noIterations = 1e9;
        for (int i = 0; i < noIterations; i++) {
            // Compute the orthogonal projection
            flag = igaPatch2D->computePointProjectionOnPatch(u, v, vertex);

            // Re-set the surface parameters
            u = .5;
            v = .2;

            // Reset the point to be projected on the NURBS surface
            vertex[0] = 4.609773274121440;
            vertex[1] = -1.500088888888889;
            vertex[2] = 9.093680000000003;
        }
    }

// Make the tests
CPPUNIT_TEST_SUITE(TestIGAMortarTypeProjectionTube);
    CPPUNIT_TEST(testProjectionOnIGAPatch);

// Make the tests for leakage
    // CPPUNIT_TEST(testProjectionOnIGAPatch4Leakage);

    CPPUNIT_TEST_SUITE_END()
    ;

};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestIGAMortarTypeProjectionTube);
