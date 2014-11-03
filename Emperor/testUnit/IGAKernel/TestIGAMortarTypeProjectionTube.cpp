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
        IGAControlPoint** controlPointNet =
                new IGAControlPoint*[uNoControlPoints * vNoControlPoints];
        double* controlPointWeights = new double[uNoControlPoints * vNoControlPoints];

        // Control Points for a NURBS patch

        controlPointNet[0] = new IGAControlPoint(0, 0.00000000000000, -1.25000000000000,
                10.00000000000000, 1.00000000000000);
        controlPointNet[1] = new IGAControlPoint(1, 0.30199894576979, -1.25000000000000,
                10.00000000000000, 0.97559223176555);
        controlPointNet[2] = new IGAControlPoint(2, 0.98892995170408, -1.25000000000000,
                10.00000000000000, 0.93491261804145);
        controlPointNet[3] = new IGAControlPoint(3, 1.76776695296637, -1.25000000000000,
                10.00000000000000, 0.90236892706218);
        controlPointNet[4] = new IGAControlPoint(4, 2.62707801776215, -1.25000000000000,
                10.00000000000000, 0.87796115882773);
        controlPointNet[5] = new IGAControlPoint(5, 3.54936114368457, -1.25000000000000,
                10.00000000000000, 0.86168931333809);
        controlPointNet[6] = new IGAControlPoint(6, 4.51184463531091, -1.25000000000000,
                10.00000000000000, 0.85355339059327);
        controlPointNet[7] = new IGAControlPoint(7, 5.48815536468909, -1.25000000000000,
                10.00000000000000, 0.85355339059327);
        controlPointNet[8] = new IGAControlPoint(8, 6.45063885631543, -1.25000000000000,
                10.00000000000000, 0.86168931333809);
        controlPointNet[9] = new IGAControlPoint(9, 7.37292198223785, -1.25000000000000,
                10.00000000000000, 0.87796115882773);
        controlPointNet[10] = new IGAControlPoint(10, 8.23223304703363, -1.25000000000000,
                10.00000000000000, 0.90236892706218);
        controlPointNet[11] = new IGAControlPoint(11, 9.01107004829591, -1.25000000000000,
                10.00000000000000, 0.93491261804146);
        controlPointNet[12] = new IGAControlPoint(12, 9.69800105423021, -1.25000000000000,
                10.00000000000000, 0.97559223176555);
        controlPointNet[13] = new IGAControlPoint(13, 10.00000000000000, -1.25000000000000,
                10.00000000000000, 1.00000000000000);
        controlPointNet[14] = new IGAControlPoint(14, 0.00000000000000, -1.50000000000000,
                9.60000000000000, 1.00000000000000);
        controlPointNet[15] = new IGAControlPoint(15, 0.30199894576979, -1.50000000000000,
                9.60000000000000, 0.97559223176555);
        controlPointNet[16] = new IGAControlPoint(16, 0.98892995170408, -1.50000000000000,
                9.60000000000000, 0.93491261804145);
        controlPointNet[17] = new IGAControlPoint(17, 1.76776695296637, -1.50000000000000,
                9.60000000000000, 0.90236892706218);
        controlPointNet[18] = new IGAControlPoint(18, 2.62707801776215, -1.50000000000000,
                9.60000000000000, 0.87796115882773);
        controlPointNet[19] = new IGAControlPoint(19, 3.54936114368457, -1.50000000000000,
                9.60000000000000, 0.86168931333809);
        controlPointNet[20] = new IGAControlPoint(20, 4.51184463531091, -1.50000000000000,
                9.60000000000000, 0.85355339059327);
        controlPointNet[21] = new IGAControlPoint(21, 5.48815536468909, -1.50000000000000,
                9.60000000000000, 0.85355339059327);
        controlPointNet[22] = new IGAControlPoint(22, 6.45063885631543, -1.50000000000000,
                9.60000000000000, 0.86168931333809);
        controlPointNet[23] = new IGAControlPoint(23, 7.37292198223785, -1.50000000000000,
                9.60000000000000, 0.87796115882773);
        controlPointNet[24] = new IGAControlPoint(24, 8.23223304703363, -1.50000000000000,
                9.60000000000000, 0.90236892706218);
        controlPointNet[25] = new IGAControlPoint(25, 9.01107004829591, -1.50000000000000,
                9.60000000000000, 0.93491261804146);
        controlPointNet[26] = new IGAControlPoint(26, 9.69800105423021, -1.50000000000000,
                9.60000000000000, 0.97559223176555);
        controlPointNet[27] = new IGAControlPoint(27, 10.00000000000000, -1.50000000000000,
                9.60000000000000, 1.00000000000000);
        controlPointNet[28] = new IGAControlPoint(28, 0.00000000000000, -1.70000000000000,
                8.80000000000000, 1.00000000000000);
        controlPointNet[29] = new IGAControlPoint(29, 0.30199894576979, -1.70000000000000,
                8.80000000000000, 0.97559223176555);
        controlPointNet[30] = new IGAControlPoint(30, 0.98892995170408, -1.70000000000000,
                8.80000000000000, 0.93491261804145);
        controlPointNet[31] = new IGAControlPoint(31, 1.76776695296637, -1.70000000000000,
                8.80000000000000, 0.90236892706218);
        controlPointNet[32] = new IGAControlPoint(32, 2.62707801776215, -1.70000000000000,
                8.80000000000000, 0.87796115882773);
        controlPointNet[33] = new IGAControlPoint(33, 3.54936114368457, -1.70000000000000,
                8.80000000000000, 0.86168931333809);
        controlPointNet[34] = new IGAControlPoint(34, 4.51184463531091, -1.70000000000000,
                8.80000000000000, 0.85355339059327);
        controlPointNet[35] = new IGAControlPoint(35, 5.48815536468909, -1.70000000000000,
                8.80000000000000, 0.85355339059327);
        controlPointNet[36] = new IGAControlPoint(36, 6.45063885631543, -1.70000000000000,
                8.80000000000000, 0.86168931333809);
        controlPointNet[37] = new IGAControlPoint(37, 7.37292198223785, -1.70000000000000,
                8.80000000000000, 0.87796115882773);
        controlPointNet[38] = new IGAControlPoint(38, 8.23223304703363, -1.70000000000000,
                8.80000000000000, 0.90236892706218);
        controlPointNet[39] = new IGAControlPoint(39, 9.01107004829591, -1.70000000000000,
                8.80000000000000, 0.93491261804146);
        controlPointNet[40] = new IGAControlPoint(40, 9.69800105423021, -1.70000000000000,
                8.80000000000000, 0.97559223176555);
        controlPointNet[41] = new IGAControlPoint(41, 10.00000000000000, -1.70000000000000,
                8.80000000000000, 1.00000000000000);
        controlPointNet[42] = new IGAControlPoint(42, 0.00000000000000, -1.28000000000000,
                7.79200000000000, 1.00000000000000);
        controlPointNet[43] = new IGAControlPoint(43, 0.30199894576979, -1.28000000000000,
                7.79200000000000, 0.97559223176555);
        controlPointNet[44] = new IGAControlPoint(44, 0.98892995170408, -1.28000000000000,
                7.79200000000000, 0.93491261804145);
        controlPointNet[45] = new IGAControlPoint(45, 1.76776695296637, -1.28000000000000,
                7.79200000000000, 0.90236892706218);
        controlPointNet[46] = new IGAControlPoint(46, 2.62707801776215, -1.28000000000000,
                7.79200000000000, 0.87796115882773);
        controlPointNet[47] = new IGAControlPoint(47, 3.54936114368457, -1.28000000000000,
                7.79200000000000, 0.86168931333809);
        controlPointNet[48] = new IGAControlPoint(48, 4.51184463531091, -1.28000000000000,
                7.79200000000000, 0.85355339059327);
        controlPointNet[49] = new IGAControlPoint(49, 5.48815536468909, -1.28000000000000,
                7.79200000000000, 0.85355339059327);
        controlPointNet[50] = new IGAControlPoint(50, 6.45063885631543, -1.28000000000000,
                7.79200000000000, 0.86168931333809);
        controlPointNet[51] = new IGAControlPoint(51, 7.37292198223785, -1.28000000000000,
                7.79200000000000, 0.87796115882773);
        controlPointNet[52] = new IGAControlPoint(52, 8.23223304703363, -1.28000000000000,
                7.79200000000000, 0.90236892706218);
        controlPointNet[53] = new IGAControlPoint(53, 9.01107004829591, -1.28000000000000,
                7.79200000000000, 0.93491261804146);
        controlPointNet[54] = new IGAControlPoint(54, 9.69800105423021, -1.28000000000000,
                7.79200000000000, 0.97559223176555);
        controlPointNet[55] = new IGAControlPoint(55, 10.00000000000000, -1.28000000000000,
                7.79200000000000, 1.00000000000000);
        controlPointNet[56] = new IGAControlPoint(56, 0.00000000000000, 0.00000000000000,
                7.29280000000000, 1.00000000000000);
        controlPointNet[57] = new IGAControlPoint(57, 0.30199894576979, 0.00000000000000,
                7.29280000000000, 0.97559223176555);
        controlPointNet[58] = new IGAControlPoint(58, 0.98892995170408, 0.00000000000000,
                7.29280000000000, 0.93491261804145);
        controlPointNet[59] = new IGAControlPoint(59, 1.76776695296637, 0.00000000000000,
                7.29280000000000, 0.90236892706218);
        controlPointNet[60] = new IGAControlPoint(60, 2.62707801776215, -0.00000000000000,
                7.29280000000000, 0.87796115882773);
        controlPointNet[61] = new IGAControlPoint(61, 3.54936114368457, -0.00000000000000,
                7.29280000000000, 0.86168931333809);
        controlPointNet[62] = new IGAControlPoint(62, 4.51184463531091, 0.00000000000000,
                7.29280000000000, 0.85355339059327);
        controlPointNet[63] = new IGAControlPoint(63, 5.48815536468909, 0.00000000000000,
                7.29280000000000, 0.85355339059327);
        controlPointNet[64] = new IGAControlPoint(64, 6.45063885631543, -0.00000000000000,
                7.29280000000000, 0.86168931333809);
        controlPointNet[65] = new IGAControlPoint(65, 7.37292198223785, -0.00000000000000,
                7.29280000000000, 0.87796115882773);
        controlPointNet[66] = new IGAControlPoint(66, 8.23223304703363, 0.00000000000000,
                7.29280000000000, 0.90236892706218);
        controlPointNet[67] = new IGAControlPoint(67, 9.01107004829591, 0.00000000000000,
                7.29280000000000, 0.93491261804146);
        controlPointNet[68] = new IGAControlPoint(68, 9.69800105423021, -0.00000000000000,
                7.29280000000000, 0.97559223176555);
        controlPointNet[69] = new IGAControlPoint(69, 10.00000000000000, 0.00000000000000,
                7.29280000000000, 1.00000000000000);
        controlPointNet[70] = new IGAControlPoint(70, 0.00000000000000, 1.28000000000000,
                7.79200000000000, 1.00000000000000);
        controlPointNet[71] = new IGAControlPoint(71, 0.30199894576979, 1.28000000000000,
                7.79200000000000, 0.97559223176555);
        controlPointNet[72] = new IGAControlPoint(72, 0.98892995170408, 1.28000000000000,
                7.79200000000000, 0.93491261804145);
        controlPointNet[73] = new IGAControlPoint(73, 1.76776695296637, 1.28000000000000,
                7.79200000000000, 0.90236892706218);
        controlPointNet[74] = new IGAControlPoint(74, 2.62707801776215, 1.28000000000000,
                7.79200000000000, 0.87796115882773);
        controlPointNet[75] = new IGAControlPoint(75, 3.54936114368457, 1.28000000000000,
                7.79200000000000, 0.86168931333809);
        controlPointNet[76] = new IGAControlPoint(76, 4.51184463531091, 1.28000000000000,
                7.79200000000000, 0.85355339059327);
        controlPointNet[77] = new IGAControlPoint(77, 5.48815536468909, 1.28000000000000,
                7.79200000000000, 0.85355339059327);
        controlPointNet[78] = new IGAControlPoint(78, 6.45063885631543, 1.28000000000000,
                7.79200000000000, 0.86168931333809);
        controlPointNet[79] = new IGAControlPoint(79, 7.37292198223785, 1.28000000000000,
                7.79200000000000, 0.87796115882773);
        controlPointNet[80] = new IGAControlPoint(80, 8.23223304703363, 1.28000000000000,
                7.79200000000000, 0.90236892706218);
        controlPointNet[81] = new IGAControlPoint(81, 9.01107004829591, 1.28000000000000,
                7.79200000000000, 0.93491261804146);
        controlPointNet[82] = new IGAControlPoint(82, 9.69800105423021, 1.28000000000000,
                7.79200000000000, 0.97559223176555);
        controlPointNet[83] = new IGAControlPoint(83, 10.00000000000000, 1.28000000000000,
                7.79200000000000, 1.00000000000000);
        controlPointNet[84] = new IGAControlPoint(84, 0.00000000000000, 1.70000000000000,
                8.80000000000000, 1.00000000000000);
        controlPointNet[85] = new IGAControlPoint(85, 0.30199894576979, 1.70000000000000,
                8.80000000000000, 0.97559223176555);
        controlPointNet[86] = new IGAControlPoint(86, 0.98892995170408, 1.70000000000000,
                8.80000000000000, 0.93491261804145);
        controlPointNet[87] = new IGAControlPoint(87, 1.76776695296637, 1.70000000000000,
                8.80000000000000, 0.90236892706218);
        controlPointNet[88] = new IGAControlPoint(88, 2.62707801776215, 1.70000000000000,
                8.80000000000000, 0.87796115882773);
        controlPointNet[89] = new IGAControlPoint(89, 3.54936114368457, 1.70000000000000,
                8.80000000000000, 0.86168931333809);
        controlPointNet[90] = new IGAControlPoint(90, 4.51184463531091, 1.70000000000000,
                8.80000000000000, 0.85355339059327);
        controlPointNet[91] = new IGAControlPoint(91, 5.48815536468909, 1.70000000000000,
                8.80000000000000, 0.85355339059327);
        controlPointNet[92] = new IGAControlPoint(92, 6.45063885631543, 1.70000000000000,
                8.80000000000000, 0.86168931333809);
        controlPointNet[93] = new IGAControlPoint(93, 7.37292198223785, 1.70000000000000,
                8.80000000000000, 0.87796115882773);
        controlPointNet[94] = new IGAControlPoint(94, 8.23223304703363, 1.70000000000000,
                8.80000000000000, 0.90236892706218);
        controlPointNet[95] = new IGAControlPoint(95, 9.01107004829591, 1.70000000000000,
                8.80000000000000, 0.93491261804146);
        controlPointNet[96] = new IGAControlPoint(96, 9.69800105423021, 1.70000000000000,
                8.80000000000000, 0.97559223176555);
        controlPointNet[97] = new IGAControlPoint(97, 10.00000000000000, 1.70000000000000,
                8.80000000000000, 1.00000000000000);
        controlPointNet[98] = new IGAControlPoint(98, 0.00000000000000, 1.50000000000000,
                9.60000000000000, 1.00000000000000);
        controlPointNet[99] = new IGAControlPoint(99, 0.30199894576979, 1.50000000000000,
                9.60000000000000, 0.97559223176555);
        controlPointNet[100] = new IGAControlPoint(100, 0.98892995170408, 1.50000000000000,
                9.60000000000000, 0.93491261804145);
        controlPointNet[101] = new IGAControlPoint(101, 1.76776695296637, 1.50000000000000,
                9.60000000000000, 0.90236892706218);
        controlPointNet[102] = new IGAControlPoint(102, 2.62707801776215, 1.50000000000000,
                9.60000000000000, 0.87796115882773);
        controlPointNet[103] = new IGAControlPoint(103, 3.54936114368457, 1.50000000000000,
                9.60000000000000, 0.86168931333809);
        controlPointNet[104] = new IGAControlPoint(104, 4.51184463531091, 1.50000000000000,
                9.60000000000000, 0.85355339059327);
        controlPointNet[105] = new IGAControlPoint(105, 5.48815536468909, 1.50000000000000,
                9.60000000000000, 0.85355339059327);
        controlPointNet[106] = new IGAControlPoint(106, 6.45063885631543, 1.50000000000000,
                9.60000000000000, 0.86168931333809);
        controlPointNet[107] = new IGAControlPoint(107, 7.37292198223785, 1.50000000000000,
                9.60000000000000, 0.87796115882773);
        controlPointNet[108] = new IGAControlPoint(108, 8.23223304703363, 1.50000000000000,
                9.60000000000000, 0.90236892706218);
        controlPointNet[109] = new IGAControlPoint(109, 9.01107004829591, 1.50000000000000,
                9.60000000000000, 0.93491261804146);
        controlPointNet[110] = new IGAControlPoint(110, 9.69800105423021, 1.50000000000000,
                9.60000000000000, 0.97559223176555);
        controlPointNet[111] = new IGAControlPoint(111, 10.00000000000000, 1.50000000000000,
                9.60000000000000, 1.00000000000000);
        controlPointNet[112] = new IGAControlPoint(112, 0.00000000000000, 1.25000000000000,
                10.00000000000000, 1.00000000000000);
        controlPointNet[113] = new IGAControlPoint(113, 0.30199894576979, 1.25000000000000,
                10.00000000000000, 0.97559223176555);
        controlPointNet[114] = new IGAControlPoint(114, 0.98892995170408, 1.25000000000000,
                10.00000000000000, 0.93491261804145);
        controlPointNet[115] = new IGAControlPoint(115, 1.76776695296637, 1.25000000000000,
                10.00000000000000, 0.90236892706218);
        controlPointNet[116] = new IGAControlPoint(116, 2.62707801776215, 1.25000000000000,
                10.00000000000000, 0.87796115882773);
        controlPointNet[117] = new IGAControlPoint(117, 3.54936114368457, 1.25000000000000,
                10.00000000000000, 0.86168931333809);
        controlPointNet[118] = new IGAControlPoint(118, 4.51184463531091, 1.25000000000000,
                10.00000000000000, 0.85355339059327);
        controlPointNet[119] = new IGAControlPoint(119, 5.48815536468909, 1.25000000000000,
                10.00000000000000, 0.85355339059327);
        controlPointNet[120] = new IGAControlPoint(120, 6.45063885631543, 1.25000000000000,
                10.00000000000000, 0.86168931333809);
        controlPointNet[121] = new IGAControlPoint(121, 7.37292198223785, 1.25000000000000,
                10.00000000000000, 0.87796115882773);
        controlPointNet[122] = new IGAControlPoint(122, 8.23223304703363, 1.25000000000000,
                10.00000000000000, 0.90236892706218);
        controlPointNet[123] = new IGAControlPoint(123, 9.01107004829591, 1.25000000000000,
                10.00000000000000, 0.93491261804146);
        controlPointNet[124] = new IGAControlPoint(124, 9.69800105423021, 1.25000000000000,
                10.00000000000000, 0.97559223176555);
        controlPointNet[125] = new IGAControlPoint(125, 10.00000000000000, 1.25000000000000,
                10.00000000000000, 1.00000000000000);

        // Test just one object of the class (that works pretty also)
        igaPatch2D = new IGAPatchSurface(id_basis, p, uNoKnots, uKnotVector, q, vNoKnots,
                vKnotVector, uNoControlPoints, vNoControlPoints, controlPointNet);

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
        CPPUNIT_ASSERT(flag == 1);

        // On the Cartesian components of the orthogonal projection
        double orthogonalProjection[] = { 4.609773274121439, -1.550611048956046, 9.096336045960706 };

        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(fabs(orthogonalProjection[i] - vertex[i]) <= relTol);

        // On the surface parametric values of the orthogonal projection
        CPPUNIT_ASSERT(fabs(u - 0.466666666666667) <= relTol);
        CPPUNIT_ASSERT(fabs(v - 0.115890578620562) <= relTol);

    }

    /***********************************************************************************************
     * \brief Test case: Test the projection of an line segement on the 2D IGA patch boundary by minimize distance
     ***********/
    void testIGAPatchSurfaceLineMinDistanceToSurfaceBoundary() {
        double P1[3] = { 1, 2, 8 };
        double P2[3] = { -1, 3, 5 };
        double t = 0.8;
        double dis = 0;
        igaPatch2D->computeLineMinimumDistanceToPatchBoundaryOnGivenEdge(t, dis, P1, P2, 2);
        CPPUNIT_ASSERT(fabs(t - 0.759200451513055) < Tol);

        double u = 0.2;
        double v = 0.8;
        igaPatch2D->computeLineMinimumDistanceToPatchBoundary(u, v, dis, P1, P2);
        CPPUNIT_ASSERT(fabs(u - 0.265254116918079) < Tol);
        CPPUNIT_ASSERT(fabs(v - 1.0) < Tol);
    }

    /***********************************************************************************************
     * \brief Test case: Test the projection of an line segement on the 2D IGA patch boundary by perpendicular mapping
     ***********/
    void testIGAPatchSurfaceLineOnSurfaceBoundary() {
        double P1[3] = { 1, 2, 8 };
        double P2[3] = { -1, 3, 5 };

        double t = 0.8;
        double div;
        double dis;
//        igaPatch2D->computePointProjectionOnPatchBoundaryOnGivenEdge(t, div,
//                dis, P1, P2, 2);
//        CPPUNIT_ASSERT(fabs(t - 0.693955750406506) < Tol);
//        CPPUNIT_ASSERT(fabs(div - 0.5) < Tol);

        double u = 0.2;
        double v = 0.8;
//        igaPatch2D->computePointProjectionOnPatchBoundary(u, v, div, dis,
//                P1, P2);
//        CPPUNIT_ASSERT(fabs(u - 0.0) < 1e-12);
//        CPPUNIT_ASSERT(fabs(v - 0.693955750406506) < Tol);
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
    CPPUNIT_TEST_SUITE (TestIGAMortarTypeProjectionTube);
    CPPUNIT_TEST (testProjectionOnIGAPatch);
    CPPUNIT_TEST (testIGAPatchSurfaceLineMinDistanceToSurfaceBoundary);
    CPPUNIT_TEST (testIGAPatchSurfaceLineOnSurfaceBoundary);

// Make the tests for leakage
    // CPPUNIT_TEST(testProjectionOnIGAPatch4Leakage);

    CPPUNIT_TEST_SUITE_END()
    ;

};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestIGAMortarTypeProjectionTube);
