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
/***********************************************************************************************//**
 * ComaAitken is a modified version of the old Aitken algorithm in CoMA, which is used to compare
 * with the Aitken algorithm in EMPIRE.
 **************************************************************************************************/
#ifndef COMAAITKEN_H_
#define COMAAITKEN_H_

#include "DataField.h"
#include <iostream>
#include <stdlib.h>

namespace EMPIRE {
using namespace std;
class ComaAitken {
private:

    int dimension;

    double startRelaxationFactor;
    double relaxationFactor;
    double aitkenFactor;
    double aitkenFactorLastSubiteration;

    int numAllDofs;

    double* disp_Current;
    double* disp_Innerloop_Old;
    double* delta_d_i;
    double* delta_d_iplus1;

    DataField *dfCurrent;
    DataField *dfOld;

    friend class TestAitken;

public:
    ComaAitken(DataField *_dfCurrent, DataField *_dfOld, double _initialAitkenFactor) :
            dfCurrent(_dfCurrent), dfOld(_dfOld), startRelaxationFactor(_initialAitkenFactor) {

        numAllDofs = dfCurrent->numLocations * dfCurrent->dimension;

        relaxationFactor = startRelaxationFactor;
        aitkenFactor = 0.0;
        aitkenFactorLastSubiteration = 0.0;

        disp_Current = new double[numAllDofs];
        for (int i = 0; i < numAllDofs; i++)
            disp_Current[i] = 0.0;
        disp_Innerloop_Old = new double[numAllDofs];
        for (int i = 0; i < numAllDofs; i++)
            disp_Innerloop_Old[i] = 0.0;
        delta_d_i = new double[numAllDofs];
        for (int i = 0; i < numAllDofs; i++)
            delta_d_i[i] = 0.0;
        delta_d_iplus1 = new double[numAllDofs];
        for (int i = 0; i < numAllDofs; i++)
            delta_d_iplus1[i] = 0.0;
    }

    virtual ~ComaAitken() {
        delete[] disp_Current;
        delete[] disp_Innerloop_Old;
        delete[] delta_d_i;
        delete[] delta_d_iplus1;
    }

    void calcNewValues(int innerLoopStepNumber) {
        //at Innerloop_Current: values just received
        //at Innerloop_old: values from last subiteration
        //relaxation: values to be send are written to innerloop_current
        //           values to be send = relaxation factor *  values just received + (1-relaxation factor) * values from last subiteration

        if (innerLoopStepNumber == 1) {
            cerr
                    << "ERROR in CoMA:\n\t constRelaxation::calcNewValues - Called in first subiteration with innerLoopStepNumber=1! \n"
                    << "Doing Relaxation here doesn't make much sense..." << endl;
            exit(1);
        }

        //write old and current values from mesh adapter to own variable management
        //loop over tagged Meshes
        for (int i = 0; i < numAllDofs; i++) {
            disp_Current[i] = dfCurrent->data[i];
            disp_Innerloop_Old[i] = dfOld->data[i];
        }

        //in first underrelaxtion use a constant value
        if (innerLoopStepNumber == 2) {
            for (int i = 0; i < numAllDofs; i++) {
                //setup delta d_i+1
                delta_d_iplus1[i] = disp_Innerloop_Old[i] - disp_Current[i]; // written in CoMA
                //setup delta d_i - delta d_i+1
                delta_d_i[i] = delta_d_i[i] - delta_d_iplus1[i]; // written in CoMA

                //loop over all nodes and apply constant relaxation
                disp_Current[i] = (1.0 - startRelaxationFactor) * disp_Innerloop_Old[i]
                        + startRelaxationFactor * disp_Current[i];

                for (int i = 0; i < numAllDofs; i++) { // added by T. Wang
                    //copy delta d_i+1 to delta d_i
                    delta_d_i[i] = delta_d_iplus1[i];
                }
            }
            relaxationFactor = startRelaxationFactor;
            aitkenFactorLastSubiteration = 1.0 - startRelaxationFactor;
        }

        //for higher underrelaxtion iterations calculate Aitken factor
        else if (innerLoopStepNumber > 2) {
            //aitken relaxation
            //based on PHD-Thesis Mok, page 131

            double numerator = 0.0; //Z�hler
            double denominator = -1.0; //Nenner
            double fraction = 0.0;

            //1. calculate Aitken Factor

            //loop over all nodes and calculate new values
            for (int i = 0; i < numAllDofs; i++) {
                //setup delta d_i+1
                delta_d_iplus1[i] = disp_Innerloop_Old[i] - disp_Current[i];

                //setup delta d_i - delta d_i+1
                delta_d_i[i] = delta_d_i[i] - delta_d_iplus1[i];
            }
            numerator = this->calcVectorProduct(delta_d_i, delta_d_iplus1, numAllDofs);
            denominator = this->calcVectorProduct(delta_d_i, delta_d_i, numAllDofs);

            for (int i = 0; i < numAllDofs; i++) {
                //copy delta d_i+1 to delta d_i
                delta_d_i[i] = delta_d_iplus1[i];
            }
            /*if (fabs(denominator) <= 1.0e-14) {
                cerr << "Warning in CoMA:\n\t Denominator in Aitken factor is nearly zero!" << endl;
            }*/
            if (denominator != 0.0) {
                fraction = numerator / denominator;

                aitkenFactor = aitkenFactorLastSubiteration
                        + (aitkenFactorLastSubiteration - 1.0) * fraction;
            } else {
                aitkenFactor = 0.0;
            }

            //2. calculate relaxation Factor
            relaxationFactor = 1.0 - aitkenFactor;

            aitkenFactorLastSubiteration = aitkenFactor;

            //loop over all nodes and calculate new values
            for (int i = 0; i < numAllDofs; i++) {
                disp_Current[i] = (1.0 - relaxationFactor) * disp_Innerloop_Old[i]
                        + relaxationFactor * disp_Current[i];
            }
        }

        //write from own variable management to values from mesh adapter
        //loop over tagged Meshes
        for (int i = 0; i < numAllDofs; i++) {
            dfCurrent->data[i] = disp_Current[i];
        }
    }

    double calcVectorProduct(double* vec1, double* vec2, int size) {

        double product = 0.0;
        for (int i = 0; i < size; i++) {
            product += vec1[i] * vec2[i];
        }
        return product;
    }
};

} /* namespace EMPIRE */
#endif /* COMAAITKEN_H_ */
