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
#include <assert.h>
#include <sstream>
#include <math.h>

#include "Aitken.h"
#include "DataField.h"
#include "Message.h"
#include "Signal.h"
#include "ConnectionIO.h"
#include "EMPEROR_Enum.h"

using namespace std;

namespace EMPIRE {

const double Aitken::LIMIT = 1E6;

Aitken::Aitken(std::string _name, double _initialAitkenFactor) :
        AbstractCouplingAlgorithm(_name), INIT_AITKEN_FACTOR(_initialAitkenFactor) {
    assert(INIT_AITKEN_FACTOR >= 0.0);
    X_out = NULL; // segmentation fault happens if it is not null during destructing (unit test)
    R_0 = NULL;
    debugMe = true;
}

Aitken::~Aitken() {
    delete X_out;
    delete R_0;
}

void Aitken::setInputAndOutput(const ConnectionIO *_input, ConnectionIO *_output) {
    assert(input==NULL);
    assert(output==NULL);
    assert(_input!=NULL);
    assert(_output!=NULL);
    input = _input;
    output = _output;
    assert(input->type == output->type);
    EMPIRE_ConnectionIO_Type IOType = input->type;
    assert(input->type == output->type);
    if (IOType == EMPIRE_ConnectionIO_DataField) {
        assert(input->dataField->dimension == output->dataField->dimension);
        assert(input->dataField->numLocations == output->dataField->numLocations);
        int numLocations = input->dataField->numLocations;
        int dimension = input->dataField->dimension;
        SIZE = numLocations * dimension;
        X_out = new double[SIZE];
        R_0 = new double[SIZE];
    } else if (IOType == EMPIRE_ConnectionIO_Signal) {
        assert(input->signal->size == output->signal->size);
        assert(input->signal->dimension == output->signal->dimension);
        SIZE = input->signal->size;
        X_out = new double[SIZE];
        R_0 = new double[SIZE];
    } else {
        assert(false);
    }
}

void Aitken::calcNewValue() {
///TODO change syntax
    /* --------------------------------------------------------------
     * FORMULAR (J. Degroote's handout P.101):
     *
     * [x_k+1] = (1-[w_k]) * [x_k] + [w_k] * [x_k_~]
     *
     * with [w_k] = - [w_k-1] * ( [r_k-1]^T * ([r_k]-[r_k-1])  /   ([r_k]-[r_k-1])^T * ([r_k]-[r_k-1]) )
     *
     * --------------------------------------------------------------
     * DEFINE:
     * X_out --- the output memory which stores [x_k], and is updated to [x_k+1] at the end
     * X_in  --- the input memory which stores ([x_k_~])
     * R_1    --- residual at this step ([r_k], defined by [r_k] = [x_k_~] - [x_k])
     * R_0    --- residual of last step ([r_k-1])
     * w_1    --- aikten factor of this step ([w_k])
     * w_0    --- aikten factor of last step ([w_k-1])
     * --------------------------------------------------------------
     */

    if (newLoop) {
        newLoop = false;
        assert(step == 1);
    }

    double *dataIn;
    double *dataOut;
    EMPIRE_ConnectionIO_Type IOType = input->type;
    if (IOType == EMPIRE_ConnectionIO_DataField) {
        dataIn = input->dataField->data; // rename
        dataOut = output->dataField->data; // rename
    } else if (IOType == EMPIRE_ConnectionIO_Signal) {
        dataIn = input->signal->array; // rename
        dataOut = output->signal->array; // rename
    } else {
        assert(false);
    }

    // X_out, R_0 and w_0 are from the last step, X_in is the input at this step
    if (step == 1) {
        const double *X_in = dataIn; // rename
        // Set X_out for next step
        for (int i = 0; i < SIZE; i++)
            X_out[i] = X_in[i];
    } else {
        // 1. initialize
        const double *X_in = dataIn; // rename
        double *tmpMemory = new double[SIZE];

        // 2. calculate w_1
        double w_1;
        if (step == 2) {
            w_1 = INIT_AITKEN_FACTOR;
        } else {
            vecMinus(X_in, X_out, tmpMemory, SIZE); // now tmpMemory stores R_1
            vecMinusEqual(tmpMemory, R_0, SIZE); // now tmpMemory stores deltaR
            double numerator = vecDotProduct(R_0, tmpMemory, SIZE);
            double denominator = vecDotProduct(tmpMemory, tmpMemory, SIZE);
            if (denominator != 0.0) {
                assert(fabs(numerator) < LIMIT * denominator);
                w_1 = -w_0 * numerator / denominator;
            } else { // this happens due to convergence is already obtained, so [r_k]=[r_k-1]=0
                // convergence is obtained, set w_1 to 1.0
                w_1 = 1.0;
            }
        }

        if (debugMe) { // write aitken factor to shell
            stringstream ss;
            ss << "Aitken relaxation factor: " << w_1;
            INDENT_OUT(1, ss.str(), infoOut);
        }

        // 3. update R_0
        vecMinus(X_in, X_out, R_0, SIZE);

        { // compute residuals
            if (step == 2) {
                initialResidual = vecL2Norm(R_0, SIZE);
                if (debugMe) { // write aitken factor to shell
                    /*stringstream ss;
                     ss << "Aitken initial residual: " << initialResidual;
                     INDENT_OUT(1, ss.str(), infoOut);*/
                }
            }
            currentResidual = vecL2Norm(R_0, SIZE);
            if (debugMe) { // write aitken factor to shell
                /*stringstream ss;
                 ss << "Aitken relative residual: " << currentResidual / initialResidual;
                 INDENT_OUT(1, ss.str(), infoOut);*/
            }
        }

        // 4. calculate new X_out
        double tmp = 1.0 - w_1;

        assert(tmp != 1.0);
        // when tmp == 1.0, w_1<1e-16, then x_k+1=x_k, does not make sense

        vecCopy(X_in, tmpMemory, SIZE);
        vecScalarMultiply(tmpMemory, w_1, SIZE);
        vecScalarMultiply(X_out, tmp, SIZE);
        vecPlusEqual(X_out, tmpMemory, SIZE);

        // 5. set w_0 for next step
        w_0 = w_1;
        delete tmpMemory;
    }

    // updata output
    for (int i = 0; i < SIZE; i++)
        dataOut[i] = X_out[i];
    step++;
}

void Aitken::setZeroState() {

}

} /* namespace EMPIRE */
