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
#include "ConstantRelaxation.h"
#include "DataField.h"
#include "ConnectionIO.h"
#include "Signal.h"
#include <assert.h>
#include <math.h>
#include <sstream>

using namespace std;

namespace EMPIRE {

ConstantRelaxation::ConstantRelaxation(std::string _name, double _relaxationFactor) :
        AbstractCouplingAlgorithm(_name), RELAXATION_FACTOR(_relaxationFactor) {
    assert(RELAXATION_FACTOR>0.0);
    debugMe = false;
}
ConstantRelaxation::~ConstantRelaxation() {

}
void ConstantRelaxation::setInputAndOutput(const ConnectionIO *_input, ConnectionIO *_output) {
    assert(input==NULL);
    assert(output==NULL);
    assert(_input!=NULL);
    assert(_output!=NULL);
    input = _input;
    output = _output;
    assert(input->type == output->type);
    EMPIRE_ConnectionIO_Type IOType = input->type;
    if (IOType == EMPIRE_ConnectionIO_DataField) {
        assert(input->dataField->dimension == output->dataField->dimension);
        assert(input->dataField->numLocations == output->dataField->numLocations);
        int numLocations = input->dataField->numLocations;
        int dimension = input->dataField->dimension;
        SIZE = numLocations * dimension;
        X_out = new double[SIZE];
    } else if (IOType == EMPIRE_ConnectionIO_Signal) {
        assert(input->signal->size == output->signal->size);
        assert(input->signal->dimension == output->signal->dimension);
        SIZE = input->signal->size;
        X_out = new double[SIZE];
    } else {
        assert(false);
    }
}
void ConstantRelaxation::calcNewValue() {
    /* --------------------------------------------------------------
     * [x_k+1] = (1-[w]) * [x_k] + [w] * [x_k_~]
     *         = [x_k] + [w] * ([x_k_~] - [x_k])
     * --------------------------------------------------------------
     * DEFINE:
     * X_out --- the output memory which stores [x_k], and is updated to [x_k+1] at the end
     * X_in  --- the input memory which stores ([x_k_~])
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

    // X_out are from the last step, X_in is the input at this step
    if (step == 1) {
        const double *X_in = dataIn; // rename
        // Set X_out for next step
        for (int i = 0; i < SIZE; i++)
            X_out[i] = X_in[i];
    } else {
        const double *X_in = dataIn; // rename
        double *residualVec = new double[SIZE];
        vecCopy(X_in, residualVec, SIZE);
        vecMinusEqual(residualVec, X_out, SIZE);
        vecScalarMultiply(residualVec, RELAXATION_FACTOR, SIZE);
        vecPlusEqual(X_out, residualVec, SIZE);

        { // compute residuals
            if (step == 2) {
                initialResidual = vecDotProduct(residualVec, residualVec, SIZE);
                initialResidual = sqrt(initialResidual);
                initialResidual /= (double) SIZE;
                if (debugMe) {
                    stringstream ss;
                    ss << "ConstantRelaxation initial residual: " << initialResidual;
                    Message::writeTextWithIndent(1, ss.str(), infoOut);
                }
            }
            currentResidual = vecDotProduct(residualVec, residualVec, SIZE);
            currentResidual = sqrt(currentResidual);
            currentResidual /= (double) SIZE;
            if (debugMe) {
                stringstream ss;
                ss << "ConstantRelaxation relative residual: " << currentResidual / initialResidual;
                Message::writeTextWithIndent(1, ss.str(), infoOut);
            }
        }

        delete residualVec;
    }

    // updata output
    for (int i = 0; i < SIZE; i++)
        dataOut[i] = X_out[i];
    step++;
}
} /* namespace EMPIRE */
