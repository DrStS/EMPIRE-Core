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
#include <assert.h>
#include <math.h>

#include "AbstractCouplingAlgorithm.h"
#include "DataField.h"
#include "Signal.h"
#include "ConnectionIO.h"
#include <iostream>

using namespace std;

namespace EMPIRE {

AbstractCouplingAlgorithm::AbstractCouplingAlgorithm(std::string _name) :
        name(_name), newLoop(false), step(0), input(NULL), output(NULL) {
}

AbstractCouplingAlgorithm::~AbstractCouplingAlgorithm() {
}

void AbstractCouplingAlgorithm::setInputAndOutput(const ConnectionIO *_input,
        ConnectionIO *_output) {
    assert(input==NULL);
    assert(output==NULL);
    assert(_input!=NULL);
    assert(_output!=NULL);
    input = _input;
    output = _output;
    assert(input->type == output->type);
}

void AbstractCouplingAlgorithm::setNewLoop() {
    newLoop = true;
    step = 1;
}

double AbstractCouplingAlgorithm::getCurrentResidual() {
    return currentResidual;
}

double AbstractCouplingAlgorithm::getInitialResidual() {
    return initialResidual;
}

void AbstractCouplingAlgorithm::vecCopy(const double *from, double *to, int size) {
    for (int i = 0; i < size; i++)
        to[i] = from[i];
}

double AbstractCouplingAlgorithm::vecDotProduct(const double *vec1, const double *vec2,
        int size) {
    double product = 0.0;
    for (int i = 0; i < size; i++)
        product += vec1[i] * vec2[i];
    return product;
}

void AbstractCouplingAlgorithm::vecScalarMultiply(double *vec, const double SCALAR,
        int size) {
    for (int i = 0; i < size; i++)
        vec[i] *= SCALAR;
}

void AbstractCouplingAlgorithm::vecPlusEqual(double *vec1, const double *vec2, int size) {
    for (int i = 0; i < size; i++)
        vec1[i] += vec2[i];
}

void AbstractCouplingAlgorithm::vecMinusEqual(double *vec1, const double *vec2, int size) {
    for (int i = 0; i < size; i++)
        vec1[i] -= vec2[i];
}

void AbstractCouplingAlgorithm::vecMinus(const double *vec1, const double *vec2, double *vec1minus2,
        int size) {
    for (int i = 0; i < size; i++)
        vec1minus2[i] = vec1[i] - vec2[i];
}

double AbstractCouplingAlgorithm::vecL2Norm(const double *vec, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += vec[i] * vec[i];
    }
    sum /= size;
    sum = sqrt(sum);
    return sum;
}

} /* namespace EMPIRE */
