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
#include "IJCSA.h"
#include "DataField.h"
#include "ConnectionIO.h"
#include "Signal.h"
#include "Residual.h"
#include "MathLibrary.h"

#include <assert.h>
#include <math.h>
#include <sstream>
#include <string.h>

using namespace std;

namespace EMPIRE {

IJCSA::IJCSA(std::string _name) :
        AbstractCouplingAlgorithm(_name) {
    debugMe = false;
    globalResidual = NULL;
    tmpVec = NULL;
}

IJCSA::~IJCSA() {
    delete[] globalResidual;
    delete[] tmpVec;
}

void IJCSA::calcNewValue() {
    /// compute the current residuals
    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++) {
        it->second->computeCurrentResidual();
    }
    /// assemble global residual vector
    int oldSize =0;
    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++) {
        Residual *residual = it->second;
    	for (int i=0; i<it->second->size; i++) {
    		globalResidual[i+oldSize]=residual->residualVector[i];
    	}
    	oldSize+=it->second->size;;
    }
    /// compute Aitken update

    /// apply the new output
    assert(outputs.size() == residuals.size());
    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++) {
        Residual *residual = it->second;
        assert(outputs.find(it->first) != outputs.end());
        CouplingAlgorithmOutput *output = outputs.find(it->first)->second;
        assert(residual->size == output->size);
        double *newOuput = new double[residual->size];
        // U_i_n+1 = U_i_n + alpha R_i_n
        for (int i=0; i<residual->size; i++) {
            newOuput[i] = output->outputCopyAtIterationBeginning[i]+ 0.1*residual->residualVector[i] ;
        }
        output->overwrite(newOuput);
        delete[] newOuput;
    }

    /// save old values
    //MathLibrary::copyDenseVector(globalResidualOld,globalResidual,globalResidualSize);

}

void IJCSA::assembleInterfaceJSystem() {
	/// tmpVec holds now globalResidualOld - globalResidual
    stringstream toOutput;
    toOutput << scientific;
    toOutput << "Aitken relaxation factor: ";
    INDENT_OUT(1, toOutput.str(), infoOut);
    cout.unsetf(ios_base::floatfield);
}


void IJCSA::init() {
    // determine global residual vector size
    globalResidualSize =0;
    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++) {
    	globalResidualSize +=it->second->size;
    }
    globalResidual    = new double [globalResidualSize];
    tmpVec            = new double [globalResidualSize];
}

} /* namespace EMPIRE */
