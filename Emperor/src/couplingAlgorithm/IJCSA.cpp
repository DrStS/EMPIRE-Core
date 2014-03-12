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
#include <string.h>


#include "IJCSA.h"
#include "DataField.h"
#include "ConnectionIO.h"
#include "Signal.h"
#include "Message.h"
#include "Residual.h"
#include "AbstractCouplingAlgorithm.h"
#include "MathLibrary.h"



using namespace std;

namespace EMPIRE {

IJCSA::IJCSA(std::string _name) :
						AbstractCouplingAlgorithm(_name) {
	debugMe = false;
	globalResidual = NULL;
	correctorVec = NULL;
}

IJCSA::~IJCSA() {
	delete[] globalResidual;
	delete[] correctorVec;
	(*interfaceJacGlobal).cleanPardiso();

}

void IJCSA::calcNewValue() {

	/// compute the current residuals
	for (map<int, Residual*>::iterator it = residuals.begin();
			it != residuals.end(); it++) {
		it->second->computeCurrentResidual();
	}
	/// assemble global residual vector
	int oldSize = 0;
	for (map<int, Residual*>::iterator it = residuals.begin();
			it != residuals.end(); it++) {
		Residual *residual = it->second;
		for (int i = 0; i < it->second->size; i++) {
			globalResidual[i + oldSize] = residual->residualVector[i];
		}
		oldSize += it->second->size;
	}
	/// get updated values for interface Jacobian
	calcInterfaceJacobian();
	assembleInterfaceJacobian();

	/// compute IJCSA update
	//(*interfaceJacGlobal).print();
	(*interfaceJacGlobal).factorize();
	(*interfaceJacGlobal).solve(correctorVec,globalResidual);
	(*interfaceJacGlobal).resetPardiso();

	/// tmpVec holds -corrector_global
	DEBUG_OUT() << std::endl;
	DEBUG_OUT() << "Corrector is:"<<std::endl;
	for (int i=0;i<globalResidualSize; i++) {
		DEBUG_OUT() << correctorVec[i] << std::endl;
	}
	DEBUG_OUT() << std::endl;
	/// apply the new output
	assert(outputs.size() == residuals.size());
	int oldResidualSize=0;
	for (map<int, Residual*>::iterator it = residuals.begin();
			it != residuals.end(); it++) {
		Residual *residual = it->second;
		assert(outputs.find(it->first) != outputs.end());
		CouplingAlgorithmOutput *output = outputs.find(it->first)->second;
		assert(residual->size == output->size);
		double *newOuput = new double[residual->size];
		// U_i_n+1 = U_i_n + alpha R_i_n

		for (int i = 0; i < residual->size; i++) {
			newOuput[i] = output->outputCopyAtIterationBeginning[i]
			                                                     - correctorVec[i+oldResidualSize];
		}
		oldResidualSize=residual->size;
		output->overwrite(newOuput);
		delete[] newOuput;
	}

}

void IJCSA::calcInterfaceJacobian() {

	/// For Automatic Differencing
	double input;
	double output;
	double denominator;
	bool localnewTimeStep;

	if (newTimeStep)
	{
		localnewTimeStep=true;
		newTimeStep=false;
	}else
	{
		localnewTimeStep=false;
	}

	for(int i=0;i<interfaceJacobianEntrys.size();i++){
		if(interfaceJacobianEntrys[i].isAutoDiff==true)
		{
			// store old function input
			assert(interfaceJacobianEntrys[i].functionInput->signal->size ==1);
			assert(interfaceJacobianEntrys[i].functionOutput->signal->size==1);

			input =interfaceJacobianEntrys[i].functionInput->signal->array[0];
			output=interfaceJacobianEntrys[i].functionOutput->signal->array[0];
			/// reset if new time step is started
			if (localnewTimeStep)
			{
				interfaceJacobianEntrys[i].value=0.0;
			}
			if (!localnewTimeStep)
			{
				denominator=input-interfaceJacobianEntrys[i].oldInput;

				if(::abs(denominator>1e-10)){
					interfaceJacobianEntrys[i].value=((output-interfaceJacobianEntrys[i].oldOutput)/denominator)*interfaceJacobianEntrys[i].coefficient;
					debugOut() << "Estimated interfaceJacobian[" << interfaceJacobianEntrys[i].indexRow<< " , "<< interfaceJacobianEntrys[i].indexColumn<< "] ="<< interfaceJacobianEntrys[i].value<< endl;
				}
				else{
					warningOut() << "Avoiding instability in automatic differencing  "<<endl;
				}
			}

			interfaceJacobianEntrys[i].oldInput=input;
			interfaceJacobianEntrys[i].oldOutput=output;

		}
	}

	// write autoDiff to file
	autoDiffFile << currentTimeStep << '\t' << currentIteration;
	autoDiffFile << scientific;
	autoDiffFile.precision(10);
	for(int i=0;i<interfaceJacobianEntrys.size();i++){
		if(interfaceJacobianEntrys[i].isAutoDiff==true)
		{
			autoDiffFile << '\t' << "[" << interfaceJacobianEntrys[i].indexRow<< " , "<< interfaceJacobianEntrys[i].indexColumn<< "]=\t"<< interfaceJacobianEntrys[i].value << "\t \t";
		}
	}
	autoDiffFile << endl;

}

void IJCSA::assembleInterfaceJacobian(){
	for(int i=0;i<interfaceJacobianEntrys.size();i++){
		(*interfaceJacGlobal)(interfaceJacobianEntrys[i].indexRow-1,
				interfaceJacobianEntrys[i].indexColumn-1)=interfaceJacobianEntrys[i].value;
	}
}


void IJCSA::init() {
	// determine global residual vector size
	globalResidualSize = 0;
	for (map<int, Residual*>::iterator it = residuals.begin();
			it != residuals.end(); it++) {
		globalResidualSize += it->second->size;
	}
	globalResidual = new double[globalResidualSize];
	correctorVec   = new double[globalResidualSize];
	interfaceJacGlobal = new MathLibrary::SparseMatrix<double>(
			globalResidualSize, false);

	assembleInterfaceJacobian();
	std::string autoDiffFileName = "automaticDifferentiation";
	autoDiffFileName.append(".log");
	autoDiffFile.open(autoDiffFileName.c_str(), ios_base::out);
	autoDiffFile << scientific;
	autoDiffFile << "timeStep" << '\t' << "iteration" << '\t' << "interfaceJacobian" << '\t'<< endl;
	assert(!autoDiffFile.fail());
}
void IJCSA::addInterfaceJacobianEntry(unsigned int _indexRow,
		unsigned int _indexColumn, double _value) {
	interfaceJacobianEntry tmp ={_indexRow,_indexColumn,true,_value,false,NULL,NULL,false,NULL,0,0,0};
	interfaceJacobianEntrys.push_back(tmp);

}

void IJCSA::addInterfaceJacobianEntry(unsigned int _indexRow,
		unsigned int _indexColumn, ConnectionIO* _functionInput, ConnectionIO* _functionOutput, double _coefficient) {

	assert( _functionInput->type==EMPIRE_ConnectionIO_Signal);
	assert(_functionOutput->type==EMPIRE_ConnectionIO_Signal);

	interfaceJacobianEntry tmp ={_indexRow,_indexColumn,false,0.0,true,_functionInput,_functionOutput,false,NULL,_coefficient,0,0};
	interfaceJacobianEntrys.push_back(tmp);

}

} /* namespace EMPIRE */
