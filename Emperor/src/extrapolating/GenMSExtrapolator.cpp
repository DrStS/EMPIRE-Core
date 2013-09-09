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
#include "GenMSExtrapolator.h"
#include "ConnectionIO.h"
#include "DataField.h"
#include "Signal.h"
#include <assert.h>
#include <cmath>

using namespace std;

namespace EMPIRE {
  GenMSExtrapolator::GenMSExtrapolator(std::string _name, 
				       int _numInput,
				       bool _sumOutput,
				       int _seqLen,
				       double _deltaTime,
				       const vector<double>* coefDot0,
				       const vector<double>* coefDot1,
				       const vector<double>* coefDot2,
				       const vector<double>* coefOut
				       )
    : AbstractExtrapolator(_name), isDataIOSet(false), numInput(_numInput), deltaTime(_deltaTime), 
      dataSize(0), out(NULL), inDot0(NULL), inDot1(NULL), inDot2(NULL), seqLen(_seqLen), coefrDot0(NULL), 
      coefrDot1(NULL), coefrDot2(NULL), coefrOut(NULL), inSeqDot0(NULL), inSeqDot1(NULL), inSeqDot2(NULL), 
      sumOutput(_sumOutput), outSeq(NULL), head(0), ip(NULL)
  {
    assert(seqLen > 0);

    // set dimensionally consistent coefficient arrays
    switch (numInput) {
    case 3:
      assert(coefDot2 != NULL);
      assert(coefDot2->size() == seqLen);
      coefrDot2 = new double[seqLen];
      initCoefficientArrays(coefDot2, deltaTime, 2, coefrDot2, seqLen);
    case 2:
      assert(coefDot1 != NULL);
      assert(coefDot1->size() == seqLen);
      coefrDot1 = new double[seqLen];
      initCoefficientArrays(coefDot1, deltaTime, 1, coefrDot1, seqLen);
    case 1:
      assert(coefDot0 != NULL);
      assert(coefDot0->size() == seqLen);
      coefrDot0 = new double[seqLen];
      initCoefficientArrays(coefDot0, deltaTime, 0, coefrDot0, seqLen);
      break;
    default:
      assert(false);
    }

    if (sumOutput) {
      assert(coefOut != NULL);
      assert(coefOut->size() == seqLen);
      coefrOut = new double[seqLen];
      initCoefficientArrays(coefOut, 1., 1, coefrOut, seqLen);
    }

    // cyclic increment
    ip = new int[seqLen];
    for (int i=0; i<seqLen; i++) {
      ip[i] = i+1;
    }
    ip[seqLen-1] = 0;

  }

  GenMSExtrapolator::~GenMSExtrapolator()
  {
    delete coefrDot0;
    delete coefrDot1;
    delete coefrDot2;
    delete coefrOut;
    delete ip;
    
    clearDataIO();
  }

  void GenMSExtrapolator::setDataIO(int _dataSize, double* _out, const double* _inDot0, 
				    const double* _inDot1, const double* _inDot2) {
    assert(_dataSize > 0);
    clearDataIO();

    dataSize = _dataSize;
    assert(_out != NULL);
    out = _out;

    switch (numInput) {
    case 3:
      assert(_inDot2 != NULL);
      inDot2 = _inDot2;
      allocMemory(inSeqDot2, seqLen, dataSize, true);
    case 2:
      assert(_inDot1 != NULL);
      inDot1 = _inDot1;
      allocMemory(inSeqDot1, seqLen, dataSize, true);
    case 1:
      assert(_inDot0 != NULL);
      inDot0 = _inDot0;
      allocMemory(inSeqDot0, seqLen, dataSize, true);
    }

    if (sumOutput) {
      allocMemory(outSeq, seqLen, dataSize, true);
    }
    isDataIOSet = true;
  }

  void GenMSExtrapolator::setInputAndOutput(const ConnectionIO *_input, ConnectionIO *_output) {
    assert(numInput == 1);
    AbstractExtrapolator::setInputAndOutput(_input, _output);
    
    setDataIOAdapter();
  }
  
  void GenMSExtrapolator::addInput(const ConnectionIO *_input) {
    AbstractExtrapolator::addInput(_input);

    if (inputVec.size() == numInput && outputVec.size() == 1) {
      assert(inputVec[0]->type == outputVec[0]->type);
      setDataIOAdapter();
    }
    else if (inputVec.size() > numInput) { // no bad stuff
      assert(false);
    }
  }

  void GenMSExtrapolator::addOutput(ConnectionIO *_output) {
    AbstractExtrapolator::addOutput(_output);

    if (inputVec.size() == numInput && outputVec.size() == 1) {
      assert(inputVec[0]->type == outputVec[0]->type);
      setDataIOAdapter();      
    }
    else if (outputVec.size() > 1) { // no bad stuff
      assert(false);
    }

  }

  void GenMSExtrapolator::extrapolate() {
    int i,j,index;
    double *tmp;

    assert(isDataIOSet);
    if (!doExtrapolate) {
      // do nothing
      return;
    }

    tmp = new double[dataSize];

    for (i=0; i<dataSize; i++) {
      tmp[i] = 0.;
    }

    // load new data for zeroth derivative input
    for (i=0; i<dataSize; i++) {
      inSeqDot0[head][i] = inDot0[i];
    }

    // sum zeroth derivative history
    index = ip[head];
    for (i=0; i<seqLen; i++) {
      for (j=0; j<dataSize; j++) {
	tmp[j] = tmp[j] + coefrDot0[i]*inSeqDot0[index][j];
      }
      index = ip[index];
    }
    
    if (numInput > 1) {
      // load new data for first derivative input
      for (i=0; i<dataSize; i++) {
	inSeqDot1[head][i] = inDot1[i];
      }

      // sum first derivative history
      index = ip[head];
      for (i=0; i<seqLen; i++) {
	for (j=0; j<dataSize; j++) {
	  tmp[j] = tmp[j] + coefrDot1[i]*inSeqDot1[index][j];
	}
	index = ip[index];
      }
      
      if (numInput > 2) {
	// load new data for second derivative input
	for (i=0; i<dataSize; i++) {
	  inSeqDot2[head][i] = inDot2[i];
	}

	// sum second derivative history
	index = ip[head];
	for (i=0; i<seqLen; i++) {
	  for (j=0; j<dataSize; j++) {
	    tmp[j] = tmp[j] + coefrDot2[i]*inSeqDot2[index][j];
	  }
	  index = ip[index];
	}
      }
    }

    if (sumOutput) {
      // sum output history
      index = head;
      for (i=0; i<seqLen; i++) {
	for (j=0; j<dataSize; j++) {
	  tmp[j] = tmp[j] + coefrOut[i]*outSeq[index][j];
	}
	index = ip[index];
      }

      // store output
      for (i=0; i<dataSize; i++) {
	outSeq[head][i] = tmp[i];
      }
    }

    // set output
    for (i=0; i<dataSize; i++) {
      out[i] = tmp[i];
    }

    head = ip[head];
    delete tmp;
    doExtrapolate = false;
  }

  /**************************************************************************************************
   * Private
   **************************************************************************************************/

  void GenMSExtrapolator::setDataIOAdapter() {
    if (inputVec[0]->type == EMPIRE_ConnectionIO_DataField) { // data field
      // check for consistent data size
      dataSize = inputVec[0]->dataField->numLocations * inputVec[0]->dataField->dimension;
      assert(outputVec[0]->dataField->numLocations * outputVec[0]->dataField->dimension == dataSize);
      for (int i=1; i<inputVec.size(); i++) {
	assert(inputVec[i]->dataField->numLocations * inputVec[i]->dataField->dimension == dataSize);
      }
      // set IO arrays
      switch (numInput) {
      case 3:
	setDataIO(dataSize, outputVec[0]->dataField->data, inputVec[0]->dataField->data,
		  inputVec[1]->dataField->data, inputVec[2]->dataField->data);
	break;
      case 2:
	setDataIO(dataSize, outputVec[0]->dataField->data, inputVec[0]->dataField->data,
		  inputVec[1]->dataField->data);
	break;
      case 1:
	setDataIO(dataSize, outputVec[0]->dataField->data, inputVec[0]->dataField->data);
	break;
      }
    }
    else if (inputVec[0]->type == EMPIRE_ConnectionIO_Signal) { // signal
      // check for consistent data size
      dataSize = inputVec[0]->signal->size;
      assert(outputVec[0]->signal->size == dataSize);
      for (int i=1; i<inputVec.size(); i++) {
	assert(inputVec[i]->signal->size == dataSize);
      }
      // set IO arrays
      switch (numInput) {
      case 3:
	setDataIO(dataSize, outputVec[0]->signal->array, inputVec[0]->signal->array,
		  inputVec[1]->signal->array, inputVec[2]->signal->array);
	break;
      case 2:
	setDataIO(dataSize, outputVec[0]->signal->array, inputVec[0]->signal->array,
		  inputVec[1]->signal->array);
	break;
      case 1:
	setDataIO(dataSize, outputVec[0]->signal->array, inputVec[0]->signal->array);
	break;
      }	
    }
    else {
      assert(false);
    }
  }

  void GenMSExtrapolator::initCoefficientArrays(const vector<double>* coef, double deltaTime, 
						int numDot, double* coefr, int len) {

    const double pre = pow(deltaTime, numDot);

    int i;
    int i_r = len-1;

    for (i=0; i<len; i++) {
      coefr[i] = pre*coef->at(i_r);
      i_r--;
    }
  }

  void GenMSExtrapolator::allocMemory(double*& arr, int n, bool zero) {
    arr = new double[n];

    if (zero) {
      for (int i=0; i<n; i++) {
	arr[i] = 0.;
      }
    }
  }

  void GenMSExtrapolator::allocMemory(double**& arr, int nr, int nc, bool zero) {
    int i, j;

    arr = new double*[nr];
    for (i=0; i<nr; i++) {
      arr[i] = new double[nc];
    }

    if (zero) {
      for (i=0; i<nr; i++) {
	for (j=0; j<nc; j++) {
	  arr[i][j] = 0.;
	}
      }
    }
  }

  void GenMSExtrapolator::freeMemory(double** arr, int nr, int nc) {
    for (int i=0; i<nr; i++) {
      delete arr[i];
    }
    delete arr;
  }

  void GenMSExtrapolator::clearDataIO() {
    if (isDataIOSet) {
      switch (numInput) {
      case 3:
	freeMemory(inSeqDot2, seqLen, dataSize);
	inDot2 = NULL;
      case 2:
	freeMemory(inSeqDot1, seqLen, dataSize);
	inDot1 = NULL;
      case 1:
	freeMemory(inSeqDot0, seqLen, dataSize);
	inDot0 = NULL;
      }

      if (sumOutput) {
	freeMemory(outSeq, seqLen, dataSize);
      }

      out = NULL;
      dataSize = 0;
      head = 0;
      isDataIOSet = false;
    }
  }

} /* namespace EMPIRE */
