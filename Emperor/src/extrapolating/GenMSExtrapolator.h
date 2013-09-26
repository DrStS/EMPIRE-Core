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
/**********************************************************************************************//**
 * \file GenMSExtrapolator.h
 * \brief This file contains the class GenMSExtrapolator.
 * \date 30/11/2012
 **************************************************************************************************/
#ifndef GENMSEXTRAPOLATOR_H_
#define GENMSEXTRAPOLATOR_H_

#include "AbstractExtrapolator.h"
#include <string.h>

namespace EMPIRE {
/**********************************************************************************************//**
 * \brief Class GenMSExtrapolator is a generic linear multistep extrapolator.
 *
 * \author Michael Andre
 **************************************************************************************************/
  class GenMSExtrapolator : public AbstractExtrapolator {
  public:
	/***********************************************************************************************
	* \brief Constructor, initialize member variables
	* \param[in] _name      Extrapolator name.
	* \param[in] _numInput  Number of dots plus one.
	* \param[in] _sumOutput Sum old extrapolated data.
	* \param[in] _seqLen    Number of old time stations to sum.
	* \param[in] _deltaTime Time step size.
	* \param[in] coefDot0 Nondimensional coefficients.
	* \param[in] coefDot1 Nondimensional coefficients.
	* \param[in] coefDot2 Nondimensional coefficients.
	* \param[in] coefOut  Nondimensional coefficients.
	* \author Michael Andre
	***********/
    GenMSExtrapolator(std::string _name,
		      int _numInput,
		      bool _sumOutput,
		      int _seqLen,
		      double _deltaTime,
		      const std::vector<double>* coefDot0,
		      const std::vector<double>* coefDot1=NULL,
		      const std::vector<double>* coefDot2=NULL,
		      const std::vector<double>* coefOut=NULL
		      );

    ~GenMSExtrapolator();
	/***********************************************************************************************
	* \brief setDataIO
	* \param[in] _dataSize  Array size.
	* \param[in] _out  Extrapolated data.
	* \param[in] _inDot0    Zeroth derivative input data.
	* \param[in] _inDot1    First derivative input data.
	* \param[in] _inDot2    Second derivative input data.
	* \author Michael Andre
	***********/
    void setDataIO(int _dataSize,
		   double* _out,
		   const double* _inDot0,
		   const double* _inDot1=NULL,
		   const double* _inDot2=NULL
		   );

    void setInputAndOutput(const ConnectionIO *_input, ConnectionIO *_output);

    void addInput(const ConnectionIO *_input);

    void addOutput(ConnectionIO *_output);

    void extrapolate();

  private:
    void setDataIOAdapter();
    void initCoefficientArrays(const std::vector<double>* coef, double deltaTime, int numDot, double* coefr, int len);

    void allocMemory(double*& arr, int n, bool zero=false);
    void allocMemory(double**& arr, int nr, int nc, bool zero=false);
    void freeMemory(double** arr, int nr, int nc);
    void clearDataIO();

    bool isDataIOSet;
    const int numInput;
    const double deltaTime;

    int dataSize;
    double* out;
    const double* inDot0;
    const double* inDot1;
    const double* inDot2;

    const int seqLen;

    double* coefrDot0;
    double* coefrDot1;
    double* coefrDot2;
    double* coefrOut;

    double** inSeqDot0;
    double** inSeqDot1;
    double** inSeqDot2;

    const bool sumOutput;
    double** outSeq;

    int head;
    int* ip;
  };
  
} /* namespace EMPIRE */
#endif /* GENMSEXTRAPOLATOR_H_ */
