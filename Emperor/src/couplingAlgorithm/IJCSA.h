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
 * \file Aitken.h
 * This file holds the class IJCSA
 * \date 6/2/2014
 **************************************************************************************************/

#ifndef IJCSA_H_
#define IJCSA_H_

#include <vector>
#include "AbstractCouplingAlgorithm.h"

namespace EMPIRE {

namespace MathLibrary{
template<typename T>
class SparseMatrix;
}

class ConnectionIO;
/********//**
 * \brief Class IJCSA does Interface-Jacobian based Co-Simulation
 ***********/
class IJCSA: public AbstractCouplingAlgorithm {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name the name of the coupling algorithm
     * \author Stefan Sicklinger
     ***********/
	IJCSA(std::string _name);
    /***********************************************************************************************
     * \brief Destructor
     * \author Stefan Sicklinger
     ***********/
    virtual ~IJCSA();
    /***********************************************************************************************
     * \brief Calculate the new value of the output (improved value for next iteration/time step)
     * \author Stefan Sicklinger
     ***********/
    void calcNewValue();
    /***********************************************************************************************
     * \brief Init IJCSA
     * \author Stefan Sicklinger
     ***********/
    void init();
    /***********************************************************************************************
     * \brief add value to interface Jacobian matrix
     * \param[in] _indexRow    row index of value
     * \param[in] _indexColumn column index of value
     * \param[in] _value
     * \author Stefan Sicklinger
     ***********/
    void addInterfaceJacobianEntry(unsigned int _indexRow, unsigned int _indexColumn, double _value);
    /***********************************************************************************************
     * \brief add value to interface Jacobian matrix
     * \param[in] _indexRow    row index of value
     * \param[in] _indexColumn column index of value
     * \param[in] _functionInput pointer to ConnectionIO  (only signal)
     * \param[in] _functionOutput pointer to ConnectionIO (only signal)
     * \param[in] _coefficient prefactor of FD computed interface jacobian
     * \author Stefan Sicklinger
     ***********/
    void addInterfaceJacobianEntry(unsigned int _indexRow, unsigned int _indexColumn, ConnectionIO* _functionInput, ConnectionIO* _functionOutput, double _coefficient);
    /***********************************************************************************************
     * \brief add value to interface Jacobian matrix
     * \param[in] _indexRow    row index of value
     * \param[in] _indexColumn column index of value
     * \param[in] _jacobianSignal pointer to ConnectionIO  (only signal)
     * \author Stefan Sicklinger
     ***********/
    void addInterfaceJacobianEntry(unsigned int _indexRow, unsigned int _indexColumn, ConnectionIO* _jacobianSignal);
private:
    /***********************************************************************************************
     * \brief Calculates interface Jacobian using FD
     * \author Stefan Sicklinger
     ***********/
    void calcInterfaceJacobian();
    /***********************************************************************************************
     * \brief Recomputes interface Jacobian
     * \author Stefan Sicklinger
     ***********/
    void assembleInterfaceJacobian();

    struct interfaceJacobianEntry{
    	unsigned int indexRow;
    	unsigned int indexColumn;
        bool isConstant;
    	double value;
        bool isAutoDiff;
        ConnectionIO * functionInput;
        ConnectionIO * functionOutput;
    	bool isSignal;
    	ConnectionIO * jacobianSignal;
    	double coefficient;
    	double oldInput;
    	double oldOutput;
    };
    /// vector of all entries of the global interface jacobian matrix
    std::vector<interfaceJacobianEntry> interfaceJacobianEntrys;
    /// whether output numbers or not
    bool debugMe;
    /// size of global residual vector
    int globalResidualSize;
    /// current global residual vector
    double *globalResidual;
    /// global corrector vector
    double *correctorVec;
    /// global interface Jacobian matrix
    MathLibrary::SparseMatrix<double> *interfaceJacGlobal;
    /// friend class in unit test
    friend class TestIJCSA;
    /// old function input
    double functionInputold;
    /// new function input
    double functionInput;
    /// old function output
    double functionOutputold;
    /// new function output
    double functionOutput;
};
}/* namespace EMPIRE */
#endif /* IJCSA_H_ */
