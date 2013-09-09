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
 * \file AbstractCouplingAlgorithm.h
 * This file holds the class AbstractCouplingAlgorithm
 * \date 5/15/2012
 **************************************************************************************************/

#ifndef ABSTRACTCOUPLINGALGORITHM_H_
#define ABSTRACTCOUPLINGALGORITHM_H_

#include <string>

namespace EMPIRE {

class ConnectionIO;
/********//**
 * \brief Class AbstractCouplingAlgorithm is the mother class of all coupling algorithms. A
 *              coupling algorithm must be attached to an iterative coupling loop
 ***********/
class AbstractCouplingAlgorithm {
public:
    /***********************************************************************************************
     * \brief Constructor, set input and output (the input and output could be the same memory)
     * \param[in] _name name of the coupling algorithm
     * \author Tianyang Wang
     ***********/
    AbstractCouplingAlgorithm(std::string _name);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~AbstractCouplingAlgorithm();
    /***********************************************************************************************
     * \brief Set the input and output
     * \param[in] _input
     * \param[in] _output
     * \author Tianyang Wang
     ***********/
    virtual void setInputAndOutput(const ConnectionIO *_input, ConnectionIO *_output);
    /***********************************************************************************************
     * \brief Calculate the output by relaxation on the input
     * \author Tianyang Wang
     ***********/
    virtual void calcNewValue() = 0;
    /***********************************************************************************************
     * \brief Set that the new implicit coupling loop starts
     * \author Tianyang Wang
     ***********/
    void setNewLoop();
    /***********************************************************************************************
     * \brief return the current residual
     * \return the current residual
     * \author Tianyang Wang
     ***********/
    double getCurrentResidual();
    /***********************************************************************************************
     * \brief return the initial residual
     * \return the initial residual
     * \author Tianyang Wang
     ***********/
    double getInitialResidual();

protected:
    /// name of the coupling algorithm
    std::string name;
    /// input
    const ConnectionIO* input;
    /// output
    ConnectionIO* output;
    /// whether it is the new loop
    bool newLoop;
    /// the implicit loop number
    int step;
    /// initial residual in fixed point iteration
    double initialResidual;
    /// current residual in fixed point iteration
    double currentResidual;

    /***********************************************************************************************
     * \brief Copy a vector
     * \param[in] from input vector (a double array)
     * \param[in] size size of the array
     * \param[out] to output vector
     * \author Tianyang Wang
     ***********/
    static void vecCopy(const double *from, double *to, int size);
    /***********************************************************************************************
     * \brief Dot product of two vectors
     * \param[in] vec1 the 1st vector
     * \param[in] vec2 the 2nd vector
     * \param[in] size size of the array
     * \return the dot product
     * \author Tianyang Wang
     ***********/
    static double vecDotProduct(const double *vec1, const double *vec2, int size);
    /***********************************************************************************************
     * \brief vector scalar multiplification
     * \param[in/out] vec input vector (a double array)
     * \param[in] size size of the arsignal    * \param[out] to output vector
     * \author Tianyang Wang
     ***********/
    static void vecScalarMultiply(double *vec, const double SCALAR, int size);
    /***********************************************************************************************
     * \brief vec1 = vec1 + vec2
     * \param[in/out] vec1 the 1st vector
     * \param[in] vec2 the 2nd vector
     * \param[in] size size of the array
     * \author Tianyang Wang
     ***********/
    static void vecPlusEqual(double *vec1, const double *vec2, int size);
    /***********************************************************************************************
     * \brief vec1 = vec1 - vec2
     * \param[in/out] vec1 the 1st vector
     * \param[in] vec2 the 2nd vector
     * \param[in] size size of the array
     * \author Tianyang Wang
     ***********/
    static void vecMinusEqual(double *vec1, const double *vec2, int size);
    /***********************************************************************************************
     * \brief vec1minus2 = vec1 - vec2
     * \param[in] vec1 the 1st vector
     * \param[in] vec2 the 2nd vector
     * \param[out] vec1minus2 the 1st vector
     * \param[in] size size of the array
     * \author Tianyang Wang
     ***********/
    static void vecMinus(const double *vec1, const double *vec2, double *vec1minus2,
            int size);
    /***********************************************************************************************
     * \brief compute L2 norm of the vector
     * \param[in] vec the 2nd vector
     * \param[out] L2 norm of the vector
     * \param[in] size size of the array
     * \author Tianyang Wang
     ***********/
    static double vecL2Norm(const double *vec, int size);
};

} /* namespace EMPIRE */
#endif /* ABSTRACTCOUPLINGALGORITHM_H_ */
