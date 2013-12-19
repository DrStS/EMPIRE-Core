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
 * \file Residual.h
 * This file holds the class Residual
 * \date 12/18/2013
 **************************************************************************************************/
#ifndef RESIDUAL_H_
#define RESIDUAL_H_

#include <string>
#include <vector>

namespace EMPIRE {

class ConnectionIO;
/********//**
 * \brief Class Residual is the residual in the coupling algorithm
 ***********/
class Residual {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _index the index of the residual
     * \author Tianyang Wang
     ***********/
    Residual(int _index);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~Residual();
    /***********************************************************************************************
     * \brief Add a component to the components
     * \param[in] _coefficient the coefficient
     * \param[in] _timeToUpdate "iterationBeginning" or "iterationEnd"
     * \param[in] _reference the io
     * \author Tianyang Wang
     ***********/
    void addComponent(double _coefficient, std::string _timeToUpdate, ConnectionIO *_reference);
    /***********************************************************************************************
     * \brief Init residualVector, must be called after all components are added
     * \author Tianyang Wang
     ***********/
    void init();
    /***********************************************************************************************
     * \brief Update components at the iteration beginning
     * \author Tianyang Wang
     ***********/
    void updateAtIterationBeginning();
    /***********************************************************************************************
     * \brief Update components at the iteration end
     * \author Tianyang Wang
     ***********/
    void updateAtIterationEnd();
    /***********************************************************************************************
     * \brief Compute the current residual (both residualVector and residualVectorL2Norm)
     * \author Tianyang Wang
     ***********/
    void computeCurrentResidual();
    /// size of the array residualVector
    int size;
    /// the residual vector
    double *residualVector;
    /// the L2 norm of the residual vector
    double residualVectorL2Norm;
    /// the unit test class
    friend class TestResidual;
private:
    /********//**
     * \brief Class Component is the component in the residual
     ***********/
    class Component {
    public:
        /***********************************************************************************************
         * \brief Constructor
         * \param[in] _coefficient the coefficient
         * \param[in] _timeToUpdate "iterationBeginning" or "iterationEnd"
         * \param[in] _reference the io
         * \author Tianyang Wang
         ***********/
        Component(double _coefficient, std::string _timeToUpdate, ConnectionIO *_reference);
        /***********************************************************************************************
         * \brief Destructor
         * \author Tianyang Wang
         ***********/
        virtual ~Component();
        /***********************************************************************************************
         * \brief Update itself at the iteration beginning
         * \author Tianyang Wang
         ***********/
        void updateAtIterationBeginning();
        /***********************************************************************************************
         * \brief Update itself at the iteration end
         * \author Tianyang Wang
         ***********/
        void updateAtIterationEnd();
        /// copy of the array in the reference
        double *dataCopy;
        /// size of the array
        int size;
        /// the coefficient of itself in the residual
        double coefficient;
    private:
        /// the reference to the io object
        ConnectionIO *reference;
        /// "iterationBeginning" or "iterationEnd"
        std::string timeToUpdate;
    };
    /// all the components
    std::vector<Component*> components;
    /// the index of the residual
    int index;
};

} /* namespace EMPIRE */
#endif /* RESIDUAL_H_ */
