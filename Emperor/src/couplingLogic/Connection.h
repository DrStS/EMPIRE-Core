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
 * \file Connection.h
 * This file holds the class Connection
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef CONNECTION_H_
#define CONNECTION_H_

#include <string>
#include <vector>
#include "AbstractCouplingLogic.h"
#include "EMPEROR_Enum.h"

namespace EMPIRE {

class ClientCode;
class AbstractFilter;
class AbstractExtrapolator;
class ConnectionIO;

/********//**
 * \brief Class Connection sends and receives data between two clients and does data processing
 ***********/
class Connection: public AbstractCouplingLogic {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name name of the Connection
     * \author Tianyang Wang
     ***********/
    Connection(std::string _name);
    /***********************************************************************************************
     * \brief Destructor, delete the filters and the extrapolator.
     * \author Tianyang Wang
     ***********/
    virtual ~Connection();
    /***********************************************************************************************
     * \brief Do coupling, which simply calls the function transferData()
     * \author Tianyang Wang
     ***********/
    void doCoupling();
    /***********************************************************************************************
     * \brief Main job of this class, i.e. receive data field, filter it by a sequence of filters, and send it
     * \author Tianyang Wang
     ***********/
    void transferData();

    void addInput(ConnectionIO *input);
    void addOutput(ConnectionIO *output);

    /***********************************************************************************************
     * \brief Add a filter into the filter sequence
     * \param[in] filter pointer to the filter
     * \author Tianyang Wang
     ***********/
    void addFilter(AbstractFilter *filter);
protected:
    /// name of the connection
    std::string name;
    /// inputs
    std::vector<ConnectionIO*> inputVec;
    /// outputs
    std::vector<ConnectionIO*> outputVec;
    /// sequence of filters
    std::vector<AbstractFilter*> filterVec;
    /// the unit test class
    friend class TestConnection;
};

} /* namespace EMPIRE */
#endif /* CONNECTION_H_ */
