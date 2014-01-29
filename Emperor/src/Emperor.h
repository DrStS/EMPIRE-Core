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
 * \file Emperor.h
 * This file holds the class of Emperor.
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef EMPEROR_H_
#define EMPEROR_H_
#include <map>
#include <string>
#include <vector>

namespace EMPIRE {

class ServerCommunication;
class ClientCode;
class DataOutput;
class Connection;
class MetaDatabase;
class MapperAdapter;
class AbstractCouplingAlgorithm;
class AbstractExtrapolator;
class AbstractCouplingLogic;
class ConnectionIO;
struct structCouplingLogic;
struct structDataFieldRef;
struct structConnectionIO;
/********//**
 * \brief This class manages the program, provides the interface functions of the program
 *        which are called in the main function
 **************************************************************************************************/
class Emperor {
public:
    /***********************************************************************************************
     * \brief Constructor, initialize member variables
     * \author Tianyang Wang
     ***********/
    Emperor();
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~Emperor();
    /***********************************************************************************************
     * \brief Initializes Meta database (parsing done) and ServerCommunication
     * \param[in] argc pointer to number of command line arguments
     * \param[in] argv pointer to char array which holds the command line arguments
     * \author Tianyang Wang, Stefan Sicklinger
     ***********/
    void initEnvironment(int *argc, char ***argv);
    /***********************************************************************************************
     * \brief Start server listening (listen to new clients)
     * \author Stefan Sicklinger
     ***********/
    void startServerListening();
    /***********************************************************************************************
     * \brief Executes the complete EMPEROR Co-Simulation.
     * \author Tianyang Wang
     ***********/
    void startServerCoupling();

private:
    /***********************************************************************************************
     * \brief Wait until all clients specified in the XML input file are connected
     * \author Stefan Sicklinger
     ***********/
    void connectAllClients();
    /***********************************************************************************************
     * \brief Disconnect all clients
     * \author Stefan Sicklinger
     ***********/
    void disconnectAllClients();
    /***********************************************************************************************
     * \brief Initialize ClientCode instances by the info sent by the real clients
     * \author Tianyang Wang
     ***********/
    void initClientCodes();
    /***********************************************************************************************
     * \brief Initialize DataOutput instances
     * \author Tianyang Wang
     ***********/
    void initDataOutputs();
    /***********************************************************************************************
     * \brief Initialize mappers based on interface meshes of different clients
     * \author Tianyang Wang
     ***********/
    void initMappers();
    /***********************************************************************************************
     * \brief Initialize coupling algorithms
     * \author Tianyang Wang
     ***********/
    void initCouplingAlgorithms();
    /***********************************************************************************************
     * \brief Initialize extrapolators
     * \author Tianyang Wang
     ***********/
    void initExtrapolators();
    /***********************************************************************************************
     * \brief Initialize connections
     * \author Tianyang Wang
     ***********/
    void initConnections();
    /***********************************************************************************************
     * \brief Initialize globalCouplingLogic
     * \author Tianyang Wang
     ***********/
    void initGlobalCouplingLogic();
    /***********************************************************************************************
     * \brief Parse StructCouplingLogic which is called by initCouplingLogicSequence()
     * \param[in] couplingLogicStruct pointer to structCouplingLogic coming from MetaDatabase
     * \return AbstractCouplingLogic which is made from information coming from couplingLogicStruct
     * \author Tianyang Wang
     ***********/
    AbstractCouplingLogic *parseStructCouplingLogic(structCouplingLogic &couplingLogicStruct);
    /***********************************************************************************************
     * \brief This does the CoSimulation management
     * \author Stefan Sicklinger
     ***********/
    void doCoSimulation();
    /***********************************************************************************************
     * \brief Construct ConnectionIO, avoid copying same code in different places
     * \author Tianyang Wang
     ***********/
    ConnectionIO *constructConnectionIO(const structConnectionIO &settingConnectionIO);

    /// Get ClientCode instance by its name
    std::map<std::string, ClientCode*> nameToClientCodeMap;
    /// Get DataOutput instance by its name
    std::map<std::string, DataOutput*> nameToDataOutputMap;
    /// Get Connection instance by its name
    std::map<std::string, Connection*> nameToConnetionMap;
    /// Get Mapper instance by its name
    std::map<std::string, MapperAdapter*> nameToMapperMap;
    /// Get coupling algorithm instance by its name
    std::map<std::string, AbstractCouplingAlgorithm*> nameToCouplingAlgorithmMap;
    /// Get extrapolator instance by its name
    std::map<std::string, AbstractExtrapolator*> nameToExtrapolatorMap;
    /// The global coupling logic (corresponds to the coSimulation block in the XML file)
    AbstractCouplingLogic *globalCouplingLogic;
    /// Collection of coupling logics which is only used for destruction
    std::vector<AbstractCouplingLogic*> couplingLogicVec;
    /// the unit test class
    friend class TestEmperor;
    /// the unit test class
    friend class TestMetaDatabase;
};

} /* namespace EMPIRE */

#endif /* EMPEROR_H_ */
