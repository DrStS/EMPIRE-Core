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
 * \file PseudoCodeOutput.h
 * This file holds the class PseudoCodeOutput
 * \date 9/3/2012
 **************************************************************************************************/
#ifndef PSEUDOCODEOUTPUT_H_
#define PSEUDOCODEOUTPUT_H_

#include <string>
#include <fstream>

namespace EMPIRE {

class MetaDatabase;
struct structConnection;
struct structCouplingLogic;
struct structClientCode;

/********//**
 * \brief This class can output pseudo codes of the server and the clients according to the input
 *        file. Not everything, but the communication part.
 ***********/
class PseudoCodeOutput {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _metaDatabase metaDatabase which parses the input file
     * \param[in] fileName name of the pseudo code file
     * \author Tianyang Wang
     ***********/
    PseudoCodeOutput(MetaDatabase *_metaDatabase, std::string fileName);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~PseudoCodeOutput();
    /***********************************************************************************************
     * \brief Write the pseudo code to the file
     * \author Tianyang Wang
     ***********/
    void writePseudoCode();
private:
    /// output file stream
    std::fstream outputStream;
    /// metadatabase (TODO use MetaDatabase::getsingleton() to get it instead of defining a member)
    MetaDatabase *metaDatabase;
    /// indents of the current line (number of "\t" )
    std::string indents;
    /***********************************************************************************************
     * \brief Write the code of the server
     * \author Tianyang Wang
     ***********/
    void writeServerCode();
    /***********************************************************************************************
     * \brief Add receive meshes to the server code
     * \author Tianyang Wang
     ***********/
    void addRecvMeshToServerCode();
    /***********************************************************************************************
     * \brief Add coupling logics to the server code
     * \author Tianyang Wang
     ***********/
    void addCouplingLogicToServerCode();
    /***********************************************************************************************
     * \brief Add coupling logic to the server code (recursive function)
     * \param[in] settingCouplingLogic the coupling logic to write
     * \author Tianyang Wang
     ***********/
    void addCouplingLogicToServerCode(structCouplingLogic &settingCouplingLogic);
    /***********************************************************************************************
     * \brief Add connections to the server code (recursive function)
     * \param[in] settingconnection the connection to write
     * \author Tianyang Wang
     ***********/
    void addConnectionToServerCode(structConnection &settingconnection);
    /***********************************************************************************************
     * \brief Write the code of the clients
     * \author Tianyang Wang
     ***********/
    void writeClientCodes();
    /***********************************************************************************************
     * \brief Write the code of a client
     * \param[in] settingClientCode the setting of the client code
     * \author Tianyang Wang
     ***********/
    void writeClientCode(structClientCode &settingClientCode);
    /***********************************************************************************************
     * \brief Write coupling logic to the client code
     * \param[in] clientCodeName name of the client code
     * \author Tianyang Wang
     ***********/
    void addCouplingLogicToClientCode(std::string clientCodeName);
    /***********************************************************************************************
     * \brief Write coupling logic to the client code (recursive function)
     * \param[in] clientCodeName name of the client code
     * \param[in] settingCouplingLogic the coupling logic to write
     * \author Tianyang Wang
     ***********/
    void addCouplingLogicToClientCode(std::string clientCodeName,
            structCouplingLogic &settingCouplingLogic);
    /***********************************************************************************************
     * \brief Write connection to the client code
     * \param[in] clientCodeName name of the client code
     * \param[in] settingConnection the connection to write
     * \author Tianyang Wang
     ***********/
    void addConnectionToClientCode(std::string clientCodeName, structConnection &settingConnection);
    /***********************************************************************************************
     * \brief Is client taking part in the coupling logic
     * \param[in] clientCodeName name of the client code
     * \param[in] settingCouplingLogic the setting of the couplingLogic
     * \author Tianyang Wang
     ***********/
    bool isClientInCouplingLogic(std::string clientCodeName,
            structCouplingLogic &settingCouplingLogic);
    /***********************************************************************************************
     * \brief Get a structConnection by its name
     * \param[in] connectionName name of the connection
     * \author Tianyang Wang
     ***********/
    structConnection &getStructConnectionByName(std::string connectionName);
    /***********************************************************************************************
     * \brief Set number of indents
     * \param[in] n number of indents
     * \author Tianyang Wang
     ***********/
    void setIndents(int n);
    /***********************************************************************************************
     * \brief Increment number of indents
     * \author Tianyang Wang
     ***********/
    void incrementIndents();
    /***********************************************************************************************
     * \brief Decrement number of indents
     * \author Tianyang Wang
     ***********/
    void decrementIndents();
    /***********************************************************************************************
     * \brief Add a block start
     * \param[in] str block information
     * \author Tianyang Wang
     ***********/
    void addBlockStart(std::string str);
    /***********************************************************************************************
     * \brief Add a few empty lines
     * \param[in] n number of empty lines
     * \author Tianyang Wang
     ***********/
    void addEmptyLines(int n);
    /***********************************************************************************************
     * \brief Add a line of description
     * \param[in] str the description
     * \author Tianyang Wang
     ***********/
    void addDescriptionLine(std::string str);
    /***********************************************************************************************
     * \brief Add quotation around the string
     * \param[in] str string without quotation
     * \return string with quotation
     * \author Tianyang Wang
     ***********/
    std::string addQuatation(std::string str);
};

} /* namespace EMPIRE */
#endif /* PSEUDOCODEOUTPUT_H_ */
