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
 * \file ClientMetaDatabase.cpp
 * This file holds the class implementation of the client meta database
 * \date 3/18/2012
 **************************************************************************************************/

#ifndef CLIENTMETADATABASE_H_
#define CLIENTMETADATABASE_H_
#include <vector>
#include <string>
#include <map>
namespace ticpp {
class Document;
}

namespace EMPIRE {
class ClientMetaDatabase {
public:
    /***********************************************************************************************
     * \brief Initialize the client metadatabase
     * \param[in] inputFileName Input file name
     * \author Tianyang Wang
     ***********/
    static void init(char *inputFileName);
    /***********************************************************************************************
     * \brief Delete the singleton
     * \author Tianyang Wang
     ***********/
    static void free();
    /***********************************************************************************************
     * \brief return the singleton
     *
     * \param[in] Input file name
     * \return the singleton
     * \author Stefan Sicklinger
     ***********/
    static ClientMetaDatabase *getSingleton();
	/***********************************************************************************************
	 * \brief getServerPortFile return name of server port file
	 *
	 * \return string serverPortFile
	 * \author Stefan Sicklinger
	 ***********/
	std::string getServerPortFile();
	/***********************************************************************************************
	 * \brief getClientName return name of this client
	 *
	 * \return string clientName
	 * \author Stefan Sicklinger
	 ***********/
	std::string getClientName();

	std::string getUserDefinedText(std::string elementName);
private:
	/***********************************************************************************************
	 * \brief disallow public constructor
	 *
	 * \author Stefan Sicklinger
	 ***********/
	ClientMetaDatabase(char *inputFileName);
	/***********************************************************************************************
	 * \brief disallow public copy constructor
	 *
	 * \author Stefan Sicklinger
	 ***********/
	ClientMetaDatabase(const ClientMetaDatabase&);
	/***********************************************************************************************
	 * \brief disallow assignment operator
	 *
	 * \author Stefan Sicklinger
	 ***********/
	ClientMetaDatabase& operator=(const ClientMetaDatabase&);
	/***********************************************************************************************
	 * \brief Destructor
	 *
	 * \author Stefan Sicklinger
	 ***********/
	virtual ~ClientMetaDatabase();
	/***********************************************************************************************
	 * \brief Fills the string clientName
	 *
	 * \author Stefan Sicklinger
	 ***********/
	void fillClientName();
	/***********************************************************************************************
	 * \brief Fills the map fillServerPortFile
	 *
	 * \author Stefan Sicklinger
	 ***********/
	void fillServerPortFile();
    /***********************************************************************************************
     * \brief Fill the map fillVerbosity
     * \author Stefan Sicklinger
     ***********/
    void fillVerbosity();
    /***********************************************************************************************
     * \brief Fill the map fillServerPortFile
     * \author Stefan Sicklinger
     ***********/
    void fillrunTimeModifiableUserDefinedText();
    /***********************************************************************************************
     * \brief Compare two string case insensitive
     * \return true or false
     * \author Stefan Sicklinger
     ***********/
    bool CompareStringInsensitive(std::string strFirst, std::string strSecond);
	/// The singleton of this class
	static ClientMetaDatabase* clientMetaDatabase;
	///Object of the tinyxml++ parser
	ticpp::Document *inputFile;
	///String for ServerPort file
	std::string serverPortFile;
	///String for ClientName (Name of this client)
	std::string clientName;
	///boolean ifgetUserDefinedText is runTimeModifiable
	bool runTimeModifiableUserDefinedText;
    /// verbosity
    std::string verbosity;
};

} /* namespace EMPIRE */
#endif /* CLIENTMETADATABASE_H_ */
