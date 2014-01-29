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
 * \file MetaDatabase.cpp
 * This file holds the class implementation of the meta database
 * \date 3/12/2012
 **************************************************************************************************/

#ifndef METADATABASE_H_
#define METADATABASE_H_
#include <vector>
#include <string>
#include <map>

#include "EMPEROR_Enum.h"
#include "MetaDataStructures.h"

namespace ticpp {
class Document;
class Element;
}

namespace EMPIRE {
/********//**
 * \brief Class MetaDatabase parses the XML input file, holds the settings of the co-simulation
 ***********/
class MetaDatabase {
public:
    /***********************************************************************************************
     * \brief Initialize the metadatabase singleton by parsing the input file
     * \param[in] inputFileName name of the input file
     * \author Tianyang Wang
     ***********/
    static void init(char *inputFileName);
    /***********************************************************************************************
     * \brief return the singleton
     * \return the singleton
     * \author Stefan Sicklinger
     ***********/
    static MetaDatabase *getSingleton();
    /***********************************************************************************************
     * \brief Destructor
     * \author Stefan Sicklinger
     ***********/
    virtual ~MetaDatabase();
    /***********************************************************************************************
     * \brief Check if client code is in database
     * \param[in] ClientCode name which should be tested
     * \return true if client code is defined in the input file (otherwise false)
     * \author Stefan Sicklinger
     ***********/
    bool checkForClientCodeName(std::string clientName);

    ///String for ServerPort file
    std::string serverPortFile;
    /// verbosity
    std::string verbosity;
    /// setting of client codes in XML input file
    std::vector<structClientCode> settingClientCodeVec;
    /// setting of data outputs in XML input file
    std::vector<structDataOutput> settingDataOutputVec;
    /// setting of mappers in XML input file
    std::vector<structMapper> settingMapperVec;
    /// setting of coupling algorithms in XML input file
    std::vector<structCouplingAlgorithm> settingCouplingAlgorithmVec;
    /// setting of extrapolators in XML input file
    std::vector<structExtrapolator> settingExtrapolatorVec;
    /// setting of connections in XML input file
    std::vector<structConnection> settingConnectionVec;
    // correspond to coSimulation block in XML input file
    structCouplingLogic settingGlobalCouplingLogic;

private:
    /***********************************************************************************************
     * \brief do nothing, only used by unit test class
     * \author Tianyang Wang
     ***********/
    MetaDatabase();
    /***********************************************************************************************
     * \brief disallow public constructor
     * \author Stefan Sicklinger
     ***********/
    MetaDatabase(char *inputFileName);
    /***********************************************************************************************
     * \brief disallow public copy constructor
     * \author Stefan Sicklinger
     ***********/
    MetaDatabase(const MetaDatabase&);
    /***********************************************************************************************
     * \brief disallow assignment operator
     * \author Stefan Sicklinger
     ***********/
    MetaDatabase& operator=(const MetaDatabase&);
    /***********************************************************************************************
     * \brief Fill the map fillServerPortFile
     * \author Stefan Sicklinger
     ***********/
    void fillServerPortFile();
    /***********************************************************************************************
     * \brief Fill the map fillVerbosity
     * \author Stefan Sicklinger
     ***********/
    void fillVerbosity();
    /***********************************************************************************************
     * \brief Fill client code setting by parsing XML input file
     * \author Tianyang Wang
     ***********/
    void fillSettingClientCodesVec();
    /***********************************************************************************************
     * \brief Fill data output setting by parsing XML input file
     * \author Tianyang Wang
     ***********/
    void fillSettingDataOutputVec();
    /***********************************************************************************************
     * \brief Fill mapper setting by parsing XML input file
     * \author Tianyang Wang
     ***********/
    void fillSettingMapperVec();
    /***********************************************************************************************
     * \brief Fill coupling algorithm setting by parsing XML input file
     * \author Tianyang Wang
     ***********/
    void fillSettingCouplingAlgorithmVec();
    /***********************************************************************************************
     * \brief Fill extrapolator setting by parsing XML input file
     * \author Tianyang Wang
     ***********/
    void fillSettingExtrapolatorVec();
    /***********************************************************************************************
     * \brief Fill connection setting by parsing XML input file
     * \author Tianyang Wang
     ***********/
    void fillSettingConnectionVec();
    /***********************************************************************************************
     * \brief Fill coupling logic setting by parsing XML input file
     * \author Tianyang Wang
     ***********/
    void fillSettingCouplingLogic();
    /***********************************************************************************************
     * \brief Parse DataFieldRef or SignalRef
     * \param[in] xmlElement xml element
     * \return the setting after parsing
     * \author Tianyang Wang
     ***********/
    structConnectionIO parseConnectionIORef(ticpp::Element *xmlElement);
    /***********************************************************************************************
     * \brief Parse a vector of DataFieldRef or SignalRef
     * \param[in] xmlElement xml element
     * \return the setting after parsing
     * \author Tianyang Wang
     ***********/
    std::vector<structConnectionIO> parseConnectionIORefs(ticpp::Element *xmlElement);

    /// The singleton of this class
    static MetaDatabase* metaDatabase;
    /// input xml file
    ticpp::Document *inputFile;
    /// the unit test class
    friend class TestEmperor;
    /// the unit test class
    friend class TestMetaDatabase;
};

} /* namespace EMPIRE */
#endif /* METADATABASE_H_ */
