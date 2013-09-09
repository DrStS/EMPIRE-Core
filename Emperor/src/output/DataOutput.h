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
 * \file DataOutput.h
 * This file holds the class DataOutput
 * \date 8/22/2012
 **************************************************************************************************/

#ifndef DATAOUTPUT_H_
#define DATAOUTPUT_H_

#include <string>
#include <vector>
#include <map>

namespace EMPIRE {

struct structDataOutput;
class ClientCode;
/********//**
 * \brief This class can output meshes and dataFields in GiD format
 ***********/
class DataOutput {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _settingDataOutput setting of this dataOutput (it is set up in MetaDatabase)
     * \param[in] _nameToClientCodeMap a reference to the nameToClientCodeMap (in class Emperor)
     * \author Tianyang Wang
     ***********/
    DataOutput(const structDataOutput &_settingDataOutput,
            std::map<std::string, ClientCode*> &_nameToClientCodeMap);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~DataOutput();
    /***********************************************************************************************
     * \brief Initialize mesh files, data field files, signal files
     * \param[in] rearPart of the name of the files, has non-empty value only in iterativeCouplingLogic
     * \author Tianyang Wang
     ***********/
    void init(std::string rearPart);
    /***********************************************************************************************
     * \brief Write data of current step
     * \param[in] step the step number
     * \author Tianyang Wang
     ***********/
    void writeCurrentStep(int step);
private:
    /// dataOutputName, starting part of outputting files
    std::string dataOutputName;
    /// setting of this DataOutput
    const structDataOutput &settingDataOutput;
    /// a reference to the nameToClientCodeMap of class Emperor
    std::map<std::string, ClientCode*> &nameToClientCodeMap;

    /***********************************************************************************************
     * \brief Write meshes to mesh files
     * \author Tianyang Wang
     ***********/
    void writeMeshes();
    /***********************************************************************************************
     * \brief Initialize all data field files
     * \author Tianyang Wang
     ***********/
    void initDataFieldFiles();
    /***********************************************************************************************
     * \brief Write data field of current step to data files
     * \param[in] step the step number
     * \author Tianyang Wang
     ***********/
    void writeDataFields(int step);
    /***********************************************************************************************
     * \brief Initialize all signal files
     * \author Tianyang Wang
     ***********/
    void initSignalFiles();
    /***********************************************************************************************
     * \brief Write signals of current step to signal files
     * \param[in] step the step number
     * \author Tianyang Wang
     ***********/
    void writeSignals(int step);
};

} /* namespace EMPIRE */
#endif /* DATAOUTPUT_H_ */
