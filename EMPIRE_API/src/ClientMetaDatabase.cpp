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
#include <iostream>
#include <assert.h>
#include "ClientMetaDatabase.h"
#include "ticpp.h"
using namespace std;
using namespace ticpp;
namespace EMPIRE {
ClientMetaDatabase *ClientMetaDatabase::clientMetaDatabase = NULL;

void ClientMetaDatabase::init(char *inputFileName) {
    assert(clientMetaDatabase == NULL);
    clientMetaDatabase = new ClientMetaDatabase(inputFileName);
}

void ClientMetaDatabase::free() {
    delete clientMetaDatabase;
    clientMetaDatabase = NULL;
}

ClientMetaDatabase *ClientMetaDatabase::getSingleton() {
    assert(clientMetaDatabase != NULL);
    return clientMetaDatabase;
}

ClientMetaDatabase::ClientMetaDatabase(char *inputFileName) {

    try {
        inputFile = new Document(inputFileName);
        inputFile->LoadFile();
        /// Fill up data base
        fillClientName();
        fillServerPortFile();
    } catch (ticpp::Exception& ex) {
        cerr << "ERROR Parser: " << ex.what() << endl;
        exit(EXIT_FAILURE);
    }

}

ClientMetaDatabase::~ClientMetaDatabase() {
    delete inputFile;
}

std::string ClientMetaDatabase::getClientName() {
    return (clientName);
}

string ClientMetaDatabase::getUserDefinedText(string elementName) {
    Element *userDefinedBlock = inputFile->FirstChildElement()->FirstChildElement("userDefined");
    Element *userDefined = userDefinedBlock->FirstChildElement(elementName, false);
    if (userDefined == NULL) {
        return "";
    }
    return userDefined->GetText();
}

void ClientMetaDatabase::fillClientName() {
    Element *pXMLElement = inputFile->FirstChildElement()->FirstChildElement("code");
    clientName = pXMLElement->GetAttribute("name");
}

void ClientMetaDatabase::fillServerPortFile() {
    Element *pXMLElement =
            inputFile->FirstChildElement()->FirstChildElement("general")->FirstChildElement(
                    "portFile");
    serverPortFile = pXMLElement->GetText();

}

string ClientMetaDatabase::getServerPortFile() {
    return (serverPortFile);
}


} /* namespace EMPIRE */
