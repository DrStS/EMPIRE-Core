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
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "cppunit/extensions/TestFactoryRegistry.h"
#include "cppunit/ui/text/TestRunner.h"
#include "cppunit/XmlOutputter.h"
#include "cppunit/TestResult.h"

#include "ProgressOutputter.h"
#include "Message.h"
#include "EMPEROR_Enum.h"
#include "ConnectionIOSetup.h"

using namespace std;

string pathToFolderOfFiles;

int main(int argc, char **argv) {
    if (argc != 2) {
        cerr << "ERROR: Path to unit test folder is not given!" << endl;
        cerr << "Exit the unit test!" << endl;
        exit(EXIT_FAILURE);
    } else {
        pathToFolderOfFiles = argv[1];
        pathToFolderOfFiles.append("/files/");
    }
    std::ofstream file("unitTestResult.xml");
    EMPIRE::Message::userSetOutputLevel = EMPIRE::Message::DEBUG; // may output some data for checking

    // 1. add registered tests to testRunner
    CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(registry.makeTest());

    // 2. set listener which outputs the progress to shell
    EMPIRE::ProgressOutputter testProgressOutputter;
    runner.eventManager().addListener(&testProgressOutputter);

    // 3. run
    runner.run();

    // 4. output to xml file
    CppUnit::XmlOutputter xmlOutputter(&runner.result(), file);
    xmlOutputter.write();

    return 0;
}


