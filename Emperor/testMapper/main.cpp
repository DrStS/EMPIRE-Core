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
#include "TestMapper.h"

using namespace EMPIRE;

int main(int argc, char** argv) {
    assert(argc == 2);
    TestMapper *testMapper = new TestMapper();
    testMapper->parseInputFile(argv[1]);
    testMapper->outputSettingToShell();
    testMapper->doMapping();
    testMapper->doConservationAnalysis();
    delete testMapper;
    cout << endl << "================================================" << endl;
    cout << "Feierabend!" << endl;
}
