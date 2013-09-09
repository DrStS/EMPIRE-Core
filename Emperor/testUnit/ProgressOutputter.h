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
#ifndef TESTPROGRESSOUTPUTTER_H_
#define TESTPROGRESSOUTPUTTER_H_

#include <iostream>

#include "cppunit/TestListener.h"
#include "cppunit/Test.h"
#include "cppunit/TextTestProgressListener.h"

namespace EMPIRE {
// This is an adaptor of the stupid TextTestProgressListener
// It is especially used to look for the source of segmentation fault
using namespace CppUnit;
using namespace std;
class ProgressOutputter: public TestListener {
public:
    ProgressOutputter() :
            progressListener() {

    }
    virtual ~ProgressOutputter() {
    }

    // Print the name of the test instead of a dot
    void startTest(Test *test) {
        //progressListener.startTest(test);
        cout << endl;
        cout << "testing:      " << test->getName();
        cout << endl;
    }

    void addFailure(const TestFailure &failure) {
        progressListener.addFailure(failure);
    }

    void endTest(Test * test) {
        progressListener.endTest(test);
    }

    void startSuite(Test * suite) {
        progressListener.startSuite(suite);
    }

    void endSuite(Test * suite) {
        progressListener.endSuite(suite);
    }

    void startTestRun(Test *test, TestResult * eventManager) {
        progressListener.startTestRun(test, eventManager);
    }

    void endTestRun(Test *test, TestResult *eventManager) {
        progressListener.endTestRun(test, eventManager);
    }
private:
    CppUnit::TextTestProgressListener progressListener;
};

} /* namespace EMPIRE */

#endif /* TESTPROGRESSOUTPUTTER_H_ */
