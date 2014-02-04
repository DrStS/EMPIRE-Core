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
#include <assert.h>

#include "PseudoCodeOutput.h"
#include "MetaDatabase.h"

using namespace std;

namespace EMPIRE {

PseudoCodeOutput::PseudoCodeOutput(MetaDatabase *_metaDatabase, std::string fileName) :
        metaDatabase(_metaDatabase), outputStream(fileName.c_str(), fstream::out) {
}

PseudoCodeOutput::~PseudoCodeOutput() {
}

void PseudoCodeOutput::writePseudoCode() {
    writeServerCode();
    writeClientCodes();
    outputStream << endl << endl;
    outputStream.close();
}

void PseudoCodeOutput::writeServerCode() {
    setIndents(0);
    addBlockStart("Pseudo code of the server (EMPEROR)");

    outputStream << indents << "int main() {" << endl;
    incrementIndents();
    {
        //addDescriptionLine("Parsing input file");
        //addDescriptionLine("Initialization");
        addDescriptionLine("Stage 1: receive meshes");
        addRecvMeshToServerCode();
        addDescriptionLine("Stage 2: do co-simulation");
        addCouplingLogicToServerCode();
        addDescriptionLine("");
        outputStream << indents << "return;" << endl;
    }
    decrementIndents();
    outputStream << indents << "}" << endl;
}

void PseudoCodeOutput::addRecvMeshToServerCode() {
    vector<structClientCode> &settingClientCodesVec = metaDatabase->settingClientCodeVec;
    for (int i = 0; i < settingClientCodesVec.size(); i++) {
        string clientCodeName = settingClientCodesVec[i].name;
        for (int j = 0; j < settingClientCodesVec[i].meshes.size(); j++) {
            string meshName = settingClientCodesVec[i].meshes[j].name;
            outputStream << indents << "receiveMeshFrom(" << addQuatation(clientCodeName) << ", "
                    << addQuatation(meshName) << ");" << endl;
        }
    }
}

void PseudoCodeOutput::addCouplingLogicToServerCode() {
    addCouplingLogicToServerCode(metaDatabase->settingGlobalCouplingLogic);
}

void PseudoCodeOutput::addCouplingLogicToServerCode(structCouplingLogic &settingCouplingLogic) {
    if (settingCouplingLogic.type == EMPIRE_connection) {
        string connectionName = settingCouplingLogic.connectionRef.connectionName;
        structConnection &settingConnection = getStructConnectionByName(connectionName);
        addDescriptionLine("");
        addConnectionToServerCode(settingConnection);
        addDescriptionLine("");
    } else if (settingCouplingLogic.type == EMPIRE_CouplingLogicSequence) {
        for (int i = 0; i < settingCouplingLogic.sequence.size(); i++) {
            addCouplingLogicToServerCode(settingCouplingLogic.sequence[i]);
        }
    } else if (settingCouplingLogic.type == EMPIRE_TimeStepLoop) {
        outputStream << indents << "for (int i=0; i<numberOfTimeSteps; i++) {"
                << "// time step loop" << endl;
        incrementIndents();
        {
            addDescriptionLine("");
            for (int i = 0; i < settingCouplingLogic.sequence.size(); i++) {
                addCouplingLogicToServerCode(settingCouplingLogic.sequence[i]);
            }
        }
        decrementIndents();
        outputStream << indents << "}" << endl;
    } else if (settingCouplingLogic.type == EMPIRE_IterativeCouplingLoop) {
        outputStream << indents << "while (!isConvergent) { // iterative coupling loop" << endl;
        incrementIndents();
        {
            for (int i = 0; i < settingCouplingLogic.sequence.size(); i++) {
                addCouplingLogicToServerCode(settingCouplingLogic.sequence[i]);
            }
            vector<string> &convergenceObservers =
                    settingCouplingLogic.iterativeCouplingLoop.convergenceObservers;
            for (int i = 0; i < convergenceObservers.size(); i++) {
                string clientCodeName = convergenceObservers[i];
                outputStream << indents << "sendConvergenceSignalTo(" << clientCodeName
                        << ", isConvergent" << ");" << endl;
            }
        }
        decrementIndents();
        outputStream << indents << "}" << endl;
    } else if (settingCouplingLogic.type == EMPIRE_OptimizationLoop) {
        outputStream << indents << "while (!isConvergent) { // optimization loop" << endl;
        incrementIndents();
        {
            for (int i = 0; i < settingCouplingLogic.sequence.size(); i++) {
                addCouplingLogicToServerCode(settingCouplingLogic.sequence[i]);
            }
            string &convergenceSignalSender =
                    settingCouplingLogic.optimizationLoop.convergenceSignalSender;
            outputStream << indents << "receiveConvergenceSignalFrom("
                    << convergenceSignalSender << ", isConvergent" << ");" << endl;

            vector<string> &convergenceSignalReceivers =
                    settingCouplingLogic.optimizationLoop.convergenceSignalReceivers;
            for (int i = 0; i < convergenceSignalReceivers.size(); i++) {
                string clientCodeName = convergenceSignalReceivers[i];
                outputStream << indents << "sendConvergenceSignalTo(" << clientCodeName
                        << ", isConvergent" << ");" << endl;
            }
        }
        decrementIndents();
        outputStream << indents << "}" << endl;
    } else {
        assert(false);
    }
}

void PseudoCodeOutput::addConnectionToServerCode(structConnection &settingConnection) {
    for (int i = 0; i < settingConnection.inputs.size(); i++) {
        if (settingConnection.inputs[i].type == EMPIRE_ConnectionIO_DataField) {
            const structDataFieldRef &input = settingConnection.inputs[i].dataFieldRef;
            outputStream << indents << "receiveDataFieldFrom(" << addQuatation(input.clientCodeName)
                    << ", " << addQuatation(input.meshName) << ", "
                    << addQuatation(input.dataFieldName) << ");" << endl;
        } else if (settingConnection.inputs[i].type == EMPIRE_ConnectionIO_Signal) {
            const structSignalRef &input = settingConnection.inputs[i].signalRef;
            outputStream << indents << "receiveSignalFrom(" << addQuatation(input.clientCodeName)
                    << ", " << addQuatation(input.signalName) << ");" << endl;
        } else {
            assert(false);
        }
    }
    for (int i = 0; i < settingConnection.outputs.size(); i++) {
        if (settingConnection.outputs[i].type == EMPIRE_ConnectionIO_DataField) {
            const structDataFieldRef &output = settingConnection.outputs[i].dataFieldRef;
            outputStream << indents << "sendDataFieldTo(" << addQuatation(output.clientCodeName)
                    << ", " << addQuatation(output.meshName) << ", "
                    << addQuatation(output.dataFieldName) << ");" << endl;
        } else if (settingConnection.outputs[i].type == EMPIRE_ConnectionIO_Signal) {
            const structSignalRef &output = settingConnection.outputs[i].signalRef;
            outputStream << indents << "sendSignalTo(" << addQuatation(output.clientCodeName)
                    << ", " << addQuatation(output.signalName) << ");" << endl;
        } else {
            assert(false);
        }
    }
}

void PseudoCodeOutput::writeClientCodes() {
    addEmptyLines(3);
    vector<structClientCode> &settingClientCodesVec = metaDatabase->settingClientCodeVec;
    for (int i = 0; i < settingClientCodesVec.size(); i++) {
        structClientCode &settingClientCode = settingClientCodesVec[i];
        writeClientCode(settingClientCode);
        addEmptyLines(3);
    }
}

void PseudoCodeOutput::writeClientCode(structClientCode &settingClientCode) {
    string clientCodeName = settingClientCode.name;
    setIndents(0);
    addBlockStart("Pseudo code of the client code (" + clientCodeName + ")");

    outputStream << indents << "int main() {" << endl;
    incrementIndents();
    {
        addDescriptionLine("Stage 1: send meshes");
        for (int i = 0; i < settingClientCode.meshes.size(); i++) {
            outputStream << indents << "sendMeshToServer("
                    << addQuatation(settingClientCode.meshes[i].name) << ")" << endl;
        }

        addDescriptionLine("Stage 2: do co-simulation");
        addCouplingLogicToClientCode(clientCodeName);
        addDescriptionLine("");
        outputStream << indents << "return;" << endl;
    }
    decrementIndents();
    outputStream << indents << "}" << endl;
}

void PseudoCodeOutput::addCouplingLogicToClientCode(std::string clientCodeName) {
    addCouplingLogicToClientCode(clientCodeName, metaDatabase->settingGlobalCouplingLogic);
}

void PseudoCodeOutput::addCouplingLogicToClientCode(std::string clientCodeName,
        structCouplingLogic &settingCouplingLogic) {
    if (isClientInCouplingLogic(clientCodeName, settingCouplingLogic)) {
        if (settingCouplingLogic.type == EMPIRE_connection) {
            string connectionName = settingCouplingLogic.connectionRef.connectionName;
            structConnection &settingConnection = getStructConnectionByName(connectionName);
            addDescriptionLine("");
            addConnectionToClientCode(clientCodeName, settingConnection);
            addDescriptionLine("");
        } else if (settingCouplingLogic.type == EMPIRE_CouplingLogicSequence) {
            for (int i = 0; i < settingCouplingLogic.sequence.size(); i++) {
                addCouplingLogicToClientCode(clientCodeName, settingCouplingLogic.sequence[i]);
            }
        } else if (settingCouplingLogic.type == EMPIRE_TimeStepLoop) {
            outputStream << indents << "for (int i=0; i<numberOfTimeSteps; i++) {"
                    << "// time step loop" << endl;
            incrementIndents();
            {
                addDescriptionLine("");
                for (int i = 0; i < settingCouplingLogic.sequence.size(); i++) {
                    addCouplingLogicToClientCode(clientCodeName, settingCouplingLogic.sequence[i]);
                }
            }
            decrementIndents();
            outputStream << indents << "}" << endl;
        } else if (settingCouplingLogic.type == EMPIRE_IterativeCouplingLoop) {
            outputStream << indents << "while (!isConvergent) { // iterative coupling loop" << endl;
            incrementIndents();
            {
                for (int i = 0; i < settingCouplingLogic.sequence.size(); i++) {
                    addCouplingLogicToClientCode(clientCodeName, settingCouplingLogic.sequence[i]);
                }
                vector<string> &convergenceObservers =
                        settingCouplingLogic.iterativeCouplingLoop.convergenceObservers;
                for (int i = 0; i < convergenceObservers.size(); i++) {
                    if (clientCodeName == convergenceObservers[i]) {
                        outputStream << indents
                                << "isConvergent = receiveConvergenceSignalFromServer();" << endl;
                    }
                }
            }
            decrementIndents();
            outputStream << indents << "}" << endl;
        } else if (settingCouplingLogic.type == EMPIRE_OptimizationLoop) {
            outputStream << indents << "while (!isConvergent) { // optimization loop" << endl;
            incrementIndents();
            {
                for (int i = 0; i < settingCouplingLogic.sequence.size(); i++) {
                    addCouplingLogicToClientCode(clientCodeName, settingCouplingLogic.sequence[i]);
                }
                vector<string> &convergenceSignalReceivers =
                        settingCouplingLogic.optimizationLoop.convergenceSignalReceivers;
                for (int i = 0; i < convergenceSignalReceivers.size(); i++) {
                    if (clientCodeName == convergenceSignalReceivers[i]) {
                        outputStream << indents
                                << "isConvergent = receiveConvergenceSignalFromServer();" << endl;
                    }
                }
                string convergenceSignalSender =
                        settingCouplingLogic.optimizationLoop.convergenceSignalSender;
                if (clientCodeName == convergenceSignalSender) {
                    outputStream << indents
                            << "sendConvergenceSignalToServer(isConvergent);" << endl;
                }
            }
            decrementIndents();
            outputStream << indents << "}" << endl;
        } else {
            assert(false);
        }
    }
}

void PseudoCodeOutput::addConnectionToClientCode(std::string clientCodeName,
        structConnection &settingConnection) {
    for (int i = 0; i < settingConnection.inputs.size(); i++) {
        if (settingConnection.inputs[i].type == EMPIRE_ConnectionIO_DataField) {
            const structDataFieldRef &input = settingConnection.inputs[i].dataFieldRef;
            if (input.clientCodeName == clientCodeName) {
                outputStream << indents << "sendDataFieldToServer(" << addQuatation(input.meshName)
                        << ", " << addQuatation(input.dataFieldName) << ");" << endl;
            }
        } else if (settingConnection.inputs[i].type == EMPIRE_ConnectionIO_Signal) {
            const structSignalRef &input = settingConnection.inputs[i].signalRef;
            if (input.clientCodeName == clientCodeName) {
                outputStream << indents << "sendSignalToServer(" << addQuatation(input.signalName)
                        << ");" << endl;
            }
        } else {
            assert(false);
        }
    }
    for (int i = 0; i < settingConnection.outputs.size(); i++) {
        if (settingConnection.outputs[i].type == EMPIRE_ConnectionIO_DataField) {
            const structDataFieldRef &output = settingConnection.outputs[i].dataFieldRef;
            if (output.clientCodeName == clientCodeName) {
                outputStream << indents << "receiveDataFieldFromServer("
                        << addQuatation(output.meshName) << ", "
                        << addQuatation(output.dataFieldName) << ");" << endl;
            }
        } else if (settingConnection.outputs[i].type == EMPIRE_ConnectionIO_Signal) {
            const structSignalRef &output = settingConnection.outputs[i].signalRef;
            if (output.clientCodeName == clientCodeName) {
                outputStream << indents << "receiveSignalFromServer("
                        << addQuatation(output.signalName) << ");" << endl;
            }
        } else {
            assert(false);
        }
    }
}

bool PseudoCodeOutput::isClientInCouplingLogic(std::string clientCodeName,
        structCouplingLogic &settingCouplingLogic) {
    // TODO: algorithm is a little bit complex, to be fully tested!!!!
    if (settingCouplingLogic.type == EMPIRE_connection) {
        string connectionName = settingCouplingLogic.connectionRef.connectionName;
        structConnection &settingConnection = getStructConnectionByName(connectionName);
        for (int i = 0; i < settingConnection.inputs.size(); i++) {
            if (settingConnection.inputs[i].type == EMPIRE_ConnectionIO_DataField) {
                const structDataFieldRef &input = settingConnection.inputs[i].dataFieldRef;
                if (input.clientCodeName == clientCodeName) {
                    return true;
                }
            } else if (settingConnection.inputs[i].type == EMPIRE_ConnectionIO_Signal) {
                const structSignalRef &input = settingConnection.inputs[i].signalRef;
                if (input.clientCodeName == clientCodeName) {
                    return true;
                }
            } else {
                assert(false);
            }
        }
        for (int i = 0; i < settingConnection.outputs.size(); i++) {
            if (settingConnection.outputs[i].type == EMPIRE_ConnectionIO_DataField) {
                const structDataFieldRef &output = settingConnection.outputs[i].dataFieldRef;
                if (output.clientCodeName == clientCodeName) {
                    return true;
                }
            } else if (settingConnection.outputs[i].type == EMPIRE_ConnectionIO_Signal) {
                const structSignalRef &output = settingConnection.outputs[i].signalRef;
                if (output.clientCodeName == clientCodeName) {
                    return true;
                }
            } else {
                assert(false);
            }
        }
        return false;
    } else {
        bool isIn = false;
        for (int i = 0; i < settingCouplingLogic.sequence.size(); i++) { // true if one is true
            isIn = isIn
                    || isClientInCouplingLogic(clientCodeName, settingCouplingLogic.sequence[i]);
        }
        return isIn;
    }
}

structConnection &PseudoCodeOutput::getStructConnectionByName(std::string connectionName) {
    std::vector<structConnection> &settingConnectionVec = metaDatabase->settingConnectionVec;
    for (int i = 0; i < settingConnectionVec.size(); i++) {
        if (connectionName == settingConnectionVec[i].name)
            return settingConnectionVec[i];
    }assert(false);
    // it means the connection cannot be found
    return settingConnectionVec[100000000];
}

void PseudoCodeOutput::setIndents(int n) {
    indents = "";
    for (int i = 0; i < n; i++)
        indents.append("\t");
}

void PseudoCodeOutput::incrementIndents() {
    indents.append("\t");
}

void PseudoCodeOutput::decrementIndents() {
    indents.erase(indents.end() - 1);
}

void PseudoCodeOutput::addBlockStart(std::string str) {
    outputStream << indents << "/*" << endl;
    outputStream << indents << " * ================================================" << endl;
    outputStream << indents << " * " << str << endl;
    outputStream << indents << " * ================================================" << endl;
    outputStream << indents << " */" << endl;
}

void PseudoCodeOutput::addEmptyLines(int n) {
    for (int i = 0; i < n; i++)
        outputStream << endl;
}

void PseudoCodeOutput::addDescriptionLine(std::string str) {
    outputStream << indents << "// " << str << " ..." << endl;
}

string PseudoCodeOutput::addQuatation(string str) {
    string toReturn("\"");
    toReturn.append(str);
    toReturn.append("\"");
    return toReturn;
}

} /* namespace EMPIRE */

