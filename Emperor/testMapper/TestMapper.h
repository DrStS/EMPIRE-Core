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
#ifndef TESTMAPPER_H_
#define TESTMAPPER_H_

#include <iostream>
#include <string>
#include <assert.h>
#include <vector>
#include <math.h>

#include "ticpp.h"

#include "GiDFileIO.h"
#include "FieldCreator.h"
#include "MortarMapper.h"
#include "AbstractMapper.h"
#include "BarycentricInterpolationMapper.h"
#include "NearestNeighborMapper.h"
#include "NearestElementMapper.h"
#include "DataFieldIntegration.h"

using namespace std;
using namespace ticpp;

namespace EMPIRE {

class TestMapper {
public:
    TestMapper();
    virtual ~TestMapper();
    void parseInputFile(char *fileName);
    void outputSettingToShell();
    void doMapping();
    void doConservationAnalysis();
    struct StructDataField {
        string format;
        string function; // if format=="function"
        string file; // if format=="GiDResult"
        string resultName;
        string analysisName;
        int stepNum;
    };
    struct StructMapper {
        string name;
        string type;
        string meshAFormat;
        string meshFileA;
        string meshBFormat;
        string meshFileB;
        struct StructMortarMapper {
            bool oppositeSurfaceNormal;
            bool dual;
            bool enforceConsistency;
        } settingMortarMapper;
        struct StructConsistentMapping {
            bool doConsistentMapping;
            string dataFieldATypeOfQuantity;
            string dataFieldBTypeOfQuantity;
            StructDataField dataFieldA;
            bool doErrorCalculation;
            StructDataField dataFieldBRef;
        } settingConsistentMapping;
        struct StructConservativeMapping {
            bool doConservativeMapping;
            string dataFieldATypeOfQuantity;
            string dataFieldBTypeOfQuantity;
            StructDataField dataFieldB;
            bool doErrorCalculation;
            StructDataField dataFieldARef;
        } settingConservativeMapping;
    };
    struct StructConservationAnalysis {
        string name;
        string displacementsATypeOfQuantity;
        string forcesATypeOfQuantity;
        string displacementsBTypeOfQuantity;
        string forcesBTypeOfQuantity;
        string meshAFormat;
        string meshFileA;
        string meshBFormat;
        string meshFileB;
        StructDataField displacementsA;
        StructDataField forcesA;
        StructDataField displacementsB;
        StructDataField forcesB;
    };
private:
    vector<StructMapper> settingMappers;
    vector<StructConservationAnalysis> settingConservationAnalysises;
    vector<AbstractMapper*> mappers;
    void parseXMLDataField(ticpp::Element *xmlDataField, StructDataField &structDataField);
    void outputXMLDataFieldToShell(StructDataField &structDataField, string dataFieldName,
            string indent);
    void initDataField(StructDataField &structDataField, int numNodes, int numElems,
            double *nodeCoors, int *nodeIDs, int* numNodesPerElem, int *elemTable, int *elemIDs,
            double *dataField);
};

} /* namespace EMPIRE */
#endif /* TESTMAPPER_H_ */
