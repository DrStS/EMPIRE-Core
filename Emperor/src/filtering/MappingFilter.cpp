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
#include "MappingFilter.h"
#include "MapperAdapter.h"
#include "ConnectionIO.h"
#include "DataField.h"
#include "DataFieldIntegrationFilter.h"

#include <assert.h>
#include <iostream>

using namespace std;

namespace EMPIRE {

MappingFilter::MappingFilter(MapperAdapter *_mapper) :
        AbstractFilter(), mapper(_mapper) {
    assert(mapper!=NULL);
    inputModifier = NULL;
    outputModifier = NULL;
    tmpInput = NULL;
    tmpOutput = NULL;
}

MappingFilter::~MappingFilter() {
    if (inputModifier != NULL) {
        assert(tmpInput != NULL);
        delete tmpInput;
        delete inputModifier; // this will also delete the ConnectionIOs or it
    }
    if (outputModifier != NULL) {
        assert(tmpOutput != NULL);
        delete tmpOutput;
        delete outputModifier; // this will also delete the ConnectionIOs or it
    }
}

void MappingFilter::filtering() {
    // 1. input -> tmpInput
    if (inputModifier != NULL) {
        inputModifier->filtering();
    }
    // 2. do mapping
    DataField *inDataField;
    DataField *outDataField;
    if (inputModifier != NULL) {
        inDataField = tmpInput;
    } else {
        inDataField = inputVec[0]->dataField;
    }
    if (outputModifier != NULL) {
        outDataField = tmpOutput;
    } else {
        outDataField = outputVec[0]->dataField;
    }
    if (consistentMapping) {
        mapper->consistentMapping(inDataField, outDataField);
    } else {
        mapper->conservativeMapping(inDataField, outDataField);
    }
    // 3. tmpOutput -> output
    if (outputModifier != NULL) {
        outputModifier->filtering();
    }
}

void MappingFilter::init() {
    assert(inputVec.size() == 1);
    assert(outputVec.size() == 1);
    assert(inputVec[0]->type == outputVec[0]->type);
    assert(inputVec[0]->type == EMPIRE_ConnectionIO_DataField);
    DataField *inDataField = inputVec[0]->dataField;
    DataField *outDataField = outputVec[0]->dataField;
    assert(inDataField->dimension == outDataField->dimension);
    assert(inDataField->location == outDataField->location);
    assert(inDataField->location == EMPIRE_DataField_atNode);
    AbstractMesh *inMesh = inputVec[0]->mesh;
    AbstractMesh *outMesh = outputVec[0]->mesh;
    assert(inMesh != outMesh); // otherwise, cannot determine whether to do consistent mapping or conservative mapping
    // consistent mapping A->B, conservative mapping B->A
    if (mapper->isMeshA(inMesh) && mapper->isMeshB(outMesh))
        consistentMapping = true;
    else if (mapper->isMeshB(inMesh) && mapper->isMeshA(outMesh))
        consistentMapping = false;
    else
        assert(false);

    if (consistentMapping) {
        // mapper does consistent mapping between fields
        if (inDataField->typeOfQuantity == EMPIRE_DataField_fieldIntegral) {
            inputModifier = new DataFieldIntegrationFilter(inMesh);
            tmpInput = new DataField("", EMPIRE_DataField_atNode, inDataField->numLocations,
                    inDataField->dimension, EMPIRE_DataField_field);
            inputModifier->addInput(
                    new ConnectionIO(inputVec[0]->clientCode, inputVec[0]->mesh,
                            inputVec[0]->dataField));
            inputModifier->addOutput(
                    new ConnectionIO(inputVec[0]->clientCode, inputVec[0]->mesh, tmpInput));
            inputModifier->init();
        }
        if (outDataField->typeOfQuantity == EMPIRE_DataField_fieldIntegral) {
            outputModifier = new DataFieldIntegrationFilter(outMesh);
            tmpOutput = new DataField("", EMPIRE_DataField_atNode, outDataField->numLocations,
                    outDataField->dimension, EMPIRE_DataField_field);
            outputModifier->addInput(
                    new ConnectionIO(outputVec[0]->clientCode, outputVec[0]->mesh, tmpOutput));
            outputModifier->addOutput(
                    new ConnectionIO(outputVec[0]->clientCode, outputVec[0]->mesh,
                            outputVec[0]->dataField));
            outputModifier->init();

        }
    } else {
        // mapper does conservative mapping between fieldIntegrals
        if (inDataField->typeOfQuantity == EMPIRE_DataField_field) {
            inputModifier = new DataFieldIntegrationFilter(inMesh);
            tmpInput = new DataField("", EMPIRE_DataField_atNode, inDataField->numLocations,
                    inDataField->dimension, EMPIRE_DataField_fieldIntegral);
            inputModifier->addInput(
                    new ConnectionIO(inputVec[0]->clientCode, inputVec[0]->mesh,
                            inputVec[0]->dataField));
            inputModifier->addOutput(
                    new ConnectionIO(inputVec[0]->clientCode, inputVec[0]->mesh, tmpInput));
            inputModifier->init();
        }
        if (outDataField->typeOfQuantity == EMPIRE_DataField_field) {
            outputModifier = new DataFieldIntegrationFilter(outMesh);
            tmpOutput = new DataField("", EMPIRE_DataField_atNode, outDataField->numLocations,
                    outDataField->dimension, EMPIRE_DataField_fieldIntegral);
            outputModifier->addInput(
                    new ConnectionIO(outputVec[0]->clientCode, outputVec[0]->mesh, tmpOutput));
            outputModifier->addOutput(
                    new ConnectionIO(outputVec[0]->clientCode, outputVec[0]->mesh,
                            outputVec[0]->dataField));
            outputModifier->init();
        }
    }
}

} /* namespace EMPIRE */
