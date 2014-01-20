#include "LinearExtrapolator.h"
#include "Signal.h"
#include "DataField.h"
#include "ConnectionIO.h"
#include "Message.h"
#include <assert.h>

using namespace std;

namespace EMPIRE {

LinearExtrapolator::LinearExtrapolator(std::string _name) :
        AbstractExtrapolator(_name) {
}

LinearExtrapolator::~LinearExtrapolator() {
    for (int i = 0; i < connectionIOs.size(); i++) {
        delete connectionIOs[i];
    }
    for (int i = 0; i < data00.size(); i++) {
        delete[] data00[i];
    }
    for (int i = 0; i < data0.size(); i++) {
        delete[] data0[i];
    }
}

void LinearExtrapolator::init() {
    assert(connectionIOs.size() > 0);
    for (int i = 0; i < connectionIOs.size(); i++) {
        if (connectionIOs[i]->type == EMPIRE_ConnectionIO_DataField) {
            DataField *dataField = connectionIOs[i]->dataField;
            double *tmp00 = new double[dataField->dimension * dataField->numLocations];
            data00.push_back(tmp00);
            double *tmp0 = new double[dataField->dimension * dataField->numLocations];
            data0.push_back(tmp0);
        } else if (connectionIOs[i]->type == EMPIRE_ConnectionIO_Signal) {
            Signal *signal = connectionIOs[i]->signal;
            double *tmp00 = new double[signal->size];
            data00.push_back(tmp00);
            double *tmp0 = new double[signal->size];
            data0.push_back(tmp0);
        } else {
            assert(false);
        }
    }
}

void LinearExtrapolator::extrapolate() {
    HEADING_OUT(4, "LinearExtrapolator", "doing linear extrapolation ...", infoOut);

    assert(connectionIOs.size() == data0.size());
    assert(data00.size() == data0.size());
    currentTimeStepNumber++;
    if (currentTimeStepNumber == 1) {
        return;
    } else if (currentTimeStepNumber == 2) {
        // store the data at time step 1 to data0
        for (int i = 0; i < connectionIOs.size(); i++) {
            if (connectionIOs[i]->type == EMPIRE_ConnectionIO_DataField) {
                DataField *dataField = connectionIOs[i]->dataField;
                for (int j = 0; j < dataField->dimension * dataField->numLocations; j++) {
                    data0[i][j] = dataField->data[j];
                }
            } else if (connectionIOs[i]->type == EMPIRE_ConnectionIO_Signal) {
                Signal *signal = connectionIOs[i]->signal;
                for (int j = 0; j < signal->size; j++) {
                    data0[i][j] = signal->array[j];
                }
            } else {
                assert(false);
            }
        }
    } else {
        // store data0 to data00
        for (int i = 0; i < connectionIOs.size(); i++) {
            if (connectionIOs[i]->type == EMPIRE_ConnectionIO_DataField) {
                DataField *dataField = connectionIOs[i]->dataField;
                for (int j = 0; j < dataField->dimension * dataField->numLocations; j++) {
                    data00[i][j] = data0[i][j];
                }
            } else if (connectionIOs[i]->type == EMPIRE_ConnectionIO_Signal) {
                Signal *signal = connectionIOs[i]->signal;
                for (int j = 0; j < signal->size; j++) {
                    data00[i][j] = data0[i][j];
                }
            } else {
                assert(false);
            }
        }

        // store data of last time step to data0
        for (int i = 0; i < connectionIOs.size(); i++) {
            if (connectionIOs[i]->type == EMPIRE_ConnectionIO_DataField) {
                DataField *dataField = connectionIOs[i]->dataField;
                for (int j = 0; j < dataField->dimension * dataField->numLocations; j++) {
                    data0[i][j] = dataField->data[j];
                }
            } else if (connectionIOs[i]->type == EMPIRE_ConnectionIO_Signal) {
                Signal *signal = connectionIOs[i]->signal;
                for (int j = 0; j < signal->size; j++) {
                    data0[i][j] = signal->array[j];
                }
            } else {
                assert(false);
            }
        }

        // do linear extrapolation with formula "new = 2 x data0 - data00"
        for (int i = 0; i < connectionIOs.size(); i++) {
            if (connectionIOs[i]->type == EMPIRE_ConnectionIO_DataField) {
                DataField *dataField = connectionIOs[i]->dataField;
                for (int j = 0; j < dataField->dimension * dataField->numLocations; j++) {
                    dataField->data[j] = 2.0 * data0[i][j] - data00[i][j];
                }
            } else if (connectionIOs[i]->type == EMPIRE_ConnectionIO_Signal) {
                Signal *signal = connectionIOs[i]->signal;
                for (int j = 0; j < signal->size; j++) {
                    signal->array[j] = 2.0 * data0[i][j] - data00[i][j];;
                }
            } else {
                assert(false);
            }
        }
    }
}

} /* namespace EMPIRE */
