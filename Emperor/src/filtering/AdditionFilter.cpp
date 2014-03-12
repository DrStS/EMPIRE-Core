#include "AdditionFilter.h"
#include "DataField.h"
#include "Signal.h"
#include "ConnectionIO.h"
#include "Message.h"
#include <assert.h>

namespace EMPIRE {

AdditionFilter::AdditionFilter(double _a, double _b) :
        a(_a), b(_b) {
}

AdditionFilter::~AdditionFilter() {
}

void AdditionFilter::init() {
    assert(inputVec.size() == 2);
    assert(outputVec.size() == 1);
    EMPIRE_ConnectionIO_Type type = outputVec[0]->type;
    assert(inputVec[0]->type == type);
    assert(inputVec[1]->type == type);
}

void AdditionFilter::filtering() {
    double *x, *y, *z;
    int n = -1;
    EMPIRE_ConnectionIO_Type type = outputVec[0]->type;
    if (type == EMPIRE_ConnectionIO_DataField) {
        x = inputVec[0]->dataField->data;
        y = inputVec[1]->dataField->data;
        z = outputVec[0]->dataField->data;

        n = outputVec[0]->dataField->dimension * outputVec[0]->dataField->numLocations;
        assert(inputVec[0]->dataField->dimension * inputVec[0]->dataField->numLocations == n);
        assert(inputVec[1]->dataField->dimension * inputVec[1]->dataField->numLocations == n);
    } else if (type == EMPIRE_ConnectionIO_Signal) {
        x = inputVec[0]->signal->array;
        y = inputVec[1]->signal->array;
        z = outputVec[0]->signal->array;

        n = outputVec[0]->signal->size;
        assert(inputVec[0]->signal->size == n);
        assert(inputVec[1]->signal->size == n);
    } else {
        assert(false);
    }

    // do z = a*x + b*y, use mathlibrary in the future
    for (int i = 0; i < n; i++) {
        z[i] = a * x[i] + b * y[i];
    }
}

} /* namespace EMPIRE */
