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
#ifndef CONSTANTRELAXATION_H_
#define CONSTANTRELAXATION_H_

#include "AbstractCouplingAlgorithm.h"

namespace EMPIRE {

class ConstantRelaxation: public AbstractCouplingAlgorithm {
public:
    ConstantRelaxation(std::string _name, double _relaxationFactor);
    virtual ~ConstantRelaxation();
    virtual void setInputAndOutput(const ConnectionIO *_input, ConnectionIO *_output);
    void calcNewValue();
private:
    const double RELAXATION_FACTOR;
    double *X_out;
    int SIZE;
    friend class TestConstantRelaxation;
    /// whether output numbers or not
    bool debugMe;
};

} /* namespace EMPIRE */
#endif /* CONSTANTRELAXATION_H_ */
