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
/***********************************************************************************************//**
 * \file DataFieldIntegrationFilter.h
 * This file holds the class DataFieldIntegrationFilter
 * \date 2/7/2013
 **************************************************************************************************/
#ifndef DATAFIELDINTEGRATIONFILTER_H_
#define DATAFIELDINTEGRATIONFILTER_H_

#include "AbstractFilter.h"

namespace EMPIRE {

class AbstractMesh;
class DataFieldIntegration;
/********//**
 * \brief Class DataFieldIntegrationFilter is an operator from traction to force or vice versa
 * \author Tianyang Wang
 ***********/
class DataFieldIntegrationFilter: public AbstractFilter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \author Tianyang Wang
     ***********/
    DataFieldIntegrationFilter(AbstractMesh *mesh);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~DataFieldIntegrationFilter();
    /***********************************************************************************************
     * \brief Copy the input to the output, if not enough input data, fill by 0
     * \author Tianyang Wang
     ***********/
    void filtering();
    /***********************************************************************************************
     * \brief Initialize data according to the inputs and outputs
     * \author Tianyang Wang
     ***********/
    void init();

private:
    /// the mesh
    AbstractMesh *mesh;
    /// whether to do integration or
    bool doIntegration;
    /// the actual one who performs integration
    DataFieldIntegration *dataFieldIntegration;
    /***********************************************************************************************
     * \brief Do integration (massMatrix*input=output).
     * \author Tianyang Wang
     ***********/
    void integrate();
    /***********************************************************************************************
     * \brief Do deintegration (massMatrix^(-1)*input=output).
     * \author Tianyang Wang
     ***********/
    void deIntegrate();
};
} /* namespace EMPIRE */
#endif /* DATAFIELDINTEGRATIONFILTER_H_ */
