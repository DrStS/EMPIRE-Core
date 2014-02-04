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
 * \file EMPEROR_Enum.h
 * This file holds all enumerations (specifying types) of the program
 * \date 9/3/2012
 **************************************************************************************************/
#ifndef EMPEROR_ENUM_H_
#define EMPEROR_ENUM_H_

/********//**
 * The rule of naming 1:
 * EMPIRE_[className]_[variableName] and
 * EMPIRE_[className]_[variableValue].
 * For example, in class "DataField", there is a member variable "dimension", it could
 * have value "scalar" or "vector", the enum below is named according to this rule.
 ***********/
enum EMPIRE_DataField_dimension {
    EMPIRE_DataField_scalar = 1, EMPIRE_DataField_vector = 3
};
enum EMPIRE_DataField_location {
    EMPIRE_DataField_atNode, EMPIRE_DataField_atElemCentroid
};
enum EMPIRE_DataField_typeOfQuantity {
    EMPIRE_DataField_field, EMPIRE_DataField_fieldIntegral
};
enum EMPIRE_ConnectionIO_Type {
    EMPIRE_ConnectionIO_Signal, EMPIRE_ConnectionIO_DataField
};
enum EMPIRE_Signal_dimension {
    EMPIRE_Signal_0D, EMPIRE_Signal_1D, EMPIRE_Signal_2D, EMPIRE_Signal_3D
};

/********//**
 * The rule of naming 2:
 * EMPIRE_[motherClassName]_type and
 * EMPIRE_[childClassName].
 * But sometimes it can be simpler, e.g. use "Mapper" instead of "AbstractMapper",
 * use "MortarMapper" instead of "MortarMapperAdapter"
 ***********/
enum EMPIRE_Mapper_type {
    EMPIRE_IGAMortarMapper, EMPIRE_MortarMapper, EMPIRE_NearestNeighborMapper, EMPIRE_BarycentricInterpolationMapper
};

enum EMPIRE_Mesh_type {
    EMPIRE_Mesh_FEMesh, EMPIRE_Mesh_IGAMesh
};

enum EMPIRE_DataFieldFilter_type {
    EMPIRE_CopyFilter,
    EMPIRE_MappingFilter,
    EMPIRE_LocationFilter,
    EMPIRE_ScalingFilter,
    EMPIRE_SetFilter,
    EMPIRE_DataFieldIntegrationFilter
};

enum EMPIRE_CouplingLogic_type {
    EMPIRE_IterativeCouplingLoop,
    EMPIRE_TimeStepLoop,
    EMPIRE_connection,
    EMPIRE_CouplingLogicSequence,
    EMPIRE_OptimizationLoop
};

enum EMPIRE_Extrapolator_type {
    EMPIRE_LinearExtrapolator
};

enum EMPIRE_CouplingAlgorithm_type {
    EMPIRE_Aitken, EMPIRE_ConstantRelaxation
};

/********//**
 * The rule of naming 3:
 * When rule 1 and rule 2 are not applicable, do whatever you want
 ***********/

#endif /* EMPEROR_ENUM_H_ */
