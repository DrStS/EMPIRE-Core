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
 * \file MetaDataStructures.h
 * This file holds the the data structure of setting (metadata)
 * \author Tianyang Wang
 * \date 6/5/2012
 **************************************************************************************************/

#ifndef METADATASTRUCTURES_H_
#define METADATASTRUCTURES_H_

#include "EMPEROR_Enum.h"

namespace EMPIRE {
/*
 * The structs defined here should correspond to the counterpart in the XML inputDataField file.
 * Please mark the difference: sometimes a struct owns another struct, sometimes it refers to another struct.
 * There is such a difference also in the XML inputDataField file.
 */
struct structMeshRef {
    std::string clientCodeName;
    std::string meshName;
};

struct structDataFieldRef {
    std::string clientCodeName;
    std::string meshName;
    std::string dataFieldName;
};

struct structSignalRef {
    std::string clientCodeName;
    std::string signalName;
};

struct structConnectionIO {
    EMPIRE_ConnectionIO_Type type;
    structSignalRef signalRef;
    structDataFieldRef dataFieldRef;
};

struct structClientCode {
    struct structMesh {
        struct structDataField {
            std::string name;
            EMPIRE_DataField_location location;
            EMPIRE_DataField_dimension dimension;
            EMPIRE_DataField_typeOfQuantity typeOfQuantity;
        };
        std::string name;
        EMPIRE_Mesh_type type;
        std::vector<structDataField> dataFields;
    };
    struct structSignal {
        std::string name;
        int size3D[3];
    };
    std::string name;
    std::vector<structMesh> meshes;
    std::vector<structSignal> signals;
};

struct structDataOutput {
    std::string name;
    int interval;
    std::vector<structDataFieldRef> dataFieldRefs;
    std::vector<structSignalRef> signalRefs;
};

struct structMapper {
    struct structMortarMapper {
        bool oppositeSurfaceNormal;
        bool dual;
        bool enforceConsistency;
    };
    struct structIGAMortarMapper{
        double tolProjectionDistance;
        int numGPsTriangle;
        int numGPsQuad;
    };
    std::string name;
    structMeshRef meshRefA;
    structMeshRef meshRefB;
    EMPIRE_Mapper_type type;
    structMortarMapper mortarMapper;
    structIGAMortarMapper igaMortarMapper;
};

struct structCouplingAlgorithm {
    struct structAitken {
        double initialAitkenFactor;
    };
    struct structConstantRelaxation {
        double relaxationFactor;
    };
    std::string name;
    EMPIRE_CouplingAlgorithm_type type;
    structAitken aitken;
    structConstantRelaxation constantRelaxation;
};

struct structFilter {
    struct structMappingFilter {
        std::string mapperName;
    };
    struct structExtrapolatingFilter {
        std::string extrapolatorName;
    };
    struct structCouplingAlgorithmFilter {
        std::string couplingAlgorithmName;
    };
    struct structScalingFilter {
        double factor;
    };
    EMPIRE_DataFieldFilter_type type;
    structMappingFilter mappingFilter;
    structCouplingAlgorithmFilter couplingAlgorithmFilter;
    structExtrapolatingFilter extrapolatingFilter;
    structScalingFilter scalingFilter;
    std::vector<structConnectionIO> inputs;
    std::vector<structConnectionIO> outputs;
};

struct structExtrapolator {
  struct structGenMSExtrapolator {
    int numInput;
    int seqLen;
    bool sumOutput;
    double deltaTime;
    std::vector<double> coefficientDot0;
    std::vector<double> coefficientDot1;
    std::vector<double> coefficientDot2;
    std::vector<double> coefficientOut;
  };
  std::string name;
  EMPIRE_Extrapolator_type type;
  structGenMSExtrapolator genMSExtrapolator;
};

struct structConnection {
    std::string name;
    std::vector<structFilter> filterSequence;
    std::vector<structConnectionIO> inputs;
    std::vector<structConnectionIO> outputs;
};

struct structCouplingLogic {
    struct structConnectionRef {
        std::string connectionName;
    };
    struct structTimeStepLoop {
        int numTimeSteps;
        std::vector<std::string> extrapolatorRefs;
        std::vector<std::string> dataOutputRefs;
    };
    struct structIterativeCouplingLoop {
        struct structConvergenceChecker {
            EMPIRE_ConvergenceChecker_whichRef whichRef;
            structDataFieldRef dataFieldRef;
            structSignalRef signalRef;
            std::string couplingAlgorithmRef;
            double absoluteTolerance;
            double relativeTolerance;
            double maxNumOfIterations;
            bool hasAbsTol;
            bool hasRelTol;
            bool hasMaxNumOfIters;
        };
        structConvergenceChecker convergenceChecker;
        std::vector<std::string> convergenceObservers;
        std::vector<std::string> couplingAlgorithmRefs;
        std::vector<std::string> dataOutputRefs;
    };
    EMPIRE_CouplingLogic_type type;
    std::vector<structCouplingLogic> sequence;
    structConnectionRef connectionRef;
    structTimeStepLoop timeStepLoop;
    structIterativeCouplingLoop iterativeCouplingLoop;
};

} /* namespace EMPIRE */
#endif /* METADATASTRUCTURES_H_ */
