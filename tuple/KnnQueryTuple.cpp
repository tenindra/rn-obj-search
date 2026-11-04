// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "KnnQueryTuple.h"

#include "../utility/utility.h"
#include <sstream>

KnnQueryTuple::KnnQueryTuple(std::string networkName, int vertices, int edges, VertexID queryVertexID, std::string queryMethod, 
                             int k, double queryTimeMs, std::string objSetType, double objSetDensity, int objSetVariable, int objSetSize, 
                             std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances): 
    Tuple(networkName,vertices,edges), queryVertexID(queryVertexID), queryMethod(queryMethod), k(k), queryTimeMs(queryTimeMs), 
    objSetType(objSetType), objSetDensity(objSetDensity), objSetVariable(objSetVariable), objSetSize(objSetSize) {

    // Create string to represented kNN result pairs (object,distance)
    this->kNNsDistancesString = utility::implodeVerticesAndDistances(kNNs,kNNDistances);
}

std::string KnnQueryTuple::getTupleString()
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;
    ss << Tuple::getTupleString() << "\t";
    ss << this->queryMethod << "\t";
    ss << this->k << "\t";
    ss << this->queryTimeMs << "\t";
    ss << this->objSetType << "\t";
    ss.precision(5);
    ss << this->objSetDensity << "\t";
    ss << this->objSetVariable << "\t";
    ss << this->objSetSize << "\t";
    ss << this->queryVertexID << "\t";
    ss << this->getAdditionalFields();
    return ss.str();
}

std::string KnnQueryTuple::getMultilineTupleString()
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;
    ss << "\nkNN Query Results: \n";
    ss << Tuple::getMultilineTupleString();
    ss << "Query Method: " << this->queryMethod << ", k: " << this->k << ", Query Vertex: " << this->queryVertexID << "\n"; 
    ss << "Query Time: " << this->queryTimeMs << "ms\n";
    ss.precision(5);
    ss << "Object Set: Type=" << this->objSetType << ", Density=" << this->objSetDensity << ", Variable=" << this->objSetVariable << ", Size=" << this->objSetSize << "\n";   
    ss << "kNNs and Distances: " + this->kNNsDistancesString + "\n";
    return ss.str();
}
