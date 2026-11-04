// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "ObjectIndexTuple.h"

#include <sstream>

ObjectIndexTuple::ObjectIndexTuple(std::string networkName, int vertices, int edges, std::string indexMethod, 
                             double constructionTimeMs, double memoryUsageKB, std::string objSetType,
                             double objSetDensity,int objSetVariable, int objSetSize): 
    Tuple(networkName,vertices,edges), indexMethod(indexMethod), constructionTimeMs(constructionTimeMs), 
    memoryUsageKB(memoryUsageKB), objSetType(objSetType), objSetDensity(objSetDensity), objSetVariable(objSetVariable), objSetSize(objSetSize) {}
    
std::string ObjectIndexTuple::getTupleString()
{
    std::stringstream ss;
    ss.precision(4);
    ss << std::fixed;
    ss << Tuple::getTupleString() << "\t";
    ss << this->indexMethod << "\t";
    ss << this->constructionTimeMs << "\t";
    ss << this->memoryUsageKB << "\t";
    ss << this->objSetType << "\t";
    ss.precision(5);
    ss << this->objSetDensity << "\t";
    ss << this->objSetVariable << "\t";
    ss << this->objSetSize << "\t";
    ss << this->getAdditionalFields();
    return ss.str();
}

std::string ObjectIndexTuple::getMultilineTupleString()
{
    std::stringstream ss;
    ss.precision(4);
    ss << std::fixed;
    ss << "\nObject Index Construction Results: \n";
    ss << Tuple::getMultilineTupleString();
    ss << "Object Index Method: " << this->indexMethod << "\n";
    ss << "Construction Time: " << this->constructionTimeMs << "ms\n";
    ss << "Memory Usage: " << this->memoryUsageKB << "KB\n";
    ss.precision(5);
    ss << "Object Set: Type=" << this->objSetType << ", Density=" << this->objSetDensity << ", Variable=" << this->objSetVariable << ", Size=" << this->objSetSize << "\n";
    return ss.str();
}
