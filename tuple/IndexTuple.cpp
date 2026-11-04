// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "IndexTuple.h"

#include <sstream>

IndexTuple::IndexTuple(std::string networkName, int vertices, int edges, std::string indexMethod, double constructionTimeMs, double memoryUsageMB): 
    Tuple(networkName,vertices,edges), indexMethod(indexMethod), constructionTimeMs(constructionTimeMs), memoryUsageMB(memoryUsageMB) {}
    
std::string IndexTuple::getTupleString()
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;
    ss << Tuple::getTupleString() << "\t";
    ss << this->indexMethod << "\t";
    ss << this->constructionTimeMs << "\t";
    ss << this->memoryUsageMB << "\t";
    ss << this->getAdditionalFields();
    return ss.str();
}

std::string IndexTuple::getMultilineTupleString()
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;
    ss << "\nIndex Construction Results: \n";
    ss << Tuple::getMultilineTupleString();
    ss << "Index Method: " << this->indexMethod << "\n";
    ss << "Construction Time: " << this->constructionTimeMs << "ms\n";
    ss << "Memory Usage: " << this->memoryUsageMB << "MB\n";
    return ss.str();
}
