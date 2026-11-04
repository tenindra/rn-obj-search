// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _OBJECTINDEXTUPLE_H
#define _OBJECTINDEXTUPLE_H

#include "Tuple.h"

#include <string>

class ObjectIndexTuple: public Tuple {
    
public:
    std::string getTupleString();
    std::string getMultilineTupleString();
    ObjectIndexTuple(std::string networkName, int vertices, int edges, std::string indexMethod, 
                     double constructionTimeMs, double memoryUsageKB, std::string objSetType, 
                     double objSetDensity, int objSetVariable, int objSetSize);
    
protected:
    std::string indexMethod;
    double constructionTimeMs;
    double memoryUsageKB;
    std::string objSetType;
    double objSetDensity;
    int objSetVariable;
    int objSetSize;
};

#endif // _OBJECTINDEXTUPLE_H
