// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _KNNQUERYTUPLE_H
#define _KNNQUERYTUPLE_H

#include "Tuple.h"
#include "../common.h"

#include <string>
#include <vector>

class KnnQueryTuple: public Tuple {
    
public:
    std::string getTupleString();
    std::string getMultilineTupleString();
    KnnQueryTuple(std::string networkName, int vertices, int edges, VertexID queryVertexID, std::string queryMethod, int k,
                  double queryTimeMs, std::string objSetType, double objSetDensity, int objSetVariable, int objSetSize, 
                  std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances);
    
private:
    VertexID queryVertexID;
    std::string queryMethod;
    int k;
    double queryTimeMs;
    std::string objSetType;
    double objSetDensity;
    int objSetVariable, objSetSize;
    std::string kNNsDistancesString;
};

#endif // _KNNQUERYTUPLE_H
