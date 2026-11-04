// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _AGGREGATEFUNCTION_H
#define _AGGREGATEFUNCTION_H

#include "../common.h"
#include <vector>

template <typename T>
class AggregateFunction {

    public:
        virtual T aggregate(std::vector<T>& distances) = 0;
        virtual T getThreshold(T distance, int numQueryVertices) = 0;
        virtual T getLowerBound(T distance, int numQueryVertices) = 0;
        virtual std::string getPrintString() = 0;
        // Assumes both vectors are sorted on edge-weights
        virtual int getMinObjListIndex(std::vector<VertexEdgeWeightPair>& objectList, std::vector<EdgeWeight>& queryVertexLmkDists) = 0;

};

#endif // _AGGREGATEFUNCTION_H