// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _INE_H
#define _INE_H

#include "Graph.h"
#include "../common.h"
#include "../utility/Statistics.h"
#include "../utility/AggregateFunction.h"

#include <vector>

class INE {

    public:
        void getKNNs(Graph& graph, std::vector<bool>& vertexObjectStatus, unsigned int k, VertexID queryVertexID, std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances);
        void getAkNNs(Graph& graph, std::vector<bool>& vertexObjectStatus, unsigned int k, std::vector<VertexID>& queryVertices, 
                      std::vector<VertexID>& AkNNs, std::vector<EdgeWeight>& AkNNDistances, AggregateFunction<EdgeWeight>& agg);
        void getFkNs(Graph& graph, std::vector<bool>& vertexObjectStatus, int numObjects, int k, VertexID queryVertexID, 
                     std::vector<VertexID>& FkNs, std::vector<EdgeWeight>& FkNDistances);
        void rangeSearch(Graph& graph, std::vector<bool>& vertexObjectStatus, int numObjects, EdgeWeight range, VertexID queryVertexID, std::vector<VertexID>& results, std::vector<EdgeWeight>& resultDistances);
        std::vector<bool> buildVertexObjectStatusVector(Graph& graph, std::vector<VertexID>& objectVertices);

};

#endif // _INE_H
