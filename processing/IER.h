// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _IER_H
#define _IER_H

#include "Graph.h"
#include "Gtree.h"
#include "pruned_highway_labeling.h"
#include "PackedRtree.h"
#include "../common.h"
#include "../utility/Statistics.h"
#include "../utility/AggregateFunction.h"

#include <vector>

class IER {

    public:
        // kNN
        void getKNNsByPHL(PackedRtree& rtree, unsigned int k, VertexID queryVertexID, 
                          std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                          Graph& graph, PrunedHighwayLabeling& phl);
        void getKNNsByPHLTravelTimes(PackedRtree& rtree, unsigned int k, VertexID queryVertexID, 
                                     std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                                     Graph& graph, PrunedHighwayLabeling& phl);

        // AkNN
        void getAkNNsByPHL(PackedRtree& rtree, unsigned int k, std::vector<VertexID>& queryVertices, 
                           std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                           Graph& graph, PrunedHighwayLabeling& phl, AggregateFunction<EdgeWeight>& agg, 
                           AggregateFunction<EuclideanDistance>& euclidAgg);
        void getAkNNsByPHLTravelTimes(PackedRtree& rtree, unsigned int k, std::vector<VertexID>& queryVertices, 
                                      std::vector<VertexID>& AkNNs, std::vector<EdgeWeight>& AkNNDistances, 
                                      Graph& graph, PrunedHighwayLabeling& phl, AggregateFunction<EdgeWeight>& agg, 
                                      AggregateFunction<EuclideanDistance>& euclidAgg);

        // Range
        void rangeSearchByPHL(PackedRtree& rtree, EdgeWeight range, VertexID queryVertexID, 
                              std::vector<VertexID>& results, std::vector<EdgeWeight>& resultDistances, 
                              Graph& graph, PrunedHighwayLabeling& phl);
        void rangeSearchByPHLTravelTimes(PackedRtree& rtree, EdgeWeight range, VertexID queryVertexID, 
                                         std::vector<VertexID>& results, std::vector<EdgeWeight>& resultDistances, 
                                         Graph& graph, PrunedHighwayLabeling& phl);

        Statistics stats;

};

#endif // _IER_H
