// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _ILBR_H
#define _ILBR_H

#include "Graph.h"
#include "Gtree.h"
#include "ALT.h"
#include "NVDQuadtree.h"
#include "pruned_highway_labeling.h"
#include "SULTree.h"
#include "../common.h"
#include "../utility/Statistics.h"
#include "../utility/AggregateFunction.h"

#include <vector>

class ILBR {

    public:
        // kNN
        void getKNNsByPHLNaive(ALT& alt, std::vector<VertexID>& objects, unsigned int k, VertexID queryVertexID, 
                               std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances, PrunedHighwayLabeling& phl);
        void getKNNsByPHLImprovedByQueryLandmark(ALT& alt, ObjectList& ol, unsigned int k, VertexID queryVertexID, std::vector<VertexID>& kNNs, 
                                                 std::vector<EdgeWeight>& kNNDistances, PrunedHighwayLabeling& phl);
        void getKNNByNVDQuadtreePHL(ALT& alt, NVDQuadtree& nvd, unsigned int k, VertexID queryVertexID, 
                                    std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                                    PrunedHighwayLabeling& phl, Graph& graph);
        void getKNNByCOLTPHL(ALT& alt, COLTree& colt, SULTree& sltree, unsigned int k, VertexID queryVertexID, std::vector<VertexID>& kNNs, 
                             std::vector<EdgeWeight>& kNNDistances, PrunedHighwayLabeling& phl);

        // AkNN
        void getAkNNByNVDQuadtreePHLByConcExps(ALT& alt, NVDQuadtree& nvd, unsigned int k, std::vector<VertexID>& queryVertices, 
                                     std::vector<VertexID>& AkNNs, std::vector<EdgeWeight>& AkNNDistances, 
                                     PrunedHighwayLabeling& phl, Graph& graph, AggregateFunction<EdgeWeight>& agg);
        void getAkNNByCOLTPHL(ALT& alt, COLTree& colt, SULTree& sltree, unsigned int k, std::vector<VertexID>& queryVertices, 
                              std::vector<VertexID>& AkNNs, std::vector<EdgeWeight>& AkNNDistances, PrunedHighwayLabeling& phl, 
                              AggregateFunction<EdgeWeight>& agg);
        
        // Range Queries
        void rangeSearchByCOLTPHL(ALT& alt, COLTree& colt, SULTree& sltree, EdgeWeight range, VertexID queryVertexID, std::vector<VertexID>& results, 
                                  std::vector<EdgeWeight>& resultDistances, PrunedHighwayLabeling& phl);
        void rangeSearchyNVDQuadtreePHL(ALT& alt, NVDQuadtree& nvd, EdgeWeight range, VertexID queryVertexID, 
                                        std::vector<VertexID>& results, std::vector<EdgeWeight>& resultDistances, 
                                        PrunedHighwayLabeling& phl, Graph& graph);
        Statistics stats;

};

#endif // _ILBR_H
