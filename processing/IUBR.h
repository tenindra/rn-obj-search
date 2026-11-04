// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _IUBR_H
#define _IUBR_H

#include "Graph.h"
#include "ALT.h"
#include "pruned_highway_labeling.h"
#include "SULTree.h"
#include "PackedRtree.h"
#include "../common.h"
#include "../utility/Statistics.h"

#include <vector>

class IUBR {

    public:
        void getFkNsByPHLNaive(ALT& alt, std::vector<VertexID>& objects, unsigned int k, VertexID queryVertexID, 
                               std::vector<VertexID>& FkNs, std::vector<EdgeWeight>& FkNDistances, PrunedHighwayLabeling& phl);
        void getFkNsByCOLTPHL(ALT& alt, COLTree& colt, SULTree& sltree, unsigned int k, VertexID queryVertexID, std::vector<VertexID>& FkNs, 
                              std::vector<EdgeWeight>& FkNDistances, PrunedHighwayLabeling& phl);
        
        Statistics stats;

};

#endif // _IUBR_H
