// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _DIJKSTRASEARCH_H
#define _DIJKSTRASEARCH_H

#include "../common.h"
#include "DynamicGraph.h"
#include "Graph.h"
#include "../queue/BinaryMinHeap.h"

#include <vector>
#include <unordered_set>
#include <unordered_map>

class DijkstraSearch {

    public:

        // Single Source Shortest Paths
        void findSSSPDistances(Graph& graph, VertexID source, std::vector<EdgeWeight>& targetDistances);
        void findSSSPDistances(Graph& graph, VertexID source, std::vector<EdgeWeight>& targetDistances, 
                               BinaryMinHeap<EdgeWeight,VertexID> *pqueue);

        // Single Source Multiple Destination Shortest Paths
        void findSSMTDistances(DynamicGraph& graph, VertexID source,
                               std::unordered_set<VertexID>& targetSet, 
                               std::unordered_map<VertexID,EdgeWeight>& results,
                               BinaryMinHeap<EdgeWeight,VertexID> *pqueue);

        void findSSMTDistances(Graph& graph, VertexID source,
                               std::unordered_set<VertexID>& targetSet, 
                               std::unordered_map<VertexID,EdgeWeight>& results,
                               BinaryMinHeap<EdgeWeight,VertexID> *pqueue);

        int findSSMTDistances(Graph& graph, VertexID source, std::unordered_set<VertexID>& subgraphVertices, 
                              std::unordered_map<VertexID,EdgeWeight>& results);

        // Subgraph-Restricted Shortest Path
        EdgeWeight findShortestPathDistanceSubgraph(Graph& graph, VertexID source, VertexID target, std::vector<bool>& edgeInSubgraph);

        // Functions to Support ALT Landmark Selection
        VertexID findFarthestVertex(Graph& graph, VertexID source, EdgeWeight& farthestDist);
        VertexID findFarthestVertexByBFS(Graph& graph, VertexID source);

        // Graph Property Computation
        EdgeWeight getLowerboundGraphDiameterByDoubleSweep(Graph& graph, int numSweeps);

    private:
        static const unsigned int DOUBLE_SWEEP_SEED = 0; // Seed of the double sweeps' start vertices

};

#endif // _DIJKSTRASEARCH_H
