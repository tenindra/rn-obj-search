// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _GRAPHPARTITIONER_H
#define _GRAPHPARTITIONER_H

#include "Graph.h"
#include "../common.h"

#include <metis.h>
#include <string>
#include <unordered_set>
#include <vector>

// Partitioning method used to split a subgraph into numParts parts (HierarchicalPartition splits its tree nodes, Gtree its
// nodes and SetGenerator the whole graph). A subgraph is given as a range of a vertex order: its vertices are
// reorderedToOriginal[start,start+numVertices) and a vertex u belongs to it iff originalToReordered[u] lies in that range.
// On return parts[i] is the part (0 to numParts-1) of the range vertex reorderedToOriginal[start+i].
class GraphPartitioner {

    public:
        virtual ~GraphPartitioner() {}
        virtual std::string getName() = 0;
        virtual void partition(Graph& graph, const std::vector<VertexID>& reorderedToOriginal, VertexID start, int numVertices,
                               const std::vector<VertexID>& originalToReordered, int numParts, std::vector<int>& parts) = 0;
        // Subgraph given as a set of vertices: the set is copied into vertexOrder (in its iteration order) and parts[i] is the part of vertexOrder[i]
        void partition(Graph& graph, const std::unordered_set<VertexID>& vertices, int numParts, std::vector<VertexID>& vertexOrder, std::vector<int>& parts);

    private:
        std::vector<VertexID> positions; // Position of each vertex of the set in vertexOrder (UNUSED_VERTEX_ID otherwise), re-used between calls

};

// METIS k-way partitioning (the only method so far)
class METISGraphPartitioner: public GraphPartitioner {

    public:
        METISGraphPartitioner(int maxNumVertices, int maxNumEdges);
        ~METISGraphPartitioner();
        METISGraphPartitioner(const METISGraphPartitioner&) = delete;
        METISGraphPartitioner& operator=(const METISGraphPartitioner&) = delete;
        using GraphPartitioner::partition;
        std::string getName();
        void partition(Graph& graph, const std::vector<VertexID>& reorderedToOriginal, VertexID start, int numVertices,
                       const std::vector<VertexID>& originalToReordered, int numParts, std::vector<int>& parts);

    private:
        void populateMETISArrays(Graph& graph, const std::vector<VertexID>& reorderedToOriginal, VertexID rangeStart, idx_t numVertices,
                                 const std::vector<VertexID>& originalToReordered);

        // METIS options and arrays in Compressed Storage Format (sized for the whole graph and re-used for every subgraph)
        idx_t options[METIS_NOPTIONS];
        idx_t ncons;
        idx_t *xadj;
        idx_t *adjncy;
        idx_t *parts;

};

#endif // _GRAPHPARTITIONER_H
