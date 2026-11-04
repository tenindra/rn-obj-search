// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "GraphPartitioner.h"

#include <cstdlib>
#include <iostream>

void GraphPartitioner::partition(Graph& graph, const std::unordered_set<VertexID>& vertices, int numParts, std::vector<VertexID>& vertexOrder, std::vector<int>& parts)
{
    // The set becomes the range [0,|vertices|) of its own iteration order, with positions marked only for its vertices
    if (this->positions.size() != static_cast<std::size_t>(graph.getNumVertices())) {
        this->positions.assign(graph.getNumVertices(),constants::UNUSED_VERTEX_ID);
    }
    vertexOrder.assign(vertices.begin(),vertices.end());
    for (std::size_t i = 0; i < vertexOrder.size(); ++i) {
        this->positions[vertexOrder[i]] = i;
    }
    this->partition(graph,vertexOrder,0,vertexOrder.size(),this->positions,numParts,parts);
    for (std::size_t i = 0; i < vertexOrder.size(); ++i) {
        this->positions[vertexOrder[i]] = constants::UNUSED_VERTEX_ID; // Reset for the next call
    }
}

METISGraphPartitioner::METISGraphPartitioner(int maxNumVertices, int maxNumEdges)
{
    // Configure METIS Partitioning Options
    METIS_SetDefaultOptions(this->options);
    this->options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // Partition graph k-ways (this will be contiguous)
    this->options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; // Minimise actual number of borders (slower than cut)
    this->options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
    this->options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM;
    this->options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
    this->options[METIS_OPTION_UFACTOR] = 500;
    this->options[METIS_OPTION_CONTIG] = 1;
    this->options[METIS_OPTION_NUMBERING] = 0;
    this->ncons = 1; // Number of balancing constraints (at least 1)

    this->xadj = new idx_t[maxNumVertices+1];
    this->adjncy = new idx_t[maxNumEdges];
    this->parts = new idx_t[maxNumVertices];
}

METISGraphPartitioner::~METISGraphPartitioner()
{
    delete[] xadj;
    delete[] adjncy;
    delete[] parts;
}

std::string METISGraphPartitioner::getName()
{
    return "metis";
}

void METISGraphPartitioner::partition(Graph& graph, const std::vector<VertexID>& reorderedToOriginal, VertexID start, int numVertices,
                                      const std::vector<VertexID>& originalToReordered, int numParts, std::vector<int>& parts)
{
    // Assuming original graph is undirected graph (if not it must made undirected
    // before calling this program in order for METIS to function)
    idx_t nvtxs = numVertices, nparts = numParts, objVal;
    this->populateMETISArrays(graph,reorderedToOriginal,start,nvtxs,originalToReordered);

    // Unit edge weights (no adjwgt) and no vertex weights or sizes
    int status = METIS_PartGraphKway(&nvtxs,&this->ncons,this->xadj,this->adjncy,NULL,NULL,NULL,&nparts,NULL,NULL,this->options,&objVal,this->parts);
    if (status != METIS_OK) {
        std::cerr << "METIS failed to partition " << numVertices << " vertices into " << numParts << " parts (status " << status << ")" << std::endl;
        std::exit(1);
    }

    parts.assign(this->parts,this->parts+nvtxs);
}

void METISGraphPartitioner::populateMETISArrays(Graph& graph, const std::vector<VertexID>& reorderedToOriginal, VertexID rangeStart, idx_t numVertices,
                                                const std::vector<VertexID>& originalToReordered)
{
    /*
     * Convert the range to Compressed Storage Format (required by METIS) using the vertex order as the METIS numbering:
     * METIS index i is reorderedToOriginal[rangeStart+i] and a neighbour u is in the subgraph iff originalToReordered[u]
     * is in [rangeStart,rangeStart+numVertices), its METIS index being that position minus rangeStart
     */
    int adjncy_pos = 0;
    int adjListStart, nextAdjListStart;
    VertexID rangeEnd = rangeStart + numVertices; // Exclusive

    for (idx_t i = 0; i < numVertices; ++i) {
        this->xadj[i] = adjncy_pos;
        VertexID vertex = reorderedToOriginal[rangeStart+i];
        adjListStart = graph.getEdgeListStartIndex(vertex);
        nextAdjListStart = graph.getEdgeListSize(vertex);
        for (int j = adjListStart; j < nextAdjListStart; ++j) {
            VertexID neighbourPos = originalToReordered[graph.edges[j].first];
            if (neighbourPos >= rangeStart && neighbourPos < rangeEnd) {
                // We only consider neighbours that are in the subgraph we are partitiong
                this->adjncy[adjncy_pos] = neighbourPos - rangeStart;
                ++adjncy_pos;
            }
        }
    }

    // Must set final position of xadj as end of adjncy
    // I.e. the index after the end of neighbours for the last vertex
    this->xadj[numVertices] = adjncy_pos;
}
