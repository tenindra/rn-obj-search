// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "DynamicGraph.h"
#include "../utility/utility.h"

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>

// Create an updateable graph from non-updateable graph
DynamicGraph::DynamicGraph(Graph& obj)
{
    this->numVertices = obj.getNumVertices();
    for (std::size_t i = 0; i < obj.coordinates.size(); ++i) {
        this->insertVertex(i,obj.coordinates[i].first,obj.coordinates[i].second);
    }

    int adjListStart, adjListSize;
    VertexID numVertices = static_cast<VertexID>(this->numVertices);
    for (VertexID i = 0; i < numVertices; ++i) {
        adjListStart = obj.getEdgeListStartIndex(i);
        adjListSize = obj.getEdgeListSize(i);
        
        for (int j = adjListStart; j < adjListSize; ++j) {
            this->insertEdge(i,obj.edges[j].first,obj.edges[j].second);
        }        
    }    
}

void DynamicGraph::insertVertex(VertexID vertexID, Coordinate x, Coordinate y)
{
    // The checks here and in insertEdge are debug-only, as the G-tree builds this graph from a Graph, whose parser checks the
    // same vertices and edges (a graph read by parseGraphFile would need explicit checks, like Graph's)
    RN_DEBUG_ASSERT(vertexID == this->vertices.size() && "Vertex not inserted in order");
    
    Vertex newVertex(vertexID,x,y);
    this->vertices.push_back(newVertex);
}

void DynamicGraph::insertEdge(VertexID sourceID, VertexID targetID, EdgeWeight edgeWeight) {
    // Note: We assume edgeWeight is valid (e.g., against the minimum Euclidean distance)
    RN_DEBUG_ASSERT(!this->vertices[sourceID].edgeExists(targetID) && "Parallel edge detected");
    RN_DEBUG_ASSERT(sourceID != targetID && "Self-loop detected");
    RN_DEBUG_ASSERT(edgeWeight != 0 && "Zero-edge weight detected");
    
    this->vertices[sourceID].adjVertices.push_back(targetID);
    this->vertices[sourceID].adjVertexWgts.push_back(edgeWeight);
}

void DynamicGraph::insertImaginaryNonInvertibleEdge(VertexID sourceID, VertexID targetID, EdgeWeight edgeWeight) {
    // This function is for the benefit of G-tree (and any other method that
    // relies on adding new edges to make Dijkstra computation faster). In this
    // we do not need inverted edges.
    this->vertices[sourceID].adjVertices.push_back(targetID);
    this->vertices[sourceID].adjVertexWgts.push_back(edgeWeight);
}

int DynamicGraph::getNumVertices() {
    return this->numVertices;
}

/*
 * DynamicGraph::Vertex
 */

bool Vertex::edgeExists(VertexID neighbour)
{
    for (std::size_t i = 0; i < this->adjVertices.size(); ++i) {
        if (this->adjVertices[i] == neighbour) {
            return true;
        }
    }
    // Edge not found so we return false
    return false;
}
