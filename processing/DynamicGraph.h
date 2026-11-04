// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _DYNAMICGRAPH_H
#define _DYNAMICGRAPH_H

#include "../common.h"
#include "Graph.h"

#include <vector>

/*
 * Updateable graph built from input Graph object when things like graph reduction is needed
 */

class Vertex {
    public:
        Vertex(VertexID id, Coordinate x, Coordinate y): id(id), x(x), y(y) {}
        Vertex() {}
        bool edgeExists(VertexID neighbour);
        
        VertexID id;
        Coordinate x, y;
        std::vector<VertexID> adjVertices;
        std::vector<EdgeWeight> adjVertexWgts;
};

class DynamicGraph {
    
    public:
        DynamicGraph() {}
        DynamicGraph(Graph &obj);
        void insertVertex(VertexID vertexID, Coordinate x, Coordinate y);
        void insertEdge (VertexID sourceID, VertexID targetID, EdgeWeight edgeWeight);
        void insertImaginaryNonInvertibleEdge (VertexID sourceID, VertexID targetID, EdgeWeight edgeWeight);
        int getNumVertices();

        std::vector<Vertex> vertices;
         
    private:
        friend class boost::serialization::access;
        
        int numVertices;
};

#endif // _DYNAMICGRAPH_H
