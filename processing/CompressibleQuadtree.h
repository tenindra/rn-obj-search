// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _COMPRESSIBLEQUADTREE_H
#define _COMPRESSIBLEQUADTREE_H

#include "Graph.h"
#include "MortonCode.h"
#include "../common.h"

#include <vector>
#include <unordered_set>

class CompressibleQuadtreeNode {

    public:
        CompressibleQuadtreeNode(MortonCode& parentCode, Direction dir, CompressibleQuadtreeNode* prevLeaf = NULL, CompressibleQuadtreeNode* nextLeaf = NULL);
        CompressibleQuadtreeNode(bool isRoot, bool isLeaf);
        ~CompressibleQuadtreeNode();
        void setNextLeaf(CompressibleQuadtreeNode* nextLeaf);
        bool isRootNode();
        bool isLeafNode();
        void setLeaf(bool isLeaf);
        int getNumVertices(); // This wil get number > 0 even for non-leaf nodes
        int getNumColours(); // This wil get number > 0 even for non-leaf nodes
        void addVertex(VertexID vertex, VertexID colour);
        void pushDownVertices(std::vector<MortonCode>& allMortonCodes, std::vector<VertexID>& colourMap);
        void pushDownVertex(VertexID vertex, MortonCode& mortonCode, VertexID colour);
        bool contains(const MortonCode& point);
        double computeIndexSizeBytes();
        
        // Child quadtree nodes (if non leaf)
        CompressibleQuadtreeNode *sw, *se, *nw, *ne;
        
        // Doubly-linked list of leaf nodes
        CompressibleQuadtreeNode *prevLeaf, *nextLeaf;
        
        MortonCode mortonCode;
        std::vector<VertexID> vertices; // Only full for leaf nodes
        std::unordered_set<VertexID> colourIDs;

    private:
        bool isRoot;
        bool isLeaf;
        int numVertices;
        int numColours;
};

class CompressibleQuadtree {

    public:
        CompressibleQuadtree(int maxLeafColours);
        ~CompressibleQuadtree();
        void buildFromGraph(Graph& graph, std::vector<MortonCode>& mortonCodes, std::vector<VertexID>& colourMap);
        void insert(VertexID vertex, std::vector<MortonCode>& mortonCodes, std::vector<VertexID>& colourMap);
        void split(CompressibleQuadtreeNode* treeNode, std::vector<MortonCode>& mortonCodes, std::vector<VertexID>& colourMap);
        int countLeafNodes(int& totalVertices, int& nonEmptyNodes);
        double computeIndexSize();
        
        CompressibleQuadtreeNode *rootNode;
        CompressibleQuadtreeNode *firstLeaf, *lastLeaf;
        
    private:
        int maxLeafColours;
};

#endif // _COMPRESSIBLEQUADTREE_H