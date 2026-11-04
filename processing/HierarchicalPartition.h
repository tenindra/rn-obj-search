// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _HIERARCHICALPARTITION_H
#define _HIERARCHICALPARTITION_H

#include "Graph.h"
#include "GraphPartitioner.h"
#include "../common.h"

#include <string>
#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

// A node of the hierarchical partition: the subgraph with vertices reorderedToOriginal[startVertexID..endVertexID]
struct HierarchicalPartitionNode {
    int parentIdx;
    int numVertices;
    bool isLeaf;
    VertexID startVertexID; // This is by reordered IDs
    VertexID endVertexID; // This is by reordered IDs (endVertexID is inclusive)
    std::vector<int> childNodeIdxs;
    std::vector<VertexID> borderVertexIDs; // Vertices with a neighbour outside the node (original IDs), none for the root
    Coordinate minX, minY, maxX, maxY; // Bounding box of the node's vertices

    HierarchicalPartitionNode(int parentIdx, int numVertices);
    HierarchicalPartitionNode();

    private:
        friend class boost::serialization::access;

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->parentIdx;
            ar & this->numVertices;
            ar & this->isLeaf;
            ar & this->startVertexID;
            ar & this->endVertexID;
            ar & this->childNodeIdxs;
            ar & this->borderVertexIDs;
            ar & this->minX;
            ar & this->minY;
            ar & this->maxX;
            ar & this->maxY;
        }
};

// Hierarchical partitioning of a road network into a tree of subgraphs, as used by the SUL-Tree: each node is split
// into fanout parts until a node has at most maxLeafSize vertices. It depends only on the graph, the partitioning
// method and these two parameters, so it is built once (lmk_exp -e partitions) and saved, and SUL-Trees with different
// landmark parameters are built from it. Vertices are reordered so that the vertices of every tree node are contiguous:
// reorderedToOriginal[startVertexID..endVertexID] are the original IDs of a node's vertices and originalToReordered is
// the inverse map (this is the reordering the SUL-Tree and COL-Tree use, the graph itself is not renumbered).
class HierarchicalPartition {

    public:
        HierarchicalPartition();
        void build(Graph& graph, GraphPartitioner& partitioner, int fanout, std::size_t maxLeafSize);
        std::string getNetworkName() const;
        int getNumVertices() const;
        int getNumEdges() const;
        int getFanout() const;
        std::size_t getMaxLeafSize() const;
        std::string getMethodName() const;
        double getBuildTimeMs() const;
        int getNumLevels() const;
        double computeIndexSize() const;
        double computeVertexNodePositionsSize() const;

        // Serialized Members
        std::vector<HierarchicalPartitionNode> nodes;
        std::vector<VertexID> reorderedToOriginal;
        std::vector<VertexID> originalToReordered;

    private:
        friend class boost::serialization::access;

        void findBordersAndBoundingBox(Graph& graph, int treeIdx);

        std::string networkName;
        int numVertices;
        int numEdges;
        int fanout;
        std::size_t maxLeafSize;
        std::string methodName;
        double buildTimeMs; // Partitioning time (incl. border detection), part of the construction time of a SUL-Tree built from it

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->networkName;
            ar & this->numVertices;
            ar & this->numEdges;
            ar & this->fanout;
            ar & this->maxLeafSize;
            ar & this->methodName;
            ar & this->buildTimeMs;
            ar & this->nodes;
            ar & this->reorderedToOriginal;
            ar & this->originalToReordered;
        }
};

#endif // _HIERARCHICALPARTITION_H
