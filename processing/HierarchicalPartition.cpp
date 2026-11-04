// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "HierarchicalPartition.h"
#include "../utility/StopWatch.h"

#include <algorithm>
#include <iostream>
#include <stack>

HierarchicalPartitionNode::HierarchicalPartitionNode(int parentIdx, int numVertices):
    parentIdx(parentIdx), numVertices(numVertices), isLeaf(false), startVertexID(0), endVertexID(0),
    minX(0), minY(0), maxX(0), maxY(0) {}

HierarchicalPartitionNode::HierarchicalPartitionNode() {}

HierarchicalPartition::HierarchicalPartition(): numVertices(0), numEdges(0), fanout(0), maxLeafSize(0), buildTimeMs(0) {}

std::string HierarchicalPartition::getNetworkName() const
{
    return this->networkName;
}

int HierarchicalPartition::getNumVertices() const
{
    return this->numVertices;
}

int HierarchicalPartition::getNumEdges() const
{
    return this->numEdges;
}

int HierarchicalPartition::getFanout() const
{
    return this->fanout;
}

std::size_t HierarchicalPartition::getMaxLeafSize() const
{
    return this->maxLeafSize;
}

std::string HierarchicalPartition::getMethodName() const
{
    return this->methodName;
}

double HierarchicalPartition::getBuildTimeMs() const
{
    return this->buildTimeMs;
}

void HierarchicalPartition::build(Graph& graph, GraphPartitioner& partitioner, int fanout, std::size_t maxLeafSize)
{
    if (fanout < 2 || static_cast<std::size_t>(fanout) > maxLeafSize) {
        std::cerr << "Fanout must be at least 2 and at most the maximum leaf size" << std::endl;
        std::exit(1);
    }
    this->networkName = graph.getNetworkName();
    this->numVertices = graph.getNumVertices();
    this->numEdges = graph.getNumEdges();
    this->fanout = fanout;
    this->maxLeafSize = maxLeafSize;
    this->methodName = partitioner.getName();

    StopWatch sw, swPhase;
    double methodTime = 0, reorderTime = 0, borderTime = 0; // Breakdown of the partitioning time for the console
    sw.start();

    // Initially the vertex order is the original order; partitioning a node permutes the vertices within its range only
    this->reorderedToOriginal = graph.getVertexIDsVector();
    this->originalToReordered = this->reorderedToOriginal;

    RN_DEBUG_ASSERT(this->nodes.empty() && "A HierarchicalPartition is built once, after it is constructed");
    this->nodes.push_back(HierarchicalPartitionNode(constants::ROOT_PARENT_INDEX,graph.getNumVertices()));
    this->nodes[0].startVertexID = 0;
    this->nodes[0].endVertexID = graph.getNumVertices()-1;
    // The root has no borders but min/max X/Y are needed for some landmark selection methods (e.g., slice)
    graph.getMinMaxCoordinates(this->nodes[0].minX,this->nodes[0].maxX,this->nodes[0].minY,this->nodes[0].maxY);

    std::vector<int> parts, partSizes(fanout), partStarts(fanout), nextPos;
    std::vector<VertexID> sortedRange;

    std::stack<int> pending;
    pending.push(0);
    while (!pending.empty()) {
        int treeIdx = pending.top();
        pending.pop();
        if (static_cast<std::size_t>(this->nodes[treeIdx].numVertices) <= maxLeafSize) {
            // With at most maxLeafSize vertices this is a leaf node and we no longer partition (its vertex order is final)
            this->nodes[treeIdx].isLeaf = true;
            continue;
        }

        VertexID start = this->nodes[treeIdx].startVertexID;
        int numNodeVertices = this->nodes[treeIdx].numVertices;
        swPhase.reset();
        swPhase.start();
        partitioner.partition(graph,this->reorderedToOriginal,start,numNodeVertices,this->originalToReordered,fanout,parts);
        swPhase.stop();
        methodTime += swPhase.getTimeMs();

        swPhase.reset();
        swPhase.start();

        // Stable counting sort of the node's vertices by part so that each child's vertices are contiguous (child 0 first).
        // The reordering of the whole node is updated before any child is processed, so that border detection can test
        // the membership of a neighbour by its position alone
        std::fill(partSizes.begin(),partSizes.end(),0);
        for (int i = 0; i < numNodeVertices; ++i) {
            ++partSizes[parts[i]];
        }
        for (int p = 1; p < fanout; ++p) { // partStarts[0] stays 0
            partStarts[p] = partStarts[p-1] + partSizes[p-1];
        }
        nextPos = partStarts;
        sortedRange.resize(numNodeVertices);
        for (int i = 0; i < numNodeVertices; ++i) {
            sortedRange[nextPos[parts[i]]++] = this->reorderedToOriginal[start+i];
        }
        for (int i = 0; i < numNodeVertices; ++i) {
            this->reorderedToOriginal[start+i] = sortedRange[i];
            this->originalToReordered[sortedRange[i]] = start+i;
        }
        swPhase.stop();
        reorderTime += swPhase.getTimeMs();

        swPhase.reset();
        swPhase.start();

        // Create a child node for each part and add it to the stack (the last child is processed
        // in the next loop iteration so tree node creation follows DFS order)
        for (int p = 0; p < fanout; ++p) {
            if (partSizes[p] == 0) {
                std::cerr << "Partitioning produced an empty part for tree node " << treeIdx << " (" << numNodeVertices
                          << " vertices, fanout " << fanout << ")" << std::endl;
                std::exit(1);
            }
            int childIdx = this->nodes.size();
            this->nodes.push_back(HierarchicalPartitionNode(treeIdx,partSizes[p]));
            this->nodes[treeIdx].childNodeIdxs.push_back(childIdx);
            this->nodes[childIdx].startVertexID = start + partStarts[p];
            this->nodes[childIdx].endVertexID = start + partStarts[p] + partSizes[p] - 1;
            this->findBordersAndBoundingBox(graph,childIdx);
            pending.push(childIdx);
        }
        swPhase.stop();
        borderTime += swPhase.getTimeMs();
    }

    sw.stop();
    this->buildTimeMs = sw.getTimeMs();
    std::cout << "Partitioning time = " << this->buildTimeMs << "ms (" << this->methodName << " " << methodTime << "ms, vertex reordering "
              << reorderTime << "ms, borders and bounding boxes " << borderTime << "ms)" << std::endl;
}

// The border vertices of a node are those with a neighbour outside the node, i.e., whose position in the vertex order is
// outside the node's range (vertices of other nodes never move into this range later). The bounding box is computed in
// the same pass over the node's vertices.
void HierarchicalPartition::findBordersAndBoundingBox(Graph& graph, int treeIdx)
{
    VertexID start = this->nodes[treeIdx].startVertexID, end = this->nodes[treeIdx].endVertexID;
    Coordinate x, y;
    int adjListStart, nextAdjListStart;
    for (VertexID pos = start; pos <= end; ++pos) {
        VertexID vertex = this->reorderedToOriginal[pos];
        adjListStart = graph.getEdgeListStartIndex(vertex);
        nextAdjListStart = graph.getEdgeListSize(vertex);
        for (int i = adjListStart; i < nextAdjListStart; ++i) {
            VertexID neighbourPos = this->originalToReordered[graph.edges[i].first];
            if (neighbourPos < start || neighbourPos > end) {
                // Break so we only add this vertex as border once (it may have several neighbours outside)
                this->nodes[treeIdx].borderVertexIDs.push_back(vertex);
                break;
            }
        }

        graph.getCoordinates(vertex,x,y);
        if (pos == start) {
            // This is the first vertex, so we need to initialize min/max coordinates
            this->nodes[treeIdx].minX = x;
            this->nodes[treeIdx].maxX = x;
            this->nodes[treeIdx].minY = y;
            this->nodes[treeIdx].maxY = y;
        } else {
            if (x < this->nodes[treeIdx].minX) {
                this->nodes[treeIdx].minX = x;
            }
            if (x > this->nodes[treeIdx].maxX) {
                this->nodes[treeIdx].maxX = x;
            }
            if (y < this->nodes[treeIdx].minY) {
                this->nodes[treeIdx].minY = y;
            }
            if (y > this->nodes[treeIdx].maxY) {
                this->nodes[treeIdx].maxY = y;
            }
        }
    }
}

int HierarchicalPartition::getNumLevels() const
{
    // Traverse level by level, counting levels without storing them
    int numLevels = 0;
    std::vector<int> currentLevel, nextLevel;
    nextLevel.push_back(0);
    while (nextLevel.size() != 0) {
        ++numLevels;
        currentLevel.swap(nextLevel);
        nextLevel.clear();
        for (int i: currentLevel) {
            nextLevel.insert(nextLevel.end(),this->nodes[i].childNodeIdxs.begin(),this->nodes[i].childNodeIdxs.end());
        }
    }
    return numLevels;
}

double HierarchicalPartition::computeIndexSize() const
{
    double memoryUsage = 0;
    for (std::size_t i = 0; i < this->nodes.size(); ++i) {
        memoryUsage += sizeof(HierarchicalPartitionNode);
        memoryUsage += sizeof(int)*this->nodes[i].childNodeIdxs.size();
        memoryUsage += sizeof(VertexID)*this->nodes[i].borderVertexIDs.size();
    }
    memoryUsage += sizeof(VertexID)*(this->reorderedToOriginal.size()+this->originalToReordered.size());
    return memoryUsage/(1024*1024);
}

// Memory the per-vertex node positions of a SUL-Tree would need without the subgraph vertex ordering, i.e. for every
// vertex its position in the vertex list of each tree node containing it (root included, as the SUL-Tree construction
// stored them before the ordering). It is computed from the hierarchy without storing them:
// a vertex of a leaf at depth d has d positions, plus the vector per vertex.
double HierarchicalPartition::computeVertexNodePositionsSize() const
{
    double memoryUsage = sizeof(std::vector<int>)*static_cast<double>(this->numVertices);
    std::vector<int> currentLevel, nextLevel;
    nextLevel.push_back(0);
    for (int depth = 1; nextLevel.size() != 0; ++depth) {
        currentLevel.swap(nextLevel);
        nextLevel.clear();
        for (int i: currentLevel) {
            if (this->nodes[i].isLeaf) {
                memoryUsage += sizeof(int)*static_cast<double>(depth)*this->nodes[i].numVertices;
            }
            nextLevel.insert(nextLevel.end(),this->nodes[i].childNodeIdxs.begin(),this->nodes[i].childNodeIdxs.end());
        }
    }
    return memoryUsage/(1024*1024);
}
