// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _GTREE_H
#define _GTREE_H

#include "Graph.h"
#include "../queue/BinaryMinHeap.h"
#include "GraphPartitioner.h"
#include "../utility/Statistics.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/unordered_map.hpp>

class OccurenceList {

    public:
        OccurenceList();
        OccurenceList(std::string setType, double setDensity, int setVariable, int setSize);
        double getObjSetDensity();
        int getObjSetSize();
        int getObjSetVariable();
        std::string getObjSetType();
        double computeIndexSize();
        
        void addLeafOccurence(int leafIdx, VertexID object);
        void addParentOccurence(int parentIdx, int childIdx);

        std::unordered_map<int,std::vector<VertexID>> leafOccurenceList;
        std::unordered_map<int,std::vector<int>> nonLeafOccurenceList;
        std::unordered_map<int,std::unordered_set<VertexID>> leafOccurenceSet;
        // Note: Using vector or unordered_map above doesn't make difference
        // as either case will lead to indirection
        
        // Non-Serialized
        std::unordered_map<int,std::unordered_set<int>> nonLeafOccurenceSet;

private:
        friend class boost::serialization::access;
        
        std::string objSetType;
        double objSetDensity;
        int objSetVariable, objSetSize;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->leafOccurenceList;
            ar & this->leafOccurenceSet; // Serialized to help Dijkstra search (could easily be reconstructed from leafOccurenceList)
            ar & this->nonLeafOccurenceList;
            ar & this->objSetType;
            ar & this->objSetDensity;
            ar & this->objSetVariable;
            ar & this->objSetSize;
        }    
    
};

class GtreeNode {

    public:
        GtreeNode(int parentIdx, int numVertices);
        GtreeNode();
        bool isLeafNode();
        void setLeafNode();
        int getParentIdx();
        int getNumVertices();
        int getParentChildIdx();
        void setParentChildIdx(int idx);
        int addChild(int treeIdx);
        int addBorder(VertexID vertex);
        std::vector<VertexID>& getBorders();
        int getNumBorders();
        int addLeafVertex(VertexID vertex);
        std::vector<VertexID>& getLeafVertices();
        std::unordered_set<VertexID>& getLeafVerticesUset();
        void addChildBorder(VertexID vertex);
        std::vector<VertexID>& getChildBorders();
        std::unordered_set<VertexID>& getChildBordersUset();
        void addGtreePathNode(int nodeIdx);
        void addGtreePathNodes(const std::vector<int>& parentPath);
        std::vector<int>& getGtreePathFromRoot();
        int getBorderIdxInChildBorderVec(int borderIdx);
        int getChildOffsetInChildBorderVec(int childIdx);
        void addChildOffsetInChildBorderVec();
        double computeIndexSizeBytes();
        
        std::vector<int> children;
        std::vector<VertexID> bordersVec;
        std::vector<VertexID> leafVerticesVec;
        std::vector<int> gtreePath;
        std::vector<EdgeWeight> distanceMatrix;
        // Corresponding childBordersVec index for each of this nodes borders (will have same size as bordersVec)
        // If it is the leaf node then it has the corresponding leafVerticesVec index
        // Note: These both cases correlate to the other order in the distanceMatrix vector
        std::vector<int> borderOffsetsInChildBorderVec;
        // The index of the first border of each of this node's children in childBordersVec (will have same size as fanout)
        std::vector<int> childOffsetsInChildBorderVec;
        int matrixRowLength;
        
        // Non-Serialized Members
        std::vector<VertexID> childBordersVec; // This doesn't need to serialized because we retrieve them from child nodes
        std::unordered_set<VertexID> leafVerticesUset;
        std::unordered_set<VertexID> childBordersUset;
        std::unordered_map<VertexID,int> childBorderToChildBorderVecIdx;
        std::unordered_map<VertexID,int> leafVerticeToLeafVerticeVecIdx;
        // Note: We only need these during construction (e.g. to build child offsets)
        
    private:
        friend class boost::serialization::access;

        int parentIdx;
        int numVertices;
        bool isLeaf;
        int parentChildIdx; // Index in parent child vector

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->parentIdx;
            ar & this->numVertices;
            ar & this->isLeaf;
            ar & this->parentChildIdx;
            ar & this->children;
            ar & this->bordersVec;
            ar & this->leafVerticesVec;
            ar & this->gtreePath;
            ar & this->distanceMatrix;
            ar & this->borderOffsetsInChildBorderVec;
            ar & this->childOffsetsInChildBorderVec;
            ar & this->matrixRowLength;
        }
};

class Gtree {

    public:
        Gtree(std::string networkName, int numVertices, int numEdges, int fanout, std::size_t maxLeafSize);
        Gtree();
        std::string getNetworkName();
        int getNumVertices();
        int getNumEdges();
        int getFanout();
        int getTreeSize();
        int getLeafIndex(VertexID vertexID);
        int getParentIndex(int treeIdx);
        std::size_t getMaxLeafSize();
        void buildGtree(Graph& graph);
        void buildTreeHierarchy(Graph& graph);
        void computeDistanceMatrix(Graph& graph);
        int getNumLevels();
        int getNumBorders();
        int getBorderToBorderRelationships();
        int getAvgPathCost(int treeIdx = 0); // By default starts at root
        std::vector<std::vector<int>> getTreeNodesByLevel();
        void addNode(int parentTreeIdx, std::unordered_set<VertexID>& subgraph, Graph& originalGraph, GraphPartitioner& partitioner);
        void addChildren(int parentTreeIdx, std::unordered_set<VertexID>& parentSubgraph, Graph& originalGraph, GraphPartitioner& partitioner);
        EdgeWeight getShortestPathDistance(Graph& graph, VertexID u, VertexID v);
        EdgeWeight SPDistLeaf(VertexID u, VertexID v, int leafNode, Graph& graph);
        EdgeWeight BorderDist(VertexID u, VertexID v, int leafNode);
        EdgeWeight DijkstraDist(VertexID u, VertexID v, int leafNode, Graph& graph);
        std::vector<int> getGtreePath(VertexID u, VertexID v, int& firstLCAChild);
        EdgeWeight SPDistToSourceLeafNode(VertexID u, int treeIdx);
        EdgeWeight SPDistToParentNode(int childTreeIdx, int parentTreeIdx, bool computeSPDist = true);
        EdgeWeight SPDistToSiblingNode(int firstLCAChildIdx, int targetLCAChildIdx, int LCAIdx, bool computeSPDist = true);
        EdgeWeight SPDistToChildNode(int childTreeIdx, int parentTreeIdx, bool computeSPDist = true);
        EdgeWeight SPDistToLeafTarget(VertexID target, int leafIdx);
        void getKNNs(OccurenceList& occList, unsigned int k, VertexID queryVertexID, std::vector<VertexID>& kNNs, 
                     std::vector<EdgeWeight>& kNNDistances, Graph& graph);
        bool getSourceLeafkNNsByINE(VertexID queryVertex, unsigned int k, std::unordered_set<VertexID>& targets, int leafNode, Graph& graph, std::vector<VertexID>& kNNs, 
                                    std::vector<EdgeWeight>& kNNDistances, BinaryMinHeap<EdgeWeight,VertexID>& pqueue);
        void initialiseGtreeQueryStructure();
        int getIdxInLeafVerticesVec(VertexID u);
        void setLeafVerticesVecIdx(VertexID u, int leafVerticesVecIdx);
        int getIdxInLeafBordersVec(VertexID u);
        void setLeafBordersVecIdx(VertexID u, int leafBordersVecIdx);
        bool isEdgeInLeafSubgraph(VertexID edge);
        void setEdgeNotInSubgraph(VertexID edge);
        double computeIndexSize();
        int getComputations(int leafIdx, int currIdx);
        
        // getRepeatedShortestPathDistance is used to find the shortest path distance to multiple
        // target from the same query vertex. This assumes the G-tree index is not being used 
        // concurrently for difference query vertices (currently it does not need to be). The visited
        // vector indicates the G-tree nodes (by ID) that have been visited in previous calls
        // to this function. To start query for a new vertex, we need only assign all elements
        // in this vector to be false.
        EdgeWeight getRepeatedShortestPathDistance(Graph& graph, VertexID u, VertexID v, std::vector<bool>& visited, std::vector<EdgeWeight>& leafVertexDistances, 
                                                   BinaryMinHeap<EdgeWeight,VertexID>& pqueue, std::unordered_set<VertexID>& leafVertexVisited);
        EdgeWeight getRepeatedSourceLeafShortestPathDistance(Graph& graph, int leafNode, VertexID v, std::vector<EdgeWeight>& leafVertexDistances, 
                                                             BinaryMinHeap<EdgeWeight,VertexID>& pqueue, std::unordered_set<VertexID>& leafVertexVisited);

        std::vector<GtreeNode> treeNodes;
        std::vector<int> leafIdxs;
        // Note: The sourceToTreeNodeBorderDist is used to store the intermediary results during
        // querying using the G-tree index. Firstly this structure is included in the index size.
        // Secondly this structure can be re-used without needing to cleared for different 
        // query vertices because G-tree queries traverse the hierarchy incrementally (source 
        // to source leaf, child to parent, LCA sibling to LCA sibling, parent to child, and
        // finally target leaf to leaf). This means that when we want to calculate the 
        // distances to a new G-tree node, we have already computed the distances to it's preceding node
        // (child, parent or sibling). That is it will always be overwritten for each new query vertex. 
        // There may be incorrect distances in it (to a previous query vertex) but this data will 
        // never be used until after it is overwritten for the current query vertex. Currently this is
        // true for getShortestPathDistance and getKNNs.
        std::vector<std::vector<EdgeWeight>> sourceToTreeNodeBorderDist;
        std::vector<int> vertexIDLeafVerticesVecIdx; // For each graph vertex we store the index occupied by that graph vertex in it's leaf Gtree node leafVerticesVec
        std::vector<int> vertexIDLeafBordersVecIdx; // For each graph vertex we store the index occupied by that graph vertex in it's leaf Gtree node bordersVec
        std::vector<bool> edgeInLeafSubgraph; // For each graph edge we store whether it is an edge that leads to another vertex in the subgraph to which the originating vertex belongs (use later to prune DijkstraDistMultiTarget)
        std::vector<int> vertexLeafIdxs; // For each graph vertex we store it's leaf Gtree node's index in the tree
        
        // Non-Serialized Members
        Statistics stats;
        
    private:
        friend class boost::serialization::access;

        std::string networkName;
        int numVertices;
        int numEdges;
        int fanout;
        std::size_t maxLeafSize;

        // Non-Serialized Members
        static const int ROOT_PARENT_INDEX = -1;
#
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->networkName;
            ar & this->numVertices;
            ar & this->numEdges;
            ar & this->fanout;
            ar & this->maxLeafSize;
            ar & this->treeNodes;
            ar & this->leafIdxs;
            ar & this->sourceToTreeNodeBorderDist;
            ar & this->vertexIDLeafVerticesVecIdx;
            ar & this->vertexIDLeafBordersVecIdx;
            ar & this->edgeInLeafSubgraph;
            ar & this->vertexLeafIdxs;
        }
        
};

#endif // _GTREE_H