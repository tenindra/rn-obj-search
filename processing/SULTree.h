// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _SULTREE_H
#define _SULTREE_H

#include "Graph.h"
#include "ALT.h"
#include "pruned_highway_labeling.h"
#include "../queue/BinaryMaxHeap.h"
#include "HierarchicalPartition.h"
#include "../utility/Statistics.h"
#include "../utility/AggregateFunction.h"

#include <random>
#include <unordered_set>
#include <boost/archive/binary_oarchive.hpp>

enum COLTREE_LANDMARK_TYPE
{
    SUBGRAPH_RANDOM,
    SUBGRAPH_SLICE_FURTHEST,
    SUBGRAPH_FURTHEST_BORDER,
    SUBGRAPH_BORDER_MINMAX
};

// Methods to compute the subgraph distance lists (SDL) of a node, i.e. the distances from each of its landmarks to its vertices.
// SDL_BORDER_SET is the method of the index; the other two are its predecessors, kept for the SDL experiment (lmk_exp -e sdl_exp)
enum SDL_METHOD
{
    SDL_HASH_TABLE, // Dijkstra search with hash-table based subgraph membership and distances
    SDL_DIJKSTRA, // Dijkstra search with vector based membership and distances, stopping when all subgraph vertices are settled
    SDL_BORDER_SET // SDL_DIJKSTRA pruned with bounds through the border set and the root landmarks
};

class SULTreeNode {

    public:
        SULTreeNode(int parentIdx, int numVertices);
        SULTreeNode();
        bool isLeafNode();
        int getParentIdx();
        unsigned int getNumVertices();
        double computeIndexSizeBytes();

        int parentIdx;
        unsigned int numVertices;
        bool isLeaf;
        VertexID startVertexID; // This is by reordered IDs
        VertexID endVertexID; // This is by reordered IDs (endVertexID is inclusive)

        // Serialized Members
        std::vector<int> childNodeIdxs;
        std::vector<VertexID> landmarks;
        
        // Non-Serialized Members
        std::vector<double> spilloverFactors;
        std::vector<bool> isDistanceListComputed; // For each landmark, whether the landmark selection already stored its distance list
        Coordinate minX, minY, maxX, maxY;
        // Note: Below members are not defined for root node
        std::vector<VertexID> borderVertexIDs;
        std::vector<EdgeWeight> minBorderDist; // For each root landmark, store the min distance to any border
        std::vector<EdgeWeight> maxBorderDist; // For each root landmark, store the max distance to any border

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
            ar & this->landmarks;
        }
};

class SULTree {

    public:
        SULTree(std::string networkName, int numVertices, int numEdges, int numNodeLandmarks);
        SULTree();
        std::string getNetworkName();
        int getNumVertices();
        int getNumEdges();
        int getFanout();
        std::size_t getMaxLeafPartitionSize();
        int getNumLandmarksPerNode();

        // Build the tree (landmarks and vertex-to-landmark distance lists) on a hierarchical partition of the graph
        void buildSULTree(Graph& graph, const HierarchicalPartition& partition, COLTREE_LANDMARK_TYPE landmarkType, SDL_METHOD sdlMethod = SDL_METHOD::SDL_BORDER_SET);
        // Phase times of the last build
        double getLandmarkTimeMs();
        double getDistanceListTimeMs();
        double getAllocationTimeMs();
        void allocateNodeVertexLmkDists();
        void computeNodeDistanceLists(Graph& graph, TreeNodeID treeIdx);
        void computeRootDistanceLists(Graph& graph);
        void computeBorderDistanceBounds(TreeNodeID treeIdx);
        void selectNodeLandmarks(Graph& graph, TreeNodeID treeIdx);

        // Distance List Computation Variants
        int findSSMTDistances(Graph& graph, int treeIdx, VertexID source, std::vector<EdgeWeight>& results);
        int findSSMTDistancesViaBorderSet(Graph& graph, int treeIdx, VertexID source, std::vector<EdgeWeight>& results);

        // Landmark Selection Methods
        COLTREE_LANDMARK_TYPE getLandmarkTypeFromString(std::string landmarkTypeString, bool& success);
        std::mt19937 getLandmarkEngine(int treeIdx);
        void selectSlicedFurthestLandmarks(Graph& originalGraph, int treeIdx, Coordinate minX, Coordinate minY, Coordinate maxX, Coordinate maxY);
        void selectRandomLandmarks(int treeIdx);
        void selectRandomFurthestBorderLandmarks(Graph& originalGraph, int treeIdx);
        void findFurthestBorderSubset(Graph& graph, int treeIdx, VertexID seedBorderID);
        void selectMinmaxLandmarksFromBorders(Graph& originalGraph, int treeIdx);
        void findMinmaxVerticesFromSourceSet(Graph& graph, int treeIdx, std::vector<VertexID>& sourceSet, std::size_t numMinmaxVertices,
                                             std::vector<int>& minmaxVertexPositions);

        // Other/Misc
        double computeIndexSize();
        std::vector<std::vector<int>> getTreeNodesByLevel();
        int getNumLevels();
        EdgeWeight computeRootLowerBound(VertexID s, VertexID t);

        // Serialized Members
        std::vector<SULTreeNode> treeNodes;
        
        // Store vertex to landmark distances for all tree nodes grouped together by tree node
        // where nodeStartIndex marks the beginning of the list of distances. The distances
        // to all landmarks are further grouped together for each vertex.
        std::vector<long> nodeStartIndex;
        std::vector<EdgeWeight> nodeVertexLmkDists; 

        // SULTree to COLTree conversion optimization relies on reordering old VertexIDs
        // such that VertexIDs for every SULTree are contiguous. As it's complicated 
        // to implement this due to change VertexIDs for all methods, we implement as
        // a mapping for SULTree/COLTree. Note that these vertexes are not include in index size
        // computations because as it is not need in standalone implementation of SULTree,
        // in which case you can just permanently reorder vertexIDs.
        std::vector<VertexID> reorderedVertexIDtoOriginalVertexID; // During construction a node's vertices are the range [startVertexID,endVertexID] of this (from the HierarchicalPartition)
        std::vector<VertexID> originalVertexIDtoReorderedVertexID;

        // Non-Serialized Members
        static const int ROOT_PARENT_INDEX = -1;
        std::vector<std::vector<double>> nodeSpilloverFactors; // Dijkstra search spillover factor for each node and each
                                                               // of its landmarks (vertices settled vs vertices subgraph)

    private:
        friend class boost::serialization::access;

        std::string networkName;
        int numVertices;
        int numEdges;
        int fanout;
        std::size_t maxLeafSize;
        unsigned int numNodeLandmarks;
        COLTREE_LANDMARK_TYPE landmarkType;

        // Seed of the landmark selection, combined with the tree node (getLandmarkEngine); any other value selects other landmarks
        static const unsigned int LANDMARK_SEED = 0;

        // Non-Serialized Members
        std::size_t nodeVertexLmkDistsSize = 0; // Size of nodeVertexLmkDists (could be larger than int)
        SDL_METHOD sdlMethod = SDL_METHOD::SDL_BORDER_SET; // Subgraph distance list method of the current build
        double landmarkTimeMs = 0, distanceListTimeMs = 0, allocationTimeMs = 0; // Phase times of the last build

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->networkName;
            ar & this->numVertices;
            ar & this->numEdges;
            ar & this->fanout;
            ar & this->maxLeafSize;
            ar & this->numNodeLandmarks;
            ar & this->landmarkType;
            ar & this->treeNodes;
            ar & this->nodeStartIndex;
            ar & this->nodeVertexLmkDists;
            ar & this->reorderedVertexIDtoOriginalVertexID;
            ar & this->originalVertexIDtoReorderedVertexID;
        }
        
};

struct COLTreeNodeQueryData {
    int leftIdx = -1, rightIdx = -1, lmIdx = -1;
    bool intialized = false;
    EdgeWeight lmkDist = 0;
    // Only intialized needs resetting: each query sets lmIdx before it reads it, and the other fields when it initialises the node
    void reset() {
        intialized = false;
    }
        
    private:
        friend class boost::serialization::access;

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->leftIdx;
            ar & this->rightIdx;
            ar & this->lmIdx;
            ar & this->lmkDist;
            ar & this->intialized;
        }
};

class COLTreeNode {

    public:
        COLTreeNode(int parentIdx, int numObjects);
        COLTreeNode();
        bool hasObjectList();
        void setObjectListStatus(bool status);
        int getParentIdx();
        unsigned int getNumObjects();
        void addChildNode(int childIdx);
        double computeIndexSizeBytes();

        // Shared Node Bound Functions
        int getIndexClosestToDistance(EdgeWeight matchDist, int listIdx);
        std::size_t getNumObjectsBeforeDistance(EdgeWeight dist, int listIdx, std::size_t numElements);
        std::size_t getNumObjectsUpToDistance(EdgeWeight dist, int listIdx, std::size_t numElements);

        // kNN Node Bound Functions
        EdgeWeight getBestNodeLowerbound(ALT& alt, VertexID queryVertex);
        EdgeWeight getBestNodeLowerboundAllCases(ALT& alt, VertexID queryVertex);
        EdgeWeight getBestNodeLowerboundAndFurthestLmk(ALT& alt, VertexID queryVertex, int& furthestLandmarkIdx);
        int getFurthestLmkIdxByLB(ALT& alt, VertexID queryVertex);

        // AkNN Node Bound Functions
        void getBestNodeLowerbounds(ALT& alt, std::vector<VertexID>& queryVertices, std::vector<EdgeWeight>& queryDistances);
        void getBestNodeLowerboundsAndFurthestSumLmk(ALT& alt, std::vector<VertexID>& queryVertices, std::vector<EdgeWeight>& queryDistances,
                                                             int& furthestSumLandmarkIdx);
        int getFurthestLmkIdxByLB(ALT& alt, std::vector<VertexID>& queryVertices);

        // FkN Node Bound Functions
        EdgeWeight getBestNodeUpperbound(ALT& alt, VertexID queryVertex);
        int getNearestLmkIdxByUB(ALT& alt, VertexID queryVertex);
        EdgeWeight getBestNodeUpperboundAndNearestLmk(ALT& alt, VertexID queryVertex, int& furthestLandmarkIdx);

        // Range Search Bound Functions
        void getBestNodeLowerAndUpperbound(ALT& alt, VertexID queryVertex, EdgeWeight& lb, EdgeWeight& ub);
        void getRangeSearchLeafBounds(ALT& alt, VertexID queryVertex, EdgeWeight& lb, EdgeWeight& ub, std::vector<EdgeWeight>& queryLmkLBs, std::vector<EdgeWeight>& queryLmkUBs);
        int getRangeSearchListIdx(EdgeWeight range, std::vector<EdgeWeight>& queryLmkLBs, std::vector<EdgeWeight>& queryLmkUBs);

#if defined(COLLECT_STATISTICS)
        // Support kNN Search
        void updateQueueWithMinObj(int listIdx, int& leftPos, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, 
                                   BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, Statistics& stats);
        bool updateQueue(int listIdx, int& leftPos, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, 
                         BinaryMinHeap<EdgeWeight,VertexID>& pqueue,ALT& alt, EdgeWeight& listLB, Statistics& stats);

        // Support AkNN Search
        void updateQueueWithMinObjAkNN(int listIdx, int& leftPos, int& rightPos, std::vector<EdgeWeight>& queryVertexLmkDists, std::vector<VertexID>& queryVertices, 
                                       std::vector<EdgeWeight>& queryDistances, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, AggregateFunction<EdgeWeight>& agg, Statistics& stats);
        bool updateQueueAkNN(int listIdx, int& leftPos, int& rightPos, std::vector<EdgeWeight>& queryVertexLmkDists, std::vector<VertexID>& queryVertices, 
                             std::vector<EdgeWeight>& queryDistances, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, EdgeWeight& listLB, AggregateFunction<EdgeWeight>& agg, Statistics& stats);

        // Support FkN Search
        void updateQueueWithMaxObjFkN(int listIdx, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, BinaryMaxHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, Statistics& stats);
        bool updateQueueFkN(int listIdx, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, BinaryMaxHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, EdgeWeight& listUB, Statistics& stats);

        // Support Range Search
        void rangeSearchObjectList(int listIdx, VertexID queryVertexID, EdgeWeight range, ALT& alt, PrunedHighwayLabeling& phl, std::vector<VertexID>& results, std::vector<VertexID>& resultDistances, Statistics& stats);
#else
        // Support kNN Search
        void updateQueueWithMinObj(int listIdx, int& leftPos, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt);
        bool updateQueue(int listIdx, int& leftPos, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, EdgeWeight& listLB);

        // Support AkNN Search
        void updateQueueWithMinObjAkNN(int listIdx, int& leftPos, int& rightPos, std::vector<EdgeWeight>& queryVertexLmkDists, std::vector<VertexID>& queryVertices, 
                                       std::vector<EdgeWeight>& queryDistances, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, AggregateFunction<EdgeWeight>& agg);
        bool updateQueueAkNN(int listIdx, int& leftPos, int& rightPos, std::vector<EdgeWeight>& queryVertexLmkDists, std::vector<VertexID>& queryVertices, 
                             std::vector<EdgeWeight>& queryDistances, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, EdgeWeight& listLB, AggregateFunction<EdgeWeight>& agg);

        // Support FkN Search
        void updateQueueWithMaxObjFkN(int listIdx, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, BinaryMaxHeap<EdgeWeight,VertexID>& pqueue, ALT& alt);
        bool updateQueueFkN(int listIdx, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, BinaryMaxHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, EdgeWeight& listUB);

        // Support Range Search
        void rangeSearchObjectList(int listIdx, VertexID queryVertexID, EdgeWeight range, ALT& alt, PrunedHighwayLabeling& phl, std::vector<VertexID>& results, std::vector<VertexID>& resultDistances);
#endif
        
        // Serialized Members
        std::vector<int> childNodeIdxs; // Only child indexes from original tree that contain objects
        // Data structures below are not filled for root node
        std::vector<VertexID> landmarks;
        std::vector<std::vector<VertexEdgeWeightPair>> lmkObjectLists; // Leaf node only
        std::vector<EdgeWeight> lmkOLMaxDist; // Max distance for an object from each landmark
        std::vector<EdgeWeight> lmkOLMinDist; // Min distance for an object from each landmark
        std::vector<EdgeWeight> globalLmkOLMaxDist; // Max distance for an object from each landmark using root level landmarks
        std::vector<EdgeWeight> globalLmkOLMinDist; // Min distance for an object from each landmark using root level landmarks

        // Non-Serialized Members
        // Temporary variables used during construction to track the first/last indicides of objects that 
        // that fall witin the range of this node's (reordered) subgraph vertex IDs in the list of all sorted object IDs
        int firstObjectIDIndex, lastObjectIDIndex; 
        
    private:
        friend class boost::serialization::access;

        // Serialized Members
        int parentIdx;
        bool hasOL;
        unsigned int numObjects;

        // Non-Serialized Members

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->parentIdx;
            ar & this->childNodeIdxs;
            ar & this->landmarks;
            ar & this->lmkObjectLists;
            ar & this->lmkOLMaxDist;
            ar & this->lmkOLMinDist;
            ar & this->globalLmkOLMaxDist;
            ar & this->globalLmkOLMinDist;
            ar & this->hasOL;
            ar & this->numObjects;
        }
};

class COLTree {

    public:
        COLTree(std::string setType, double setDensity, int setVariable, int setSize);
        COLTree();
        void buildCOLTTopDown(SULTree& sltree, std::vector<VertexID>& objectIDs, unsigned int maxOLSize);
        double getObjSetDensity();
        int getObjSetSize();
        int getObjSetVariable();
        std::string getObjSetType();
        double computeIndexSize();
        int getMaxOLSize();
        int getNumNonEmptyNodes();
        int getNodeDepth(int treeNodeIdx);
        int getTreeDepth();
        EdgeWeight computeGlobalLowerbound(int coltNodeIdx, VertexID queryVertex, SULTree& sltree);
        EdgeWeight computeGlobalUpperbound(int coltNodeIdx, VertexID queryVertex, SULTree& sltree);
        void computeGlobalLowerAndUpperbound(int coltNodeIdx, VertexID queryVertex, SULTree& sltree, EdgeWeight& lb, EdgeWeight& ub);
        void computeGlobalAggregateLowerbound(int coltNodeIdx, std::vector<VertexID>& queryVertices, std::vector<EdgeWeight>& queryDistances, SULTree& sltree, bool& used);

        std::vector<COLTreeNode> treeNodes;
        // Structure for keeping track of object list evaluation in COL-Trees.
        // This is reset in every query (and included in query time cost). For
        // a real-world implementation, several vectors would be allocated in
        // a shared pool that can be used by each concurrent query thread.
        std::vector<COLTreeNodeQueryData> queryData;
        // Bounds on the query's distance to each landmark of the leaf node last bounded by range search.
        // Like queryData, they hold the state of one query at a time (not serialized)
        std::vector<EdgeWeight> rangeQueryLmkLBs, rangeQueryLmkUBs;

    private:
        friend class boost::serialization::access;

        std::string objSetType;
        double objSetDensity;
        int objSetVariable, objSetSize;
        unsigned int maxOLSize;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->objSetType;
            ar & this->objSetDensity;
            ar & this->objSetVariable;
            ar & this->objSetSize;
            ar & this->treeNodes;
            ar & this->maxOLSize;
            ar & this->queryData;
        }
        
};

#endif // _SULTREE_H
