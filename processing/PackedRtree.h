// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _PACKEDRTREE_H
#define _PACKEDRTREE_H

#include "../queue/BinaryMinHeap.h"
#include "../queue/BinaryMaxHeap.h"
#include "../common.h"
#include "../utility/AggregateFunction.h"

#include <vector>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

/* 
 * 2D In-Memory Static Rtree Implementation for Point Data
 */

// Rectangles of the R-tree, and the pointer-linked nodes and entries that PackedRtree::bulkLoad builds the tree
// from before copying it into its flat node vector

class PtrRtreeNode;

class Rectangle {

    public:
        Rectangle() {}
        void setDimensions(Coordinate left, Coordinate bottom, Coordinate right, Coordinate top);
        void setDimensions(Rectangle& rect);
        void setCentroid();
        void expand(Rectangle& rect);
        EuclideanDistanceSqrd getMinDistSqrd(Coordinate x, Coordinate y);
        EuclideanDistanceSqrd getMinMaxDistSqrd(Coordinate x, Coordinate y);
        EuclideanDistanceSqrd getMaxDistSqrd(Coordinate x, Coordinate y);
        EuclideanDistanceSqrd getObjectDistSqrd(Coordinate x, Coordinate y);
        double computeIndexSizeBytes();
        bool intersects(Rectangle& rect);
        
        Coordinate left;
        Coordinate bottom;
        Coordinate right;
        Coordinate top;
        Coordinate centreX;
        Coordinate centreY;

    private:
        friend class boost::serialization::access;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->left;
            ar & this->bottom;
            ar & this->right;
            ar & this->top;
            ar & this->centreX;
            ar & this->centreY;
        }
};

class PtrRtreeEntry {

    public:
        PtrRtreeEntry() {}
        void setEntryData(CoordinatePair& ptCoords, PtrRtreeNode* rtreeNode, VertexID id = 0);
        void setEntryData(Rectangle rect, PtrRtreeNode* rtreeNode, VertexID id = 0);
        static bool compareX(const PtrRtreeEntry& entry1, const PtrRtreeEntry& entry2);
        static bool compareY(const PtrRtreeEntry& entry1, const PtrRtreeEntry& entry2);
        
        Rectangle mbr;
        PtrRtreeNode* rtreeNode;
        VertexID id; // Object ID for entries of leaf node

};


class PtrRtreeNode {

    public:
        PtrRtreeNode(int level);
        void addEntry(PtrRtreeEntry entry);
        bool isLeaf();
        
        std::vector<PtrRtreeEntry> entries;
        int level; // 0 for leaf nodes

};

struct PackedRtreeDataTuple {
    bool isNode;
    int nodeIdx;
    VertexID objectID;
    void setNodeIdx(int nodeIdx) {
        this->isNode = true;
        this->nodeIdx = nodeIdx;
    }
    void setObjectID(VertexID objectID) {
        this->isNode = false;
        this->objectID = objectID;
    }
};

class PackedRtreeEntry {

    public:
        PackedRtreeEntry() {}
        void setEntryData(Rectangle rect, int nodeIdx, VertexID id = 0);
        double computeIndexSizeBytes();
        
        Rectangle mbr;
        int nodeIdx;
        VertexID id; // Object ID for entries of leaf node
        
    private:
        friend class boost::serialization::access;

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->mbr;
            ar & this->nodeIdx;
            ar & this->id;
        }
};

class PackedRtreeNode {

    public:
        PackedRtreeNode() {}
        PackedRtreeNode(int level);
        void addEntry(PackedRtreeEntry entry);
        bool isLeaf();
        double computeIndexSizeBytes();
        
        std::vector<PackedRtreeEntry> entries;
        int level; // 0 for leaf nodes

    private:
        friend class boost::serialization::access;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->entries;
            ar & this->level;
        }
};


class PackedRtree {

    public:
        PackedRtree() {}
        PackedRtree(int branchFactor, std::string setType, double setDensity, int setVariable, int setSize);
        void bulkLoad(std::vector<VertexID>& objectIDs, std::vector<CoordinatePair>& objectCoords);
        PtrRtreeNode* buildLevel(int level, std::vector<PtrRtreeEntry>& levelEntries);
        BinaryMinHeap<EuclideanDistanceSqrd,PackedRtreeDataTuple> getKNNs(unsigned int k, Coordinate queryPointX, Coordinate queryPointY, 
                                                                    std::vector<VertexID>& kNNs, std::vector<EuclideanDistanceSqrd>& kNNDistancesSqrd);
        bool getNextNearestNeighbour(BinaryMinHeap<EuclideanDistanceSqrd,PackedRtreeDataTuple>& heap, Coordinate queryPointX, 
                                     Coordinate queryPointY, VertexID& nn, EuclideanDistanceSqrd& nnDistSqrd);
        std::vector<PtrRtreeNode*> getAllNodePtrs(PtrRtreeNode* root);
        double getObjSetDensity();
        int getObjSetSize();
        int getObjSetVariable();
        std::string getObjSetType();
        int getBranchFactor();
        double computeIndexSize();
        void vectorise(PtrRtreeNode* root);
        std::vector<VertexID> windowQuery(Rectangle& rect);
        
        // Support for Aggregate kNNs: the heap keys are the aggregate of the query points' Euclidean distances, a lower bound on
        // the aggregate network distance (aggregating squared distances and taking the square root would be a much weaker bound for
        // a sum)
        BinaryMinHeap<EuclideanDistance,PackedRtreeDataTuple> getAkNNs(unsigned int k, std::vector<Coordinate>& queryPointsX, std::vector<Coordinate>& queryPointsY, 
                                                                   std::vector<VertexID>& AkNNs, std::vector<EuclideanDistance>& AkNNDistances, 
                                                                   AggregateFunction<EuclideanDistance>& agg, std::vector<EuclideanDistance>& queryDistances);
        bool getNextAggregateNearestNeighbour(BinaryMinHeap<EuclideanDistance,PackedRtreeDataTuple>& heap, std::vector<Coordinate>& queryPointsX, 
                                              std::vector<Coordinate>& queryPointsY, VertexID& ann, EuclideanDistance& annDist, 
                                              AggregateFunction<EuclideanDistance>& agg, std::vector<EuclideanDistance>& queryDistances);

    private:
        friend class boost::serialization::access;
        
        std::vector<PackedRtreeNode> nodes;
        int branchFactor;
        std::string objSetType;
        double objSetDensity;
        int objSetVariable,  objSetSize;

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->nodes;
            ar & this->branchFactor;
            ar & this->objSetType;
            ar & this->objSetDensity;
            ar & this->objSetVariable;
            ar & this->objSetSize;
        }
};


#endif // _PACKEDRTREE_H