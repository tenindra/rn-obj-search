// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _COMPRESSIBLEMORTONLIST_H
#define _COMPRESSIBLEMORTONLIST_H

#include "Graph.h"
#include "CompressibleQuadtree.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

class CompressibleMortonBlock {

    public:
        CompressibleMortonBlock() {};
        CompressibleMortonBlock(MortonCode code, VertexID objectID);
        VertexID getColourID();
        void printBlock();
        MortonNumber getMortonNumber();
        MortonLength getMortonLength();
        int compareMortonCode(const MortonCode& code) const;
        MortonCode getMortonCode();
        void decompose(int maxRangeExp, Coordinate& x, Coordinate& y, int& width);
        bool operator==(const MortonCode& code) const;
        double computeIndexSizeBytes();
        
    private:
        friend class boost::serialization::access;
        
        MortonNumber zNumber;
        MortonLength len; // zNumber length (i.e. used bits)
        // Note: len actually refers to the number of useful bits
        // so when len is 0 the block contains all points
        VertexID colourID;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->zNumber;
            ar & this->len;
            ar & this->colourID;
        }        

};

class CompressibleMortonList {

    public:
        CompressibleMortonList() {};
        void buildCompressibleMortonList(CompressibleQuadtree& qt);
        VertexID getColourID(const MortonCode& targetCode);
        std::size_t getBlockIdxContainingTargetByBST(const MortonCode& targetCode);
        std::size_t getLeftChild(std::size_t parent);
        std::size_t getRightChild(std::size_t parent);
        bool isTreeSorted(int pos);
        bool areSubtreesValid(int pos);
        bool isLeftSubtreeValid(int leftChild, MortonCode& parentCode);
        bool isRightSubtreeValid(int rightChild, MortonCode& parentCode);
        int getNumBlocks();
        void printBST();
        double computeIndexSizeBytes();
        void clear();
    private:
        friend class boost::serialization::access;
        
        std::vector<CompressibleMortonBlock> binarySearchTreeCompressibleMortonList;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->binarySearchTreeCompressibleMortonList;
        }

};


#endif // _COMPRESSIBLEMORTONLIST_H