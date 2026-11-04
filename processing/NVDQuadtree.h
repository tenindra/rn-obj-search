// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _NVDQUADTREE_H
#define _NVDQUADTREE_H

#include "Graph.h"
#include "ALT.h"
#include "CompressibleMortonList.h"
#include "pruned_highway_labeling.h"
#include "../common.h"
#include "../queue/BinaryMinHeap.h"
#include "../utility/AggregateFunction.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <google/dense_hash_set>
#include <google/dense_hash_map>

class NVDQuadtree {

    public:
        NVDQuadtree(std::string setType, double setDensity, int setVariable, int setSize);
        NVDQuadtree() {};
        void buildNVD(Graph& graph, std::vector<VertexID>& objectIDs, std::vector<MortonCode>& mortonCodes, int maxLeafColours);
        double getObjSetDensity();
        int getObjSetSize();
        int getObjSetVariable();
        std::string getObjSetType();
        double computeIndexSize();

        VertexID getNearestNeighbourHashSet(MortonCode& queryVertexMC, google::dense_hash_set<VertexID>& neighbourVnsIDsAdded, VertexID& vnsID);
        int generateNeighbourCandidates(ALT& alt, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, VertexID queryVertexID, VertexID vnsID, 
                                         EdgeWeight spdist, google::dense_hash_set<VertexID>& neighbourVnsIDsAdded);

        int generateNeighbourCandidatesAkNN(ALT& alt, BinaryMinHeap<EdgeWeight,VertexPair>& pqueue, VertexID queryVertexID, VertexID vnsID, 
                                             EdgeWeight spdist, google::dense_hash_set<VertexID>& neighbourVnsIDsAdded, unsigned int queryVertexIdx);

        int generateNeighbourCandidatesRangeSearch(ALT& alt, BinaryMinHeap<EdgeWeight,VertexEdgeWeightPair>& pqueue, VertexID queryVertexID, VertexID vnsID, 
                                                   EdgeWeight spdist, google::dense_hash_set<VertexID>& neighbourVnsIDsAdded);

        std::vector<VertexID> vnsObjectIDs; // Object VertexID associated with dominance region
        std::vector<VertexDistancePair> neighbourIDsAndDist; // For each dominance region, list of neighbour dominance regions and dist between objects
        std::vector<VertexID> firstNeighbourIdx;
        CompressibleMortonList mortonList;
        
    private:
        friend class boost::serialization::access;
        
        std::string objSetType;
        double objSetDensity;
        int objSetVariable, objSetSize;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->objSetType;
            ar & this->objSetDensity;
            ar & this->objSetVariable;
            ar & this->objSetSize;
            ar & this->vnsObjectIDs;
            ar & this->neighbourIDsAndDist;
            ar & this->firstNeighbourIdx;
            ar & this->mortonList;
        }

};

#endif // _NVDQUADTREE_H