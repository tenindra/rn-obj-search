// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _ALT_H
#define _ALT_H

#include "Graph.h"
#include "../queue/BinaryMinHeap.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <google/dense_hash_set>
#include <random>

enum LANDMARK_TYPE
{
    RANDOM,
    FARTHEST,
    FARTHEST_BFS
};

class ALT {

    public:
        ALT(std::string networkName, int numVertices, int numEdges);
        ALT() {};
        void buildALT(Graph& graph, LANDMARK_TYPE landmarkType, unsigned int numLandmarks);
        std::vector<VertexID> selectFarthestLandmarks(Graph& graph, unsigned int numLandmarks, std::mt19937& engine);
        std::vector<VertexID> selectFarthestBFSLandmarks(Graph& graph, unsigned int numLandmarks, std::mt19937& engine);
        std::string getNetworkName();
        unsigned int getNumVertices();
        unsigned int getNumEdges();
        unsigned int getNumLandmarks();
        double computeIndexSize();
        EdgeWeight getLowerBound(VertexID s, VertexID t);
        EdgeWeight getUpperBound(VertexID s, VertexID t);
        VertexID getNearestLandmarkIndex(VertexID s, EdgeWeight& dist);
        void getLowerAndUpperBound(VertexID s, VertexID t, EdgeWeight& lb, EdgeWeight& ub);
        LANDMARK_TYPE getLandmarkType(std::string selectionMethod, bool& success);

        std::vector<VertexID> landmarks;
        std::vector<EdgeWeight> vertexFromLandmarkDistances;

        // Seed of the landmark selection (buildALT); any other value selects other landmarks
        static const unsigned int LANDMARK_SEED = 0;
        
    private:
        friend class boost::serialization::access;
        
        std::string networkName;
        unsigned int numVertices;
        unsigned int numEdges;
        unsigned int numLandmarks;
        
        // Non-Serialized Members
    
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & networkName;
            ar & numVertices;
            ar & numEdges;
            ar & numLandmarks;
            ar & landmarks;
            ar & vertexFromLandmarkDistances;
        }    
};

class ObjectList {

    public:
        ObjectList(std::string setType, double setDensity, int setVariable, int setSize);
        ObjectList() {};
        void buildObjectList(ALT& alt, std::vector<VertexID>& objectIDs);
        double getObjSetDensity();
        int getObjSetSize();
        int getObjSetVariable();
        std::string getObjSetType();
        double computeIndexSize();
        int getClosestDistanceIndex(std::vector<VertexDistancePair>& distances, EdgeWeight matchDist);
        BinaryMinHeap<EdgeWeight,VertexID> getNearestCandidateByQueryLandmark(ALT& alt, VertexID queryVertex, EdgeWeight queryVertexLMDistance, 
                                                                            int& leftPos, int& rightPos, VertexID& candidate, VertexID& candidateLMIndex);
        bool getNextNearestCandidateByQueryLandmark(ALT& alt, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, VertexID queryVertex, 
                                                    EdgeWeight queryVertexLMDistance, int& leftPos, int& rightPos, VertexID& candidate, 
                                                    VertexID& candidateLMIndex, EdgeWeight& candidateLBDist);
        void fullyExploreSingleObjectList(VertexID landmarkIndex, ALT& alt, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, VertexID queryVertex, 
                                          EdgeWeight queryVertexLMDistance, int& leftPos, int& rightPos);

        std::vector<std::vector<VertexDistancePair>> lmObjDistances;
        
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
            ar & this->lmObjDistances;
        }

};

#endif // _ALT_H
