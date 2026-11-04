// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "ILBR.h"

#include "../queue/BinaryMaxHeap.h"
#include "../queue/BinaryMinHeap.h"

#include <algorithm>
#include <set>
#include <google/dense_hash_map>
#include <iostream>
#include <stack>

void ILBR::getKNNsByPHLNaive(ALT& alt, std::vector<VertexID>& objects, unsigned int k, VertexID queryVertexID, 
                        std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances, PrunedHighwayLabeling& phl)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
#endif

    // Compute lowerbounds for all objects and create sorted list
    BinaryMinHeap<EdgeWeight,VertexID> candidateQueue;
    VertexID candidate;
    EdgeWeight lowerBound, spDist, Dk;
    for (std::size_t i = 0; i < objects.size(); ++i) {
        lowerBound = alt.getLowerBound(queryVertexID,objects[i]);
        candidateQueue.insert(objects[i],lowerBound);
    }
    // We compute the network distances to k best candidates
    BinaryMaxHeap<EdgeWeight,VertexID> knnCandidates;
    for (std::size_t i = 0; i < k && candidateQueue.size() > 0; ++i) {
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("spdist_calls",1);
#endif
        candidate = candidateQueue.extractMinElement();
        spDist = phl.Query(queryVertexID,candidate);
        knnCandidates.insert(candidate,spDist);
    }
    
    // While the lowerbound distance to the next closest candidate
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    while (candidateQueue.size() > 0) {
        lowerBound = candidateQueue.getMinKey();
        if (lowerBound < Dk) {
#if defined(COLLECT_STATISTICS)
            // This means some candidate (existing or new) was a poor choice
            this->stats.incrementStatistic("false_hit",1);
            this->stats.incrementStatistic("spdist_calls",1);
#endif
            candidate = candidateQueue.extractMinElement();
            spDist = phl.Query(queryVertexID,candidate);
            if (spDist < Dk) {
#if defined(COLLECT_STATISTICS)
                // This means some old candidate was a poor choice
                this->stats.incrementStatistic("false_old_candidate",1);
#endif
                // Only insert if it is a better kNN candidate
                knnCandidates.insert(candidate,spDist);
                knnCandidates.extractMaxElement();
                Dk = knnCandidates.getMaxKey();
            }
#if defined(COLLECT_STATISTICS)
            // This means the new candidate was the bad choice
              else {
                this->stats.incrementStatistic("false_new_candidate",1);
            }
#endif
        } else {
            break;
        }
    }

    knnCandidates.populateKNNs(kNNs,kNNDistances);
}

void ILBR::getKNNsByPHLImprovedByQueryLandmark(ALT& alt, ObjectList& ol, unsigned int k, VertexID queryVertexID, std::vector<VertexID>& kNNs, 
                                               std::vector<EdgeWeight>& kNNDistances, PrunedHighwayLabeling& phl)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
    this->stats.initialiseStatistic("lbs_computed",0);
#endif

    // Compute lowerbounds for all objects and create sorted list
    EdgeWeight queryLMDistance;
    int leftPos, rightPos;
    VertexID candidate, candidateLMIndex;
    EdgeWeight candidateLBDist;
    
    candidateLMIndex = alt.getNearestLandmarkIndex(queryVertexID,queryLMDistance);

    // Retrieve k candidates by landmark lower bounds and compute distance to each of them
    BinaryMaxHeap<EdgeWeight,VertexID> knnCandidates;
    EdgeWeight spDist, Dk;
    bool runAlgorithm = true;
    
    BinaryMinHeap<EdgeWeight,VertexID> candidateQueue = ol.getNearestCandidateByQueryLandmark(alt,queryVertexID,queryLMDistance,leftPos,rightPos,candidate,candidateLMIndex);
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("spdist_calls",1);
    this->stats.incrementStatistic("lbs_computed",1);
#endif
    spDist = phl.Query(queryVertexID,candidate);
    knnCandidates.insert(candidate,spDist);
    for (unsigned int i = 1; i < k; ++i) {
        // Since we have retrieved first candidate i starts at 1 abgove
        if (ol.getNextNearestCandidateByQueryLandmark(alt,candidateQueue,queryVertexID,queryLMDistance,leftPos,rightPos,candidate,candidateLMIndex,candidateLBDist)) {
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("spdist_calls",1);
            this->stats.incrementStatistic("lbs_computed",1);
#endif
            spDist = phl.Query(queryVertexID,candidate);
            knnCandidates.insert(candidate,spDist);
        } else {
            runAlgorithm = false;
        }
    }
    
    // While the lowerbound distance to the next nearest candidate
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    while (runAlgorithm) {
        if (ol.getNextNearestCandidateByQueryLandmark(alt,candidateQueue,queryVertexID,queryLMDistance,leftPos,rightPos,candidate,candidateLMIndex,candidateLBDist)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we do not reach here
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("lbs_computed",1);
#endif
            if (candidateLBDist < Dk) {
#if defined(COLLECT_STATISTICS)
                // This means some candidate (existing or new) was a poor choice
                this->stats.incrementStatistic("false_hit",1);
                this->stats.incrementStatistic("spdist_calls",1);
#endif
                spDist = phl.Query(queryVertexID,candidate);
                if (spDist < Dk) {
#if defined(COLLECT_STATISTICS)
                    // This means some old candidate was a poor choice
                    this->stats.incrementStatistic("false_old_candidate",1);
#endif
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(candidate,spDist);
                    knnCandidates.extractMaxElement();
                    Dk = knnCandidates.getMaxKey();
                }
#if defined(COLLECT_STATISTICS)
                // This means the new candidate was the bad choice
                  else {
                    this->stats.incrementStatistic("false_new_candidate",1);
                }
#endif
            } else {
                break;
            }
        } else {
            // This mean no nearest neighbours were found (we have reported
            // all objects) so we can stop the search
            break;
        }
    }

#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("lbs_computed",static_cast<int>(candidateQueue.size()));
#endif
    knnCandidates.populateKNNs(kNNs,kNNDistances);
}

void ILBR::getKNNByNVDQuadtreePHL(ALT& alt, NVDQuadtree& nvd, unsigned int k, VertexID queryVertexID, std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances, PrunedHighwayLabeling& phl, Graph& graph)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("cand_unpruned",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
    this->stats.initialiseStatistic("lbs_computed",0);
#endif

    // Retrieve k candidates by landmark lower bounds and compute distance to each of them
    VertexID candidate;
    google::dense_hash_set<VertexID> neighboursAdded;
    neighboursAdded.set_empty_key(constants::UNUSED_VERTEX_ID);
    BinaryMaxHeap<EdgeWeight,VertexID> knnCandidates;
    EdgeWeight spDist, Dk = 0, candidateLBDist;
    bool DkInfinity = true;
    [[maybe_unused]] int numCandidates;
    
    BinaryMinHeap<EdgeWeight,VertexID> nvdQueue; // VertexID element is vnsID (NOT objectID)
    VertexID vnsID;    
    candidate = nvd.getNearestNeighbourHashSet(graph.mortonCodes[queryVertexID],neighboursAdded,vnsID);
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("cand_spdist",1); // Note: This is not considered unpruned because we know it is a result (k >= 1)
    this->stats.incrementStatistic("spdist_calls",1);
#endif
    spDist = phl.Query(queryVertexID,candidate);
    knnCandidates.insert(candidate,spDist);
    if (knnCandidates.size() != k) {
        numCandidates = nvd.generateNeighbourCandidates(alt,nvdQueue,queryVertexID,vnsID,spDist,neighboursAdded);
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("cand_unpruned",numCandidates);
        this->stats.incrementStatistic("lbs_computed",numCandidates);
#endif
    }
    // If k=1 then we don't need to generate anymore candidates as
    // as the NVD always return correct 1NN (following loop will not 
    // run as nvdQueue.size() is 0
    
    // While the lowerbound distance to the next nearest candidate
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    while (nvdQueue.size() > 0 && (DkInfinity || nvdQueue.getMinKey() < Dk)) {
        candidateLBDist = nvdQueue.getMinKey();
        vnsID = nvdQueue.extractMinElement();
        if (DkInfinity) {
#if defined(COLLECT_STATISTICS)
            // This means some candidate (existing or new) was a poor choice
            this->stats.incrementStatistic("cand_spdist",1);
            this->stats.incrementStatistic("spdist_calls",1);
#endif
            candidate = nvd.vnsObjectIDs[vnsID];
            spDist = phl.Query(queryVertexID,candidate);
            knnCandidates.insert(candidate,spDist);
            numCandidates = nvd.generateNeighbourCandidates(alt,nvdQueue,queryVertexID,vnsID,spDist,neighboursAdded);
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("cand_unpruned",numCandidates);
            this->stats.incrementStatistic("lbs_computed",numCandidates);
#endif
            if (knnCandidates.size() == k) {
                DkInfinity = false;
                Dk = knnCandidates.getMaxKey();
            }            
        } else if (candidateLBDist < Dk) {
#if defined(COLLECT_STATISTICS)
            // This means some candidate (existing or new) was a poor choice
            this->stats.incrementStatistic("cand_spdist",1);
            this->stats.incrementStatistic("false_hit",1);
            this->stats.incrementStatistic("spdist_calls",1);
#endif
            candidate = nvd.vnsObjectIDs[vnsID];
            spDist = phl.Query(queryVertexID,candidate);
            if (spDist < Dk) {
#if defined(COLLECT_STATISTICS)
                // This means some old candidate was a poor choice
                this->stats.incrementStatistic("false_old_candidate",1);
#endif
                // Only insert if it is a better kNN candidate
                knnCandidates.insert(candidate,spDist);
                numCandidates = nvd.generateNeighbourCandidates(alt,nvdQueue,queryVertexID,vnsID,spDist,neighboursAdded);
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("cand_unpruned",numCandidates);
                this->stats.incrementStatistic("lbs_computed",numCandidates);
#endif
                // Note: We only generate candidates from candidates that can be kNNs
                // any neighbour objects from this must necessarily be further from 
                // the query vertex so cannot be objects
                knnCandidates.extractMaxElement();
                Dk = knnCandidates.getMaxKey();
            }
#if defined(COLLECT_STATISTICS)
            // This means the new candidate was the bad choice
                else {
                this->stats.incrementStatistic("false_new_candidate",1);
            }
#endif
        }
    }

    knnCandidates.populateKNNs(kNNs,kNNDistances);
}

void ILBR::getAkNNByNVDQuadtreePHLByConcExps(ALT& alt, NVDQuadtree& nvd, unsigned int k, std::vector<VertexID>& queryVertices, 
                                   std::vector<VertexID>& AkNNs, std::vector<EdgeWeight>& AkNNDistances, 
                                   PrunedHighwayLabeling& phl, Graph& graph, AggregateFunction<EdgeWeight>& agg)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("cand_unpruned",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
    this->stats.initialiseStatistic("lbs_computed",0);
#endif

    // Retrieve k candidates by landmark lower bounds and compute distance to each of them
    VertexID candidate, queryVertexIdx;
    std::vector<google::dense_hash_set<VertexID>> neighboursAdded(queryVertices.size());
    BinaryMaxHeap<EdgeWeight,VertexID> knnCandidates;
    EdgeWeight Dk = 0, candidateLBDist, aggScore, spDist;
    bool DkInfinity = true;
    std::vector<EdgeWeight> queryDistances(queryVertices.size());
    [[maybe_unused]] int numCandidates;
    
    BinaryMinHeap<EdgeWeight,VertexPair> nvdQueue; // VertexID element is vnsID (NOT objectID)
    VertexID vnsID;
    for (std::size_t i = 0; i < queryVertices.size(); ++i) {
        neighboursAdded[i].set_empty_key(constants::UNUSED_VERTEX_ID);
        candidate = nvd.getNearestNeighbourHashSet(graph.mortonCodes[queryVertices[i]],neighboursAdded[i],vnsID);
        candidateLBDist = alt.getLowerBound(queryVertices[i],candidate);
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("cand_unpruned",1);
        this->stats.incrementStatistic("lbs_computed",1);
#endif
        nvdQueue.insert(VertexPair(vnsID,i),candidateLBDist);
    }

    google::dense_hash_map<VertexID,VertexID> settledObjects;
    settledObjects.set_empty_key(constants::UNUSED_VERTEX_ID);
    settledObjects.set_deleted_key(constants::UNUSED_VERTEX_ID+1); 
    
    int settledCount = 0;
    // Note: When an object is settled by a search (i.e. popped off the queue), we increment
    // settledCount. When that object has been settled from all searches from queries query vertices
    // we decrease settledCount by the number of queryVertices. This way, we can capture any objects
    // that are settled by a search and may still be results, stopping us from terminating
    while (nvdQueue.size() > 0 && (DkInfinity || agg.getLowerBound(nvdQueue.getMinKey(),queryVertices.size()) < Dk)) {
        // Note: We terminate when no object in nvdQueue can improve on Dk and there are not objects settled by
        // one search, but not all the searches. This handles the case where an object that is yet to be settled from
        // one search and is a result, but the lower-bounds from other searches are greater than Dk.
        
        spDist = nvdQueue.getMinKey(); // Initialize as input to NVD expansion (in case we don't compute true spDist)
        VertexPair minElement = nvdQueue.extractMinElement();
        vnsID = minElement.first;
        queryVertexIdx = minElement.second;
        
        if (settledObjects.find(vnsID) == settledObjects.end()) {
            settledObjects[vnsID] = 1;
        } else {
            ++settledObjects[vnsID];
        }
        ++settledCount;
        
        if (settledObjects[vnsID] == queryVertices.size()) {
            // Candidate has been settled from all searches from Vns associated with 
            candidate = nvd.vnsObjectIDs[vnsID];
            settledObjects.erase(vnsID);
            settledCount = settledCount-queryVertices.size();
            
            if (DkInfinity) {
                for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                    queryDistances[i] = phl.Query(queryVertices[i],candidate);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("cand_spdist",1);
                    this->stats.incrementStatistic("spdist_calls",1);
#endif
                }
                aggScore = agg.aggregate(queryDistances);
                spDist = queryDistances[queryVertexIdx];
                knnCandidates.insert(candidate,aggScore);
                if (knnCandidates.size() == k) {
                    DkInfinity = false;
                    Dk = knnCandidates.getMaxKey();
                }                
            } else {
                // Compute lower-bound score first to filter bad candidates
                for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                    queryDistances[i] = alt.getLowerBound(queryVertices[i],candidate);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("lbs_computed",1);
#endif
                }
                aggScore = agg.aggregate(queryDistances);                
                if (aggScore < Dk) {
                    for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                        queryDistances[i] = phl.Query(queryVertices[i],candidate);
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("cand_spdist",1);
                        this->stats.incrementStatistic("spdist_calls",1);
#endif
                    }
                    aggScore = agg.aggregate(queryDistances);
                    spDist = queryDistances[queryVertexIdx];
                    
                    if (aggScore < Dk) {
#if defined(COLLECT_STATISTICS)
                        // This means some old candidate was a poor choice
                        this->stats.incrementStatistic("false_hit",1);
                        this->stats.incrementStatistic("false_old_candidate",1);
#endif
                        knnCandidates.insert(candidate,aggScore);
                        knnCandidates.extractMaxElement();
                        Dk = knnCandidates.getMaxKey();
                    }
#if defined(COLLECT_STATISTICS)
                        // This means the new candidate was the bad choice
                        else {
                        this->stats.incrementStatistic("false_hit",1);
                        this->stats.incrementStatistic("false_new_candidate",1);
                    }
#endif
                }
            }
        }

        numCandidates = nvd.generateNeighbourCandidatesAkNN(alt,nvdQueue,queryVertices[queryVertexIdx],vnsID,spDist,neighboursAdded[queryVertexIdx],queryVertexIdx);
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("cand_unpruned",numCandidates);
        this->stats.incrementStatistic("lbs_computed",numCandidates);
#endif
    }

    if (settledCount != 0) {
        // Note: nvdQueue.size() != 0 here, otherwise settledCount would also be zero (because all searches visited all Vns')
        // Note: DkInfinity must be false here, if DkInfinity was true then nvdQueue.size() == 0 (which it cannot be)
        for (auto it = settledObjects.begin(); it != settledObjects.end(); ++it) {
            candidate = nvd.vnsObjectIDs[it->first];
            // Compute lower-bound score first to filter bad candidates
            for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                queryDistances[i] = alt.getLowerBound(queryVertices[i],candidate);
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("lbs_computed",1);
#endif
            }
            aggScore = agg.aggregate(queryDistances);            
            if (aggScore < Dk) {
                for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                    queryDistances[i] = phl.Query(queryVertices[i],candidate);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("cand_spdist",1);
                    this->stats.incrementStatistic("spdist_calls",1);
#endif
                }
                aggScore = agg.aggregate(queryDistances);
                if (aggScore < Dk) {
#if defined(COLLECT_STATISTICS)
                    // This means some old candidate was a poor choice
                    this->stats.incrementStatistic("false_hit",1);
                    this->stats.incrementStatistic("false_old_candidate",1);
#endif
                    knnCandidates.insert(candidate,aggScore);
                    knnCandidates.extractMaxElement();
                    Dk = knnCandidates.getMaxKey();
                }
#if defined(COLLECT_STATISTICS)
                // This means the new candidate was the bad choice
                    else {
                    this->stats.incrementStatistic("false_hit",1);
                    this->stats.incrementStatistic("false_new_candidate",1);
                }
#endif
            }
        }
    }


#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("lbs_computed",static_cast<int>(nvdQueue.size()));
#endif
    
    knnCandidates.populateKNNs(AkNNs,AkNNDistances);
}

// Version where object candidates are inserted into pqueue by their ALT lower-bound
void ILBR::getKNNByCOLTPHL(ALT& alt, COLTree& colt, SULTree& sltree, unsigned int k, VertexID queryVertexID, 
                           std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances, PrunedHighwayLabeling& phl) 
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("cand_unpruned",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
    this->stats.initialiseStatistic("lbs_computed",0);
    this->stats.initialiseStatistic("lbs_computed_node",0);
    this->stats.initialiseStatistic("lbs_computed_obj",0);
    this->stats.initialiseStatistic("ols_started",0);
    this->stats.initialiseStatistic("ols_finished",0);
    this->stats.initialiseStatistic("nodes_total",colt.getNumNonEmptyNodes());
    this->stats.initialiseStatistic("nodes_inserted",0);
    this->stats.initialiseStatistic("nodes_nonleaf_dequeued",0);
    this->stats.initialiseStatistic("nodes_leaf_dequeued_first",0);
    this->stats.initialiseStatistic("nodes_leaf_dequeued_repeat",0);
    this->stats.initialiseStatistic("nodes_lbzero",0);
    this->stats.initialiseStatistic("obj_inserted",0);
    this->stats.initialiseStatistic("obj_dequeued",0);
    this->stats.initialiseStatistic("alt_used",0);
    this->stats.initialiseStatistic("alt_tried",0);
    this->stats.initialiseStatistic("global_landmark_used",0);
    this->stats.initialiseStatistic("global_landmark_attempted",0);
    bool collectLevelStats = false;
    if (collectLevelStats) {
        int treeDepth = colt.getTreeDepth();
        this->stats.initialiseStatistic("y_total_levels",treeDepth);
        for (int i = 0; i < treeDepth; ++i) {
            this->stats.initialiseStatistic("z_nodes_inserted_level" + std::to_string(i),0);
            this->stats.initialiseStatistic("z_nodes_nonleaf_dequeued_level" + std::to_string(i),0);
            this->stats.initialiseStatistic("z_nodes_leaf_dequeued_first_level" + std::to_string(i),0);
            this->stats.initialiseStatistic("z_global_landmark_used_level" + std::to_string(i),0);
        }
    }
#endif

    BinaryMaxHeap<EdgeWeight,VertexID> knnCandidates;
    BinaryMinHeap<EdgeWeight,VertexID> pqueue;
    
    VertexID objectID;
    TreeNodeID nodeIdx, childIdx;
    EdgeWeight Dk = 0, spDist, lbDist;
    bool DkInfinity = true;
    
    // Note: In order to use a priority queue with primitive types we
    // map the tree indexes as VertexID greater than the largest VertexID
    // for the current graph. Whenever we deque we can tell if the element
    // is a tree node or graph vertex by check if it's larger than the
    // number of vertices in the graph (since they are number 0 to n-1)
    VertexID firstMappedTreeNodeID = static_cast<VertexID>(alt.getNumVertices());

    pqueue.insert(0+firstMappedTreeNodeID,0); // Start with root node in COLTree (with index 0)
    colt.queryData[0].reset();
    if (colt.treeNodes[0].hasObjectList()) {
        // The first a OL node is dequeued, it assumes that lmIdx has been set earlier
        // so if root has an OL (because it is the only node), we need to set it here
        colt.queryData[0].lmIdx = colt.treeNodes[0].getFurthestLmkIdxByLB(alt,queryVertexID);
#if defined(COLLECT_STATISTICS)
        // getFurthestLmkIdxByLB computes a lower-bound from query vertex to each landmark to get best landmark
        this->stats.incrementStatistic("lbs_computed",colt.treeNodes[0].landmarks.size());
        this->stats.incrementStatistic("lbs_computed_node",colt.treeNodes[0].landmarks.size());
#endif
    }

    while (pqueue.size() > 0 && (DkInfinity || pqueue.getMinKey() < Dk)) {
        objectID = pqueue.extractMinElement();
        if (objectID < firstMappedTreeNodeID) {
            // Element is an object
            spDist = phl.Query(queryVertexID,objectID);
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("cand_spdist",1);
            this->stats.incrementStatistic("obj_dequeued",1);
            this->stats.incrementStatistic("spdist_calls",1);
#endif
            if (DkInfinity) {
                knnCandidates.insert(objectID,spDist);
                if (knnCandidates.size() == k) {
                    DkInfinity = false;
                    Dk = knnCandidates.getMaxKey();
                }            
            } else if (spDist < Dk) {
                knnCandidates.insert(objectID,spDist);
                knnCandidates.extractMaxElement();
                Dk = knnCandidates.getMaxKey();
#if defined(COLLECT_STATISTICS)
                // This means some old candidate was a poor choice
                this->stats.incrementStatistic("false_hit",1);
                this->stats.incrementStatistic("false_old_candidate",1);
#endif
            }
#if defined(COLLECT_STATISTICS)
            // This means the new candidate was the bad choice
            else {
                this->stats.incrementStatistic("false_hit",1);
                this->stats.incrementStatistic("false_new_candidate",1);
            }
#endif
        } else {
            nodeIdx = objectID-firstMappedTreeNodeID;
            if (!colt.treeNodes[nodeIdx].hasObjectList()) {
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("nodes_nonleaf_dequeued",1);
                if (collectLevelStats) {
                    this->stats.incrementStatistic("z_nodes_nonleaf_dequeued_level" + std::to_string(colt.getNodeDepth(nodeIdx)),1);
                }
#endif
                // No OL, so we expand search towards children using COLTree
                // Note: These will not include children with no objects
                for (std::size_t i = 0; i < colt.treeNodes[nodeIdx].childNodeIdxs.size(); ++i) {
                    childIdx = colt.treeNodes[nodeIdx].childNodeIdxs[i];
                    RN_DEBUG_ASSERT(colt.treeNodes[childIdx].getNumObjects() != 0 && "A COL-Tree child node holds at least one object");
                    colt.queryData[childIdx].reset();
                    if (!colt.treeNodes[childIdx].hasObjectList()) {
                        lbDist = colt.treeNodes[childIdx].getBestNodeLowerbound(alt,queryVertexID);
                    } else {
                        // If this is a node with OL, save furthest landmark simultaneously.
                        // (This will reduce the number of duplicate ALT calls when it is searched)
                        int lmIdx;
                        lbDist = colt.treeNodes[childIdx].getBestNodeLowerboundAndFurthestLmk(alt,queryVertexID,lmIdx);
                        colt.queryData[childIdx].lmIdx = lmIdx;
                    }
                    EdgeWeight altLBDist = colt.computeGlobalLowerbound(childIdx,queryVertexID,sltree);
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("global_landmark_attempted",1);
#endif
                    if (altLBDist > lbDist) {
                        lbDist = altLBDist;
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("global_landmark_used",1);
                        if (collectLevelStats) {
                            this->stats.incrementStatistic("z_global_landmark_used_level" + std::to_string(colt.getNodeDepth(childIdx)),1);
                        }
#endif
                    } 
                    pqueue.insert(childIdx+firstMappedTreeNodeID,lbDist);
#if defined(COLLECT_STATISTICS)
                    // getBestNodeLowerbound computes a lower-bound from query vertex to each landmark to get best lower-bound for node
                    this->stats.incrementStatistic("lbs_computed",colt.treeNodes[childIdx].landmarks.size());
                    this->stats.incrementStatistic("lbs_computed_node",colt.treeNodes[childIdx].landmarks.size());
                    this->stats.incrementStatistic("nodes_inserted",1);
                    if (collectLevelStats) {
                        this->stats.incrementStatistic("z_nodes_inserted_level" + std::to_string(colt.getNodeDepth(childIdx)),1);
                    }
                    if (lbDist == 0) {
                        this->stats.incrementStatistic("nodes_lbzero",1);
                    }
#endif
                }
            } else {
                // Handle depending on whether OL for this node has been searched before
                if (!colt.queryData[nodeIdx].intialized) {
                    colt.queryData[nodeIdx].intialized = true;
                    colt.queryData[nodeIdx].lmkDist = phl.Query(queryVertexID,colt.treeNodes[nodeIdx].landmarks[colt.queryData[nodeIdx].lmIdx]);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("nodes_leaf_dequeued_first",1);
                    if (collectLevelStats) {
                        this->stats.incrementStatistic("z_nodes_leaf_dequeued_first_level" + std::to_string(colt.getNodeDepth(nodeIdx)),1);
                    }
                    this->stats.incrementStatistic("spdist_calls",1);
                    this->stats.incrementStatistic("ols_started",1);
#endif
                    
                    // Find the object in this landmark's object list and udpate result set with it
                    // and also initializes left and right indexes from this object
#if defined(COLLECT_STATISTICS)
                    colt.treeNodes[nodeIdx].updateQueueWithMinObj(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].leftIdx,colt.queryData[nodeIdx].rightIdx,colt.queryData[nodeIdx].lmkDist,queryVertexID,pqueue,alt,this->stats);
#else
                    colt.treeNodes[nodeIdx].updateQueueWithMinObj(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].leftIdx,colt.queryData[nodeIdx].rightIdx,colt.queryData[nodeIdx].lmkDist,queryVertexID,pqueue,alt);
#endif
                }
#if defined(COLLECT_STATISTICS)
                else {
                    stats.incrementStatistic("nodes_leaf_dequeued_repeat",1);
                }
#endif

                // Update queue by searching left and right from current left and right indexes 
                // in object list of chosen landmark (terminate when object lower-bound is higher than
                // priority queue lower-bound - i.e. there maybe a bette object in another obj list)
#if defined(COLLECT_STATISTICS)
                if (colt.treeNodes[nodeIdx].updateQueue(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].leftIdx,colt.queryData[nodeIdx].rightIdx,colt.queryData[nodeIdx].lmkDist,queryVertexID,pqueue,alt,lbDist,this->stats)) {
#else
                if (colt.treeNodes[nodeIdx].updateQueue(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].leftIdx,colt.queryData[nodeIdx].rightIdx,colt.queryData[nodeIdx].lmkDist,queryVertexID,pqueue,alt,lbDist)) {
#endif
                    // We insert objectID because that's the mapped VertexID for nodeIdx
                    pqueue.insert(objectID,lbDist);
                }
#if defined(COLLECT_STATISTICS)
                else {
                    stats.incrementStatistic("ols_finished",1);
                }
#endif
            }
        }
    }

    knnCandidates.populateKNNs(kNNs,kNNDistances);
}

void ILBR::getAkNNByCOLTPHL(ALT& alt, COLTree& colt, SULTree& sltree, unsigned int k, std::vector<VertexID>& queryVertices, 
                            std::vector<VertexID>& AkNNs, std::vector<EdgeWeight>& AkNNDistances, PrunedHighwayLabeling& phl, 
                            AggregateFunction<EdgeWeight>& agg) 
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("cand_unpruned",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
    this->stats.initialiseStatistic("lbs_computed",0);
    this->stats.initialiseStatistic("lbs_computed_node",0);
    this->stats.initialiseStatistic("lbs_computed_obj",0);
    this->stats.initialiseStatistic("ols_started",0);
    this->stats.initialiseStatistic("ols_finished",0);
    this->stats.initialiseStatistic("nodes_total",colt.getNumNonEmptyNodes());
    this->stats.initialiseStatistic("nodes_inserted",0);
    this->stats.initialiseStatistic("nodes_nonleaf_dequeued",0);
    this->stats.initialiseStatistic("nodes_leaf_dequeued_first",0);
    this->stats.initialiseStatistic("nodes_leaf_dequeued_repeat",0);
    this->stats.initialiseStatistic("nodes_lbzero",0);
    this->stats.initialiseStatistic("obj_inserted",0);
    this->stats.initialiseStatistic("obj_dequeued",0);
    this->stats.initialiseStatistic("alt_used",0);
    this->stats.initialiseStatistic("alt_tried",0);
    this->stats.initialiseStatistic("global_landmark_used",0);
    this->stats.initialiseStatistic("global_landmark_attempted",0);
    bool collectLevelStats = false;
    if (collectLevelStats) {
        int treeDepth = colt.getTreeDepth();
        this->stats.initialiseStatistic("y_total_levels",treeDepth);
        for (int i = 0; i < treeDepth; ++i) {
            this->stats.initialiseStatistic("z_nodes_inserted_level" + std::to_string(i),0);
            this->stats.initialiseStatistic("z_nodes_nonleaf_dequeued_level" + std::to_string(i),0);
            this->stats.initialiseStatistic("z_nodes_leaf_dequeued_first_level" + std::to_string(i),0);
            this->stats.initialiseStatistic("z_global_landmark_used_level" + std::to_string(i),0);
        }
    }
#endif

    BinaryMaxHeap<EdgeWeight,VertexID> knnCandidates;
    BinaryMinHeap<EdgeWeight,VertexID> pqueue;
    std::vector<EdgeWeight> queryDistances(queryVertices.size(),0);
    std::vector<std::vector<EdgeWeight>> queryVertexLmkDists; // Contains distances from query vertices to an OL node's chosen landmark
    
    VertexID objectID;
    TreeNodeID nodeIdx, childIdx;
    EdgeWeight Dk = 0, lbScore, aggScore;
    bool DkInfinity = true, improvedByGlobalLB = false;
    
    // Note: In order to use a priority queue with primitive types we
    // map the tree indexes as VertexID greater than the largest VertexID
    // for the current graph. Whenever we deque we can tell if the element
    // is a tree node or graph vertex by check if it's larger than the
    // number of vertices in the graph (since they are number 0 to n-1)
    VertexID firstMappedTreeNodeID = static_cast<VertexID>(alt.getNumVertices());

    pqueue.insert(0+firstMappedTreeNodeID,0); // Start with root node in COLTree (with index 0)
    colt.queryData[0].reset();
    if (colt.treeNodes[0].hasObjectList()) {
        // The first a OL node is dequeued, it assumes that lmIdx has been set earlier
        // so if root has an OL (because it is the only node), we need to set it here
        colt.queryData[0].lmIdx = colt.treeNodes[0].getFurthestLmkIdxByLB(alt,queryVertices);
#if defined(COLLECT_STATISTICS)
        // getFurthestLmkIdxByLB computes a lower-bound from each query vertex to each landmark to get best landmark
        this->stats.incrementStatistic("lbs_computed",colt.treeNodes[0].landmarks.size());
        this->stats.incrementStatistic("lbs_computed_node",colt.treeNodes[0].landmarks.size());
#endif
    }

    while (pqueue.size() > 0 && (DkInfinity || pqueue.getMinKey() < Dk)) {
        objectID = pqueue.extractMinElement();
        if (objectID < firstMappedTreeNodeID) {
            // Element is an object
            for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                queryDistances[i] = phl.Query(queryVertices[i],objectID);
            }
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("cand_spdist",queryVertices.size());
            this->stats.incrementStatistic("obj_dequeued",1);
            this->stats.incrementStatistic("spdist_calls",queryVertices.size());
#endif
            aggScore = agg.aggregate(queryDistances);
            if (DkInfinity) {
                knnCandidates.insert(objectID,aggScore);
                if (knnCandidates.size() == k) {
                    DkInfinity = false;
                    Dk = knnCandidates.getMaxKey();
                }            
            } else if (aggScore < Dk) {
                knnCandidates.insert(objectID,aggScore);
                knnCandidates.extractMaxElement();
                Dk = knnCandidates.getMaxKey();
#if defined(COLLECT_STATISTICS)
                // This means some old candidate was a poor choice
                this->stats.incrementStatistic("false_hit",1);
                this->stats.incrementStatistic("false_old_candidate",1);
#endif
            }
#if defined(COLLECT_STATISTICS)
            // This means the new candidate was the bad choice
            else {
                this->stats.incrementStatistic("false_hit",1);
                this->stats.incrementStatistic("false_new_candidate",1);
            }
#endif
        } else {
            nodeIdx = objectID-firstMappedTreeNodeID;
            if (!colt.treeNodes[nodeIdx].hasObjectList()) {
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("nodes_nonleaf_dequeued",1);
                if (collectLevelStats) {
                    this->stats.incrementStatistic("z_nodes_nonleaf_dequeued_level" + std::to_string(colt.getNodeDepth(nodeIdx)),1);
                }
#endif
                // No OL, so we expand search towards children using COLTree
                // Note: These will not include children with no objects
                for (std::size_t i = 0; i < colt.treeNodes[nodeIdx].childNodeIdxs.size(); ++i) {
                    childIdx = colt.treeNodes[nodeIdx].childNodeIdxs[i];
                    RN_DEBUG_ASSERT(colt.treeNodes[childIdx].getNumObjects() != 0 && "A COL-Tree child node holds at least one object");
                    colt.queryData[childIdx].reset();
                    if (!colt.treeNodes[childIdx].hasObjectList()) {
                        colt.treeNodes[childIdx].getBestNodeLowerbounds(alt,queryVertices,queryDistances);
                    } else {
                        // If this is a node with OL, save furthest landmark simultaneously.
                        // (This will reduce the number of duplicate ALT calls when it is searched)
                        colt.treeNodes[childIdx].getBestNodeLowerboundsAndFurthestSumLmk(alt,queryVertices,queryDistances,colt.queryData[childIdx].lmIdx);
                    }
#if defined(COLLECT_STATISTICS)
                    improvedByGlobalLB = false;
#endif
                    colt.computeGlobalAggregateLowerbound(childIdx,queryVertices,queryDistances,sltree,improvedByGlobalLB);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("global_landmark_attempted",1);
                    if (improvedByGlobalLB) {
                        this->stats.incrementStatistic("global_landmark_used",1);
                        if (collectLevelStats) {
                            this->stats.incrementStatistic("z_global_landmark_used_level" + std::to_string(colt.getNodeDepth(childIdx)),1);
                        }
                    }
#endif
                    aggScore = agg.aggregate(queryDistances);
                    pqueue.insert(childIdx+firstMappedTreeNodeID,aggScore);
#if defined(COLLECT_STATISTICS)
                    // getBestNodeLowerbounds computes a lower-bound from query vertex to each landmark to get best lower-bound for node
                    this->stats.incrementStatistic("lbs_computed",colt.treeNodes[childIdx].landmarks.size()*queryVertices.size());
                    this->stats.incrementStatistic("lbs_computed_node",colt.treeNodes[childIdx].landmarks.size()*queryVertices.size());
                    this->stats.incrementStatistic("nodes_inserted",1);
                    if (collectLevelStats) {
                        this->stats.incrementStatistic("z_nodes_inserted_level" + std::to_string(colt.getNodeDepth(childIdx)),1);
                    }
                    if (aggScore == 0) {
                        this->stats.incrementStatistic("nodes_lbzero",1);
                    }
#endif
                }
            } else {
                // Handle depending on whether OL for this node has been searched before
                if (!colt.queryData[nodeIdx].intialized) {
                    colt.queryData[nodeIdx].intialized = true;
                    std::vector<EdgeWeight> lmkDists(queryVertices.size());
                    for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                        lmkDists[i] = phl.Query(queryVertices[i],colt.treeNodes[nodeIdx].landmarks[colt.queryData[nodeIdx].lmIdx]);
                    }
                    colt.queryData[nodeIdx].lmkDist = queryVertexLmkDists.size();
                    queryVertexLmkDists.push_back(std::move(lmkDists));
                    // Note: colt.queryData[nodeIdx].lmkDist marks index of the vector in queryVertexLmkDists with the
                    // distances from query vertices to landmark indicated by lmkIdx
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("nodes_leaf_dequeued_first",1);
                    if (collectLevelStats) {
                        this->stats.incrementStatistic("z_nodes_leaf_dequeued_first_level" + std::to_string(colt.getNodeDepth(nodeIdx)),1);
                    }
                    this->stats.incrementStatistic("spdist_calls",queryVertices.size());
                    this->stats.incrementStatistic("ols_started",1);
#endif
                    
                    // Find the object in this landmark's object list and udpate result set with it
                    // and also initializes left and right indexes from this object
#if defined(COLLECT_STATISTICS)
                    colt.treeNodes[nodeIdx].updateQueueWithMinObjAkNN(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].leftIdx,colt.queryData[nodeIdx].rightIdx,queryVertexLmkDists[colt.queryData[nodeIdx].lmkDist],queryVertices,queryDistances,pqueue,alt,agg,this->stats);
#else
                    colt.treeNodes[nodeIdx].updateQueueWithMinObjAkNN(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].leftIdx,colt.queryData[nodeIdx].rightIdx,queryVertexLmkDists[colt.queryData[nodeIdx].lmkDist],queryVertices,queryDistances,pqueue,alt,agg);
#endif
                }
#if defined(COLLECT_STATISTICS)
                else {
                    stats.incrementStatistic("nodes_leaf_dequeued_repeat",1);
                }
#endif
                // Update queue by searching left and right from current left and right indexes 
                // in object list of chosen landmark (terminate when object lower-bound is higher than
                // priority queue lower-bound - i.e. there maybe a bette object in another obj list)
#if defined(COLLECT_STATISTICS)
                if (colt.treeNodes[nodeIdx].updateQueueAkNN(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].leftIdx,colt.queryData[nodeIdx].rightIdx,queryVertexLmkDists[colt.queryData[nodeIdx].lmkDist],queryVertices,queryDistances,pqueue,alt,lbScore,agg,this->stats)) {
#else
                if (colt.treeNodes[nodeIdx].updateQueueAkNN(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].leftIdx,colt.queryData[nodeIdx].rightIdx,queryVertexLmkDists[colt.queryData[nodeIdx].lmkDist],queryVertices,queryDistances,pqueue,alt,lbScore,agg)) {
#endif
                    // We insert objectID because that's the mapped VertexID for nodeIdx
                    pqueue.insert(objectID,lbScore);
                }
#if defined(COLLECT_STATISTICS)
                else {
                    stats.incrementStatistic("ols_finished",1);
                }
#endif
            }
        }
    }

    knnCandidates.populateKNNs(AkNNs,AkNNDistances);
}

/*
 * Range Search
 */

void ILBR::rangeSearchByCOLTPHL(ALT& alt, COLTree& colt, SULTree& sltree, EdgeWeight range, VertexID queryVertexID, std::vector<VertexID>& results, 
                                 std::vector<EdgeWeight>& resultDistances, PrunedHighwayLabeling& phl) 
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("cand_unpruned",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("spdist_calls",0);
    this->stats.initialiseStatistic("lbs_computed",0);
    this->stats.initialiseStatistic("lbs_computed_node",0);
    this->stats.initialiseStatistic("lbs_computed_obj",0);
    this->stats.initialiseStatistic("ols_finished",0);
    this->stats.initialiseStatistic("nodes_total",colt.getNumNonEmptyNodes());
    this->stats.initialiseStatistic("nodes_inserted",0);
    this->stats.initialiseStatistic("nodes_nonleaf_dequeued",0);
    this->stats.initialiseStatistic("nodes_leaf_dequeued",0);
    this->stats.initialiseStatistic("nodes_lbzero",0);
    this->stats.initialiseStatistic("global_landmark_lb_used",0);
    this->stats.initialiseStatistic("global_landmark_lb_attempted",0);
    this->stats.initialiseStatistic("global_landmark_ub_used",0);
    this->stats.initialiseStatistic("global_landmark_ub_attempted",0);
    this->stats.initialiseStatistic("results",0);
#endif

    // Note: stack only contains tree nodes (never objects)
    std::stack<TreeNodeID> stack;
    
    TreeNodeID nodeIdx, childIdx;
    EdgeWeight lbDist, ubDist, altLBDist, altUBDist;
    
    stack.push(0); // Start with root node in COLTree (with index 0)
    if (colt.treeNodes[0].hasObjectList()) {
        // The first time an OL node is dequeued, it assumes that lmIdx has been set earlier
        // so if root has an OL (because it is the only node), we need to set it here
        colt.treeNodes[0].getRangeSearchLeafBounds(alt,queryVertexID,lbDist,ubDist,colt.rangeQueryLmkLBs,colt.rangeQueryLmkUBs);
        colt.queryData[0].lmIdx = colt.treeNodes[0].getRangeSearchListIdx(range,colt.rangeQueryLmkLBs,colt.rangeQueryLmkUBs);
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("lbs_computed",colt.treeNodes[0].landmarks.size());
        this->stats.incrementStatistic("lbs_computed_node",colt.treeNodes[0].landmarks.size());
#endif
    } 

    while (stack.size() > 0) {
        nodeIdx = stack.top();
        stack.pop();
        if (!colt.treeNodes[nodeIdx].hasObjectList()) {
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("nodes_nonleaf_dequeued",1);
#endif
            // No OL, so we expand search towards children using COLTree
            // Note: These will not include children with no objects
            for (std::size_t i = 0; i < colt.treeNodes[nodeIdx].childNodeIdxs.size(); ++i) {
                childIdx = colt.treeNodes[nodeIdx].childNodeIdxs[i];
                RN_DEBUG_ASSERT(colt.treeNodes[childIdx].getNumObjects() != 0 && "A COL-Tree child node holds at least one object");
                if (!colt.treeNodes[childIdx].hasObjectList()) {
                    colt.treeNodes[childIdx].getBestNodeLowerAndUpperbound(alt,queryVertexID,lbDist,ubDist);
                } else {
                    // If this is a leaf node (i.e., has an OL), also retrieve the bounds to its landmarks, from which
                    // its object list is chosen below if the node is searched
                    colt.treeNodes[childIdx].getRangeSearchLeafBounds(alt,queryVertexID,lbDist,ubDist,colt.rangeQueryLmkLBs,colt.rangeQueryLmkUBs);
                }
                colt.computeGlobalLowerAndUpperbound(childIdx,queryVertexID,sltree,altLBDist,altUBDist);
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("global_landmark_lb_attempted",1);
                this->stats.incrementStatistic("global_landmark_ub_attempted",1);
                // getBestNodeLowerAndUpperbound computes a lower-bound from query vertex to each landmark to get 
                // best lower-bound for node (the upper-bound computation is essentially free)
                this->stats.incrementStatistic("lbs_computed",colt.treeNodes[childIdx].landmarks.size());
                this->stats.incrementStatistic("lbs_computed_node",colt.treeNodes[childIdx].landmarks.size());
#endif
                if (altLBDist > lbDist) {
                    lbDist = altLBDist;
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("global_landmark_lb_used",1);
#endif
                } 
                if (altUBDist < ubDist) {
                    ubDist = altUBDist;
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("global_landmark_ub_used",1);
#endif
                }
                if (ubDist <= range) {
                    // All objects contained in this treeNode must be results so rapidly expand children until
                    // leaf nodes are reached then compute distance to and insert all objects into result set
                    std::stack<TreeNodeID> innerStack;
                    innerStack.push(childIdx);
                    while(innerStack.size() > 0) {
                        TreeNodeID innerNodeIdx = innerStack.top();
                        innerStack.pop();
                        if (!colt.treeNodes[innerNodeIdx].hasObjectList()) {
                            for (std::size_t j = 0; j < colt.treeNodes[innerNodeIdx].childNodeIdxs.size(); ++j) {
                                innerStack.push(colt.treeNodes[innerNodeIdx].childNodeIdxs[j]);
                            }
                        } else {
                            // Just use first object list (it must have at least one) since we just iterate voer all
                            for (std::size_t j = 0; j < colt.treeNodes[innerNodeIdx].lmkObjectLists[0].size(); ++j) {
                                VertexID innerObjectID = colt.treeNodes[innerNodeIdx].lmkObjectLists[0][j].first;
                                results.push_back(innerObjectID);
                                resultDistances.push_back(0);
                            }
                        }
                    }
                } else if (lbDist <= range) {
                    // Range intersects [lbDist,ubDist] so there could be some results objects in this tree node subgraph
                    if (colt.treeNodes[childIdx].hasObjectList()) {
                        // Search the object list with the fewest candidates
                        colt.queryData[childIdx].lmIdx = colt.treeNodes[childIdx].getRangeSearchListIdx(range,colt.rangeQueryLmkLBs,colt.rangeQueryLmkUBs);
                    }
                    stack.push(childIdx);
#if defined(COLLECT_STATISTICS)
                    stats.incrementStatistic("nodes_inserted",1);
#endif
                }
                // Note: This implicitly prunes nodes with lbDist > range (cannot contain any reuslts)
#if defined(COLLECT_STATISTICS)
                if (lbDist == 0) {
                    this->stats.incrementStatistic("nodes_lbzero",1);
                }
#endif
            }
        } else {
#if defined(COLLECT_STATISTICS)
            colt.treeNodes[nodeIdx].rangeSearchObjectList(colt.queryData[nodeIdx].lmIdx,queryVertexID,range,alt,phl,results,resultDistances, stats);
            stats.incrementStatistic("nodes_leaf_dequeued",1);
            stats.incrementStatistic("ols_finished",1);
#else
            colt.treeNodes[nodeIdx].rangeSearchObjectList(colt.queryData[nodeIdx].lmIdx,queryVertexID,range,alt,phl,results,resultDistances);
#endif
        }
    }
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("results",results.size());
#endif
}

// Version without spDist usage to set tighter lower-bound
void ILBR::rangeSearchyNVDQuadtreePHL(ALT& alt, NVDQuadtree& nvd, EdgeWeight range, VertexID queryVertexID, 
                                      std::vector<VertexID>& results, std::vector<EdgeWeight>& resultDistances, 
                                    PrunedHighwayLabeling& phl, Graph& graph)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("cand_unpruned",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("spdist_calls",0);
    this->stats.initialiseStatistic("lbs_computed",0);
    this->stats.initialiseStatistic("results",0);
#endif

    // Retrieve k candidates by landmark lower bounds and compute distance to each of them
    VertexID candidate;
    google::dense_hash_set<VertexID> neighboursAdded;
    neighboursAdded.set_empty_key(constants::UNUSED_VERTEX_ID);
    EdgeWeight spDist, lbDist, ubDist;
    [[maybe_unused]] int numLBs;
    
    BinaryMinHeap<EdgeWeight,VertexEdgeWeightPair> nvdQueue; // VertexID element is vnsID (NOT objectID)
    VertexID vnsID;
    spDist = 0; // Need to make sure this is initialized as generateNeighbourCandidatesRangeSearch will try to use it
    candidate = nvd.getNearestNeighbourHashSet(graph.mortonCodes[queryVertexID],neighboursAdded,vnsID);
    alt.getLowerAndUpperBound(queryVertexID,candidate,lbDist,ubDist);
#if defined(COLLECT_STATISTICS)
    stats.incrementStatistic("cand_unpruned",1);
    // A lower-bound must have been computed for this candidate
    this->stats.incrementStatistic("lbs_computed",1);
#endif
    if (ubDist <= range) {
        results.push_back(candidate);
        resultDistances.push_back(0);
    } else if (lbDist <= range) {
        spDist = phl.Query(queryVertexID,candidate);
#if defined(COLLECT_STATISTICS)
        stats.incrementStatistic("cand_spdist",1);
        this->stats.incrementStatistic("spdist_calls",1);
#endif
        if (spDist <= range) {
            results.push_back(candidate);
            resultDistances.push_back(0);
        } else {
            // Nearest neighbour is outside range, so we do not need to search anymore
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("false_hit",1);
#endif
            return;
        }
    } else {
        // Nearest neighbour is outside range, so we do not need to search anymore
        return;
    }
    numLBs = nvd.generateNeighbourCandidatesRangeSearch(alt,nvdQueue,queryVertexID,vnsID,spDist,neighboursAdded);
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("lbs_computed",numLBs);
    this->stats.incrementStatistic("cand_unpruned",numLBs); // numLBs is also the number candidates inserted into queue
#endif

    // Terminate when lower-bound to next nearest neighbour is greater than range or no more candidates
    while (nvdQueue.size() > 0 && nvdQueue.getMinKey() <= range) {
        VertexEdgeWeightPair element = nvdQueue.extractMinElement();
        vnsID = element.first;
        ubDist = element.second;
        spDist = 0; // Need to make sure this is initialized as generateNeighbourCandidatesRangeSearch will try to use it
        candidate = nvd.vnsObjectIDs[vnsID];
        if (ubDist <= range) {
            results.push_back(candidate);
            resultDistances.push_back(0);
        } else {
            // lbDist is actually element's key, which we know is <= range by loop definition
            spDist = phl.Query(queryVertexID,candidate);
#if defined(COLLECT_STATISTICS)
            stats.incrementStatistic("cand_spdist",1);
            this->stats.incrementStatistic("spdist_calls",1);
#endif
            if (spDist <= range) {
                results.push_back(candidate);
                resultDistances.push_back(0);
#if defined(COLLECT_STATISTICS)
            } else {
                this->stats.incrementStatistic("false_hit",1);
#endif
            }
        } 
        numLBs = nvd.generateNeighbourCandidatesRangeSearch(alt,nvdQueue,queryVertexID,vnsID,spDist,neighboursAdded);
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("lbs_computed",numLBs);
        this->stats.incrementStatistic("cand_unpruned",numLBs); // numLBs is also the number candidates inserted into queue
#endif
    }
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("results",results.size());
#endif
}

