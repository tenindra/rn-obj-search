// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "IUBR.h"

#include "../queue/BinaryMaxHeap.h"
#include "../queue/BinaryMinHeap.h"

#include <algorithm>
#include <set>
#include <google/dense_hash_map>
#include <iostream>
#include <cmath>

void IUBR::getFkNsByPHLNaive(ALT& alt, std::vector<VertexID>& objects, unsigned int k, VertexID queryVertexID, 
                             std::vector<VertexID>& FkNs, std::vector<EdgeWeight>& FkNDistances, PrunedHighwayLabeling& phl)
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

    // Quickly set Dk using the first k objects (this will be random but better than nothing)
    EdgeWeight upperBound, spDist, Dk;
    BinaryMinHeap<EdgeWeight,VertexID> fknCandidates;
    VertexID candidate;
    std::size_t i = 0;
    for (; i < k && i < objects.size(); ++i) {
        candidate = objects[i];
        spDist = phl.Query(queryVertexID,candidate);
        fknCandidates.insert(candidate,spDist);
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("cand_spdist",1);
        this->stats.incrementStatistic("cand_unpruned",1);
        this->stats.incrementStatistic("spdist_calls",1);
#endif
    }
    Dk = fknCandidates.getMinKey();

    // Reusing i so we don't visit same objects
    for (; i < objects.size(); ++i) {
        candidate = objects[i];
        upperBound = alt.getUpperBound(queryVertexID,candidate);
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("cand_unpruned",1);
        this->stats.incrementStatistic("lbs_computed",1);
#endif
        if (upperBound > Dk) {
#if defined(COLLECT_STATISTICS)
            // This means some candidate (existing or new) was a poor choice
            this->stats.incrementStatistic("cand_spdist",1);
            this->stats.incrementStatistic("false_hit",1);
            this->stats.incrementStatistic("spdist_calls",1);
#endif
            spDist = phl.Query(queryVertexID,candidate);
            if (spDist > Dk) {
#if defined(COLLECT_STATISTICS)
                // This means some old candidate was a poor choice
                this->stats.incrementStatistic("false_old_candidate",1);
#endif
                // Only insert if it is a better kNN candidate
                fknCandidates.insert(candidate,spDist);
                fknCandidates.extractMinElement();
                Dk = fknCandidates.getMinKey();
            }
#if defined(COLLECT_STATISTICS)
            // This means the new candidate was the bad choice
              else {
                this->stats.incrementStatistic("false_new_candidate",1);
            }
#endif
        }
    }
    fknCandidates.populateFkNs(FkNs,FkNDistances);
}

void IUBR::getFkNsByCOLTPHL(ALT& alt, COLTree& colt, SULTree& sltree, unsigned int k, VertexID queryVertexID, std::vector<VertexID>& FkNs, 
                            std::vector<EdgeWeight>& FkNDistances, PrunedHighwayLabeling& phl) 
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
    this->stats.initialiseStatistic("global_landmark_usage",0);
    bool collectLevelStats = false;
    if (collectLevelStats) {
        int treeDepth = colt.getTreeDepth();
        this->stats.initialiseStatistic("y_total_levels",treeDepth);
        for (int i = 0; i < treeDepth; ++i) {
            this->stats.initialiseStatistic("z_nodes_inserted_level" + std::to_string(i),0);
            this->stats.initialiseStatistic("z_nodes_nonleaf_dequeued_level" + std::to_string(i),0);
            this->stats.initialiseStatistic("z_nodes_leaf_dequeued_first_level" + std::to_string(i),0);
            this->stats.initialiseStatistic("z_global_landmark_usage_level" + std::to_string(i),0);
        }
    }
#endif

    BinaryMinHeap<EdgeWeight,VertexID> fknCandidates;
    BinaryMaxHeap<EdgeWeight,VertexID> pqueue;
    
    VertexID objectID;
    TreeNodeID nodeIdx, childIdx;
    EdgeWeight Dk = 0, spDist, ubDist;
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
        colt.queryData[0].lmIdx = colt.treeNodes[0].getNearestLmkIdxByUB(alt,queryVertexID);
#if defined(COLLECT_STATISTICS)
        // getNearestLmkIdxByUB computes a upper-bound from query vertex to each landmark to get best landmark
        this->stats.incrementStatistic("lbs_computed",colt.treeNodes[0].landmarks.size());
        this->stats.incrementStatistic("lbs_computed_node",colt.treeNodes[0].landmarks.size());
#endif
    }

    while (pqueue.size() > 0 && (DkInfinity || pqueue.getMaxKey() > Dk)) {
        objectID = pqueue.extractMaxElement();
        if (objectID < firstMappedTreeNodeID) {
            // Element is an object
            spDist = phl.Query(queryVertexID,objectID);
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("cand_spdist",1);
            this->stats.incrementStatistic("obj_dequeued",1);
            stats.incrementStatistic("spdist_calls",1);
#endif
            if (DkInfinity) {
                fknCandidates.insert(objectID,spDist);
                if (fknCandidates.size() == k) {
                    DkInfinity = false;
                    Dk = fknCandidates.getMinKey();
                }            
            } else if (spDist > Dk) {
                fknCandidates.insert(objectID,spDist);
                fknCandidates.extractMinElement();
                Dk = fknCandidates.getMinKey();
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
                    colt.queryData[childIdx].reset();
                    if (!colt.treeNodes[childIdx].hasObjectList()) {
                        ubDist = colt.treeNodes[childIdx].getBestNodeUpperbound(alt,queryVertexID);
                    } else {
                        // If this is a node with OL, save furthest landmark simultaneously.
                        // (This will reduce the number of duplicate ALT calls when it is searched)
                        int lmIdx;
                        ubDist = colt.treeNodes[childIdx].getBestNodeUpperboundAndNearestLmk(alt,queryVertexID,lmIdx);
                        colt.queryData[childIdx].lmIdx = lmIdx;
                    }
                    EdgeWeight altUBDist = colt.computeGlobalUpperbound(childIdx,queryVertexID,sltree);
                    if (altUBDist < ubDist) {
                        ubDist = altUBDist;
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("global_landmark_usage",1);
                        if (collectLevelStats) {
                            this->stats.incrementStatistic("z_global_landmark_usage_level" + std::to_string(colt.getNodeDepth(childIdx)),1);
                        }
#endif
                    } 
                    pqueue.insert(childIdx+firstMappedTreeNodeID,ubDist);
#if defined(COLLECT_STATISTICS)
                    // getBestNodeUpperboundAndNearestLmk computes a lower-bound from query vertex to each landmark to get best lower-bound for node
                    this->stats.incrementStatistic("lbs_computed",colt.treeNodes[childIdx].landmarks.size());
                    this->stats.incrementStatistic("lbs_computed_node",colt.treeNodes[childIdx].landmarks.size());
                    this->stats.incrementStatistic("nodes_inserted",1);
                    if (collectLevelStats) {
                        this->stats.incrementStatistic("z_nodes_inserted_level" + std::to_string(colt.getNodeDepth(childIdx)),1);
                    }
                    if (ubDist == 0) {
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
                    colt.treeNodes[nodeIdx].updateQueueWithMaxObjFkN(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].rightIdx,colt.queryData[nodeIdx].lmkDist,queryVertexID,pqueue,alt,this->stats);
#else
                    colt.treeNodes[nodeIdx].updateQueueWithMaxObjFkN(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].rightIdx,colt.queryData[nodeIdx].lmkDist,queryVertexID,pqueue,alt);
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
                if (colt.treeNodes[nodeIdx].updateQueueFkN(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].rightIdx,colt.queryData[nodeIdx].lmkDist,queryVertexID,pqueue,alt,ubDist,this->stats)) {
#else
                if (colt.treeNodes[nodeIdx].updateQueueFkN(colt.queryData[nodeIdx].lmIdx,colt.queryData[nodeIdx].rightIdx,colt.queryData[nodeIdx].lmkDist,queryVertexID,pqueue,alt,ubDist)) {
#endif
                    // We insert objectID because that's the mapped VertexID for nodeIdx
                    pqueue.insert(objectID,ubDist);
                }
#if defined(COLLECT_STATISTICS)
                else {
                    stats.incrementStatistic("ols_finished",1);
                }
#endif
            }
        }
    }

    fknCandidates.populateFkNs(FkNs,FkNDistances);
}

