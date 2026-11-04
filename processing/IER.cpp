// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "IER.h"

#include "DijkstraSearch.h"
#include "../queue/BinaryMaxHeap.h"

#include <cmath>

void IER::getKNNsByPHL(PackedRtree& rtree, unsigned int k, VertexID queryVertexID, std::vector<VertexID>& kNNs, 
                       std::vector<EdgeWeight>& kNNDistances, Graph& graph, PrunedHighwayLabeling& phl)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
#endif
    
    // Retrieve kNN by euclidean distance
    std::vector<VertexID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    euclideanKNNs.reserve(k);
    euclideanKNNDistances.reserve(k);
    Coordinate queryVertexX, queryVertexY;
    graph.getCoordinates(queryVertexID,queryVertexX,queryVertexY);
    BinaryMinHeap<EuclideanDistanceSqrd,PackedRtreeDataTuple> heap = rtree.getKNNs(k,queryVertexX,queryVertexY,euclideanKNNs,euclideanKNNDistances);
    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to each of these
    BinaryMaxHeap<EdgeWeight,VertexID> knnCandidates;
    EdgeWeight spDist, Dk, euclidDist;
    for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("cand_spdist",1);
        this->stats.incrementStatistic("spdist_calls",1);
#endif
        spDist = phl.Query(queryVertexID,euclideanKNNs[i]);
        knnCandidates.insert(euclideanKNNs[i],spDist);
    }

    // While the euclidean distance to the next euclidean nearest neigbour
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    VertexID nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    while (true) {
        if (rtree.getNextNearestNeighbour(heap,queryVertexX,queryVertexY,nextEuclidNN,currEuclidDistSqrd)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            euclidDist = std::floor(std::sqrt(currEuclidDistSqrd)); // Floor as it's a lowerbound
            if (euclidDist < Dk) {
#if defined(COLLECT_STATISTICS)
                // This means some candidate (existing or new) was a poor choice
                this->stats.incrementStatistic("cand_spdist",1);
                this->stats.incrementStatistic("false_hit",1);
                this->stats.incrementStatistic("spdist_calls",1);
#endif
                spDist = phl.Query(queryVertexID,nextEuclidNN);
                if (spDist < Dk) {
#if defined(COLLECT_STATISTICS)
                    // This means some old candidate was a poor choice
                    this->stats.incrementStatistic("false_old_candidate",1);
#endif
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(nextEuclidNN,spDist);
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

    knnCandidates.populateKNNs(kNNs,kNNDistances);
}


void IER::getKNNsByPHLTravelTimes(PackedRtree& rtree, unsigned int k, VertexID queryVertexID, std::vector<VertexID>& kNNs, 
                                  std::vector<EdgeWeight>& kNNDistances, Graph& graph, PrunedHighwayLabeling& phl)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
#endif

    // Retrieve kNN by euclidean distance
    std::vector<VertexID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    euclideanKNNs.reserve(k);
    euclideanKNNDistances.reserve(k);
    Coordinate queryVertexX, queryVertexY;
    graph.getCoordinates(queryVertexID,queryVertexX,queryVertexY);

    BinaryMinHeap<EuclideanDistanceSqrd,PackedRtreeDataTuple> heap = rtree.getKNNs(k,queryVertexX,queryVertexY,euclideanKNNs,euclideanKNNDistances);

    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to each of these
    BinaryMaxHeap<EdgeWeight,VertexID> knnCandidates;
    EdgeWeight spDist, Dk, minTimeByEuclid;
    for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("cand_spdist",1);
        this->stats.incrementStatistic("spdist_calls",1);
#endif
        spDist = phl.Query(queryVertexID,euclideanKNNs[i]);
        knnCandidates.insert(euclideanKNNs[i],spDist);
    }

    // While the euclidean distance to the next euclidean nearest neigbour
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    VertexID nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    while (true) {
        if (rtree.getNextNearestNeighbour(heap,queryVertexX,queryVertexY,nextEuclidNN,currEuclidDistSqrd)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            minTimeByEuclid = std::floor(std::sqrt(currEuclidDistSqrd)/graph.getMaxGraphSpeedByEdge()); // Floor as it's a lowerbound
            if (minTimeByEuclid < Dk) {
#if defined(COLLECT_STATISTICS)
                // This means some candidate (existing or new) was a poor choice
                this->stats.incrementStatistic("cand_spdist",1);
                this->stats.incrementStatistic("false_hit",1);
                this->stats.incrementStatistic("spdist_calls",1);
#endif
                spDist = phl.Query(queryVertexID,nextEuclidNN);
                if (spDist < Dk) {
#if defined(COLLECT_STATISTICS)
                    // This means some old candidate was a poor choice
                    this->stats.incrementStatistic("false_old_candidate",1);
#endif
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(nextEuclidNN,spDist);
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

    knnCandidates.populateKNNs(kNNs,kNNDistances);
}

void IER::getAkNNsByPHL(PackedRtree& rtree, unsigned int k, std::vector<VertexID>& queryVertices, 
                        std::vector<VertexID>& AkNNs, std::vector<EdgeWeight>& AkNNDistances, 
                        Graph& graph, PrunedHighwayLabeling& phl, AggregateFunction<EdgeWeight>& agg, 
                        AggregateFunction<EuclideanDistance>& euclidAgg) 
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
#endif
    // Retrieve kNN by euclidean distance
    std::vector<VertexID> euclideanAkNNs;
    std::vector<EuclideanDistance> euclideanAkNNDistances;
    euclideanAkNNs.reserve(k);
    euclideanAkNNDistances.reserve(k);
    std::vector<Coordinate> queryVerticesX(queryVertices.size()), queryVerticesY(queryVertices.size());
    std::vector<EdgeWeight> queryDistances(queryVertices.size());
    std::vector<EuclideanDistance> euclidQueryDistances(queryVertices.size());
    for (std::size_t i = 0; i < queryVertices.size(); ++i) {
        graph.getCoordinates(queryVertices[i],queryVerticesX[i],queryVerticesY[i]);
    }

    BinaryMinHeap<EuclideanDistance,PackedRtreeDataTuple> heap = rtree.getAkNNs(k,queryVerticesX,queryVerticesY,euclideanAkNNs,euclideanAkNNDistances,euclidAgg,euclidQueryDistances);

    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to each of these
    BinaryMaxHeap<EdgeWeight,VertexID> knnCandidates;
    EdgeWeight Dk, euclidDist, aggDist;
    for (std::size_t i = 0; i < euclideanAkNNs.size(); ++i) {
        for (std::size_t j = 0; j < queryVertices.size(); ++j) {
            queryDistances[j] = phl.Query(queryVertices[j],euclideanAkNNs[i]);
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("cand_spdist",1);
            this->stats.incrementStatistic("spdist_calls",1);
#endif
        }
        aggDist = agg.aggregate(queryDistances);
        knnCandidates.insert(euclideanAkNNs[i],aggDist);
    }

    // While the euclidean distance to the next euclidean nearest neigbour
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    VertexID nextEuclidNN;
    EuclideanDistance currEuclidDist;
    while (true) {
        if (rtree.getNextAggregateNearestNeighbour(heap,queryVerticesX,queryVerticesY,nextEuclidNN,currEuclidDist,euclidAgg,euclidQueryDistances)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            euclidDist = std::floor(currEuclidDist); // Floor as it's a lowerbound
            if (euclidDist < Dk) {
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("false_hit",1);
#endif
                for (std::size_t j = 0; j < queryVertices.size(); ++j) {
                    queryDistances[j] = phl.Query(queryVertices[j],nextEuclidNN);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("cand_spdist",1);
                    this->stats.incrementStatistic("spdist_calls",1);
#endif
                }
                aggDist = agg.aggregate(queryDistances);
                if (aggDist < Dk) {
#if defined(COLLECT_STATISTICS)
                    // This means some old candidate was a poor choice
                    this->stats.incrementStatistic("false_old_candidate",1);
#endif
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(nextEuclidNN,aggDist);
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

    knnCandidates.populateKNNs(AkNNs,AkNNDistances);
}

void IER::getAkNNsByPHLTravelTimes(PackedRtree& rtree, unsigned int k, std::vector<VertexID>& queryVertices, 
                                   std::vector<VertexID>& AkNNs, std::vector<EdgeWeight>& AkNNDistances, 
                                   Graph& graph, PrunedHighwayLabeling& phl, AggregateFunction<EdgeWeight>& agg, 
                                   AggregateFunction<EuclideanDistance>& euclidAgg) 
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
#endif
    // Retrieve kNN by euclidean distance
    std::vector<VertexID> euclideanAkNNs;
    std::vector<EuclideanDistance> euclideanAkNNDistances;
    euclideanAkNNs.reserve(k);
    euclideanAkNNDistances.reserve(k);
    std::vector<Coordinate> queryVerticesX(queryVertices.size()), queryVerticesY(queryVertices.size());
    std::vector<EdgeWeight> queryDistances(queryVertices.size());
    std::vector<EuclideanDistance> euclidQueryDistances(queryVertices.size());
    for (std::size_t i = 0; i < queryVertices.size(); ++i) {
        graph.getCoordinates(queryVertices[i],queryVerticesX[i],queryVerticesY[i]);
    }

    BinaryMinHeap<EuclideanDistance,PackedRtreeDataTuple> heap = rtree.getAkNNs(k,queryVerticesX,queryVerticesY,euclideanAkNNs,euclideanAkNNDistances,euclidAgg,euclidQueryDistances);

    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to each of these
    BinaryMaxHeap<EdgeWeight,VertexID> knnCandidates;
    EdgeWeight Dk, minTimeByEuclid, aggDist;
    for (std::size_t i = 0; i < euclideanAkNNs.size(); ++i) {
        for (std::size_t j = 0; j < queryVertices.size(); ++j) {
            queryDistances[j] = phl.Query(queryVertices[j],euclideanAkNNs[i]);
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("cand_spdist",1);
            this->stats.incrementStatistic("spdist_calls",1);
#endif
        }
        aggDist = agg.aggregate(queryDistances);
        knnCandidates.insert(euclideanAkNNs[i],aggDist);
    }

    // While the euclidean distance to the next euclidean nearest neigbour
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    VertexID nextEuclidNN;
    EuclideanDistance currEuclidDist;
    while (true) {
        if (rtree.getNextAggregateNearestNeighbour(heap,queryVerticesX,queryVerticesY,nextEuclidNN,currEuclidDist,euclidAgg,euclidQueryDistances)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            minTimeByEuclid = std::floor(currEuclidDist/graph.getMaxGraphSpeedByEdge()); // Floor as it's a lowerbound
            if (minTimeByEuclid < Dk) {
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("false_hit",1);
#endif
                for (std::size_t j = 0; j < queryVertices.size(); ++j) {
                    queryDistances[j] = phl.Query(queryVertices[j],nextEuclidNN);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("cand_spdist",1);
                    this->stats.incrementStatistic("spdist_calls",1);
#endif
                }
                aggDist = agg.aggregate(queryDistances);
                if (aggDist < Dk) {
#if defined(COLLECT_STATISTICS)
                    // This means some old candidate was a poor choice
                    this->stats.incrementStatistic("false_old_candidate",1);
#endif
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(nextEuclidNN,aggDist);
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

    knnCandidates.populateKNNs(AkNNs,AkNNDistances);
}

void IER::rangeSearchByPHL(PackedRtree& rtree, EdgeWeight range, VertexID queryVertexID, 
                           std::vector<VertexID>& results, std::vector<EdgeWeight>& resultDistances, 
                           Graph& graph, PrunedHighwayLabeling& phl)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("spdist_calls",0);
    this->stats.initialiseStatistic("results",0);
#endif
    
    // Retrieve kNN by euclidean distance
    std::vector<VertexID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    Coordinate queryVertexX, queryVertexY;
    graph.getCoordinates(queryVertexID,queryVertexX,queryVertexY);
    BinaryMinHeap<EuclideanDistanceSqrd,PackedRtreeDataTuple> heap = rtree.getKNNs(1,queryVertexX,queryVertexY,euclideanKNNs,euclideanKNNDistances);
    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to the first nearest neighbour by Euclidean distance
    EdgeWeight spDist, euclidDist;
    if (euclideanKNNs.size() == 1) {
        // Note: rtree.getKNNs also return Euclidean distance squared
        euclidDist = std::floor(std::sqrt(euclideanKNNDistances[0])); // Floor as it's a lowerbound
        if (euclidDist <= range) {
            spDist = phl.Query(queryVertexID,euclideanKNNs[0]);
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("cand_spdist",1);
            this->stats.incrementStatistic("spdist_calls",1);
#endif
            if (spDist <= range) {
                results.push_back(euclideanKNNs[0]);
                resultDistances.push_back(0);
            }
#if defined(COLLECT_STATISTICS)
            // This means the candidate was a bad choice
                else {
                this->stats.incrementStatistic("false_hit",1);
            }
#endif
        } else {
            // Either no results were returned or nearest neighbour is outside range
            // in which case we do not need to search anymore
            return;
        }
    }

    VertexID nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    while (true) {
        if (rtree.getNextNearestNeighbour(heap,queryVertexX,queryVertexY,nextEuclidNN,currEuclidDistSqrd)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            euclidDist = std::floor(std::sqrt(currEuclidDistSqrd)); // Floor as it's a lowerbound
            if (euclidDist <= range) {
                spDist = phl.Query(queryVertexID,nextEuclidNN);
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("cand_spdist",1);
                this->stats.incrementStatistic("spdist_calls",1);
#endif
                if (spDist <= range) {
                    results.push_back(nextEuclidNN);
                    resultDistances.push_back(0);
                }
#if defined(COLLECT_STATISTICS)
                // This means the candidate was a bad choice
                  else {
                    this->stats.incrementStatistic("false_hit",1);
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
}

void IER::rangeSearchByPHLTravelTimes(PackedRtree& rtree, EdgeWeight range, VertexID queryVertexID, 
                                      std::vector<VertexID>& results, std::vector<EdgeWeight>& resultDistances, 
                                      Graph& graph, PrunedHighwayLabeling& phl)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("cand_spdist",0);
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("spdist_calls",0);
    this->stats.initialiseStatistic("results",0);
#endif
    
    // Retrieve kNN by euclidean distance
    std::vector<VertexID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    Coordinate queryVertexX, queryVertexY;
    graph.getCoordinates(queryVertexID,queryVertexX,queryVertexY);
    BinaryMinHeap<EuclideanDistanceSqrd,PackedRtreeDataTuple> heap = rtree.getKNNs(1,queryVertexX,queryVertexY,euclideanKNNs,euclideanKNNDistances);
    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to the first nearest neighbour by Euclidean distance
    EdgeWeight spDist, minTimeByEuclid;
    if (euclideanKNNs.size() == 1) {
        // Note: rtree.getKNNs also return Euclidean distance squared
        minTimeByEuclid = std::floor(std::sqrt(euclideanKNNDistances[0])/graph.getMaxGraphSpeedByEdge()); // Floor as it's a lowerbound
        if (minTimeByEuclid <= range) {
            spDist = phl.Query(queryVertexID,euclideanKNNs[0]);
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("cand_spdist",1);
            this->stats.incrementStatistic("spdist_calls",1);
#endif
            if (spDist <= range) {
                results.push_back(euclideanKNNs[0]);
                resultDistances.push_back(0);
            }
#if defined(COLLECT_STATISTICS)
            // This means the candidate was a bad choice
                else {
                this->stats.incrementStatistic("false_hit",1);
            }
#endif
        } else {
            // Either no results were returned or nearest neighbour is outside range
            // in which case we do not need to search anymore
            return;
        }
    }

    VertexID nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    while (true) {
        if (rtree.getNextNearestNeighbour(heap,queryVertexX,queryVertexY,nextEuclidNN,currEuclidDistSqrd)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            minTimeByEuclid = std::floor(std::sqrt(currEuclidDistSqrd)/graph.getMaxGraphSpeedByEdge()); // Floor as it's a lowerbound
            if (minTimeByEuclid <= range) {
                spDist = phl.Query(queryVertexID,nextEuclidNN);
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("cand_spdist",1);
                this->stats.incrementStatistic("spdist_calls",1);
#endif
                if (spDist <= range) {
                    results.push_back(nextEuclidNN);
                    resultDistances.push_back(0);
                }
#if defined(COLLECT_STATISTICS)
                // This means the candidate was a bad choice
                  else {
                    this->stats.incrementStatistic("false_hit",1);
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
    this->stats.incrementStatistic("results",results.size());
#endif
}

