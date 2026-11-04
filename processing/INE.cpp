// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "INE.h"

#include "../queue/BinaryMinHeap.h"
#include "../queue/BinaryMaxHeap.h"

#include <iostream>

void INE::getKNNs(Graph& graph, std::vector<bool>& vertexObjectStatus, unsigned int k, VertexID queryVertexID, std::vector<VertexID>& kNNs, std::vector<EdgeWeight>& kNNDistances)
{
    // Note: We assumpe vectors passed to this method are empty
    BinaryMinHeap<EdgeWeight,VertexID> pqueue;
    EdgeWeight minDist, newDistance;
    VertexID minDistVertexID, adjVertex;
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with query vertex ID
    pqueue.insert(queryVertexID,0);
    
    while (pqueue.size() > 0) {
        // Extract and remove vertex with smallest distance from query point
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistVertexID = pqueue.extractMinElement();
        if (!isVertexSettled[minDistVertexID]) {
            isVertexSettled[minDistVertexID] = 1;

            if (vertexObjectStatus[minDistVertexID]) {
                // If the minimum is an object we have found a kNN
                kNNs.push_back(minDistVertexID);
                kNNDistances.push_back(minDist);
                if (kNNs.size() == k) {
                    // If this is the kth nearest neighbour object
                    // then we no longer need to inspect neighbours
                    break;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);
            // Note: We have to make sure we don't exceed size of graph.edges vector
            
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                if (!isVertexSettled[adjVertex]) {
                    // Only update those we haven't already settled
                    newDistance = minDist + graph.edges[i].second;
                    pqueue.insert(adjVertex,newDistance);
                }
            }
        }
    }
}

void INE::getAkNNs(Graph& graph, std::vector<bool>& vertexObjectStatus, unsigned int k, std::vector<VertexID>& queryVertices, 
                  std::vector<VertexID>& AkNNs, std::vector<EdgeWeight>& AkNNDistances, AggregateFunction<EdgeWeight>& agg)
{
    // Note: We assumpe vectors passed to this method are empty
    EdgeWeight Dk = 0, minDist, newDistance, aggDist, threshold = 0;
    VertexID minDistVertexID, adjVertex, queryVertexIdx;
    int adjListStart, nextAdjListStart;
    
    BinaryMaxHeap<EdgeWeight,VertexID> knnCandidates;
    
    BinaryMinHeap<EdgeWeight,VertexPair> pqueue;
    std::vector<std::unordered_set<VertexID>> settledVertexIDs(queryVertices.size());
    std::unordered_map<VertexID,VertexID> settledObjects;
    std::unordered_map<VertexID,std::vector<EdgeWeight>> settledObjectsDistances;

    // Initialize with priority queue with query vertex ID
    for (std::size_t i = 0; i < queryVertices.size(); ++i) {
        pqueue.insert(VertexPair(queryVertices[i],i),0);
    }
    
    bool DkInfinity = true;
    
    int settledCount = 0;
    while (pqueue.size() > 0) {
        minDist = pqueue.getMinKey();
        if (!DkInfinity) {
            threshold = agg.getThreshold(Dk,queryVertices.size());
            if (settledCount == 0 && minDist >= threshold) {
                knnCandidates.populateKNNs(AkNNs,AkNNDistances);
                return;
            }
        }
        
        // Extract and remove vertex with smallest distance from query point
        // and mark it as "settled" so we do not inspect again
        VertexPair minElement = pqueue.extractMinElement();
        minDistVertexID = minElement.first;
        queryVertexIdx = minElement.second;
        
        if (settledVertexIDs[queryVertexIdx].find(minDistVertexID) == settledVertexIDs[queryVertexIdx].end()) {
            settledVertexIDs[queryVertexIdx].insert(minDistVertexID);
            // Vertex has not been settled from the search from queryVertexIdx yet
            if (vertexObjectStatus[minDistVertexID]) {
                // Vertex is an object!
                if (settledObjects.find(minDistVertexID) == settledObjects.end()) {
                    settledObjects[minDistVertexID] = 1;
                    std::vector<EdgeWeight> distances;
                    settledObjectsDistances[minDistVertexID] = distances;
                } else {
                    ++settledObjects[minDistVertexID];
                }
                settledObjectsDistances[minDistVertexID].push_back(minDist);
                ++settledCount;
                
                if (settledObjects[minDistVertexID] == queryVertices.size()) {
                    settledCount = settledCount-queryVertices.size();
                    aggDist = agg.aggregate(settledObjectsDistances[minDistVertexID]);
                    
                    if (DkInfinity) {
                        knnCandidates.insert(minDistVertexID,aggDist);
                        if (knnCandidates.size() == k) {
                            DkInfinity = false;
                            Dk = knnCandidates.getMaxKey();
                        }            
                    } else if (aggDist < Dk) {
                        knnCandidates.insert(minDistVertexID,aggDist);
                        knnCandidates.extractMaxElement();
                        Dk = knnCandidates.getMaxKey();
                    }
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);
            // Note: We have to make sure we don't exceed size of graph.edges vector
            
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                // Only update those we haven't already settled
                newDistance = minDist + graph.edges[i].second;
                pqueue.insert(VertexPair(adjVertex,queryVertexIdx),newDistance);
            }
        }
    }
    // Reached when all vertices have been settled by every search (e.g., k is larger than the number of objects)
    knnCandidates.populateKNNs(AkNNs,AkNNDistances);
}

void INE::getFkNs(Graph& graph, std::vector<bool>& vertexObjectStatus, int numObjects, int k, VertexID queryVertexID, 
                  std::vector<VertexID>& FkNs, std::vector<EdgeWeight>& FkNDistances)
{
    // Note: We assumpe vectors passed to this method are empty
    BinaryMinHeap<EdgeWeight,VertexID> pqueue;
    EdgeWeight minDist, newDistance;
    VertexID minDistVertexID, adjVertex;
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);
    int adjListStart, nextAdjListStart, objectsRemaining = numObjects;
    
    // Initialize with priority queue with query vertex ID
    pqueue.insert(queryVertexID,0);

    while (pqueue.size() > 0 && objectsRemaining > 0) {
        // Extract and remove vertex with smallest distance from query point
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistVertexID = pqueue.extractMinElement();

        if (!isVertexSettled[minDistVertexID]) {
            isVertexSettled[minDistVertexID] = 1;

            if (vertexObjectStatus[minDistVertexID]) {
                --objectsRemaining;
                // Each of the last k objects to be discovered are results
                if (objectsRemaining < k) {
                    // Note: Less than k so we include the current object (it is kth furthest)
                    FkNs.push_back(minDistVertexID);
                    FkNDistances.push_back(minDist);
                    if (objectsRemaining == 0) {
                        // This is the last object found so must be the farthest
                        break;
                    }
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);
            // Note: We have to make sure we don't exceed size of graph.edges vector
            
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                if (!isVertexSettled[adjVertex]) {
                    // Only update those we haven't already settled
                    newDistance = minDist + graph.edges[i].second;
                    pqueue.insert(adjVertex,newDistance);
                }
            }
        }
    }
}

void INE::rangeSearch(Graph& graph, std::vector<bool>& vertexObjectStatus, int numObjects, EdgeWeight range, VertexID queryVertexID, std::vector<VertexID>& results, std::vector<EdgeWeight>& resultDistances)
{
    // Note: We assume vectors passed to this method are empty

    BinaryMinHeap<EdgeWeight,VertexID> pqueue;
    EdgeWeight minDist, newDistance;
    VertexID minDistVertexID, adjVertex, objectsRemaining = numObjects;
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with query vertex ID
    pqueue.insert(queryVertexID,0);
    
    while (pqueue.size() > 0 && pqueue.getMinKey() <= range) {
        // Extract and remove vertex with smallest distance from query point
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistVertexID = pqueue.extractMinElement();
        if (!isVertexSettled[minDistVertexID]) {
            isVertexSettled[minDistVertexID] = 1;
            if (vertexObjectStatus[minDistVertexID]) {
                --objectsRemaining;
                // minDist must be within range if reach here, so just add to results
                results.push_back(minDistVertexID);
                resultDistances.push_back(0);
                if (objectsRemaining == 0) {
                    // Terminate search early if all objects found
                    break; // Break instead of return so we can collect stats below
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);
            // Note: We have to make sure we don't exceed size of graph.edges vector
            
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                if (!isVertexSettled[adjVertex]) {
                    // Only update those we haven't already settled
                    newDistance = minDist + graph.edges[i].second;
                    pqueue.insert(adjVertex,newDistance);
                }
            }
        }
    }
}

std::vector<bool> INE::buildVertexObjectStatusVector (Graph& graph, std::vector<VertexID>& objectVertices)
{
    std::vector<bool> objectVertexStatus(graph.getNumVertices(),false);
    for(std::size_t i = 0; i < objectVertices.size(); ++i) {
        objectVertexStatus[objectVertices[i]] = true;
    }
    return objectVertexStatus;
}

