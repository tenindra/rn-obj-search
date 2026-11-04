// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "DijkstraSearch.h"

#include "../queue/RadixHeap.h"
#include "../utility/geometry.h"

#include <cmath>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <random>

// We assume priority has been cleared
void DijkstraSearch::findSSSPDistances(Graph& graph, VertexID source, std::vector<EdgeWeight>& targetDistances, 
                                       BinaryMinHeap<EdgeWeight,VertexID> *pqueue)
{
    // We assume targetDistances has been resized for the size of the graph
    // However it does not need to have zero values as we overwrite them
    
    // Note: Using vector of bools is significantly faster than unordered_set.
    // Unlike the case where Dijkstra's search is between one source and
    // and one target, this is will never be less space efficient than unordered_set
    // because they will have the same number of elements. In fact std::vector<bool>
    // is more space efficient than a standard vector which is more space
    // and time efficient than a unordered_set (array vs. hash table)
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);
    int adjListStart, nextAdjListStart;
    
    EdgeWeight minDist;
    VertexID minDistVertexID, adjVertex;
    
    // Initialize with priority queue with source vertex
    pqueue->insert(source,0);
    
    while (pqueue->size() > 0) {
        // Extract and remove vertex with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue->getMinKey();
        minDistVertexID = pqueue->extractMinElement();
        if (!isVertexSettled[minDistVertexID]) {
            isVertexSettled[minDistVertexID] = true;
            targetDistances[minDistVertexID] = minDist;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isVertexSettled[adjVertex]) {
                    pqueue->insert(adjVertex,minDist+graph.edges[i].second);
                }
            }
        }
    }
}

void DijkstraSearch::findSSSPDistances(Graph& graph, VertexID source, std::vector<EdgeWeight>& targetDistances)
{
    // We assume targetDistances has been resized for the size of the graph
    // However it does not need to have zero values as we overwrite them
    
    // Note: Using vector of bools is significantly faster than unordered_set.
    // Unlike the case where Dijkstra's search is between one source and
    // and one target, this is will never be less space efficient than unordered_set
    // because they will have the same number of elements. In fact std::vector<bool>
    // is more space efficient than a standard vector which is more space
    // and time efficient than a unordered_set (array vs. hash table)
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);
    int adjListStart, nextAdjListStart;
    
    RadixHeap<VertexID> pqueue;
    
    EdgeWeight minDist;
    VertexID minDistVertexID, adjVertex;
    
    // Initialize with priority queue with source vertex
    pqueue.push(0,source);
    
    while (!pqueue.empty()) {
        // Extract and remove vertex with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDistVertexID = pqueue.top_value();
        minDist = pqueue.top_key();
        pqueue.pop();
        if (!isVertexSettled[minDistVertexID]) {
            isVertexSettled[minDistVertexID] = true;
            targetDistances[minDistVertexID] = minDist;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isVertexSettled[adjVertex]) {
                    pqueue.push(minDist+graph.edges[i].second,adjVertex);
                }
            }
        }
    }
}

// We assume priority has been cleared
void DijkstraSearch::findSSMTDistances(DynamicGraph& graph, VertexID source, 
                                       std::unordered_set<VertexID>& targetSet, 
                                       std::unordered_map<VertexID,EdgeWeight>& results,
                                       BinaryMinHeap<EdgeWeight,VertexID> *pqueue)
{
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);

    EdgeWeight minDist;
    VertexID minDistVertexID, neighbourVertexID;
    
    // Initialize with priority queue with source vertex
    pqueue->insert(source,0);
    std::size_t targetsFound = 0;
    
    while (pqueue->size() > 0) {
        // Extract and remove vertex with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue->getMinKey();
        minDistVertexID = pqueue->extractMinElement();
        if (!isVertexSettled[minDistVertexID]) {
            isVertexSettled[minDistVertexID] = 1;

            if (targetSet.find(minDistVertexID) != targetSet.end()) {
                // We have found the SP distance to one of the target set
                // So we add it to result set and remove it from input set
                results[minDistVertexID] = minDist; // This is the final SP distance from source
                ++targetsFound;
                if (targetsFound == targetSet.size()) {
                    break;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            for (std::size_t i = 0; i < graph.vertices[minDistVertexID].adjVertices.size(); ++i) {
                neighbourVertexID = graph.vertices[minDistVertexID].adjVertices[i];
                // Only update those we haven't already settled
                if (!isVertexSettled[neighbourVertexID]) {
                    pqueue->insert(neighbourVertexID,minDist + graph.vertices[minDistVertexID].adjVertexWgts[i]);
                }
            }
        }
    }

    if (targetsFound != targetSet.size()) {
        std::cerr << "INTERNAL ERROR: the Dijkstra search from vertex " << source << " reached only " << targetsFound << " of its "
                  << targetSet.size() << " targets" << std::endl;
        std::exit(1);
    }
}

// We assume priority has been cleared
void DijkstraSearch::findSSMTDistances(Graph& graph, VertexID source, 
                                       std::unordered_set<VertexID>& targetSet, 
                                       std::unordered_map<VertexID,EdgeWeight>& results,
                                       BinaryMinHeap<EdgeWeight,VertexID> *pqueue)
{
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);
    
    EdgeWeight minDist;
    VertexID minDistVertexID, adjVertex;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with source vertex
    pqueue->insert(source,0);
    std::size_t targetsFound = 0;
    
    while (pqueue->size() > 0) {
        // Extract and remove vertex with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue->getMinKey();
        minDistVertexID = pqueue->extractMinElement();
        if (!isVertexSettled[minDistVertexID]) {
            isVertexSettled[minDistVertexID] = true;

            if (targetSet.find(minDistVertexID) != targetSet.end()) {
                // We have found the SP distance to one of the target set
                // So we add it to result set and remove it from input set
                results[minDistVertexID] = minDist; // This is the final SP distance from source
                ++targetsFound;
                if (targetsFound == targetSet.size()) {
                    break;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isVertexSettled[adjVertex]) {
                    // This is the first time we have seen this vertex (i.e. distance infinity)
                    pqueue->insert(adjVertex,minDist+graph.edges[i].second);
                }
            }
        }
    }

    if (targetsFound != targetSet.size()) {
        std::cerr << "INTERNAL ERROR: the Dijkstra search from vertex " << source << " reached only " << targetsFound << " of its "
                  << targetSet.size() << " targets" << std::endl;
        std::exit(1);
    }
}

int DijkstraSearch::findSSMTDistances(Graph& graph, VertexID source, std::unordered_set<VertexID>& targetSet, 
                                      std::unordered_map<VertexID,EdgeWeight>& results)
{
    int numSettledVertices = 0;
    
    RadixHeap<VertexID> pqueue;
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);
    
    EdgeWeight minDist;
    VertexID minDistVertexID, adjVertex;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with source vertex
    pqueue.push(0,source);
    std::size_t targetsFound = 0;
    
    while (!pqueue.empty()) {
        // Extract and remove vertex with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDistVertexID = pqueue.top_value();
        minDist = pqueue.top_key();
        pqueue.pop();
        if (!isVertexSettled[minDistVertexID]) {
            isVertexSettled[minDistVertexID] = true;
            ++numSettledVertices;

            if (targetSet.find(minDistVertexID) != targetSet.end()) {
                // We have found the SP distance to one of the target set
                // So we add it to result set and remove it from input set
                results[minDistVertexID] = minDist; // This is the final SP distance from source
                ++targetsFound;
                if (targetsFound == targetSet.size()) {
                    break;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isVertexSettled[adjVertex]) {
                    // This is the first time we have seen this vertex (i.e. distance infinity)
                    pqueue.push(minDist+graph.edges[i].second,adjVertex);
                }
            }
        }
    }

    if (targetsFound != targetSet.size()) {
        std::cerr << "INTERNAL ERROR: the Dijkstra search from vertex " << source << " reached only " << targetsFound << " of its "
                  << targetSet.size() << " targets" << std::endl;
        std::exit(1);
    }
    return numSettledVertices;
}

EdgeWeight DijkstraSearch::findShortestPathDistanceSubgraph(Graph& graph, VertexID source, VertexID target, std::vector<bool>& edgeInSubgraph)
{
    BinaryMinHeap<EdgeWeight,VertexID> pqueue;
    // We use unordered_set because we are only searching subgraph (no benefit 
    // from allocating std::vector<bool> for all vertices, just overhead)
    std::unordered_set<VertexID> settledVertexIDs;
    
    EdgeWeight minDist, distanceToTarget = 0;
    VertexID minDistVertexID, adjVertex;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with source vertex
    pqueue.insert(source,0);
    
    while (pqueue.size() > 0) {
        // Extract and remove vertex with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistVertexID = pqueue.extractMinElement();
        if (settledVertexIDs.find(minDistVertexID) == settledVertexIDs.end()) {
            settledVertexIDs.insert(minDistVertexID); // Mark it as "seen" so we can avoid later

            if (minDistVertexID == target) {
                // If the minimum is the target we have finished searching
                distanceToTarget = minDist;
                break;
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                // Only update those we haven't already settled and within subgraph
                if (settledVertexIDs.find(adjVertex) == settledVertexIDs.end() && edgeInSubgraph[i]) {
                    // This is the first time we have seen this vertex (i.e. distance infinity)
                    pqueue.insert(adjVertex,minDist+graph.edges[i].second);
                }
            }
        }
    }

    return distanceToTarget;
}

VertexID DijkstraSearch::findFarthestVertex(Graph& graph, VertexID source, EdgeWeight& farthestDist)
{
    BinaryMinHeap<EdgeWeight,VertexID> pqueue;
    
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);
    
    EdgeWeight minDist = 0;
    VertexID minDistVertexID = source, adjVertex;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with source vertex
    pqueue.insert(source,0);
    
    while (pqueue.size() > 0) {
        // Extract and remove vertex with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistVertexID = pqueue.extractMinElement();
        if (!isVertexSettled[minDistVertexID]) {
            isVertexSettled[minDistVertexID] = true;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isVertexSettled[adjVertex]) {
                    // This is the first time we have seen this vertex (i.e. distance infinity)
                    pqueue.insert(adjVertex,minDist+graph.edges[i].second);
                }
            }
        }
    }

    // When pqueue is empty, the furthest vertex is the last extracted element
    farthestDist = minDist;
    return minDistVertexID;
}

VertexID DijkstraSearch::findFarthestVertexByBFS(Graph& graph, VertexID source)
{
    std::deque<VertexID> pqueue;
    
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);
    
    VertexID minDistVertexID = source, adjVertex;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with source vertex
    pqueue.push_back(source);
    
    while (pqueue.size() > 0) {
        // Extract and remove vertex with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDistVertexID = pqueue.front();
        pqueue.pop_front();
        if (!isVertexSettled[minDistVertexID]) {
            isVertexSettled[minDistVertexID] = true;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isVertexSettled[adjVertex]) {
                    // This is the first time we have seen this vertex
                    pqueue.push_back(adjVertex);
                }
            }
        }
    }

    // When pqueue is empty, the furthest vertex is the last extracted element
    return minDistVertexID;

}

/*
 * Graph Property Computation
 */

EdgeWeight DijkstraSearch::getLowerboundGraphDiameterByDoubleSweep(Graph &graph, int numSweeps)
{
    // Choose random seed vertices in subgraph, from a fixed seed so that the estimate depends only on the graph
    std::mt19937 e1(DOUBLE_SWEEP_SEED);
    std::uniform_int_distribution<int> uniform_dist2(0,graph.getNumVertices()-1);
    std::unordered_set<VertexID> seeds;

    if (numSweeps > graph.getNumVertices()) {
        std::cerr << "Cannot perform " << numSweeps << " double-sweeps from unique seeds on a graph with " << graph.getNumVertices() << " vertices" << std::endl;
        std::exit(1);
    }

    EdgeWeight diameter = 0;

    // Repeatedly perform double-sweep up to numSeeds
    for (int i = 0; i < numSweeps; ++i) {
        VertexID seed = 0;
        bool uniqueSeedFound = false;
        while (!uniqueSeedFound) {
            seed = uniform_dist2(e1);
            if (seeds.find(seed) == seeds.end()) {
                uniqueSeedFound = true;
                seeds.insert(seed);
            }
        }
        // Perform double-sweep from seed
        EdgeWeight farthestVertexDist;
        VertexID farthestVertex = this->findFarthestVertex(graph,seed,farthestVertexDist);
        this->findFarthestVertex(graph,farthestVertex,farthestVertexDist);
        if (diameter < farthestVertexDist) {
            diameter = farthestVertexDist;
        }
    }
    return diameter;
}