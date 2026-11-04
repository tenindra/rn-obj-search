// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "NVDQuadtree.h"

#include "../queue/BinaryMinHeap.h"
#include "../queue/RadixHeap.h"
#include "../utility/utility.h"

#include <unordered_map>
#include <google/dense_hash_map>
#include <iostream>

NVDQuadtree::NVDQuadtree(std::string setType, double setDensity, int setVariable, int setSize):
                             objSetType(setType), objSetDensity(setDensity), objSetVariable(setVariable), objSetSize(setSize)
{
}

// Assumes that buildNVD has not be called before
void NVDQuadtree::buildNVD(Graph& graph, std::vector<VertexID>& objectIDs, std::vector<MortonCode>& mortonCodes, int maxLeafColours)
{
    this->vnsObjectIDs.reserve(objectIDs.size());
    
    std::vector<google::dense_hash_map<VertexID,EdgeWeight>> tempVnsNeighbourIDDistUB(objectIDs.size());
    for (std::size_t i = 0; i < objectIDs.size(); ++i) {
        tempVnsNeighbourIDDistUB[i].set_empty_key(constants::UNUSED_VERTEX_ID); // Required by google::dense_hash_map
    }
    
    RadixHeap<VertexPair> pqueue;
    std::vector<IntDistancePair> vertexVnsID(graph.getNumVertices(),IntDistancePair(-1,0));
    // A temporary copy of the starting colour, as binding the vector's argument to the constant itself alters the generated code
    // of this timed build
    std::vector<VertexID> colourMap(graph.getNumVertices(),static_cast<VertexID>(quadtree_constants::UNASSIGNED_COLOUR));;
    
    EdgeWeight minDist;
    VertexID minDistVertexID, minVnsID, adjVertex;
    int adjListStart, nextAdjListStart;
    
    // We do not need to insert the object into queue, we can directly assigned to their own cells
    unsigned int numNeighbourRelationships = 0; // Keep track of this to reserve vector size later
    for (VertexID i = 0; i < objectIDs.size(); ++i) {
        vertexVnsID[objectIDs[i]].first = i; // vertexVnsID[objectIDs[i]].second already initialized (to 0)
        this->vnsObjectIDs.push_back(objectIDs[i]);
        colourMap[objectIDs[i]] = i;
        
        // Inspect each neighbour and update pqueue using edge weights
        adjListStart = graph.getEdgeListStartIndex(objectIDs[i]);
        nextAdjListStart = graph.getEdgeListSize(objectIDs[i]);

        for (int j = adjListStart; j < nextAdjListStart; ++j) {
            adjVertex = graph.edges[j].first;
            if (vertexVnsID[adjVertex].first == -1) {
                // We only need to insert adjacent vertices who are unassigned or assigned
                // and a neighbour relationship hasn't been established
                pqueue.push(graph.edges[j].second,VertexPair(adjVertex,i));
            } else if (static_cast<VertexID>(vertexVnsID[adjVertex].first) != i) {
                // The encountered vertex is already in another dominance region so we 
                // have found a neighbouring to region that inserted this vertex
                if (tempVnsNeighbourIDDistUB[i].find(vertexVnsID[adjVertex].first) == 
                        tempVnsNeighbourIDDistUB[i].end()) {
                    tempVnsNeighbourIDDistUB[i][vertexVnsID[adjVertex].first] = graph.edges[j].second + vertexVnsID[adjVertex].second;
                    ++numNeighbourRelationships;
                
                    // Note: The above distance is only an upperbound on the distance
                    // not the exact distance to the object in the adjacent dominance
                    // region. This is because the shortest path to that object
                    // may lead through another adjacent dominance region
                } else if (tempVnsNeighbourIDDistUB[i][vertexVnsID[adjVertex].first] > graph.edges[j].second + vertexVnsID[adjVertex].second) {
                    tempVnsNeighbourIDDistUB[i][vertexVnsID[adjVertex].first] = graph.edges[j].second + vertexVnsID[adjVertex].second;
                }
                
                // Note: Reciprocal neighbourhood relationship will be added when
                // the search for that region meets this region
                
                // Note: The search is implicity stopped once a vertex from another 
                // dominance region is encountered as any further vertex encountered
                // will not be added to the current dominance region
            }
        }
    }
    
    while (!pqueue.empty()) {
        // Extract and remove vertex with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.top_key();
        VertexPair minElement = pqueue.top_value();
        minDistVertexID = minElement.first;
        minVnsID = minElement.second;
        pqueue.pop();
        if (vertexVnsID[minDistVertexID].first == -1) {
            // The encountered vertex has not be assigned to a dominance region
            // and since it has the minimum distance it must belong to the
            // region that insert it
            vertexVnsID[minDistVertexID].first = minVnsID;
            vertexVnsID[minDistVertexID].second = minDist;
            colourMap[minDistVertexID] = minVnsID;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                if (vertexVnsID[adjVertex].first == -1) {
                    // We only need to insert adjacent vertices who are unassigned or assigned
                    // and a neighbour relationship hasn't been established
                    pqueue.push(minDist+graph.edges[i].second,VertexPair(adjVertex,minVnsID));
                } else if (static_cast<VertexID>(vertexVnsID[adjVertex].first) != minVnsID) {
                    // The encountered vertex is already in another dominance region so we 
                    // have found a neighbouring to region that inserted this vertex
                    if (tempVnsNeighbourIDDistUB[minVnsID].find(vertexVnsID[adjVertex].first) == 
                            tempVnsNeighbourIDDistUB[minVnsID].end()) {
                        tempVnsNeighbourIDDistUB[minVnsID][vertexVnsID[adjVertex].first] = minDist + graph.edges[i].second + vertexVnsID[adjVertex].second;
                        ++numNeighbourRelationships;
                    
                        // Note: The above distance is only an upperbound on the distance
                        // not the exact distance to the object in the adjacent dominance
                        // region. This is because the shortest path to that object
                        // may lead through another adjacent dominance region
                    } else if (tempVnsNeighbourIDDistUB[minVnsID][vertexVnsID[adjVertex].first] > + graph.edges[i].second + vertexVnsID[adjVertex].second) {
                        tempVnsNeighbourIDDistUB[minVnsID][vertexVnsID[adjVertex].first] = minDist + graph.edges[i].second + vertexVnsID[adjVertex].second;
                    }
                    
                    // Note: Reciprocal neighbourhood relationship will be added when
                    // the search for that region meets this region
                    
                    // Note: The search is implicity stopped once a vertex from another 
                    // dominance region is encountered as any further vertex encountered
                    // will not be added to the current dominance region
                }
            }

        } else if (static_cast<VertexID>(vertexVnsID[minDistVertexID].first) != minVnsID) {
            // Note: Even though we do not insert adjacent vertices already in another VNS into priority queue
            // it is possible we have dequeued that was no assigned to another VNS at the time, but has 
            // since been assigned to one, in which case we still need to add the adjacency relationship
            
            // The encountered vertex is already in another dominance region so we 
            // have found a neighbouring to region that inserted this vertex
            if (tempVnsNeighbourIDDistUB[minVnsID].find(vertexVnsID[minDistVertexID].first) == 
                    tempVnsNeighbourIDDistUB[minVnsID].end()) {
                tempVnsNeighbourIDDistUB[minVnsID][vertexVnsID[minDistVertexID].first] = minDist + vertexVnsID[minDistVertexID].second;
                ++numNeighbourRelationships;
            
                // Note: The above distance is only an upperbound on the distance
                // not the exact distance to the object in the adjacent dominance
                // region. This is because the shortest path to that object
                // may lead through another adjacent dominance region
            } else if (tempVnsNeighbourIDDistUB[minVnsID][vertexVnsID[minDistVertexID].first] > minDist + vertexVnsID[minDistVertexID].second) {
                tempVnsNeighbourIDDistUB[minVnsID][vertexVnsID[minDistVertexID].first] = minDist + vertexVnsID[minDistVertexID].second;
            }
            
            // Note: Reciprocal neighbourhood relationship will be added when
            // the search for that region meets this region
            
            // Note: The search is implicity stopped once a vertex from another 
            // dominance region is encountered as any further vertex encountered
            // will not be added to the current dominance region
        }
    }
    
    EdgeWeight neighbourDist;
    this->firstNeighbourIdx.reserve(tempVnsNeighbourIDDistUB.size()+1);
    this->neighbourIDsAndDist.reserve(numNeighbourRelationships);
    for (VertexID i = 0; i < tempVnsNeighbourIDDistUB.size(); ++i) {
        // Mark the beginning of the list
        this->firstNeighbourIdx.push_back(this->neighbourIDsAndDist.size());
        for (auto keyValue: tempVnsNeighbourIDDistUB[i]) {
            neighbourDist = keyValue.second; // This is only upper bound as explained above
            this->neighbourIDsAndDist.push_back(VertexDistancePair(keyValue.first,neighbourDist));
        }
    }
    this->firstNeighbourIdx.push_back(this->neighbourIDsAndDist.size());

    // Convert vertexVnsID to compressed MortonList
    CompressibleQuadtree qt(maxLeafColours);
    qt.buildFromGraph(graph,mortonCodes,colourMap);
    this->mortonList.buildCompressibleMortonList(qt);
}


std::string NVDQuadtree::getObjSetType()
{
    return this->objSetType;
}

double NVDQuadtree::getObjSetDensity()
{
    return this->objSetDensity;
}

int NVDQuadtree::getObjSetSize()
{
    return this->objSetSize;
}

int NVDQuadtree::getObjSetVariable()
{
    return this->objSetVariable;
}


double NVDQuadtree::computeIndexSize()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(VertexID)*this->vnsObjectIDs.size();
    memoryUsage += sizeof(VertexDistancePair)*this->neighbourIDsAndDist.size();
    memoryUsage += sizeof(VertexID)*this->firstNeighbourIdx.size();
    memoryUsage += this->mortonList.computeIndexSizeBytes();
    return memoryUsage/(1024*1024);
}

VertexID NVDQuadtree::getNearestNeighbourHashSet(MortonCode& queryVertexMC, google::dense_hash_set<VertexID>& neighbourVnsIDsAdded, VertexID& vnsID)
{
    vnsID = this->mortonList.getColourID(queryVertexMC);
    neighbourVnsIDsAdded.insert(vnsID);
    return vnsObjectIDs[vnsID];
}

int NVDQuadtree::generateNeighbourCandidates(ALT& alt, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, VertexID queryVertexID, VertexID vnsID, 
                                              EdgeWeight spdist, google::dense_hash_set<VertexID>& neighbourVnsIDsAdded)
{
    VertexID neighbourVnsID, neighbourObjectID;
    EdgeWeight neighbourObjectLB, neighbourDist, nvdLB;
    int numCandidates = 0;
    
    VertexID neighbourListStart = firstNeighbourIdx[vnsID];
    VertexID nextNeighbourListStart = firstNeighbourIdx[vnsID+1];

    for (VertexID i = neighbourListStart; i < nextNeighbourListStart; ++i) {
        neighbourVnsID = neighbourIDsAndDist[i].first;
        if (neighbourVnsIDsAdded.find(neighbourVnsID) == neighbourVnsIDsAdded.end()) {
            neighbourObjectID = vnsObjectIDs[neighbourVnsID];
            neighbourObjectLB = alt.getLowerBound(queryVertexID,neighbourObjectID);
#if defined(COLLECT_STATISTICS)
            ++numCandidates; // Only need this for stats experiments
#endif
            neighbourDist = neighbourIDsAndDist[i].second;
            nvdLB = (spdist > neighbourDist) ? spdist-neighbourDist : 0; // The other side is zero because neighbourDist is only an upperbound
            if (nvdLB > neighbourObjectLB) {
                neighbourObjectLB = nvdLB;
            }
            pqueue.insert(neighbourVnsID,neighbourObjectLB);
            neighbourVnsIDsAdded.insert(neighbourVnsID);
        }
    }
    return numCandidates;
}

int NVDQuadtree::generateNeighbourCandidatesRangeSearch(ALT& alt, BinaryMinHeap<EdgeWeight,VertexEdgeWeightPair>& pqueue, VertexID queryVertexID, VertexID vnsID, 
                                                        EdgeWeight spdist, google::dense_hash_set<VertexID>& neighbourVnsIDsAdded)
{
    VertexID neighbourVnsID, neighbourObjectID;
    EdgeWeight neighbourObjectLB, neighbourObjectUB, neighbourDist, nvdLB;
    int numCandidates = 0;
    
    VertexID neighbourListStart = firstNeighbourIdx[vnsID];
    VertexID nextNeighbourListStart = firstNeighbourIdx[vnsID+1];

    for (VertexID i = neighbourListStart; i < nextNeighbourListStart; ++i) {
        neighbourVnsID = neighbourIDsAndDist[i].first;
        if (neighbourVnsIDsAdded.find(neighbourVnsID) == neighbourVnsIDsAdded.end()) {
            neighbourObjectID = vnsObjectIDs[neighbourVnsID];
            alt.getLowerAndUpperBound(queryVertexID,neighbourObjectID,neighbourObjectLB,neighbourObjectUB);
            // Retrieve the upper-bound at the same time so it can be saved and re-used later 
            // Otherwise an additional ALT call would be needed to filter by upper-bound
#if defined(COLLECT_STATISTICS)
            ++numCandidates; // Only need this for stats experiments
#endif
            neighbourDist = neighbourIDsAndDist[i].second;
            nvdLB = (spdist > neighbourDist) ? spdist-neighbourDist : 0; // The other side is zero because neighbourDist is only an upperbound
            if (nvdLB > neighbourObjectLB) {
                neighbourObjectLB = nvdLB;
            }
            pqueue.insert(VertexEdgeWeightPair(neighbourVnsID,neighbourObjectUB),neighbourObjectLB);
            neighbourVnsIDsAdded.insert(neighbourVnsID);
        }
    }

    return numCandidates;
}

int NVDQuadtree::generateNeighbourCandidatesAkNN(ALT& alt, BinaryMinHeap<EdgeWeight,VertexPair>& pqueue, VertexID queryVertexID, VertexID vnsID, 
                                                  EdgeWeight spdist, google::dense_hash_set<VertexID>& neighbourVnsIDsAdded, unsigned int queryVertexIdx)
{
    VertexID neighbourVnsID, neighbourObjectID;
    EdgeWeight neighbourObjectLB, neighbourDist, nvdLB;
    int numCandidates = 0;
    
    VertexID neighbourListStart = firstNeighbourIdx[vnsID];
    VertexID nextNeighbourListStart = firstNeighbourIdx[vnsID+1];

    for (VertexID i = neighbourListStart; i < nextNeighbourListStart; ++i) {
        neighbourVnsID = neighbourIDsAndDist[i].first;
        if (neighbourVnsIDsAdded.find(neighbourVnsID) == neighbourVnsIDsAdded.end()) {
            neighbourObjectID = vnsObjectIDs[neighbourVnsID];
            neighbourObjectLB = alt.getLowerBound(queryVertexID,neighbourObjectID);
#if defined(COLLECT_STATISTICS)
            ++numCandidates; // Only need this for stats experiments
#endif
            neighbourDist = neighbourIDsAndDist[i].second;
            nvdLB = (spdist > neighbourDist) ? spdist-neighbourDist : 0; // The other side is zero because neighbourDist is only an upperbound
            if (nvdLB > neighbourObjectLB) {
                neighbourObjectLB = nvdLB;
            }
            pqueue.insert(VertexPair(neighbourVnsID,queryVertexIdx),neighbourObjectLB);
            neighbourVnsIDsAdded.insert(neighbourVnsID);
        }
    }

    return numCandidates;
}