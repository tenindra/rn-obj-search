// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "ALT.h"

#include "SetGenerator.h"
#include "DijkstraSearch.h"
#include "../queue/BinaryMinHeap.h"

#include <iostream>

ALT::ALT(std::string _networkName, int _numVertices, int _numEdges): 
    networkName(_networkName), numVertices(_numVertices), numEdges(_numEdges) {
}
    
std::string ALT::getNetworkName()
{
    return networkName;
}

unsigned int ALT::getNumEdges()
{
    return numEdges;
}

unsigned int ALT::getNumVertices()
{
    return numVertices;
}

unsigned int ALT::getNumLandmarks()
{
    return numLandmarks;
}

double ALT::computeIndexSize()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(VertexID)*landmarks.size();
    memoryUsage += sizeof(EdgeWeight)*vertexFromLandmarkDistances.size();
    return memoryUsage/(1024*1024);
}

void ALT::buildALT(Graph& graph, LANDMARK_TYPE landmarkType, unsigned int _numLandmarks)
{
    unsigned int numVertices = graph.getNumVertices();
    numLandmarks = _numLandmarks;

    // Every random draw of the selection comes from one engine with a fixed seed, so a rebuild selects the same landmarks. Random
    // landmarks then depend only on the number of vertices and landmarks, so a network with other edge weights gets the same ones.
    // To select random landmarks again, seed the engine with std::random_device{}() instead; changing LANDMARK_SEED selects
    // other fixed landmarks
    std::mt19937 engine(LANDMARK_SEED);
    if (landmarkType == LANDMARK_TYPE::RANDOM) {
        SetGenerator sg;
        std::vector<VertexID> randomVertices = sg.generateRandomSampleSet(numVertices, numLandmarks, engine);
        landmarks.swap(randomVertices);
    } else if (landmarkType == LANDMARK_TYPE::FARTHEST) {
        landmarks = this->selectFarthestLandmarks(graph,numLandmarks,engine);
    } else if (landmarkType == LANDMARK_TYPE::FARTHEST_BFS) {
        landmarks = this->selectFarthestBFSLandmarks(graph,numLandmarks,engine);
    } else {
        std::cerr << "Unknown landmark type" << std::endl;
        std::exit(1);
    }
    
    // Allocated space for distance from each landmark to each graph vertex
    vertexFromLandmarkDistances.resize(numLandmarks*numVertices);
    
    // Use Dijkstra's to populate above vector for each landmark
    DijkstraSearch dijk;
    std::vector<EdgeWeight> landmarkDistances(numVertices,0);
    BinaryMinHeap<EdgeWeight,VertexID> *pqueue = new BinaryMinHeap<EdgeWeight,VertexID>();
    for (std::size_t i = 0; i < landmarks.size(); ++i) {
        pqueue->clear();
        dijk.findSSSPDistances(graph,landmarks[i],landmarkDistances,pqueue);
        for (std::size_t j = 0; j < numVertices; ++j) {
            vertexFromLandmarkDistances[j*numLandmarks+i] = landmarkDistances[j];
            // Note: This will not be cache efficient, because the next vector position
            // we write to will be numLandmarks away. But when we access vertexFromLandmarkDistances
            // vector later, it will be in cache efficient order so queries will be faster.
        }
    }
    delete pqueue;
}

std::vector<VertexID> ALT::selectFarthestLandmarks(Graph& graph, unsigned int numLandmarks, std::mt19937& engine)
{
    std::vector<VertexID> landmarks;
    std::unordered_set<VertexID> landmarksUset;
    
    // Pick a random vertex
    SetGenerator sg;
    std::vector<VertexID> randomVertices = sg.generateRandomSampleSet(numVertices,1,engine);
    VertexID randomVertex = randomVertices[0];
    
    DijkstraSearch dijk;
    EdgeWeight farthestDist;
    VertexID candidate = dijk.findFarthestVertex(graph,randomVertex,farthestDist), newCandidate;
    landmarks.push_back(candidate);
    landmarksUset.insert(candidate);
    while (landmarks.size() < numLandmarks) {
        newCandidate = dijk.findFarthestVertex(graph,candidate,farthestDist);
        if (landmarksUset.find(newCandidate) == landmarksUset.end()) {
            landmarks.push_back(newCandidate);
            landmarksUset.insert(newCandidate);
            candidate = newCandidate;
        } else {
            // If the farthest is already in the set, then we need need restart randomly
            std::cout << "Existing landmark found as farthest again\n";
            randomVertices = sg.generateRandomSampleSet(numVertices,1,engine);
            candidate = randomVertices[0];
        }
    }
    
    return landmarks;
}

std::vector<VertexID> ALT::selectFarthestBFSLandmarks(Graph& graph, unsigned int numLandmarks, std::mt19937& engine)
{
    std::vector<VertexID> landmarks;
    std::unordered_set<VertexID> landmarksUset;
    
    // Pick a random vertex
    SetGenerator sg;
    std::vector<VertexID> randomVertices = sg.generateRandomSampleSet(numVertices,1,engine);
    VertexID randomVertex = randomVertices[0];
    
    DijkstraSearch dijk;
    VertexID candidate = dijk.findFarthestVertexByBFS(graph,randomVertex), newCandidate;
    landmarks.push_back(candidate);
    landmarksUset.insert(candidate);
    while (landmarks.size() < numLandmarks) {
        newCandidate = dijk.findFarthestVertexByBFS(graph,candidate);
        if (landmarksUset.find(newCandidate) == landmarksUset.end()) {
            landmarks.push_back(newCandidate);
            landmarksUset.insert(newCandidate);
            candidate = newCandidate;
        } else {
            // If the farthest is already in the set, then we need need restart randomly
            std::cout << "Existing landmark found as farthest by BFS again\n";
            randomVertices = sg.generateRandomSampleSet(numVertices,1,engine);
            candidate = randomVertices[0];
        }
    }
    
    return landmarks;
}

EdgeWeight ALT::getLowerBound(VertexID s, VertexID t)
{
    EdgeWeight globalLB = 0, currentLB;
    EdgeWeight sourceLBDist, targetLBDist;
    // Reading the two rows through local pointers lets GCC vectorise this loop: GCC 15 does not when the loop indexes
    // vertexFromLandmarkDistances, because it then reloads numLandmarks and the vector's data pointer in every iteration
    const EdgeWeight* sourceDists = vertexFromLandmarkDistances.data()+s*numLandmarks;
    const EdgeWeight* targetDists = vertexFromLandmarkDistances.data()+t*numLandmarks;
    const std::size_t numLmks = landmarks.size();
    for (std::size_t i = 0; i < numLmks; ++i) {
        sourceLBDist = sourceDists[i];
        targetLBDist = targetDists[i];
        currentLB = (sourceLBDist > targetLBDist) ? sourceLBDist-targetLBDist : targetLBDist - sourceLBDist;
        if (currentLB > globalLB) {
            globalLB = currentLB;
        }
    }
    return globalLB;
}

EdgeWeight ALT::getUpperBound(VertexID s, VertexID t)
{
    EdgeWeight globalUB = std::numeric_limits<EdgeWeight>::max(), currentUB;
    // Rows read through local pointers to keep the loop vectorised (see getLowerBound)
    const EdgeWeight* sourceDists = vertexFromLandmarkDistances.data()+s*numLandmarks;
    const EdgeWeight* targetDists = vertexFromLandmarkDistances.data()+t*numLandmarks;
    const std::size_t numLmks = landmarks.size();
    for (std::size_t i = 0; i < numLmks; ++i) {
        currentUB = sourceDists[i]+targetDists[i];
        if (currentUB < globalUB) {
            globalUB = currentUB;
        }
    }
    return globalUB;
}

void ALT::getLowerAndUpperBound(VertexID s, VertexID t, EdgeWeight& lb, EdgeWeight& ub)
{
    EdgeWeight globalLB = 0, globalUB = std::numeric_limits<EdgeWeight>::max(), currentLB, currentUB;
    EdgeWeight sourceLBDist, targetLBDist;
    // Rows read through local pointers to keep the loop vectorised (see getLowerBound)
    const EdgeWeight* sourceDists = vertexFromLandmarkDistances.data()+s*numLandmarks;
    const EdgeWeight* targetDists = vertexFromLandmarkDistances.data()+t*numLandmarks;
    const std::size_t numLmks = landmarks.size();
    for (std::size_t i = 0; i < numLmks; ++i) {
        sourceLBDist = sourceDists[i];
        targetLBDist = targetDists[i];
        currentLB = (sourceLBDist > targetLBDist) ? sourceLBDist-targetLBDist : targetLBDist - sourceLBDist;
        currentUB = sourceLBDist+targetLBDist;
        if (currentLB > globalLB) {
            globalLB = currentLB;
        }
        if (currentUB < globalUB) {
            globalUB = currentUB;
        }
    }
    lb = globalLB;
    ub = globalUB;
}

LANDMARK_TYPE ALT::getLandmarkType(std::string selectionMethod, bool& success)
{
    if (selectionMethod == "random") {
        success = true;
        return LANDMARK_TYPE::RANDOM;
    } else if (selectionMethod == "farthest") {
        success = true;
        return LANDMARK_TYPE::FARTHEST;
    } else if (selectionMethod == "farthest_bfs") {
        success = true;
        return LANDMARK_TYPE::FARTHEST_BFS;
    } else {
        // Return default method but indicate method was not found in success
        success = false; 
        return LANDMARK_TYPE::RANDOM;
    }
}

VertexID ALT::getNearestLandmarkIndex(VertexID s, EdgeWeight& dist)
{
    EdgeWeight minLMDist = vertexFromLandmarkDistances[s*numLandmarks];
    VertexID minLMIdx = 0;
    for (std::size_t i = 1; i < landmarks.size(); ++i) {
        if (vertexFromLandmarkDistances[s*numLandmarks+i] < minLMDist) {
            minLMDist = vertexFromLandmarkDistances[s*numLandmarks+i];
            minLMIdx = i;
        }
    }
    dist = minLMDist;
    return minLMIdx;
}

/*
 * ObjectList
 */

ObjectList::ObjectList(std::string setType, double setDensity, int setVariable, int setSize):
                             objSetType(setType), objSetDensity(setDensity), objSetVariable(setVariable), objSetSize(setSize)
{
}

// Assumes that buildObjectList has not be called before
void ObjectList::buildObjectList(ALT& alt, std::vector<VertexID>& objectIDs)
{
    unsigned int numLandmarks = alt.getNumLandmarks();
    this->lmObjDistances.resize(numLandmarks);
    for (std::size_t i = 0; i < alt.landmarks.size(); ++i) {
        for (VertexID object: objectIDs) {
            this->lmObjDistances[i].push_back(VertexDistancePair(object,alt.vertexFromLandmarkDistances[object*numLandmarks+i]));
        }
        std::sort(this->lmObjDistances[i].begin(), this->lmObjDistances[i].end(), [](const VertexDistancePair &left, const VertexDistancePair &right) {
            return left.second < right.second;
        });
    }
}

std::string ObjectList::getObjSetType()
{
    return this->objSetType;
}

double ObjectList::getObjSetDensity()
{
    return this->objSetDensity;
}

int ObjectList::getObjSetSize()
{
    return this->objSetSize;
}

int ObjectList::getObjSetVariable()
{
    return this->objSetVariable;
}


double ObjectList::computeIndexSize()
{
    double memoryUsage = 0;
    for (std::size_t i = 0; i < this->lmObjDistances.size(); ++i) {
        memoryUsage += sizeof(VertexDistancePair)*this->lmObjDistances[i].size();
    }
    return memoryUsage/(1024*1024);
}

int ObjectList::getClosestDistanceIndex(std::vector<VertexDistancePair>& distances, EdgeWeight matchDist)
{
    int left = 0, right = distances.size()-1, mid;
    while(left <= right) {
        mid = (right+left)/2;
        if (distances[mid].second == matchDist) {
            return mid;
        } else if (distances[mid].second < matchDist) {
            left = mid+1;
        } else /* if (distances[mid].second > matchDist)*/ {
            right = mid-1;
        }
    }
    
    // At this point left is greater than right by 1 
    // One of left or right is the closest value
    EdgeWeight leftDiff, rightDiff;
    if (static_cast<std::size_t>(left) < distances.size() && right >= 0) {
        leftDiff = (distances[left].second > matchDist) ? distances[left].second - matchDist : matchDist - distances[left].second;
        rightDiff = (distances[right].second > matchDist) ? distances[right].second - matchDist : matchDist - distances[right].second;
        if (leftDiff < rightDiff) {
            return left;
        } else {
            return right;
        }
    } else if (static_cast<std::size_t>(left) >= distances.size()) {
        return right;
    } else {
        return left;
    }

}

// Assumes candidateLMIndex is the landmark associated with the query
BinaryMinHeap<EdgeWeight,VertexID> ObjectList::getNearestCandidateByQueryLandmark(ALT& alt, VertexID queryVertex, EdgeWeight queryVertexLMDistance, 
                                                                                int& leftPos, int& rightPos, VertexID& candidate, VertexID& candidateLMIndex)
{
    BinaryMinHeap<EdgeWeight,VertexID> pqueue;

    int minObjPos = getClosestDistanceIndex(lmObjDistances[candidateLMIndex],queryVertexLMDistance);
    VertexID minObject = lmObjDistances[candidateLMIndex][minObjPos].first;
    EdgeWeight minObjLBDist = alt.getLowerBound(queryVertex,minObject);
    pqueue.insert(minObject,minObjLBDist);

    leftPos = (minObjPos-1 >= 0) ? minObjPos-1 : -1;
    rightPos = (static_cast<std::size_t>(minObjPos+1) < lmObjDistances[candidateLMIndex].size()) ? minObjPos+1 : lmObjDistances[candidateLMIndex].size();
    
    fullyExploreSingleObjectList(candidateLMIndex,alt,pqueue,queryVertex,queryVertexLMDistance,leftPos,rightPos);
    
    // Extract minimum element as the candidate
    // Note: We don't save minKey because this is the first object
    // and must be exactly evaluated anyway
    candidate = pqueue.extractMinElement();
    
    return pqueue;
}

// Return false if there are no more candidates remaining
// Assumes that candidateLMIndex has not change since last candidate retrieval
bool ObjectList::getNextNearestCandidateByQueryLandmark(ALT& alt, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, VertexID queryVertex, 
                                                        EdgeWeight queryVertexLMDistance, int& leftPos, int& rightPos, VertexID& candidate, 
                                                        VertexID& candidateLMIndex, EdgeWeight& candidateLBDist)
{
    fullyExploreSingleObjectList(candidateLMIndex,alt,pqueue,queryVertex,queryVertexLMDistance,leftPos,rightPos);

    if (pqueue.size() > 0) {
        // Extract minimum element as the candidate
        candidateLBDist = pqueue.getMinKey();
        candidate = pqueue.extractMinElement();
        return true;
    } else {
        return false;
    }
}

void ObjectList::fullyExploreSingleObjectList(VertexID landmarkIndex, ALT& alt, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, VertexID queryVertex, 
                                              EdgeWeight queryVertexLMDistance, int& leftPos, int& rightPos)
{
    bool continueLeft, continueRight;
    VertexID minObject;
    EdgeWeight minObjLBDist, objectLMDist;

    continueLeft = (leftPos >= 0) ? true : false;
    continueRight = (static_cast<std::size_t>(rightPos) < lmObjDistances[landmarkIndex].size()) ? true : false;
    while (continueLeft || continueRight) {
        if (continueLeft) {
            minObject = lmObjDistances[landmarkIndex][leftPos].first;
            objectLMDist = lmObjDistances[landmarkIndex][leftPos].second;
            minObjLBDist = (objectLMDist > queryVertexLMDistance) ? objectLMDist - queryVertexLMDistance : queryVertexLMDistance - objectLMDist;
            if (pqueue.size() == 0 || minObjLBDist < pqueue.getMinKey()) {
                minObjLBDist = alt.getLowerBound(queryVertex,minObject);
                pqueue.insert(minObject,minObjLBDist);
                --leftPos;
                if (leftPos < 0) {
                    continueLeft = false;
                }
            } else {
                continueLeft = false;
            }
        }
        if (continueRight) {
            minObject = lmObjDistances[landmarkIndex][rightPos].first;
            objectLMDist = lmObjDistances[landmarkIndex][rightPos].second;
            minObjLBDist = (objectLMDist > queryVertexLMDistance) ? objectLMDist - queryVertexLMDistance : queryVertexLMDistance - objectLMDist;
            if (pqueue.size() == 0 || minObjLBDist < pqueue.getMinKey()) {
                minObjLBDist = alt.getLowerBound(queryVertex,minObject);
                pqueue.insert(minObject,minObjLBDist);
                ++rightPos;
                if (static_cast<std::size_t>(rightPos) == lmObjDistances[landmarkIndex].size()) {
                    continueRight = false;
                }
            } else {
                continueRight = false;
            }            
        }
    }
}
