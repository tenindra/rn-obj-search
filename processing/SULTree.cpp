// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "SULTree.h"

#include "DijkstraSearch.h"
#include "PackedRtree.h"
#include "../queue/RadixHeap.h"
#include "../utility/utility.h"
#include "../utility/StopWatch.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <limits>
#include <random>
#include <stack>
#include <unordered_map>
#include <tuple>

SULTree::SULTree(std::string networkName, int numVertices, int numEdges, int numNodeLandmarks):
    networkName(networkName), numVertices(numVertices), numEdges(numEdges), fanout(0), maxLeafSize(0), numNodeLandmarks(numNodeLandmarks) {
    // Note: fanout, maxLeafSize and the vertex reordering are taken from the HierarchicalPartition given to buildSULTree
}

SULTree::SULTree() {}

std::string SULTree::getNetworkName()
{
    return this->networkName;
}

int SULTree::getNumEdges()
{
    return this->numEdges;
}

int SULTree::getNumVertices()
{
    return this->numVertices;
}

int SULTree::getFanout()
{
    return this->fanout;
}

std::size_t SULTree::getMaxLeafPartitionSize()
{
     return this->maxLeafSize;
}

int SULTree::getNumLandmarksPerNode()
{
    return this->numNodeLandmarks;
}

double SULTree::getLandmarkTimeMs()
{
    return this->landmarkTimeMs;
}

double SULTree::getDistanceListTimeMs()
{
    return this->distanceListTimeMs;
}

double SULTree::getAllocationTimeMs()
{
    return this->allocationTimeMs;
}

// Build the tree on a hierarchical partition of the graph: landmark selection and distance list computation per node
void SULTree::buildSULTree(Graph& graph, const HierarchicalPartition& partition, COLTREE_LANDMARK_TYPE landmarkType, SDL_METHOD sdlMethod)
{
    this->landmarkType = landmarkType;
    this->sdlMethod = sdlMethod;
    if (partition.getNetworkName() != this->networkName || partition.getNumVertices() != this->numVertices) {
        std::cerr << "Partition was built for a different road network (" << partition.getNetworkName() << ")" << std::endl;
        std::exit(1);
    }

    StopWatch sw;
    double landmarkTime = 0, distanceListTime = 0;

    // Copy the hierarchy and the vertex reordering. During construction the vertices of a node are the range
    // reorderedVertexIDtoOriginalVertexID[startVertexID..endVertexID] of the partition's vertex order
    sw.start();
    this->fanout = partition.getFanout();
    this->maxLeafSize = partition.getMaxLeafSize();
    this->reorderedVertexIDtoOriginalVertexID = partition.reorderedToOriginal;
    this->originalVertexIDtoReorderedVertexID = partition.originalToReordered;
    RN_DEBUG_ASSERT(this->treeNodes.empty() && "A SULTree is built once, after it is constructed");
    this->treeNodes.reserve(partition.nodes.size());
    for (std::size_t i = 0; i < partition.nodes.size(); ++i) {
        const HierarchicalPartitionNode& pnode = partition.nodes[i];
        SULTreeNode node(pnode.parentIdx,pnode.numVertices);
        node.isLeaf = pnode.isLeaf;
        node.startVertexID = pnode.startVertexID;
        node.endVertexID = pnode.endVertexID;
        node.childNodeIdxs = pnode.childNodeIdxs;
        node.borderVertexIDs = pnode.borderVertexIDs;
        node.minX = pnode.minX;
        node.minY = pnode.minY;
        node.maxX = pnode.maxX;
        node.maxY = pnode.maxY;
        this->treeNodes.push_back(std::move(node));
    }
    sw.stop();
    std::cout << "Partitioning time = " << partition.getBuildTimeMs() << "ms (pre-computed), partition copy time = " << sw.getTimeMs() << "ms" << std::endl;

    // The distance lists of every node are written in place in nodeVertexLmkDists, so the index is held once during construction
    sw.reset();
    sw.start();
    this->allocateNodeVertexLmkDists();
    sw.stop();
    this->allocationTimeMs = sw.getTimeMs();
    std::cout << "Distance list allocation time = " << sw.getTimeMs() << "ms" << std::endl;

    // Process root node outside of loop so it can be parallelized more easily
    // (will read information from root node, e.g., to find min/max borders)
    sw.reset();
    sw.start();
    this->selectNodeLandmarks(graph,0);
    this->computeRootDistanceLists(graph);
    sw.stop();
    std::cout << "Root processing time = " << sw.getTimeMs() << "ms" << std::endl;

    sw.reset();
    sw.start();

    // Note: Each node only reads the root's lists and writes its own, so the nodes can be processed in parallel. The pragma below is
    // compiled only into rn_obj_search_par (PARALLEL_BUILDS, see common.h), which builds the parameter-test SUL-Trees. In rn_obj_search,
    // whose construction times are reported, and in rn_obj_search_stats it expands to nothing and the nodes are processed serially.
    // In rn_obj_search_par the number of threads comes from the caller: LMKExperimentsCommand::buildSULTreeParameterIndexes sets it with
    // omp_set_num_threads in each of its threads, as otherwise every SUL-Tree being built would start a team of all available threads

    // Note: With the pragma, landmarkTime and distanceListTime are summed over the threads (CPU time rather than wall time), so the
    // landmark selection and distance list times that rn_obj_search_par prints are not wall times

    PARALLEL_BUILDS_PRAGMA(omp parallel for schedule(dynamic) reduction(+:landmarkTime,distanceListTime))
    for (TreeNodeID treeIdx = 1; treeIdx < this->treeNodes.size(); ++treeIdx) {
        StopWatch sw2;
        sw2.reset();
        sw2.start();
        this->selectNodeLandmarks(graph,treeIdx);
        sw2.stop();
        landmarkTime += sw2.getTimeMs();

        sw2.reset();
        sw2.start();
        this->computeNodeDistanceLists(graph,treeIdx);
        sw2.stop();
        distanceListTime += sw2.getTimeMs();
    }

    sw.stop();
    std::cout << "Non-root node processing time = " << sw.getTimeMs() << "ms" << std::endl;

    std::cout << "Landmark selection time = " << landmarkTime << "ms" << std::endl;
    std::cout << "Distance list computation time = " << distanceListTime << "ms" << std::endl;
    this->landmarkTimeMs = landmarkTime;
    this->distanceListTimeMs = distanceListTime;

    // Spillover factors are kept by each node during construction so that nodes can be processed in parallel
    this->nodeSpilloverFactors.resize(this->treeNodes.size());
    for (TreeNodeID treeIdx = 0; treeIdx < this->treeNodes.size(); ++treeIdx) {
        this->nodeSpilloverFactors[treeIdx] = this->treeNodes[treeIdx].spilloverFactors;
    }
}

// Allocates nodeVertexLmkDists for the distance lists of all tree nodes and sets where each node's lists start (nodeStartIndex).
// Nodes are stored in tree index order, each with the distances to its numNodeLandmarks landmarks for each of its vertices.
void SULTree::allocateNodeVertexLmkDists()
{
    RN_DEBUG_ASSERT(this->nodeStartIndex.empty() && this->nodeVertexLmkDistsSize == 0 && "The distance lists are allocated once, by buildSULTree");
    this->nodeStartIndex.reserve(this->treeNodes.size());
    for (TreeNodeID treeIdx = 0; treeIdx < this->treeNodes.size(); ++treeIdx) {
        this->nodeStartIndex.push_back(this->nodeVertexLmkDistsSize);
        this->nodeVertexLmkDistsSize += static_cast<long>(this->treeNodes[treeIdx].getNumVertices())*this->numNodeLandmarks;
    }

    if (this->nodeVertexLmkDistsSize > std::numeric_limits<int>::max()) {
        std::cout << "WARNING: Required size of nodeVertexLmkDistsSize is larger than max integer" << std::endl;
        std::cout << "\t Required std::vector size: " << nodeVertexLmkDistsSize << std::endl;
    }

    this->nodeVertexLmkDists.assign(this->nodeVertexLmkDistsSize,0);
}

// Assumes 
// - 0th element in treeNodes is root
// - Root node vertex distance lists have been computed
// - treeIdx function argument is not root
// - Landmarks for this node have been selected (the lists of landmarks marked in isDistanceListComputed are already stored)
void SULTree::computeNodeDistanceLists(Graph& graph, TreeNodeID treeIdx) 
{
    [[maybe_unused]] int parentIdx = this->treeNodes[treeIdx].getParentIdx();
    RN_DEBUG_ASSERT(parentIdx != this->ROOT_PARENT_INDEX && "computeNodeDistanceLists cannot process root node");

    // For this node, pre-compute the minimum and maximum distances from the root node landmarks
    // to this node's border vertices so that Dijkstra search can be optimized later
    this->computeBorderDistanceBounds(treeIdx);

    // Compute distances from each landmark to all vertices in node (subgraph)
    // then store in node data structure for object index construction
    int numVerticesSettled;

    // The hash-table based method needs the node's vertices as a set
    std::unordered_set<VertexID> subgraph;
    if (this->sdlMethod == SDL_METHOD::SDL_HASH_TABLE) {
        subgraph.insert(this->reorderedVertexIDtoOriginalVertexID.begin()+this->treeNodes[treeIdx].startVertexID,
                        this->reorderedVertexIDtoOriginalVertexID.begin()+this->treeNodes[treeIdx].endVertexID+1);
    }

    long nodeStart = this->nodeStartIndex[treeIdx]; // Start of this node's lists in nodeVertexLmkDists
    for (std::size_t i = 0; i < this->treeNodes[treeIdx].landmarks.size(); ++i) {
        if (i < this->treeNodes[treeIdx].isDistanceListComputed.size() && this->treeNodes[treeIdx].isDistanceListComputed[i]) {
            continue; // Stored by the landmark selection, from the same search (see findFurthestBorderSubset)
        }
        VertexID landmarkVertexID = this->treeNodes[treeIdx].landmarks[i];
        std::vector<EdgeWeight> distancesVec(this->treeNodes[treeIdx].getNumVertices()); // Position j is the vertex with reordered ID startVertexID+j

        if (this->sdlMethod == SDL_METHOD::SDL_BORDER_SET) {
            numVerticesSettled = this->findSSMTDistancesViaBorderSet(graph,treeIdx,landmarkVertexID,distancesVec);
        } else if (this->sdlMethod == SDL_METHOD::SDL_DIJKSTRA) {
            numVerticesSettled = this->findSSMTDistances(graph,treeIdx,landmarkVertexID,distancesVec);
        } else {
            DijkstraSearch dijk;
            std::unordered_map<VertexID,EdgeWeight> distances;
            numVerticesSettled = dijk.findSSMTDistances(graph,landmarkVertexID,subgraph,distances);
            for (unsigned int j = 0; j < this->treeNodes[treeIdx].getNumVertices(); ++j) {
                distancesVec[j] = distances[this->reorderedVertexIDtoOriginalVertexID[this->treeNodes[treeIdx].startVertexID+j]];
            }
        }

        this->treeNodes[treeIdx].spilloverFactors.push_back(numVerticesSettled/static_cast<double>(this->treeNodes[treeIdx].getNumVertices()));

        for (unsigned int j = 0; j < this->treeNodes[treeIdx].getNumVertices(); ++j) {
            // Store landmark distance for each vertex together because that is how 
            // they will access later (for cache performance benefits)
            int vertexLmkDistIdx = j*this->treeNodes[treeIdx].landmarks.size()+i; // int is fine here, only accessing node list (which is at most |V|)
            this->nodeVertexLmkDists[nodeStart+vertexLmkDistIdx] = distancesVec[j];
        }
    }
}

// Computes the minimum and maximum distances from each root landmark to the node's border vertices (minBorderDist and
// maxBorderDist), which the border-set search (findSSMTDistancesViaBorderSet) reads. Assumes the root's distance lists have
// been computed and treeIdx is not the root
void SULTree::computeBorderDistanceBounds(TreeNodeID treeIdx)
{
    int numRootLandmarks, borderVertexDistIdx;
    // Initialise min and max border distances for this node for each root landmark
    this->treeNodes[treeIdx].maxBorderDist.assign(this->treeNodes[0].landmarks.size(),0);
    this->treeNodes[treeIdx].minBorderDist.assign(this->treeNodes[0].landmarks.size(),std::numeric_limits<EdgeWeight>::max());
    // For each border, find distance to each root landmark and update min/max distances
    for (std::size_t i = 0; i < this->treeNodes[treeIdx].borderVertexIDs.size(); ++i) {
        VertexID borderVertexID = this->treeNodes[treeIdx].borderVertexIDs[i];
        for (std::size_t j = 0; j < this->treeNodes[0].landmarks.size(); ++j) {
            numRootLandmarks = this->treeNodes[0].landmarks.size();
            // Vertex node position for root's landmark distance lists is just its ID (no reordering)
            borderVertexDistIdx = borderVertexID*numRootLandmarks+j;
            EdgeWeight landmarkDist = this->nodeVertexLmkDists[this->nodeStartIndex[0]+borderVertexDistIdx];
            if (landmarkDist > this->treeNodes[treeIdx].maxBorderDist[j]) {
                this->treeNodes[treeIdx].maxBorderDist[j] = landmarkDist;
            }
            if (landmarkDist < this->treeNodes[treeIdx].minBorderDist[j]) {
                this->treeNodes[treeIdx].minBorderDist[j] = landmarkDist;
            }
        }
    }
}

// Select landmarks and compute node distance lists for root node separately
// using some assumptions that can be applied since root is just original graph
// Assumes: Root node has not been process before and landmarks have been selected
// Note: The landmarks' searches are independent, so in rn_obj_search_par (PARALLEL_BUILDS, see common.h) they run in parallel, with
// the threads of the caller's team (as the node loop in buildSULTree). Each list is a column of the interleaved root lists, so the
// copies are made one at a time, as threads writing neighbouring columns at once would share cache lines. In rn_obj_search and
// rn_obj_search_stats both pragmas expand to nothing and the landmarks are processed serially
void SULTree::computeRootDistanceLists(Graph &graph)
{
    DijkstraSearch dijk; // Has no state, so the threads can share it
    PARALLEL_BUILDS_PRAGMA(omp parallel for schedule(dynamic))
    for (std::size_t i = 0; i < this->treeNodes[0].landmarks.size(); ++i) {
        VertexID landmarkVertexID = this->treeNodes[0].landmarks[i];
        std::vector<EdgeWeight> distances(graph.getNumVertices());
        dijk.findSSSPDistances(graph,landmarkVertexID,distances);
        PARALLEL_BUILDS_PRAGMA(omp critical(rootDistanceListCopy))
        for (VertexID j = 0; j < static_cast<VertexID>(graph.getNumVertices()); ++j) {
            int vertexLmkDistIdx = j*this->treeNodes[0].landmarks.size()+i; // int is fine here, only accessing node list (which is at most |V|)
            this->nodeVertexLmkDists[this->nodeStartIndex[0]+vertexLmkDistIdx] = distances[j];
        }
    }
    // Note: this->nodeSpilloverFactor does not include root node because is always one.
}

// Assumes minX/minY/maxX/maxY set for treeIdx (including root node) for some landmark types
void SULTree::selectNodeLandmarks(Graph& graph, TreeNodeID treeIdx) 
{
    // Choose landmarks and compute vertex distances if node is not root
    if (this->landmarkType == COLTREE_LANDMARK_TYPE::SUBGRAPH_RANDOM) {
        this->selectRandomLandmarks(treeIdx);
    } else if (this->landmarkType == COLTREE_LANDMARK_TYPE::SUBGRAPH_SLICE_FURTHEST) {
        this->selectSlicedFurthestLandmarks(graph,treeIdx,
            this->treeNodes[treeIdx].minX,this->treeNodes[treeIdx].minY,
            this->treeNodes[treeIdx].maxX,this->treeNodes[treeIdx].maxY);
    } else if (this->landmarkType == COLTREE_LANDMARK_TYPE::SUBGRAPH_FURTHEST_BORDER) {
        this->selectRandomFurthestBorderLandmarks(graph,treeIdx);
    } else if (this->landmarkType == COLTREE_LANDMARK_TYPE::SUBGRAPH_BORDER_MINMAX) {
        if (treeIdx != 0) {
            // Doesn't make sense for root because it doesn't have any borders
            this->selectMinmaxLandmarksFromBorders(graph,treeIdx);
        } else {
            this->selectRandomLandmarks(treeIdx);
        }
    } else {
        std::cerr << "Unknown landmark type" << std::endl;
        std::exit(1);
    }
}

double SULTree::computeIndexSize()
{
    double memoryUsage = 0;
    for (std::size_t i = 0; i < this->treeNodes.size(); ++i) {
        memoryUsage += this->treeNodes[i].computeIndexSizeBytes();
    }
    memoryUsage += sizeof(int)*this->nodeStartIndex.size(); // This is only long for safety, only extreme case during exhaustive parameters testing were bigger so we use int
    memoryUsage += sizeof(EdgeWeight)*this->nodeVertexLmkDists.size();

    // Note: We exclude originalVertexIDtoReorderedVertexID from index size because
    // it is not needed as you can simply reorder vertices permanently. This is fair
    // as we only use it to make implementation separate from other methods. Otherwise
    // we would have to rebuild all competitor indexes each time we run experiments for
    // SULTree, e.g., to have the same query vertex IDs. reorderedVertexIDtoOriginalVertexID
    // is excluded because it is not used except debugging.
    return memoryUsage/(1024*1024);
}

std::vector< std::vector<int>> SULTree::getTreeNodesByLevel()
{
    std::vector<std::vector<int>> treeNodeLevel;
    std::vector<int> currentLevel, nextLevel;
    nextLevel.push_back(0);

    while(nextLevel.size() != 0){
            treeNodeLevel.push_back(nextLevel);
            currentLevel.swap(nextLevel);
            nextLevel.clear();
            for (int i: currentLevel){
                for (std::size_t j = 0; j < this->treeNodes[i].childNodeIdxs.size(); j++ ){
                    nextLevel.push_back( this->treeNodes[i].childNodeIdxs[j] );
                }
            }
    }

    return treeNodeLevel;
}

int SULTree::getNumLevels()
{
    // Traverse level by level, counting levels without storing them
    int numLevels = 0;
    std::vector<int> currentLevel, nextLevel;
    nextLevel.push_back(0);
    while (nextLevel.size() != 0) {
        ++numLevels;
        currentLevel.swap(nextLevel);
        nextLevel.clear();
        for (int i: currentLevel) {
            nextLevel.insert(nextLevel.end(),this->treeNodes[i].childNodeIdxs.begin(),this->treeNodes[i].childNodeIdxs.end());
        }
    }
    return numLevels;
}

// Random number engine of a node's landmark selection. It is seeded from a fixed seed and the node, so a build from the same partition
// selects the same landmarks whatever order the nodes are processed in (the node loop is parallel in rn_obj_search_par)
std::mt19937 SULTree::getLandmarkEngine(int treeIdx)
{
    // To select random landmarks again, replace LANDMARK_SEED below with a seed drawn once per process:
    //     static const unsigned int runSeed = std::random_device{}();
    // The static is initialised once (thread-safe), so the builds of one process (e.g. the three of the SDL experiment) still select
    // the same landmarks, and runs differ; print runSeed to be able to reproduce a run. Changing LANDMARK_SEED selects other fixed landmarks
    std::seed_seq seed{LANDMARK_SEED,static_cast<unsigned int>(treeIdx)};
    return std::mt19937(seed);
}

void SULTree::selectSlicedFurthestLandmarks(Graph& originalGraph, int treeIdx, Coordinate minX, Coordinate minY, Coordinate maxX, Coordinate maxY)
{
    // The node's vertices are the range of the reordered vertex IDs (constant time random access)
    VertexID start = this->treeNodes[treeIdx].startVertexID;
    unsigned int numNodeVertices = this->treeNodes[treeIdx].getNumVertices();
    std::unordered_set<VertexID> verticesInSet;

    Coordinate centreX, centreY, x, y;
    
    std::mt19937 e1 = this->getLandmarkEngine(treeIdx);
    VertexID randomNum, randomVertexID;
    
    // Compute Euclidean centre of the subgraph
    centreX = (minX+maxX)/2;
    centreY = (minY+maxY)/2;

    // One slice for each landmark the node still needs (all of them when it has none yet)

    RN_DEBUG_ASSERT(this->treeNodes[treeIdx].landmarks.size() < this->numNodeLandmarks && "Already chosen required number of landmarks");
    std::size_t numSlices = this->numNodeLandmarks - this->treeNodes[treeIdx].landmarks.size();
    
    // Partition subgraph into slices and group subgraph vertices by slice
    std::vector<std::vector<VertexID>> sliceBorders(numSlices);
    std::vector<VertexID> sliceFurthestVertex(numSlices); // For each slice stores the furthest vertex
    std::vector<bool> sliceFurthestVertexSet(numSlices,false); // For each slice stores the furthest vertex
    std::vector<EuclideanDistanceSqrd> sliceFurthestVertexDist(numSlices); // For each slice stores the furthest vertex
    std::vector<Coordinate> sliceFurthestVertexAbsX(numSlices); // For each slice stores the furthest vertex
    std::vector<Coordinate> sliceFurthestVertexAbsY(numSlices); // For each slice stores the furthest vertex
    
    double angle;
    for (std::size_t i = 0; i < this->treeNodes[treeIdx].borderVertexIDs.size(); ++i) {
        originalGraph.getCoordinates(this->treeNodes[treeIdx].borderVertexIDs[i],x,y);
        angle = std::atan2(y-centreY,x-centreX);
        if (angle < 0) {
            angle += 2*3.1415926535897;
        }
        unsigned int sliceID = static_cast<unsigned int>(angle/(2.0*3.1415926535897)*numSlices)%numSlices;
        sliceBorders[sliceID].push_back(this->treeNodes[treeIdx].borderVertexIDs[i]);
    }
    for (unsigned int i = 0; i < numNodeVertices; ++i) {
        originalGraph.getCoordinates(this->reorderedVertexIDtoOriginalVertexID[start+i],x,y);
        angle = std::atan2(y-centreY,x-centreX);
        if (angle < 0) {
            angle += 2*3.1415926535897;
        }
        unsigned int sliceID = static_cast<unsigned int>(angle/(2.0*3.1415926535897)*numSlices)%numSlices;
        if (!sliceFurthestVertexSet[sliceID]) {
            sliceFurthestVertexSet[sliceID] = true;
            sliceFurthestVertex[sliceID] = this->reorderedVertexIDtoOriginalVertexID[start+i];
            sliceFurthestVertexAbsX[sliceID] = std::abs(x-centreX);
            sliceFurthestVertexAbsY[sliceID] = std::abs(y-centreY);
            sliceFurthestVertexDist[sliceID] = static_cast<EuclideanDistanceSqrd>(sliceFurthestVertexAbsX[sliceID])*sliceFurthestVertexAbsX[sliceID]
                                               + static_cast<EuclideanDistanceSqrd>(sliceFurthestVertexAbsY[sliceID])*sliceFurthestVertexAbsY[sliceID];
        } else {
            // For new distance to be bigger than the current furthest vertex
            // one of the absolute values of the relative difference in x and y
            // coordinates has to be greater
            Coordinate absX = std::abs(x-centreX);
            Coordinate absY = std::abs(y-centreY);
            if (absX > sliceFurthestVertexAbsX[sliceID] || absY > sliceFurthestVertexAbsY[sliceID]) {
                EuclideanDistanceSqrd newDist = static_cast<EuclideanDistanceSqrd>(absX)*absX + static_cast<EuclideanDistanceSqrd>(absY)*absY;
                if (newDist > sliceFurthestVertexDist[sliceID]) {
                    sliceFurthestVertex[sliceID] = this->reorderedVertexIDtoOriginalVertexID[start+i];
                    sliceFurthestVertexAbsX[sliceID] = absX;
                    sliceFurthestVertexAbsY[sliceID] = absY;
                    sliceFurthestVertexDist[sliceID] = newDist;
                }
            }
        }
    }        
    
    // Choose borders from each slice and choose non-border if a slice doesn't have any borders
    for (std::size_t i = 0; i < numSlices; ++i) {
        if (sliceBorders[i].size() > 0) {
            std::uniform_int_distribution<int> uniform_dist2(0, sliceBorders[i].size()-1);
            randomNum = uniform_dist2(e1);
            randomVertexID = sliceBorders[i][randomNum];
            RN_DEBUG_ASSERT(verticesInSet.find(randomVertexID) == verticesInSet.end() && "Shouldn't happen because each vertex is only in one slice");
            this->treeNodes[treeIdx].landmarks.push_back(randomVertexID);
            verticesInSet.insert(randomVertexID);

        } else if (sliceFurthestVertexSet[i]) {
            RN_DEBUG_ASSERT(verticesInSet.find(sliceFurthestVertex[i]) == verticesInSet.end() && "Shouldn't happen because each vertex is only in one slice");
            this->treeNodes[treeIdx].landmarks.push_back(sliceFurthestVertex[i]);
            verticesInSet.insert(sliceFurthestVertex[i]);
        }
    }
    
    // If all landmarks haven't been selected, pick some vertices in subgraph at random 
    if (this->treeNodes[treeIdx].landmarks.size() < this->numNodeLandmarks) {
        if (numNodeVertices < this->numNodeLandmarks) {
            std::cout << "Subgraph only has " << numNodeVertices << " vertices to select " << this->numNodeLandmarks << " from." << std::endl;
            std::cout << "Choose a larger maximum leaf size or smaller fanout to increase leaf vertices" << std::endl;
            std::exit(1);
        }
        std::uniform_int_distribution<int> uniform_dist2(0, numNodeVertices-1);
        while (this->treeNodes[treeIdx].landmarks.size() < this->numNodeLandmarks) {
            randomNum = uniform_dist2(e1);
            randomVertexID = this->reorderedVertexIDtoOriginalVertexID[start+randomNum];
            if (verticesInSet.find(randomVertexID) == verticesInSet.end()) {
                verticesInSet.insert(randomVertexID);
                this->treeNodes[treeIdx].landmarks.push_back(randomVertexID);
            }
        }
    }

    RN_DEBUG_ASSERT(this->treeNodes[treeIdx].landmarks.size() == this->numNodeLandmarks && "Did not find enough landmarks");
}

void SULTree::selectRandomLandmarks(int treeIdx)
{
    if (this->treeNodes[treeIdx].getNumVertices() < this->numNodeLandmarks) {
        std::cout << "Subgraph only has " << this->treeNodes[treeIdx].getNumVertices() << " vertices to select " << this->numNodeLandmarks << " from." << std::endl;
        std::cout << "Choose a larger maximum leaf size or smaller fanout to increase leaf vertices" << std::endl;
        std::exit(1);
    }

    std::unordered_set<VertexID> verticesInSet;
    
    // If some landmarks have already been set by another method
    if (this->treeNodes[treeIdx].landmarks.size() > 0) {
        verticesInSet.insert(this->treeNodes[treeIdx].landmarks.begin(),this->treeNodes[treeIdx].landmarks.end());
    }
    
    std::mt19937 e1 = this->getLandmarkEngine(treeIdx);
    VertexID randomNum, randomVertexID;
    
    if (this->treeNodes[treeIdx].landmarks.size() < this->numNodeLandmarks) {
        // The node's vertices are the range of the reordered vertex IDs (constant time random access)
        VertexID start = this->treeNodes[treeIdx].startVertexID;
        std::uniform_int_distribution<int> uniform_dist2(0, this->treeNodes[treeIdx].getNumVertices()-1);
        while (this->treeNodes[treeIdx].landmarks.size() < this->numNodeLandmarks) {
            randomNum = uniform_dist2(e1);
            randomVertexID = this->reorderedVertexIDtoOriginalVertexID[start+randomNum];
            if (verticesInSet.find(randomVertexID) == verticesInSet.end()) {
                verticesInSet.insert(randomVertexID);
                this->treeNodes[treeIdx].landmarks.push_back(randomVertexID);
            }
        }
    }

    RN_DEBUG_ASSERT(this->treeNodes[treeIdx].landmarks.size() == this->numNodeLandmarks && "Did not find enough landmarks");
}

// If the number of borders less than this->numNodeLandmarks, then any remaining borders are selected randomly
void SULTree::selectRandomFurthestBorderLandmarks(Graph& originalGraph, int treeIdx)
{
    if (this->treeNodes[treeIdx].borderVertexIDs.size() > this->numNodeLandmarks) {
        // Choose random seed vertex in from border vertex set
        std::mt19937 e1 = this->getLandmarkEngine(treeIdx);
        VertexID randomNum, randomVertexID;
        std::uniform_int_distribution<int> uniform_dist2(0,this->treeNodes[treeIdx].borderVertexIDs.size()-1);

        // Choose "dissimilar" borders by distances, starting with the seed
        randomNum = uniform_dist2(e1);
        randomVertexID = this->treeNodes[treeIdx].borderVertexIDs[randomNum];
        this->findFurthestBorderSubset(originalGraph,treeIdx,randomVertexID);

        RN_DEBUG_ASSERT(this->treeNodes[treeIdx].landmarks.size() == this->numNodeLandmarks && "Did not find enough landmarks");
    } else {
        this->treeNodes[treeIdx].landmarks = this->treeNodes[treeIdx].borderVertexIDs;
        this->selectRandomLandmarks(treeIdx);
    }
}

// Chooses the tree node's landmarks among its border vertices so that they lie far apart: the first is seedBorderID, and each
// next one the border vertex whose minimum distance to the landmarks already chosen is largest (ties broken by order in
// borderVertexIDs), so the second is the border vertex furthest from the seed. The distances come from the border-set search,
// one from each landmark except the last. That search is the search of the landmark's distance list, so when the distance
// lists are computed with the border-set search, they are stored here and marked in isDistanceListComputed, and
// computeNodeDistanceLists only computes the last landmark's list (the build's landmark selection time then includes the others)
// Assumes that treeIdx is not root node, it has more border vertices than landmarks and the root's distance lists have been computed
void SULTree::findFurthestBorderSubset(Graph& graph, int treeIdx, VertexID seedBorderID)
{
    std::vector<VertexID>& landmarks = this->treeNodes[treeIdx].landmarks;
    std::vector<VertexID>& borderVertexIDs = this->treeNodes[treeIdx].borderVertexIDs;
    VertexID startVertexID = this->treeNodes[treeIdx].startVertexID;
    int numNodeVertices = this->treeNodes[treeIdx].getNumVertices();
    std::vector<EdgeWeight> distances(numNodeVertices); // From the last landmark chosen, by position in the node
    std::vector<EdgeWeight> minDistances(borderVertexIDs.size(),std::numeric_limits<EdgeWeight>::max()); // For each border, to the landmarks
    bool storeDistanceLists = this->sdlMethod == SDL_METHOD::SDL_BORDER_SET;
    long nodeStart = this->nodeStartIndex[treeIdx]; // Start of this node's lists in nodeVertexLmkDists

    this->computeBorderDistanceBounds(treeIdx); // Read by the border-set search
    landmarks.assign(1,seedBorderID);
    this->treeNodes[treeIdx].isDistanceListComputed.assign(this->numNodeLandmarks,false);
    while (landmarks.size() < static_cast<std::size_t>(this->numNodeLandmarks)) {
        // Search from the last landmark chosen, and store it as that landmark's distance list
        int numVerticesSettled = this->findSSMTDistancesViaBorderSet(graph,treeIdx,landmarks.back(),distances);
        if (storeDistanceLists) {
            std::size_t landmarkIdx = landmarks.size()-1;
            this->treeNodes[treeIdx].spilloverFactors.push_back(numVerticesSettled/static_cast<double>(numNodeVertices));
            for (int j = 0; j < numNodeVertices; ++j) {
                this->nodeVertexLmkDists[nodeStart+j*this->numNodeLandmarks+landmarkIdx] = distances[j];
            }
            this->treeNodes[treeIdx].isDistanceListComputed[landmarkIdx] = true;
        }

        // The next landmark is the border with the largest minimum distance to the landmarks (a landmark is at distance zero
        // from itself, so it is never chosen again)
        EdgeWeight maxMinDist = 0;
        std::size_t nextBorderIdx = borderVertexIDs.size();
        for (std::size_t i = 0; i < borderVertexIDs.size(); ++i) {
            EdgeWeight borderDist = distances[this->originalVertexIDtoReorderedVertexID[borderVertexIDs[i]]-startVertexID];
            if (borderDist < minDistances[i]) {
                minDistances[i] = borderDist;
            }
            if (minDistances[i] > maxMinDist) {
                maxMinDist = minDistances[i];
                nextBorderIdx = i;
            }
        }
        if (nextBorderIdx == borderVertexIDs.size()) {
            // Cannot happen, as edge weights are positive and there are more borders than landmarks; stop rather than index past
            // the borders, and the caller's check of the number of landmarks reports it
            break;
        }
        landmarks.push_back(borderVertexIDs[nextBorderIdx]);
    }
}

// Assumes that treeIdx is not root node (conversely treeIdx must have at least 1 border vertex)
void SULTree::selectMinmaxLandmarksFromBorders(Graph& originalGraph, int treeIdx)
{
    // Randomly pick some borders, then find the subgraph vertices whose maximum distance to any of them is smallest

    std::mt19937 e1 = this->getLandmarkEngine(treeIdx);
    VertexID randomNum, randomBorderID;
    std::uniform_int_distribution<int> uniform_dist2(0,this->treeNodes[treeIdx].borderVertexIDs.size()-1);
    std::unordered_set<VertexID> chosenBorders;
    
    // Search from at least 5 borders but check total number of borders available
    std::size_t numBorders = (5 < this->treeNodes[treeIdx].borderVertexIDs.size()) ? 5 : this->treeNodes[treeIdx].borderVertexIDs.size();
    
    if (numBorders == this->treeNodes[treeIdx].borderVertexIDs.size()) {
        chosenBorders.insert(this->treeNodes[treeIdx].borderVertexIDs.begin(),this->treeNodes[treeIdx].borderVertexIDs.end());
    } else {
        for (std::size_t i = 0; i < numBorders; ++i) {
            randomNum = uniform_dist2(e1);
            randomBorderID = this->treeNodes[treeIdx].borderVertexIDs[randomNum];
            while (chosenBorders.find(randomBorderID) != chosenBorders.end()) {
                // Find a border we haven't chosen already 
                randomNum = uniform_dist2(e1);
                randomBorderID = this->treeNodes[treeIdx].borderVertexIDs[randomNum];
            }
            chosenBorders.insert(randomBorderID);
        }
    }
    std::vector<VertexID> chosenBordersVec(chosenBorders.begin(),chosenBorders.end());

    // Choose the subgraph vertices with minimum maximum distance to any of the chosen border vertices
    // Note: The positions are in the node (i.e. reordered vertex ID minus startVertexID), smallest maximum distance first
    std::vector<int> minmaxVertexPositions;
    this->findMinmaxVerticesFromSourceSet(originalGraph,treeIdx,chosenBordersVec,this->numNodeLandmarks,minmaxVertexPositions);

    // Check capacity even when landmarks already exist; the final count check alone cannot cover that case
    RN_DEBUG_ASSERT(this->numNodeLandmarks <= this->treeNodes[treeIdx].getNumVertices() && "Not enough subgraph vertices to complete landmarks set");
    for (std::size_t i = 0; i < minmaxVertexPositions.size() && this->treeNodes[treeIdx].landmarks.size() < this->numNodeLandmarks; ++i) {
        this->treeNodes[treeIdx].landmarks.push_back(this->reorderedVertexIDtoOriginalVertexID[this->treeNodes[treeIdx].startVertexID+minmaxVertexPositions[i]]);
    }

    if (this->treeNodes[treeIdx].landmarks.size() != this->numNodeLandmarks) {
        std::cout << "ERROR: Did not find enough landmarks" << std::endl;
        std::cout << "\t\"this->treeNodes[treeIdx].landmarks.size()\" = " << this->treeNodes[treeIdx].landmarks.size() << std::endl;
        std::cout << "\t\"minmaxVertexPositions.size() = " << minmaxVertexPositions.size() << std::endl;
        std::cout << "\t\"chosenBorders.size() = " << chosenBorders.size() << std::endl;
        std::cout << "\t\"chosenBordersVec.size() = " << chosenBordersVec.size() << std::endl;
        std::exit(1);
    }
}

// Finds the numMinmaxVertices vertices of the tree node's subgraph whose maximum shortest path distance from the vertices in sourceSet
// is smallest, and returns their positions in the node (i.e. reordered vertex ID minus startVertexID), smallest maximum distance
// first and ties broken by position. One Dijkstra search runs from each source, and the searches advance together: the vertex
// settled next is always the one with the smallest distance in any search. A subgraph vertex is complete once every search has
// settled it, which happens at its maximum distance from the sources, so subgraph vertices complete in order of that distance.
// The searches stop once numMinmaxVertices vertices are complete and no other vertex can complete at the same distance, so each
// search only settles the vertices within that distance of its source
void SULTree::findMinmaxVerticesFromSourceSet(Graph& graph, int treeIdx, std::vector<VertexID>& sourceSet, std::size_t numMinmaxVertices,
                                              std::vector<int>& minmaxVertexPositions)
{
    std::size_t numSources = sourceSet.size();
    std::vector<RadixHeap<VertexID>> pqueues(numSources);
    std::vector<std::vector<bool>> isVertexSettled(numSources,std::vector<bool>(graph.getNumVertices(),false));
    std::vector<unsigned int> numSearchesSettled(this->treeNodes[treeIdx].getNumVertices(),0); // For each subgraph vertex
    std::vector<std::pair<EdgeWeight,int>> completeVertices; // Maximum distance and position of each complete subgraph vertex
    std::size_t numSubgraphVertices = this->treeNodes[treeIdx].getNumVertices();

    EdgeWeight minDist, lastMaxDist = 0;
    VertexID minDistVertexID, adjVertex, reorderedID;
    int adjListStart, nextAdjListStart;

    // Initialize each search's priority queue with its source vertex
    for (std::size_t i = 0; i < numSources; ++i) {
        pqueues[i].push(0,sourceSet[i]);
    }

    while (completeVertices.size() < numSubgraphVertices) {
        // Choose the search whose next vertex is closest to its source
        std::size_t searchIdx = numSources;
        minDist = 0;
        for (std::size_t i = 0; i < numSources; ++i) {
            if (!pqueues[i].empty() && (searchIdx == numSources || pqueues[i].top_key() < minDist)) {
                searchIdx = i;
                minDist = pqueues[i].top_key();
            }
        }
        if (searchIdx == numSources || (completeVertices.size() >= numMinmaxVertices && minDist > lastMaxDist)) {
            // Every search is exhausted, or every vertex with the same maximum distance as the last one needed is complete
            break;
        }

        // Extract and remove vertex with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDistVertexID = pqueues[searchIdx].top_value();
        pqueues[searchIdx].pop();
        if (!isVertexSettled[searchIdx][minDistVertexID]) {
            isVertexSettled[searchIdx][minDistVertexID] = true;

            reorderedID = this->originalVertexIDtoReorderedVertexID[minDistVertexID];
            if (reorderedID >= this->treeNodes[treeIdx].startVertexID && reorderedID <= this->treeNodes[treeIdx].endVertexID) {
                int vertexNodePos = reorderedID - this->treeNodes[treeIdx].startVertexID;
                ++numSearchesSettled[vertexNodePos];
                if (numSearchesSettled[vertexNodePos] == numSources) {
                    // Every search has now settled this subgraph vertex, the last of them at its maximum distance
                    completeVertices.push_back(std::make_pair(minDist,vertexNodePos));
                    if (completeVertices.size() == numMinmaxVertices) {
                        lastMaxDist = minDist;
                    }
                }
            }

            // Inspect each neighbour and update this search's pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);

            for (int j = adjListStart; j < nextAdjListStart; ++j) {
                adjVertex = graph.edges[j].first;
                // Only update those we haven't already settled
                if (!isVertexSettled[searchIdx][adjVertex]) {
                    pqueues[searchIdx].push(minDist+graph.edges[j].second,adjVertex);
                }
            }
        }
    }

    // Order by maximum distance, then by position (vertices complete in order of maximum distance, but not of position)
    std::sort(completeVertices.begin(),completeVertices.end());
    minmaxVertexPositions.clear();
    for (std::size_t i = 0; i < completeVertices.size() && i < numMinmaxVertices; ++i) {
        minmaxVertexPositions.push_back(completeVertices[i].second);
    }
}

// Find shortest path from source to subgraph vertices. Assumptions
// - treeIdx is not be 0 (i.e. root)
int SULTree::findSSMTDistances(Graph& graph, int treeIdx, VertexID source, std::vector<EdgeWeight>& results)
{
    RN_DEBUG_ASSERT(treeIdx != 0 && "Cannot compute SSMT for root node as targets are the entire graph");
    int numSettledVertices = 0;
    
    RadixHeap<VertexID> pqueue;
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);
    
    EdgeWeight minDist;
    VertexID minDistVertexID, adjVertex, reorderedID;
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

            reorderedID = this->originalVertexIDtoReorderedVertexID[minDistVertexID];
            if (reorderedID >= this->treeNodes[treeIdx].startVertexID && reorderedID <= this->treeNodes[treeIdx].endVertexID) {
                // We have found the SP distance to one of the target set
                // So we add it to result set and remove it from input set
                int vertexNodePos = reorderedID - this->treeNodes[treeIdx].startVertexID;
                results[vertexNodePos] = minDist; // This is the final SP distance from source
                ++targetsFound;
                if (targetsFound == this->treeNodes[treeIdx].getNumVertices()) {
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

    return numSettledVertices;
}

// Find shortest path from source to subgraph vertices. Assumptions
// - Root node distance lists have been computed and stored in the root's part of this->nodeVertexLmkDists
// - Conversely, treeIdx cannot be 0 (i.e. root)
// - Root landmark min/max distances to border have been computed
// - All borders are stored in the the corresponding SULTreeNode and the follow the
//   definition that any shortest path from outside must pass through one of them
int SULTree::findSSMTDistancesViaBorderSet(Graph& graph, int treeIdx, VertexID source, std::vector<EdgeWeight>& results)
{
    RN_DEBUG_ASSERT(treeIdx != 0 && "Cannot compute SSMT for root node as targets are the entire graph");
    int numSettledVertices = 0;
    
    RadixHeap<VertexDistancePair> pqueue;
    std::vector<bool> isVertexSettled(graph.getNumVertices(),false);
    
    EdgeWeight minDist, newKey, prevKey, sourceToAdjVertexDist, minAdjVertexToBorderDist, subgraphLB;
    VertexID minDistVertexID, adjVertex, reorderedID;
    int adjListStart, nextAdjListStart;
    const EdgeWeight* rootLmkDists = this->nodeVertexLmkDists.data()+this->nodeStartIndex[0]; // Root's lists (position is the vertex ID)

    // Initialize with priority queue with source vertex
    pqueue.push(0,VertexDistancePair(source,0));
    std::size_t targetsFound = 0;

    // Once every border vertex is settled, a shortest path to a remaining target cannot leave the node: it would have to come back
    // through a border vertex, which is already settled at a distance no larger than that of any unsettled vertex outside the node
    std::vector<bool> isBorderPos(this->treeNodes[treeIdx].getNumVertices(),false);
    std::size_t numBorders = 0, bordersSettled = 0;
    for (VertexID borderVertexID : this->treeNodes[treeIdx].borderVertexIDs) {
        VertexID borderPos = this->originalVertexIDtoReorderedVertexID[borderVertexID]-this->treeNodes[treeIdx].startVertexID;
        if (!isBorderPos[borderPos]) {
            isBorderPos[borderPos] = true;
            ++numBorders;
        }
    }

    while (!pqueue.empty()) {
        // Extract and remove vertex with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        VertexDistancePair minElement = pqueue.top_value();
        prevKey = pqueue.top_key();
        minDistVertexID = minElement.first;
        minDist = minElement.second;
        pqueue.pop();
        reorderedID = this->originalVertexIDtoReorderedVertexID[minDistVertexID];
        bool isInsideNode = reorderedID >= this->treeNodes[treeIdx].startVertexID && reorderedID <= this->treeNodes[treeIdx].endVertexID;
        if (!isVertexSettled[minDistVertexID] && (isInsideNode || bordersSettled < numBorders)) {
            isVertexSettled[minDistVertexID] = true;
            ++numSettledVertices;

            if (isInsideNode) {
                // We have found the SP distance to one of the target set
                // So we add it to result set and remove it from input set
                int vertexNodePos = reorderedID - this->treeNodes[treeIdx].startVertexID;
                results[vertexNodePos] = minDist; // This is the final SP distance from source
                ++targetsFound;
                if (targetsFound == this->treeNodes[treeIdx].getNumVertices()) {
                    break;
                }
                if (isBorderPos[vertexNodePos]) {
                    ++bordersSettled;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isVertexSettled[adjVertex]) {
                    sourceToAdjVertexDist = minDist+graph.edges[i].second;
                    newKey = sourceToAdjVertexDist;
                    // Note: If the preceeding vertex was outside the box but the adjacent vertex 
                    // is inside, then the path will no longer be consider outside
                    reorderedID = this->originalVertexIDtoReorderedVertexID[adjVertex];
                    if (reorderedID < this->treeNodes[treeIdx].startVertexID || reorderedID > this->treeNodes[treeIdx].endVertexID) {
                        if (bordersSettled == numBorders) {
                            continue; // No vertex outside the node is needed any more
                        }
                        minAdjVertexToBorderDist = 0; // Initialize to zero
                        // If the adjacent vertex is outside subgraph we compute the key as the distance to it 
                        // plus the minimum distance to any border vertex
                        for (std::size_t j = 0; j < this->treeNodes[0].landmarks.size(); ++j) {
                            long vertexDistIdx = adjVertex*this->treeNodes[0].landmarks.size()+j;
                            if (rootLmkDists[vertexDistIdx] > this->treeNodes[treeIdx].maxBorderDist[j]) {
                                subgraphLB = rootLmkDists[vertexDistIdx] - this->treeNodes[treeIdx].maxBorderDist[j];
                                if (subgraphLB > minAdjVertexToBorderDist) {
                                    minAdjVertexToBorderDist = subgraphLB;
                                }
                            }
                            if (rootLmkDists[vertexDistIdx] < this->treeNodes[treeIdx].minBorderDist[j]) {
                                subgraphLB = this->treeNodes[treeIdx].minBorderDist[j] - rootLmkDists[vertexDistIdx];
                                if (subgraphLB > minAdjVertexToBorderDist) {
                                    minAdjVertexToBorderDist = subgraphLB;
                                }
                            }
                        }
                        // Handle scenario where lower-bound might be worse than lower-bound
                        // to mindist. This is not possible because the lower-bound distance
                        // to the border set through that vertex cannot become smaller through
                        // an adjacent vertex (it can only increase if it is a bad path or
                        // stay the same if it is the shortest path)
                        newKey += minAdjVertexToBorderDist;
                        if (newKey < prevKey) {
                            newKey = prevKey;
                        }
                        pqueue.push(newKey,VertexDistancePair(adjVertex,sourceToAdjVertexDist));
                    } else {
                        // For adjacent vertices inside the subgrpah, the new key is the just the
                        // distance to the settled vertex plus the edge weight to the adjacent vertex
                        pqueue.push(newKey,VertexDistancePair(adjVertex,sourceToAdjVertexDist));
                    }
                }
            }
        }
    }

    if (targetsFound != this->treeNodes[treeIdx].getNumVertices()) {
        std::cerr << "INTERNAL ERROR: the search from landmark " << source << " reached " << targetsFound << " of the "
                  << this->treeNodes[treeIdx].getNumVertices() << " vertices of SUL-Tree node " << treeIdx << std::endl;
        std::exit(1);
    }
    return numSettledVertices;
}

// Function to use root landmark distances in the same was ALT
EdgeWeight SULTree::computeRootLowerBound(VertexID s, VertexID t)
{
    EdgeWeight bestLB = 0, currentLB;
    EdgeWeight sourceLmkDist, targetLmkDist;
    for (std::size_t i = 0; i < this->treeNodes[0].landmarks.size(); ++i) {
        long sourceVertexDistIdx = s*this->treeNodes[0].landmarks.size()+i;
        long targetVertexDistIdx = t*this->treeNodes[0].landmarks.size()+i;
        sourceLmkDist = this->nodeVertexLmkDists[this->nodeStartIndex[0]+sourceVertexDistIdx];
        targetLmkDist = this->nodeVertexLmkDists[this->nodeStartIndex[0]+targetVertexDistIdx];
        currentLB = (sourceLmkDist > targetLmkDist) ? sourceLmkDist - targetLmkDist : targetLmkDist - sourceLmkDist;
        if (currentLB > bestLB) {
            bestLB = currentLB;
        }
    }
    return bestLB;
}

/*
 * SULTreeNode
 */

SULTreeNode::SULTreeNode(int parentIdx, int numVertices): 
    parentIdx(parentIdx), numVertices(numVertices), isLeaf(false) {}

SULTreeNode::SULTreeNode() {}
    
bool SULTreeNode::isLeafNode()
{
    return this->isLeaf;
}

int SULTreeNode::getParentIdx()
{
    return this->parentIdx;
}

unsigned int SULTreeNode::getNumVertices()
{
    return this->numVertices;
}

double SULTreeNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(int)*2;
    memoryUsage += sizeof(bool)*1;
    memoryUsage += sizeof(VertexID)*2;
    memoryUsage += sizeof(this->childNodeIdxs);
    memoryUsage += sizeof(int)*this->childNodeIdxs.size();
    memoryUsage += sizeof(this->landmarks);
    memoryUsage += sizeof(int)*this->landmarks.size();
    return memoryUsage;
}

COLTREE_LANDMARK_TYPE SULTree::getLandmarkTypeFromString(std::string landmarkTypeString, bool& success)
{
    if (landmarkTypeString == "random") {
        success = true;
        return COLTREE_LANDMARK_TYPE::SUBGRAPH_RANDOM;
    } else if (landmarkTypeString == "slice_furthest") {
        success = true;
        return COLTREE_LANDMARK_TYPE::SUBGRAPH_SLICE_FURTHEST;
    } else if (landmarkTypeString == "furthest_border") {
        success = true;
        return COLTREE_LANDMARK_TYPE::SUBGRAPH_FURTHEST_BORDER;
    } else if (landmarkTypeString == "border_minmax") {
        success = true;
        return COLTREE_LANDMARK_TYPE::SUBGRAPH_BORDER_MINMAX;
    } else {
        // Return default method but indicate method was not found in success
        success = false; 
        return COLTREE_LANDMARK_TYPE::SUBGRAPH_SLICE_FURTHEST;
    }
}

/*
 * Compact Object-Landmark Tree (COLTree)
 */

COLTree::COLTree(std::string setType, double setDensity, int setVariable, int setSize): 
    objSetType(setType), objSetDensity(setDensity), objSetVariable(setVariable), objSetSize(setSize) {
}

COLTree::COLTree() {}

// Note: COLTree will prune children without objects, which means point containment cannot
// be performed since part of the parents range of vertex ID may not be covered
void COLTree::buildCOLTTopDown(SULTree& sltree, std::vector<VertexID>& objectIDs, unsigned int maxOLSize)
{
    this->maxOLSize = maxOLSize;

    std::stack<TreeNodePair> stack; // Pair of new COLTree treeNode index and corresponding SULTreeNode index

    // Convert objectIDs to reordered vertex IDs
    // Note: Will not be necessary if permanently reorder vertices
    std::vector<VertexPair> reObjectIDs; // first = reordered ID, second = original ID
    reObjectIDs.reserve(objectIDs.size());
    for (std::size_t i = 0; i < objectIDs.size(); ++i) {
        reObjectIDs.push_back(VertexPair(sltree.originalVertexIDtoReorderedVertexID[objectIDs[i]],objectIDs[i]));
    }
    std::sort(reObjectIDs.begin(),reObjectIDs.end()); // Sorted by first by default

    // Create root node add it to stack for top-down construction from the SULTree
    COLTreeNode root(sltree.treeNodes[0].getParentIdx(),objectIDs.size());
    root.landmarks = sltree.treeNodes[0].landmarks;
    root.lmkObjectLists.resize(root.landmarks.size());
    root.firstObjectIDIndex = 0;
    root.lastObjectIDIndex = reObjectIDs.size()-1;
    root.lmkOLMaxDist.resize(root.landmarks.size(),0);
    root.lmkOLMinDist.resize(root.landmarks.size(),std::numeric_limits<EdgeWeight>::max());
    root.globalLmkOLMaxDist.resize(root.landmarks.size(),0);
    root.globalLmkOLMinDist.resize(root.landmarks.size(),std::numeric_limits<EdgeWeight>::max());
    this->treeNodes.push_back(std::move(root));

    if (sltree.treeNodes[0].isLeafNode()) {
        std::cout << "WARNING: SULTree node root is also leaf, can't build meaningful COLTree based such a SULTree" << std::endl;
    }

    if (objectIDs.size() <= this->maxOLSize) {
        std::cout << "INFO: Entire COLTree has one node (not necessarily a bad thing, but lambda could be too large)" << std::endl;
        this->treeNodes[0].setObjectListStatus(true);
        long vertexDistIdx;
        for (std::size_t i = 0; i < this->treeNodes[0].landmarks.size(); ++i) {
            this->treeNodes[0].lmkObjectLists[i].reserve(objectIDs.size());
            for (std::size_t j = 0; j < objectIDs.size(); ++j) {
                // Note: newNodeIdx is guaranteed not to be root here so we don't need to handle it separately
                VertexID objectID = objectIDs[j];
                int vertexNodePosition = objectID;
                vertexDistIdx = sltree.nodeStartIndex[0]+vertexNodePosition*this->treeNodes[0].landmarks.size()+i;
                EdgeWeight landmarkDist = sltree.nodeVertexLmkDists[vertexDistIdx];
                this->treeNodes[0].lmkObjectLists[i].push_back(VertexEdgeWeightPair(objectID,landmarkDist));
            }

            // Sort object lists according to landmark distance
            std::sort(this->treeNodes[0].lmkObjectLists[i].begin(),this->treeNodes[0].lmkObjectLists[i].end(), [](const VertexEdgeWeightPair &left, const VertexEdgeWeightPair &right) {
                return left.second < right.second;
            });
        }
    } else {
        stack.push(TreeNodePair(0,0));
    }

    while (stack.size() > 0) {
        TreeNodePair element = stack.top();
        stack.pop();

        TreeNodeID coltNodeIdx = element.first;
        TreeNodeID slTreeNodeIdx = element.second;

        // Note: guaranteed that this node is a not a leaf because it would not have been added to stack
        RN_DEBUG_ASSERT(this->treeNodes[coltNodeIdx].getNumObjects() > this->maxOLSize && "A COLTreeNode that should have been a leaf was added to stack");

        // Sort the childNodexIdx by their startVertexID
        std::vector<IntVertexPair> sortedChildNodeIdxs;
        for (int childIdx: sltree.treeNodes[slTreeNodeIdx].childNodeIdxs) {
            sortedChildNodeIdxs.push_back(IntVertexPair(childIdx,sltree.treeNodes[childIdx].startVertexID));
        }
        std::sort(sortedChildNodeIdxs.begin(),sortedChildNodeIdxs.end(),[](auto &left, auto &right) {
            return left.second < right.second;
        });

        RN_DEBUG_ASSERT(sltree.treeNodes[slTreeNodeIdx].startVertexID == sltree.treeNodes[sortedChildNodeIdxs[0].first].startVertexID && "startVertexID of first sorted childIdx should be same as the startVertexID of the parent");
        RN_DEBUG_ASSERT(sltree.treeNodes[slTreeNodeIdx].endVertexID == sltree.treeNodes[sortedChildNodeIdxs[sortedChildNodeIdxs.size()-1].first].endVertexID && "endVertexID of last sorted childIdx should be same as the endVertexID of the parent");

        // Start with the index of the first object belonging to the parent, iterate over 
        // the reordered object ID list and increment object count for each child until
        // we exceed the child vertex ID range (then move on to next child)
        int childNodesAdded = 0;
        int reObjectIDsIdx = this->treeNodes[coltNodeIdx].firstObjectIDIndex;
        int lastParentReObjectIDsIdx = this->treeNodes[coltNodeIdx].lastObjectIDIndex;
        // Note: reObjectIDsIdx <= lastParentReObjectIDsIdx in the for loop condition to ensure
        // save any pre-allocation for children that won't be created anyway
        for (std::size_t i = 0; i < sortedChildNodeIdxs.size() && reObjectIDsIdx <= lastParentReObjectIDsIdx; ++i) {
            int sltreeChildIdx = sortedChildNodeIdxs[i].first;

            if (reObjectIDs[reObjectIDsIdx].first > sltree.treeNodes[sltreeChildIdx].endVertexID) {
                // This SULTreeNode does contain any objects, so we skip
                continue;
            }

            // We can also compute COLTreeNode min/max distances while looping over its objects
            std::vector<EdgeWeight> lmkOLMaxDist(sltree.treeNodes[sltreeChildIdx].landmarks.size(),0);
            std::vector<EdgeWeight> lmkOLMinDist(sltree.treeNodes[sltreeChildIdx].landmarks.size(),std::numeric_limits<EdgeWeight>::max());
            std::vector<EdgeWeight> globalLmkOLMaxDist(sltree.treeNodes[0].landmarks.size(),0);
            std::vector<EdgeWeight> globalLmkOLMinDist(sltree.treeNodes[0].landmarks.size(),std::numeric_limits<EdgeWeight>::max());
            // Note: The above may not be used if this child contains no objects but the 
            // size are small as the number landmarks is a small constant

            unsigned int childObjects = 0;
            int firstChildObjectIDIndex = reObjectIDsIdx;
            // reObjectIDsIdx points to the first object of the parent
            // or one element after the last object that belong to previous child
            for (; reObjectIDsIdx <= lastParentReObjectIDsIdx; ++reObjectIDsIdx) {
                VertexID reObjectID = reObjectIDs[reObjectIDsIdx].first;
                // The object is within range of this child vertex ID range: the first one by the check at the start of this
                // child, and each later one because the check below did not end the child
                RN_DEBUG_ASSERT(reObjectID <= sltree.treeNodes[sltreeChildIdx].endVertexID && "The object belongs to this child");
                ++childObjects;

                // Also keep track of needed landmark distances
                int vertexNodePosition = reObjectID - sltree.treeNodes[sltreeChildIdx].startVertexID;
                int numLandmarks = sltree.treeNodes[sltreeChildIdx].landmarks.size();
                long vertexDistIdx;
                for (std::size_t j = 0; j < sltree.treeNodes[sltreeChildIdx].landmarks.size(); ++j) {
                    vertexDistIdx = sltree.nodeStartIndex[sltreeChildIdx]+vertexNodePosition*numLandmarks+j;
                    EdgeWeight landmarkDist = sltree.nodeVertexLmkDists[vertexDistIdx];
                    if (landmarkDist > lmkOLMaxDist[j]) {
                        lmkOLMaxDist[j] = landmarkDist;
                    }
                    if (landmarkDist < lmkOLMinDist[j]) {
                        // Existing min has already been initiazed to max possible values above
                        lmkOLMinDist[j] = landmarkDist;
                    }
                }
                // For global landmark we need to retrieve the original vertex ID (second)
                vertexNodePosition = reObjectIDs[reObjectIDsIdx].second;
                numLandmarks = sltree.treeNodes[0].landmarks.size();
                for (std::size_t j = 0; j < sltree.treeNodes[0].landmarks.size(); ++j) {
                    vertexDistIdx = sltree.nodeStartIndex[0]+vertexNodePosition*numLandmarks+j;
                    EdgeWeight landmarkDist = sltree.nodeVertexLmkDists[vertexDistIdx];
                    if (landmarkDist > globalLmkOLMaxDist[j]) {
                        globalLmkOLMaxDist[j] = landmarkDist;
                    }
                    if (landmarkDist < globalLmkOLMinDist[j]) {
                        // Existing min has already been initiazed to max possible values above
                        globalLmkOLMinDist[j] = landmarkDist;
                    }
                }
                // Note: We do not check the startVertexID because both reorderedObjectIDs
                // and sortedChildNodeIdxs are sorted on (reordered) vertex ID
                
                // If we've reached last element for the parents object IDs or if they next object ID is greater than this
                // child's range, then we have found all objects for this child, so we create the COLTreeNode and break out of 
                // the inner loop (iterating over object IDs) and continue the outer loop for the next child (if any)
                if (reObjectIDsIdx == lastParentReObjectIDsIdx ||
                    reObjectIDs[reObjectIDsIdx+1].first > sltree.treeNodes[sltreeChildIdx].endVertexID) {
                    // This object belongs in the next child, so we are done with this child node
                    ++childNodesAdded;

                    // Add index of new node to parent COLTreeNode child node list
                    int newNodeIdx = this->treeNodes.size();
                    this->treeNodes[coltNodeIdx].addChildNode(newNodeIdx);

                    COLTreeNode newNode(coltNodeIdx,childObjects);
                    newNode.landmarks = sltree.treeNodes[sltreeChildIdx].landmarks;
                    newNode.lmkObjectLists.resize(newNode.landmarks.size());
                    newNode.firstObjectIDIndex = firstChildObjectIDIndex;
                    newNode.lastObjectIDIndex = reObjectIDsIdx;
                    newNode.lmkOLMaxDist.swap(lmkOLMaxDist);
                    newNode.lmkOLMinDist.swap(lmkOLMinDist);
                    newNode.globalLmkOLMaxDist.swap(globalLmkOLMaxDist);
                    newNode.globalLmkOLMinDist.swap(globalLmkOLMinDist);
                    // Note: We don't copy childNodeIdxs because that is filled
                    // based on actual object occurences (i.e., children with no
                    // objects contained in its subgraph are not added in the COLTree)
                    this->treeNodes.push_back(std::move(newNode));

                    // We need to handle case where this node cannot be further processed
                    // because the correspond SULTreeNode is a leaf. So we must create 
                    // leaf COLTreeNode even if it has more than lambda objects
                    if (sltree.treeNodes[sltreeChildIdx].isLeafNode() && childObjects > maxOLSize) {
                        std::cout << "WARNING: SULTree node " << sltreeChildIdx << " is a leaf with " << childObjects
                            << " objects, which is more than lambda=" << maxOLSize << ". Consider decreasing alpha for SULTree?" << std::endl;
                    }

                    // We only add to stack if it is non-leaf, if it is a leaf we create ODL and set other variables
                    if (childObjects <= this->maxOLSize || sltree.treeNodes[sltreeChildIdx].isLeafNode()) {
                        this->treeNodes[newNodeIdx].setObjectListStatus(true);
                        // We use the first and last index position in reorderedObjectIds for objects
                        // assocate with this COLTreeNode to quickly create ODL
                        int numLandmarks = sltree.treeNodes[sltreeChildIdx].landmarks.size();
                        long vertexDistIdx;
                        for (std::size_t j = 0; j < this->treeNodes[newNodeIdx].landmarks.size(); ++j) {
                            this->treeNodes[newNodeIdx].lmkObjectLists[j].reserve(childObjects);
                            for (int k = this->treeNodes[newNodeIdx].firstObjectIDIndex; k <= this->treeNodes[newNodeIdx].lastObjectIDIndex; ++k) {
                                // Note: newNodeIdx is guaranteed not to be root here so we don't need to handle it separately
                                VertexID odlReObjectID = reObjectIDs[k].first;
                                VertexID odlOriginalObjectID = reObjectIDs[k].second;
                                int vertexNodePosition = odlReObjectID - sltree.treeNodes[sltreeChildIdx].startVertexID;
                                vertexDistIdx = sltree.nodeStartIndex[sltreeChildIdx]+vertexNodePosition*numLandmarks+j;
                                EdgeWeight landmarkDist = sltree.nodeVertexLmkDists[vertexDistIdx];
                                this->treeNodes[newNodeIdx].lmkObjectLists[j].push_back(VertexEdgeWeightPair(odlOriginalObjectID,landmarkDist));
                            }

                            // Sort object lists according to landmark distance
                            std::sort(this->treeNodes[newNodeIdx].lmkObjectLists[j].begin(),this->treeNodes[newNodeIdx].lmkObjectLists[j].end(), [](const VertexEdgeWeightPair &left, const VertexEdgeWeightPair &right) {
                                return left.second < right.second;
                            });
                        }

                        // Note: We iterate over objects again to create ODL instead of using existing iteration
                        // because otherwise we create many unnecessary extra ODLs for non-leaf nodes
                    } else {
                        // Insert into stack for further top-down construction
                        stack.push(TreeNodePair(newNodeIdx,sltreeChildIdx));
                    }

                    ++reObjectIDsIdx; // Increment so it points to the next child's first object or cause the outer end because reObjectIDsIdx > lastParentReObjectIDsIdx
                    break;
                }
            }
        }

        if (childNodesAdded == 1) {
            std::cout << "WARNING: COLTreeNode " << coltNodeIdx << " (SULTreeNode " << slTreeNodeIdx << ", "
                      << this->treeNodes[coltNodeIdx].getNumObjects() << " objects) spawned only one child" << std::endl;
            // If only one child is added (e.g., because all objects were clustered)
            // then we could not create that child (results waste tree search) and
            // isntead swap it with parent (keeping the child's more local landmarks)
            // Note: It also possible that this could happen recursively with a chain
            // single parent-child nodes but even one is very unlikely
        }
    }

    // Create struct for each node to be used during querying. This
    // saves the initialisation cost of the vector, and offers faster
    // cache peformance. Querying uses this struct in an on-demand
    // way, so it is simply reset when the node is first encountered
    this->queryData.resize(this->treeNodes.size());
}

std::string COLTree::getObjSetType()
{
    return this->objSetType;
}

double COLTree::getObjSetDensity()
{
    return this->objSetDensity;
}

int COLTree::getObjSetSize()
{
    return this->objSetSize;
}

int COLTree::getObjSetVariable()
{
    return this->objSetVariable;
}

double COLTree::computeIndexSize()
{
    double memoryUsage = 0;
    for (std::size_t i = 0; i < this->treeNodes.size(); ++i) {
        memoryUsage += this->treeNodes[i].computeIndexSizeBytes();
    }
    // Note: We exclude queryData because it just contains empty structs
    // to be fill during query (could be easily created during loading of the index)
    return memoryUsage/(1024*1024);
}

int COLTree::getMaxOLSize()
{
    return this->maxOLSize;
}

int COLTree::getNumNonEmptyNodes()
{
    int numNodes = 0;
    for (std::size_t i = 0; i < this->treeNodes.size(); ++i) {
        if (this->treeNodes[i].getNumObjects() > 0){
            ++numNodes;
        }
    }
    return numNodes;
}

// Get the distance from the root (root depth is 0)
int COLTree::getNodeDepth(int treeNodeIdx)
{
    int currentNodeIdx = treeNodeIdx, depth = 0;
    while (currentNodeIdx != 0) {
        currentNodeIdx = this->treeNodes[currentNodeIdx].getParentIdx();
        ++depth;
    }
    return depth;
}

// Note: Nodes in each level should be group by their parent
int COLTree::getTreeDepth()
{
    int numLevels = 0;
    std::vector<int> currentLevel, nextLevel;
    nextLevel.push_back(0);
    while(nextLevel.size() != 0){
        ++numLevels;
        currentLevel.swap(nextLevel);
        nextLevel.clear();
        // For each node group in current level
        for (std::size_t i = 0; i < currentLevel.size(); ++i) {
            int nodeIdx = currentLevel[i];
            for (std::size_t k = 0; k < this->treeNodes[nodeIdx].childNodeIdxs.size(); ++k){
                int childIdx = this->treeNodes[nodeIdx].childNodeIdxs[k];
                if (this->treeNodes[childIdx].getNumObjects() > 0) {
                    nextLevel.push_back(childIdx);
                }
            }
        }
    }
    return numLevels;
}

EdgeWeight COLTree::computeGlobalLowerbound(int coltNodeIdx, VertexID queryVertex, SULTree& sltree)
{
    EdgeWeight lb = 0, queryLmkDist, newLB;
    long vertexDistIdx;
    // Not defined for the root node, whose global landmark distances are not computed (if the COL-Tree only has a
    // root node, bounds are not needed anyway); every caller passes a child node
    RN_DEBUG_ASSERT(coltNodeIdx != 0 && "The global landmark bounds are defined for child nodes only");
    for (std::size_t i = 0; i < this->treeNodes[0].landmarks.size(); ++i) {
        vertexDistIdx = sltree.nodeStartIndex[0]+queryVertex*this->treeNodes[0].landmarks.size()+i;
        queryLmkDist = sltree.nodeVertexLmkDists[vertexDistIdx];
        if (queryLmkDist > this->treeNodes[coltNodeIdx].globalLmkOLMaxDist[i]) {
            newLB = queryLmkDist - this->treeNodes[coltNodeIdx].globalLmkOLMaxDist[i];
            if (newLB > lb) {
                lb = newLB;
            }
        }
        // Might improve lower-bound when closest object is further 
        // from the landmark than query vertex
        if (queryLmkDist < this->treeNodes[coltNodeIdx].globalLmkOLMinDist[i]) {
            newLB = this->treeNodes[coltNodeIdx].globalLmkOLMinDist[i] - queryLmkDist;
            if (newLB > lb) {
                lb = newLB;
            }
        }
    }
    return lb;
}

EdgeWeight COLTree::computeGlobalUpperbound(int coltNodeIdx, VertexID queryVertex, SULTree& sltree)
{
    EdgeWeight ub = std::numeric_limits<EdgeWeight>::max(), queryLmkDist, newUB;
    long vertexDistIdx;
    // Not defined for the root node, whose global landmark distances are not computed (if the COL-Tree only has a
    // root node, bounds are not needed anyway); every caller passes a child node
    RN_DEBUG_ASSERT(coltNodeIdx != 0 && "The global landmark bounds are defined for child nodes only");
    for (std::size_t i = 0; i < this->treeNodes[0].landmarks.size(); ++i) {
        vertexDistIdx = sltree.nodeStartIndex[0]+queryVertex*this->treeNodes[0].landmarks.size()+i;
        queryLmkDist = sltree.nodeVertexLmkDists[vertexDistIdx];
        newUB = queryLmkDist + this->treeNodes[coltNodeIdx].globalLmkOLMaxDist[i];
        if (newUB < ub) {
            ub = newUB;
        }
    }
    return ub;
}

void COLTree::computeGlobalLowerAndUpperbound(int coltNodeIdx, VertexID queryVertex, SULTree& sltree, EdgeWeight& lb, EdgeWeight& ub)
{
    lb = 0, ub = std::numeric_limits<EdgeWeight>::max();
    EdgeWeight queryLmkDist, newLB, newUB;
    long vertexDistIdx;
    // Not defined for the root node, whose global landmark distances are not computed (if the COL-Tree only has a
    // root node, bounds are not needed anyway); every caller passes a child node
    RN_DEBUG_ASSERT(coltNodeIdx != 0 && "The global landmark bounds are defined for child nodes only");
    for (std::size_t i = 0; i < this->treeNodes[0].landmarks.size(); ++i) {
        vertexDistIdx = sltree.nodeStartIndex[0]+queryVertex*this->treeNodes[0].landmarks.size()+i;
        queryLmkDist = sltree.nodeVertexLmkDists[vertexDistIdx];
        if (queryLmkDist > this->treeNodes[coltNodeIdx].globalLmkOLMaxDist[i]) {
            newLB = queryLmkDist - this->treeNodes[coltNodeIdx].globalLmkOLMaxDist[i];
            if (newLB > lb) {
                lb = newLB;
            }
        }
        // Might improve lower-bound when closest object is further 
        // from the landmark than query vertex
        if (queryLmkDist < this->treeNodes[coltNodeIdx].globalLmkOLMinDist[i]) {
            newLB = this->treeNodes[coltNodeIdx].globalLmkOLMinDist[i] - queryLmkDist;
            if (newLB > lb) {
                lb = newLB;
            }
        }
        newUB = queryLmkDist + this->treeNodes[coltNodeIdx].globalLmkOLMaxDist[i];
        if (newUB < ub) {
            ub = newUB;
        }
    }
}

// Assumes queryDistances is already initialised using COLTree node landmarks
void COLTree::computeGlobalAggregateLowerbound(int coltNodeIdx, std::vector<VertexID>& queryVertices, std::vector<EdgeWeight>& queryDistances, SULTree& sltree, bool& used)
{
    EdgeWeight queryLmkDist, newLB;
    long vertexDistIdx;
    // For each query vertex, find the tightest (maximum) lower-bound distance to an object in this node over all landmarks
    for (std::size_t j = 0; j < queryVertices.size(); ++j) {
        for (std::size_t i = 0; i < this->treeNodes[0].landmarks.size(); ++i) {
            // Assumes queryDistances is already initialised using COLTree node landmarks
            vertexDistIdx = sltree.nodeStartIndex[0]+queryVertices[j]*this->treeNodes[0].landmarks.size()+i;
            queryLmkDist = sltree.nodeVertexLmkDists[vertexDistIdx];
            if (queryLmkDist > this->treeNodes[coltNodeIdx].globalLmkOLMaxDist[i]) {
                newLB = queryLmkDist - this->treeNodes[coltNodeIdx].globalLmkOLMaxDist[i];
                if (newLB > queryDistances[j]) {
                    queryDistances[j] = newLB;
#if defined(COLLECT_STATISTICS)
                    used = true;
#endif
                }
            }
            // Might improve lower-bound when closest object is further 
            // from the landmark than query vertex
            if (queryLmkDist < this->treeNodes[coltNodeIdx].globalLmkOLMinDist[i]) {
                newLB = this->treeNodes[coltNodeIdx].globalLmkOLMinDist[i] - queryLmkDist;
                if (newLB > queryDistances[j]) {
                    queryDistances[j] = newLB;
#if defined(COLLECT_STATISTICS)
                    used = true;
#endif
                }
            }
        }
    }
}

/*
 * Object List Tree Node
 */

COLTreeNode::COLTreeNode(int parentIdx, int numObjects): 
    parentIdx(parentIdx), hasOL(false), numObjects(numObjects) {}

COLTreeNode::COLTreeNode() {}
    
bool COLTreeNode::hasObjectList()
{
    return this->hasOL;
}

void COLTreeNode::setObjectListStatus(bool status)
{
    this->hasOL = status;
}

int COLTreeNode::getParentIdx()
{
    return this->parentIdx;
}

unsigned int COLTreeNode::getNumObjects()
{
    return this->numObjects;
}

void COLTreeNode::addChildNode(int childIdx)
{
    this->childNodeIdxs.push_back(childIdx);
}

double COLTreeNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(int)*2;
    memoryUsage += sizeof(bool)*1;
    memoryUsage += sizeof(this->childNodeIdxs);
    memoryUsage += sizeof(int)*this->childNodeIdxs.size();
    memoryUsage += sizeof(this->landmarks);
    memoryUsage += sizeof(int)*this->landmarks.size();
    memoryUsage += sizeof(this->lmkObjectLists);
    for (std::size_t i = 0; i < this->lmkObjectLists.size(); ++i) {
        memoryUsage += sizeof(this->lmkObjectLists[i]);
        memoryUsage += sizeof(VertexEdgeWeightPair)*this->lmkObjectLists[i].size();
    }
    memoryUsage += sizeof(this->lmkOLMaxDist);
    memoryUsage += sizeof(int)*this->lmkOLMaxDist.size();
    memoryUsage += sizeof(this->lmkOLMinDist);
    memoryUsage += sizeof(int)*this->lmkOLMinDist.size();
    memoryUsage += sizeof(this->globalLmkOLMaxDist);
    memoryUsage += sizeof(int)*this->globalLmkOLMaxDist.size();
    memoryUsage += sizeof(this->globalLmkOLMinDist);
    memoryUsage += sizeof(int)*this->globalLmkOLMinDist.size();
    return memoryUsage;
}

/* Shared Node Bound Functions */

int COLTreeNode::getIndexClosestToDistance(EdgeWeight matchDist, int listIdx)
{
    int left = 0, right = lmkObjectLists[listIdx].size()-1, mid;
    while(left <= right) {
        mid = (right+left)/2;
        if (lmkObjectLists[listIdx][mid].second == matchDist) {
            return mid;
        } else if (lmkObjectLists[listIdx][mid].second < matchDist) {
            left = mid+1;
        } else /* if (lmkObjectLists[listIdx][mid].second > matchDist)*/ {
            right = mid-1;
        }
    }
    
    // At this point left is greater than right by 1 
    // One of left or right is the closest value
    EdgeWeight leftDiff, rightDiff;
    if (static_cast<std::size_t>(left) < lmkObjectLists[listIdx].size() && right >= 0) {
        leftDiff = (lmkObjectLists[listIdx][left].second > matchDist) ? lmkObjectLists[listIdx][left].second - matchDist : matchDist - lmkObjectLists[listIdx][left].second;
        rightDiff = (lmkObjectLists[listIdx][right].second > matchDist) ? lmkObjectLists[listIdx][right].second - matchDist : matchDist - lmkObjectLists[listIdx][right].second;
        if (leftDiff < rightDiff) {
            return left;
        } else {
            return right;
        }
    } else if (static_cast<std::size_t>(left) >= lmkObjectLists[listIdx].size()) {
        return right;
    } else {
        return left;
    }
}

// Returns the number of objects among the first numElements of the OL that are closer to the landmark than dist,
// i.e., the index of the first object at least dist from it. The search keeps the range [left,right) of objects
// that may be that object (unsigned indices, so halving needs no rounding correction for negative values)
std::size_t COLTreeNode::getNumObjectsBeforeDistance(EdgeWeight dist, int listIdx, std::size_t numElements)
{
    std::vector<VertexEdgeWeightPair>& lmkObjectList = this->lmkObjectLists[listIdx];
    std::size_t left = 0, right = numElements, mid;
    while (left < right) {
        mid = (left+right)/2;
        if (lmkObjectList[mid].second < dist) {
            left = mid+1;
        } else {
            right = mid;
        }
    }
    return left;
}

// Returns the number of objects among the first numElements of the OL that are no further from the landmark than dist,
// i.e., the index of the first object further than dist from it (searched as in getNumObjectsBeforeDistance)
std::size_t COLTreeNode::getNumObjectsUpToDistance(EdgeWeight dist, int listIdx, std::size_t numElements)
{
    std::vector<VertexEdgeWeightPair>& lmkObjectList = this->lmkObjectLists[listIdx];
    std::size_t left = 0, right = numElements, mid;
    while (left < right) {
        mid = (left+right)/2;
        if (lmkObjectList[mid].second <= dist) {
            left = mid+1;
        } else {
            right = mid;
        }
    }
    return left;
}

/* kNN Node Bound Functions */

EdgeWeight COLTreeNode::getBestNodeLowerbound(ALT& alt, VertexID queryVertex)
{
    EdgeWeight lb = 0, queryLmkLB, newLB;
    for (std::size_t i = 0; i < this->landmarks.size(); ++i) {
        queryLmkLB = alt.getLowerBound(queryVertex,landmarks[i]);
        if (queryLmkLB > this->lmkOLMaxDist[i]) {
            newLB = queryLmkLB - this->lmkOLMaxDist[i];
            if (newLB > lb) {
                lb = newLB;
            }
        }
    }
    // Note: We don't consider second case using upper-bound on query-to-landmark because
    // this requires the query to be closer to the landmark than any object within the 
    // node containing that landmark. This only happens at very high-levels of the 
    // COL-Tree or when the query is in that node's subgraph (faster/simpler just to skip)
    return lb;
}

EdgeWeight COLTreeNode::getBestNodeLowerboundAllCases(ALT& alt, VertexID queryVertex)
{
    // Lower-bound seems to only rarely be useful in reducing lower-bound
    EdgeWeight lb = 0, queryLmkLB, queryLmkUB, newLB;
    for (std::size_t i = 0; i < this->landmarks.size(); ++i) {
        alt.getLowerAndUpperBound(queryVertex,landmarks[i],queryLmkLB,queryLmkUB);
        if (queryLmkLB > this->lmkOLMaxDist[i]) {
            newLB = queryLmkLB - this->lmkOLMaxDist[i];
            if (newLB > lb) {
                lb = newLB;
            }
        }
        // Might improve lower-bound when closest object is further 
        // from the landmark than query vertex
        if (queryLmkUB < this->lmkOLMinDist[i]) {
            newLB = this->lmkOLMinDist[i] - queryLmkUB;
            if (newLB > lb) {
                lb = newLB;
            }
        }
    }
    return lb;
}

EdgeWeight COLTreeNode::getBestNodeLowerboundAndFurthestLmk(ALT& alt, VertexID queryVertex, int& furthestLandmarkIdx)
{
    EdgeWeight lb = 0, queryLmkLB, newLB, furthestLandmarkDist = 0;
    furthestLandmarkIdx = 0;
    for (std::size_t i = 0; i < this->landmarks.size(); ++i) {
        queryLmkLB = alt.getLowerBound(queryVertex,landmarks[i]);
        if (queryLmkLB > this->lmkOLMaxDist[i]) {
            newLB = queryLmkLB - this->lmkOLMaxDist[i];
            if (newLB > lb) {
                lb = newLB;
            }
        }
        if (queryLmkLB > furthestLandmarkDist) {
            furthestLandmarkIdx = i;
            furthestLandmarkDist = queryLmkLB;
        }
    }
    return lb;
}

// Returns the furthest landmark from the query vertex by their lower-bound distance
int COLTreeNode::getFurthestLmkIdxByLB(ALT& alt, VertexID queryVertex)
{
    int lmkIdx = 0;
    EdgeWeight lb = alt.getLowerBound(queryVertex,landmarks[0]), newLB;
    for (std::size_t i = 1; i < this->landmarks.size(); ++i) {
        newLB = alt.getLowerBound(queryVertex,landmarks[i]);
        if (newLB > lb) {
            lb = newLB;
            lmkIdx = i;
        }
    }
    return lmkIdx;
}

/* AkNN Node Bound Functions */

void COLTreeNode::getBestNodeLowerbounds(ALT& alt, std::vector<VertexID>& queryVertices, std::vector<EdgeWeight>& queryDistances)
{
    EdgeWeight queryLmkLB, newLB;
    // For each query vertex, find the tightest (maximum) lower-bound distance to an object in this node over all landmarks
    for (std::size_t i = 0; i < this->landmarks.size(); ++i) {
        for (std::size_t j = 0; j < queryVertices.size(); ++j) {
            if (i == 0) {
                // If this is the first iteration of outer loop (i.e. first landmark),
                // initialise the best maximum distance for this query vertex
                queryDistances[j] = 0;
            }
            queryLmkLB = alt.getLowerBound(queryVertices[j],landmarks[i]);
            if (queryLmkLB > this->lmkOLMaxDist[i]) {
                newLB = queryLmkLB - this->lmkOLMaxDist[i];
                if (newLB > queryDistances[j]) {
                    queryDistances[j] = newLB;
                }
            }
        }
    }
}

// Does not assume furthestSumLandmarkIdx is initialized
void COLTreeNode::getBestNodeLowerboundsAndFurthestSumLmk(ALT& alt, std::vector<VertexID>& queryVertices, std::vector<EdgeWeight>& queryDistances,
                                                                   int& furthestSumLandmarkIdx)
{
    furthestSumLandmarkIdx = 0; // Initialize (rare case but needed if all lower-bounds are zero)
    EdgeWeight queryLmkLB, newLB, furthestDistSum = 0;
    EdgeWeight newDistSum; // Not using long because it is highly unliky to overflow even for USA data
    // For each query vertex, find the tightest (maximum) lower-bound distance to an object in this node over all landmarks
    for (std::size_t i = 0; i < this->landmarks.size(); ++i) {
        newDistSum = 0;
        for (std::size_t j = 0; j < queryVertices.size(); ++j) {
            if (i == 0) {
                // If this is the first iteration of outer loop (i.e. first landmark),
                // initialise the best maximum distance for this query vertex
                queryDistances[j] = 0;
            }
            queryLmkLB = alt.getLowerBound(queryVertices[j],landmarks[i]);
            if (queryLmkLB > this->lmkOLMaxDist[i]) {
                newLB = queryLmkLB - this->lmkOLMaxDist[i];
                if (newLB > queryDistances[j]) {
                    queryDistances[j] = newLB;
                }
            }
            newDistSum += queryLmkLB;
        }
        if (newDistSum > furthestDistSum) {
            furthestSumLandmarkIdx = i;
            furthestDistSum = newDistSum;
        }
    }
}

// Returns the furthest landmark from the query vertices by their lower-bound distance
int COLTreeNode::getFurthestLmkIdxByLB(ALT& alt, std::vector<VertexID>& queryVertices)
{
    int lmkIdx = 0;
    EdgeWeight lbSum = 0, newLBSum;
    for (std::size_t i = 0; i < this->landmarks.size(); ++i) {
        newLBSum = 0;
        for (std::size_t j = 0; j < queryVertices.size(); ++j) {
            newLBSum += alt.getLowerBound(queryVertices[j],landmarks[i]);
        }
        if (newLBSum > lbSum) {
            lbSum = newLBSum;
            lmkIdx = i;
        }
    }
    return lmkIdx;
}

/* FkN Node Bound Functions */

EdgeWeight COLTreeNode::getBestNodeUpperbound(ALT& alt, VertexID queryVertex)
{
    EdgeWeight ub = std::numeric_limits<EdgeWeight>::max(), queryLmkUB, newUB;
    for (std::size_t i = 0; i < this->landmarks.size(); ++i) {
        queryLmkUB = alt.getUpperBound(queryVertex,landmarks[i]);
        newUB = queryLmkUB + this->lmkOLMaxDist[i];
        if (newUB < ub) {
            ub = newUB;
        }
    }
    return ub;
}

// Note: Nearest landmark is decided using upper-bound distances.
// This is likely to work better for FkN queries which need accurate upper-bounds
// so landmarks should be "in-between" the query vertex and the COLTree leaf subgraph
EdgeWeight COLTreeNode::getBestNodeUpperboundAndNearestLmk(ALT& alt, VertexID queryVertex, int& nearestLandmarkIdx)
{
    EdgeWeight ub = std::numeric_limits<EdgeWeight>::max(), queryLmkUB, newUB, nearestLandmarkDist = std::numeric_limits<EdgeWeight>::max();
    nearestLandmarkIdx = 0;
    for (std::size_t i = 0; i < this->landmarks.size(); ++i) {
        queryLmkUB = alt.getUpperBound(queryVertex,landmarks[i]);
        newUB = queryLmkUB + this->lmkOLMaxDist[i];
        if (newUB < ub) {
            ub = newUB;
        }
        if (queryLmkUB < nearestLandmarkDist) {
            nearestLandmarkIdx = i;
            nearestLandmarkDist = queryLmkUB;
        }
    }
    return ub;
}

// Returns the nearest landmark from the query vertex by their upper-bound distance
int COLTreeNode::getNearestLmkIdxByUB(ALT& alt, VertexID queryVertex)
{
    int lmkIdx = 0;
    EdgeWeight ub = alt.getUpperBound(queryVertex,landmarks[0]), newUB;
    for (std::size_t i = 1; i < this->landmarks.size(); ++i) {
        newUB = alt.getUpperBound(queryVertex,landmarks[i]);
        if (newUB < ub) {
            ub = newUB;
            lmkIdx = i;
        }
    }
    return lmkIdx;
}

void COLTreeNode::getBestNodeLowerAndUpperbound(ALT& alt, VertexID queryVertex, EdgeWeight& lb, EdgeWeight& ub)
{
    // Lower-bound seems to only rarely be useful in reducing lower-bound
    EdgeWeight queryLmkLB, queryLmkUB, newLB, newUB;

    // Initialize using 0th landmark
    alt.getLowerAndUpperBound(queryVertex,landmarks[0],queryLmkLB,queryLmkUB);
    lb = (queryLmkLB > this->lmkOLMaxDist[0]) ? queryLmkLB - this->lmkOLMaxDist[0] : 0;
    if (queryLmkUB < this->lmkOLMinDist[0]) {
        newLB = this->lmkOLMinDist[0] - queryLmkUB;
        if (newLB > lb) {
            lb = newLB;
        }
    }
    ub = this->lmkOLMaxDist[0] + queryLmkUB;

    // Iterate over remaining landmarks and update if improvement
    for (std::size_t i = 1; i < this->landmarks.size(); ++i) {
        alt.getLowerAndUpperBound(queryVertex,landmarks[i],queryLmkLB,queryLmkUB);
        if (queryLmkLB > this->lmkOLMaxDist[i]) {
            newLB = queryLmkLB - this->lmkOLMaxDist[i];
            if (newLB > lb) {
                lb = newLB;
            }
        }
        // Might improve lower-bound when closest object is further 
        // from the landmark than query vertex
        if (queryLmkUB < this->lmkOLMinDist[i]) {
            newLB = this->lmkOLMinDist[i] - queryLmkUB;
            if (newLB > lb) {
                lb = newLB;
            }
        }
        newUB = this->lmkOLMaxDist[i] + queryLmkUB;
        if (newUB < ub) {
            ub = newUB;
        }
    }
}

// Version of function that also retrieves the lower- and upper-bound on the query-to-landmark distance of each landmark,
// from which getRangeSearchListIdx chooses the OL to search if the node is neither pruned nor returned whole
void COLTreeNode::getRangeSearchLeafBounds(ALT& alt, VertexID queryVertex, EdgeWeight& lb, EdgeWeight& ub, std::vector<EdgeWeight>& queryLmkLBs, std::vector<EdgeWeight>& queryLmkUBs)
{
    EdgeWeight queryLmkLB, queryLmkUB, newLB, newUB;
    queryLmkLBs.resize(this->landmarks.size());
    queryLmkUBs.resize(this->landmarks.size());

    // Initialize using 0th landmark
    alt.getLowerAndUpperBound(queryVertex,landmarks[0],queryLmkLB,queryLmkUB);
    queryLmkLBs[0] = queryLmkLB;
    queryLmkUBs[0] = queryLmkUB;
    lb = (queryLmkLB > this->lmkOLMaxDist[0]) ? queryLmkLB - this->lmkOLMaxDist[0] : 0;
    if (queryLmkUB < this->lmkOLMinDist[0]) {
        newLB = this->lmkOLMinDist[0] - queryLmkUB;
        if (newLB > lb) {
            lb = newLB;
        }
    }
    ub = this->lmkOLMaxDist[0] + queryLmkUB;

    // Iterate over remaining landmarks and update if improvement
    for (std::size_t i = 1; i < this->landmarks.size(); ++i) {
        alt.getLowerAndUpperBound(queryVertex,landmarks[i],queryLmkLB,queryLmkUB);
        queryLmkLBs[i] = queryLmkLB;
        queryLmkUBs[i] = queryLmkUB;
        if (queryLmkLB > this->lmkOLMaxDist[i]) {
            newLB = queryLmkLB - this->lmkOLMaxDist[i];
            if (newLB > lb) {
                lb = newLB;
            }
        }
        // Might improve lower-bound when closest object is further 
        // from the landmark than query vertex
        if (queryLmkUB < this->lmkOLMinDist[i]) {
            newLB = this->lmkOLMinDist[i] - queryLmkUB;
            if (newLB > lb) {
                lb = newLB;
            }
        }
        newUB = this->lmkOLMaxDist[i] + queryLmkUB;
        if (newUB < ub) {
            ub = newUB;
        }
    }
}

// Returns the index of the landmark whose OL has the fewest candidates for a range search, given the lower- and upper-bound
// on each query-to-landmark distance (from getRangeSearchLeafBounds). An object is a result if its distance to the landmark
// plus the upper-bound is within range, and it can be pruned if its distance to the landmark is smaller than the lower-bound
// by more than range or larger than the upper-bound by more than range. The first case needs an upper-bound within range
// and the second a lower-bound greater than range, so at most one of them applies to an OL, and as the OL is sorted by
// distance to the landmark, its candidates are the objects between that prefix and the suffix pruned by the last case
int COLTreeNode::getRangeSearchListIdx(EdgeWeight range, std::vector<EdgeWeight>& queryLmkLBs, std::vector<EdgeWeight>& queryLmkUBs)
{
    int listIdx = 0;
    std::size_t numCandidates, minNumCandidates = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = 0; i < this->landmarks.size(); ++i) {
        // Objects further from the landmark than maxObjectLMDist can be pruned
        EdgeWeight maxObjectLMDist = queryLmkUBs[i] + range;
        std::size_t candidatesEnd = getNumObjectsUpToDistance(maxObjectLMDist,i,this->lmkObjectLists[i].size());
        std::size_t candidatesStart = 0;
        if (queryLmkUBs[i] <= range) {
            // Landmark is within range, so objects up to maxResultLMDist from it are results by upper-bound
            EdgeWeight maxResultLMDist = range - queryLmkUBs[i];
            candidatesStart = getNumObjectsUpToDistance(maxResultLMDist,i,candidatesEnd);
        } else if (queryLmkLBs[i] > range) {
            // Landmark is outside range, so objects closer to it than minObjectLMDist can be pruned
            EdgeWeight minObjectLMDist = queryLmkLBs[i] - range;
            candidatesStart = getNumObjectsBeforeDistance(minObjectLMDist,i,candidatesEnd);
        }
        numCandidates = candidatesEnd - candidatesStart;
        if (numCandidates < minNumCandidates) {
            minNumCandidates = numCandidates;
            listIdx = i;
            if (numCandidates == 0) {
                break;
            }
        }
    }
    return listIdx;
}

#if defined(COLLECT_STATISTICS)
void COLTreeNode::updateQueueWithMinObj(int listIdx, int& leftPos, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, 
                                        BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, Statistics& stats) 
#else
void COLTreeNode::updateQueueWithMinObj(int listIdx, int& leftPos, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, 
                                        BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt) 
#endif
{
    int minObjPos = getIndexClosestToDistance(queryLmkDist,listIdx);
    VertexID minObject = lmkObjectLists[listIdx][minObjPos].first;
    EdgeWeight objectLMDist = lmkObjectLists[listIdx][minObjPos].second;
    EdgeWeight minObjLBDist = (objectLMDist > queryLmkDist) ? objectLMDist - queryLmkDist : queryLmkDist - objectLMDist;
    // Use ALT index to improve lower-bound from Object List
    EdgeWeight altMinObjLBDist = alt.getLowerBound(queryVertexID,minObject);
    if (altMinObjLBDist > minObjLBDist) {
        pqueue.insert(minObject,altMinObjLBDist);
#if defined(COLLECT_STATISTICS)
        stats.incrementStatistic("alt_used",1);
#endif
    } else {
        pqueue.insert(minObject,minObjLBDist);
    }
#if defined(COLLECT_STATISTICS)
    stats.incrementStatistic("cand_unpruned",1);
    stats.incrementStatistic("obj_inserted",1);
    stats.incrementStatistic("lbs_computed",1);
    stats.incrementStatistic("lbs_computed_obj",1);
    stats.incrementStatistic("alt_tried",1);
#endif
    
    leftPos = minObjPos-1;
    rightPos = minObjPos+1;

}

// Explore object list and update result (if needed) upto next smallest LB in priority queue
// Also updates leftPos and rightPos with next indexes to search in subsequent search
// Return false when object list has been fully evaluated
// If returning true, then updates listLB with minimum list lower-bound (otherwise it does not)
#if defined(COLLECT_STATISTICS)
bool COLTreeNode::updateQueue(int listIdx, int& leftPos, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, 
                              BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, EdgeWeight& listLB, Statistics& stats) 
#else
bool COLTreeNode::updateQueue(int listIdx, int& leftPos, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, 
                              BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, EdgeWeight& listLB) 
#endif
{
    bool continueLeft, continueRight, leftFinished, rightFinished;
    VertexID minObject;
    EdgeWeight minObjLBDist, altMinObjLBDist, objectLMDist;

    continueLeft = (leftPos >= 0) ? true : false;
    leftFinished = !continueLeft;
    continueRight = (static_cast<std::size_t>(rightPos) < lmkObjectLists[listIdx].size()) ? true : false;
    rightFinished = !continueRight;
    while (continueLeft || continueRight) {
        if (continueLeft) {
            minObject = lmkObjectLists[listIdx][leftPos].first;
            objectLMDist = lmkObjectLists[listIdx][leftPos].second;
            // Since we're going left from mid point, objectLMDist must be smaller than queryLmkDist
            minObjLBDist = queryLmkDist - objectLMDist;
            if (pqueue.size() == 0 || minObjLBDist <= pqueue.getMinKey()) {
                // Use ALT index to improve lower-bound from Object List
                altMinObjLBDist = alt.getLowerBound(queryVertexID,minObject);
                if (altMinObjLBDist > minObjLBDist) {
                    pqueue.insert(minObject,altMinObjLBDist);
#if defined(COLLECT_STATISTICS)
                    stats.incrementStatistic("alt_used",1);
#endif
                } else {
                    pqueue.insert(minObject,minObjLBDist);
                }
#if defined(COLLECT_STATISTICS)
                stats.incrementStatistic("cand_unpruned",1);
                stats.incrementStatistic("obj_inserted",1);
                stats.incrementStatistic("lbs_computed",1);
                stats.incrementStatistic("lbs_computed_obj",1);
                stats.incrementStatistic("alt_tried",1);
#endif

                --leftPos;
                if (leftPos < 0) {
                    continueLeft = false;
                    leftFinished = true;
                }
            } else {
                continueLeft = false;
            }
        }
        if (continueRight) {
            minObject = lmkObjectLists[listIdx][rightPos].first;
            objectLMDist = lmkObjectLists[listIdx][rightPos].second;
            // Since we're going right from mid point, objectLMDist must be larger than queryLmkDist
            minObjLBDist = objectLMDist - queryLmkDist;
            if (pqueue.size() == 0 || minObjLBDist <= pqueue.getMinKey()) {
                // Use ALT index to improve lower-bound from Object List
                altMinObjLBDist = alt.getLowerBound(queryVertexID,minObject);
                if (altMinObjLBDist > minObjLBDist) {
                    pqueue.insert(minObject,altMinObjLBDist);
#if defined(COLLECT_STATISTICS)
                    stats.incrementStatistic("alt_used",1);
#endif
                } else {
                    pqueue.insert(minObject,minObjLBDist);
                }
#if defined(COLLECT_STATISTICS)
                stats.incrementStatistic("cand_unpruned",1);
                stats.incrementStatistic("obj_inserted",1);
                stats.incrementStatistic("lbs_computed",1);
                stats.incrementStatistic("lbs_computed_obj",1);
                stats.incrementStatistic("alt_tried",1);
#endif

                ++rightPos;
                if (static_cast<std::size_t>(rightPos) == lmkObjectLists[listIdx].size()) {
                    continueRight = false;
                    rightFinished = true;
                }
            } else {
                continueRight = false;
            }            
        }
    }
    
    if (leftFinished == true && rightFinished == true) {
        // This means no object list has been fully evaluated
        return false;
    } else {
        EdgeWeight objectLMDist, tempLB;
        if (!leftFinished) {
            objectLMDist = lmkObjectLists[listIdx][leftPos].second;
            // Since we're going left from mid point, objectLMDist must be smaller than queryLmkDist
            listLB = queryLmkDist - objectLMDist;
        }
        if (!rightFinished) {
            objectLMDist = lmkObjectLists[listIdx][rightPos].second;
            // Since we're going right from mid point, objectLMDist must be larger than queryLmkDist
            tempLB = objectLMDist - queryLmkDist;
            if (leftFinished || tempLB < listLB) {
                listLB = tempLB;
            }
        }
        return true;
    }
}

#if defined(COLLECT_STATISTICS)
void COLTreeNode::updateQueueWithMinObjAkNN(int listIdx, int& leftPos, int& rightPos, std::vector<EdgeWeight>& queryVertexLmkDists, std::vector<VertexID>& queryVertices, 
                                            std::vector<EdgeWeight>& queryDistances, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, AggregateFunction<EdgeWeight>& agg, Statistics& stats) 
#else
void COLTreeNode::updateQueueWithMinObjAkNN(int listIdx, int& leftPos, int& rightPos, std::vector<EdgeWeight>& queryVertexLmkDists, std::vector<VertexID>& queryVertices, 
                                            std::vector<EdgeWeight>& queryDistances, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, AggregateFunction<EdgeWeight>& agg)
#endif
{
    int minObjPos = agg.getMinObjListIndex(lmkObjectLists[listIdx],queryVertexLmkDists);
    VertexID minObject = lmkObjectLists[listIdx][minObjPos].first;
    EdgeWeight objectLMDist = lmkObjectLists[listIdx][minObjPos].second;
    EdgeWeight minObjLBDist, altMinObjLBDist;
    for (std::size_t i = 0; i < queryVertices.size(); ++i) {
        minObjLBDist = (objectLMDist > queryVertexLmkDists[i]) ? objectLMDist - queryVertexLmkDists[i] : queryVertexLmkDists[i] - objectLMDist;
        // Use ALT index to improve lower-bound from Object List
        altMinObjLBDist = alt.getLowerBound(queryVertices[i],minObject);
        if (altMinObjLBDist > minObjLBDist) {
            queryDistances[i] = altMinObjLBDist;
#if defined(COLLECT_STATISTICS)
            stats.incrementStatistic("alt_used",1);
#endif
        } else {
            queryDistances[i] = minObjLBDist;
        }
    }
#if defined(COLLECT_STATISTICS)
    stats.incrementStatistic("cand_unpruned",1);
    stats.incrementStatistic("obj_inserted",1);
    stats.incrementStatistic("lbs_computed",queryVertices.size());
    stats.incrementStatistic("lbs_computed_obj",queryVertices.size());
    stats.incrementStatistic("alt_tried",queryVertices.size());
#endif
    EdgeWeight aggLBScore = agg.aggregate(queryDistances);
    pqueue.insert(minObject,aggLBScore);

    leftPos = minObjPos-1;
    rightPos = minObjPos+1;

}

// Explore object list and update result (if needed) upto next smallest LB in priority queue
// Also updates leftPos and rightPos with next indexes to search in subsequent search
// Return false when object list has been fully evaluated
// If returning true, then updates listLB with minimum list lower-bound (otherwise it does not)
#if defined(COLLECT_STATISTICS)
bool COLTreeNode::updateQueueAkNN(int listIdx, int& leftPos, int& rightPos, std::vector<EdgeWeight>& queryVertexLmkDists, std::vector<VertexID>& queryVertices, 
                                  std::vector<EdgeWeight>& queryDistances, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, EdgeWeight& listLB, AggregateFunction<EdgeWeight>& agg, Statistics& stats) 
#else
bool COLTreeNode::updateQueueAkNN(int listIdx, int& leftPos, int& rightPos, std::vector<EdgeWeight>& queryVertexLmkDists, std::vector<VertexID>& queryVertices, 
                                  std::vector<EdgeWeight>& queryDistances, BinaryMinHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, EdgeWeight& listLB, AggregateFunction<EdgeWeight>& agg)
#endif
{
    bool continueLeft, continueRight, leftFinished, rightFinished;
    VertexID minObject;
    EdgeWeight minObjLBScore, altMinObjLBDist, objectLMDist;

    continueLeft = (leftPos >= 0) ? true : false;
    leftFinished = !continueLeft;
    continueRight = (static_cast<std::size_t>(rightPos) < lmkObjectLists[listIdx].size()) ? true : false;
    rightFinished = !continueRight;
    while (continueLeft || continueRight) {
        if (continueLeft) {
            minObject = lmkObjectLists[listIdx][leftPos].first;
            objectLMDist = lmkObjectLists[listIdx][leftPos].second;
            for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                queryDistances[i] = (objectLMDist > queryVertexLmkDists[i]) ? objectLMDist - queryVertexLmkDists[i] : queryVertexLmkDists[i] - objectLMDist;
            }
            minObjLBScore = agg.aggregate(queryDistances);
            if (pqueue.size() == 0 || minObjLBScore <= pqueue.getMinKey()) {
                // Use ALT index to improve lower-bound from Object List (if possible)
                for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                    altMinObjLBDist = alt.getLowerBound(queryVertices[i],minObject);
                    if (altMinObjLBDist > queryDistances[i]) {
                        queryDistances[i] = altMinObjLBDist;
#if defined(COLLECT_STATISTICS)
                        stats.incrementStatistic("alt_used",1);
#endif
                    }
                }
                minObjLBScore = agg.aggregate(queryDistances);
                pqueue.insert(minObject,minObjLBScore);
#if defined(COLLECT_STATISTICS)
                stats.incrementStatistic("cand_unpruned",1);
                stats.incrementStatistic("obj_inserted",1);
                stats.incrementStatistic("lbs_computed",queryVertices.size());
                stats.incrementStatistic("lbs_computed_obj",queryVertices.size());
                stats.incrementStatistic("alt_tried",queryVertices.size());
#endif

                --leftPos;
                if (leftPos < 0) {
                    continueLeft = false;
                    leftFinished = true;
                }
            } else {
                continueLeft = false;
            }
        }
        if (continueRight) {
            minObject = lmkObjectLists[listIdx][rightPos].first;
            objectLMDist = lmkObjectLists[listIdx][rightPos].second;
            for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                queryDistances[i] = (objectLMDist > queryVertexLmkDists[i]) ? objectLMDist - queryVertexLmkDists[i] : queryVertexLmkDists[i] - objectLMDist;
            }
            minObjLBScore = agg.aggregate(queryDistances);
            if (pqueue.size() == 0 || minObjLBScore <= pqueue.getMinKey()) {
                // Use ALT index to improve lower-bound from Object List (if possible)
                for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                    altMinObjLBDist = alt.getLowerBound(queryVertices[i],minObject);
                    if (altMinObjLBDist > queryDistances[i]) {
                        queryDistances[i] = altMinObjLBDist;
#if defined(COLLECT_STATISTICS)
                        stats.incrementStatistic("alt_used",1);
#endif
                    }
                }
                minObjLBScore = agg.aggregate(queryDistances);
                pqueue.insert(minObject,minObjLBScore);
#if defined(COLLECT_STATISTICS)
                stats.incrementStatistic("cand_unpruned",1);
                stats.incrementStatistic("obj_inserted",1);
                stats.incrementStatistic("lbs_computed",queryVertices.size());
                stats.incrementStatistic("lbs_computed_obj",queryVertices.size());
                stats.incrementStatistic("alt_tried",queryVertices.size());
#endif

                ++rightPos;
                if (static_cast<std::size_t>(rightPos) == lmkObjectLists[listIdx].size()) {
                    continueRight = false;
                    rightFinished = true;
                }
            } else {
                continueRight = false;
            }
        }
    }
    
    if (leftFinished == true && rightFinished == true) {
        // This means no object list has been fully evaluated
        return false;
    } else {
        EdgeWeight objectLMDist, tempLB;
        if (!leftFinished) {
            objectLMDist = lmkObjectLists[listIdx][leftPos].second;
            for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                queryDistances[i] = (objectLMDist > queryVertexLmkDists[i]) ? objectLMDist - queryVertexLmkDists[i] : queryVertexLmkDists[i] - objectLMDist;
            }
            listLB = agg.aggregate(queryDistances);
        }
        if (!rightFinished) {
            objectLMDist = lmkObjectLists[listIdx][rightPos].second;
            for (std::size_t i = 0; i < queryVertices.size(); ++i) {
                queryDistances[i] = (objectLMDist > queryVertexLmkDists[i]) ? objectLMDist - queryVertexLmkDists[i] : queryVertexLmkDists[i] - objectLMDist;
            }
            tempLB = agg.aggregate(queryDistances);
            if (leftFinished || tempLB < listLB) {
                listLB = tempLB;
            }
        }
        return true;
    }
}

#if defined(COLLECT_STATISTICS)
void COLTreeNode::updateQueueWithMaxObjFkN(int listIdx, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, BinaryMaxHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, Statistics& stats) 
#else
void COLTreeNode::updateQueueWithMaxObjFkN(int listIdx, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, BinaryMaxHeap<EdgeWeight,VertexID>& pqueue, ALT& alt) 
#endif
{
    // As object lists are sorted, the largest upper-bound will be provided by the right-most object
    rightPos = lmkObjectLists[listIdx].size()-1;
    VertexID maxUBObj = lmkObjectLists[listIdx][rightPos].first;
    EdgeWeight maxUBDist = queryLmkDist + lmkObjectLists[listIdx][rightPos].second;
    EdgeWeight altMinObjUBDist = alt.getUpperBound(queryVertexID,maxUBObj);
    if (altMinObjUBDist < maxUBDist) {
        pqueue.insert(maxUBObj,altMinObjUBDist);
#if defined(COLLECT_STATISTICS)
        stats.incrementStatistic("alt_used",1);
#endif
    } else {
        pqueue.insert(maxUBObj,maxUBDist);
    }
#if defined(COLLECT_STATISTICS)
    stats.incrementStatistic("cand_unpruned",1);
    stats.incrementStatistic("obj_inserted",1);
    stats.incrementStatistic("lbs_computed",1);
    stats.incrementStatistic("lbs_computed_obj",1);
    stats.incrementStatistic("alt_tried",1);
#endif
    --rightPos; // Move pointer to object to the left (which will have the next largest upper-bound)
}

// Explore object list and update result (if needed) upto next largest UB in priority queue
// Also updates rightPos with next indexes to search in subsequent search
// Return false when object list has been fully evaluated
// If returning true, then updates listLB with minimum list lower-bound (otherwise it does not)
#if defined(COLLECT_STATISTICS)
bool COLTreeNode::updateQueueFkN(int listIdx,  int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, 
                                 BinaryMaxHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, EdgeWeight& listUB, Statistics& stats) 
#else
bool COLTreeNode::updateQueueFkN(int listIdx, int& rightPos, EdgeWeight queryLmkDist, VertexID queryVertexID, 
                                 BinaryMaxHeap<EdgeWeight,VertexID>& pqueue, ALT& alt, EdgeWeight& listUB) 
#endif
{
    bool continueSearch = (rightPos >= 0) ? true : false;
    VertexID maxUBObj;
    EdgeWeight maxUBDist, altMinObjUBDist;

    // We assume that rightPos has been moved to the next candidate already
    while (continueSearch) {
        maxUBObj = lmkObjectLists[listIdx][rightPos].first;
        maxUBDist = lmkObjectLists[listIdx][rightPos].second + queryLmkDist;
        if (pqueue.size() == 0 || maxUBDist >= pqueue.getMaxKey()) {
            // If maxUBDist is greater than max upper-bound found so then
            // we should insert the associated object into max priority queue.
            // When we insert, we use the ALT upper-bound, which may find a
            // tighter (smaller, more accurate) upper-bound for the object.
            // By using ALT in this way, we prevent candidates that are
            // just not ordered well in the object list, from mistakenly
            // terminating the loop sooner (we would end up back here anyway)
            altMinObjUBDist = alt.getUpperBound(queryVertexID,maxUBObj);
            if (altMinObjUBDist < maxUBDist) {
                pqueue.insert(maxUBObj,altMinObjUBDist);
#if defined(COLLECT_STATISTICS)
                stats.incrementStatistic("alt_used",1);
#endif
            } else {
                pqueue.insert(maxUBObj,maxUBDist);
            }
#if defined(COLLECT_STATISTICS)
            stats.incrementStatistic("cand_unpruned",1);
            stats.incrementStatistic("obj_inserted",1);
            stats.incrementStatistic("lbs_computed",1);
            stats.incrementStatistic("lbs_computed_obj",1);
            stats.incrementStatistic("alt_tried",1);
#endif
            --rightPos;
            if (rightPos < 0) {
                continueSearch = false;
            }
        } else {
            continueSearch = false; // No further object can improve maximum upper-bound
        }
    }

    if (rightPos < 0) {
        // This means object list has been fully evaluated
        return false;
    } else {
        // Otherwise node associated with object list will be re-insert into priority
        // we update listUB to indicate the next best upper-bound for remaining objects
        listUB = lmkObjectLists[listIdx][rightPos].second + queryLmkDist;
        return true;
    }
}

// This function retrieves all result objects from the object list. This leaf node does not need to be evaluated again
// as, unlike k object queries, whether a leaf node object is a result is completely indepedent of other objects.
// The query-to-landmark distance plus an object's distance to the landmark is an upper-bound on the query-to-object
// distance, and their difference is a lower-bound. As the OL is sorted by distance to the landmark, the objects with
// upper-bound within range are a prefix, which needs the landmark to be within range, the objects pruned for being
// closer to the landmark than the query by more than range are a prefix, which needs the landmark to be outside range,
// and the objects pruned for being further from it than the query by more than range are a suffix. So the candidates
// are the objects between the suffix and whichever prefix applies
#if defined(COLLECT_STATISTICS)
void COLTreeNode::rangeSearchObjectList(int listIdx, VertexID queryVertexID, EdgeWeight range, ALT& alt, PrunedHighwayLabeling& phl, std::vector<VertexID>& results, std::vector<VertexID>& resultDistances, Statistics& stats) 
#else
void COLTreeNode::rangeSearchObjectList(int listIdx, VertexID queryVertexID, EdgeWeight range, ALT& alt, PrunedHighwayLabeling& phl, std::vector<VertexID>& results, std::vector<VertexID>& resultDistances) 
#endif
{
    VertexID objectID;
    EdgeWeight queryLmkDist = phl.Query(queryVertexID,this->landmarks[listIdx]), spDist, lb, ub;
#if defined(COLLECT_STATISTICS)
    stats.incrementStatistic("spdist_calls",1);
#endif
    std::vector<VertexEdgeWeightPair>& lmkObjectList = this->lmkObjectLists[listIdx];

    std::size_t elementIdx = 0;
    if (queryLmkDist <= range) {
        // Landmark is within range, so return objects while their upper-bound is within range
        // Note: No need to compute spDist because the upper-bound is already within range
        while (elementIdx < lmkObjectList.size() && queryLmkDist + lmkObjectList[elementIdx].second <= range) {
            results.push_back(lmkObjectList[elementIdx].first);
            resultDistances.push_back(0);
            ++elementIdx;
        }
    } else {
        // Landmark is outside range, so skip objects closer to it than minObjectLMDist (lower-bound greater than range)
        EdgeWeight minObjectLMDist = queryLmkDist - range;
        elementIdx = getNumObjectsBeforeDistance(minObjectLMDist,listIdx,lmkObjectList.size());
    }

    // Evaluate the candidates until an object is further from the landmark than maxObjectLMDist, as all
    // remaining objects then have lower-bound greater than range
    EdgeWeight maxObjectLMDist = queryLmkDist + range;
    for (; elementIdx < lmkObjectList.size() && lmkObjectList[elementIdx].second <= maxObjectLMDist; ++elementIdx) {
        objectID = lmkObjectList[elementIdx].first;
        alt.getLowerAndUpperBound(queryVertexID,objectID,lb,ub);
#if defined(COLLECT_STATISTICS)
        stats.incrementStatistic("cand_unpruned",1);
        stats.incrementStatistic("lbs_computed",1);
        stats.incrementStatistic("lbs_computed_obj",1);
#endif
        if (ub <= range) {
            results.push_back(objectID);
            resultDistances.push_back(0);
        } else if (lb <= range) {
            spDist = phl.Query(queryVertexID,objectID);
#if defined(COLLECT_STATISTICS)
            stats.incrementStatistic("cand_spdist",1);
            stats.incrementStatistic("spdist_calls",1);
#endif
            if (spDist <= range) {
                results.push_back(objectID);
                resultDistances.push_back(0);
            }
#if defined(COLLECT_STATISTICS)
            // This means the candidate was the bad choice
            else {
                stats.incrementStatistic("false_hit",1);
            }
#endif
        }
        // Note: if lb > range, we prune the object but we can't terminate
    }
}