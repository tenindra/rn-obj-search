// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "Gtree.h"
#include "DynamicGraph.h"

#include "DijkstraSearch.h"
#include "../utility/StopWatch.h"
#include "../utility/utility.h"

#include <cstdlib>
#include <iostream>

Gtree::Gtree(std::string networkName, int numVertices, int numEdges, int fanout, std::size_t maxLeafSize): 
    networkName(networkName), numVertices(numVertices), numEdges(numEdges), fanout(fanout), maxLeafSize(maxLeafSize) {
    this->vertexIDLeafVerticesVecIdx.assign(this->numVertices,-1);
    this->vertexIDLeafBordersVecIdx.assign(this->numVertices,-1);
    this->edgeInLeafSubgraph.assign(this->numEdges,true);
    this->vertexLeafIdxs.assign(this->numVertices,-1);
}

Gtree::Gtree() {}

std::string Gtree::getNetworkName()
{
    return this->networkName;
}

int Gtree::getNumEdges()
{
    return this->numEdges;
}

int Gtree::getNumVertices()
{
    return this->numVertices;
}

int Gtree::getFanout()
{
    return this->fanout;
}

int Gtree::getTreeSize()
{
    return this->treeNodes.size();
}

std::size_t Gtree::getMaxLeafSize()
{
     return this->maxLeafSize;
}

int Gtree::getLeafIndex(VertexID vertexID)
{
    return this->vertexLeafIdxs[vertexID];
}

int Gtree::getParentIndex(int treeIdx)
{
    return this->treeNodes[treeIdx].getParentIdx();
}

void Gtree::buildGtree(Graph& graph){ 
    
    this->buildTreeHierarchy(graph);
    this->computeDistanceMatrix(graph);
    this->initialiseGtreeQueryStructure();
    
}

void Gtree::buildTreeHierarchy(Graph& graph) 
{
    METISGraphPartitioner partitioner(this->numVertices,this->numEdges);

    std::unordered_set<VertexID> subgraph = graph.getVertexIDsUset();
    this->addNode(ROOT_PARENT_INDEX,subgraph,graph,partitioner);
}

void Gtree::addNode(int parentIdx, std::unordered_set<VertexID>& subgraph, Graph& originalGraph, GraphPartitioner& partitioner)
{
    // Create Gtree node and add to tree
    int treeIdx = this->treeNodes.size(), childvecIdx;
    GtreeNode treeNode(parentIdx,subgraph.size());
    this->treeNodes.push_back(treeNode);

    if (parentIdx != ROOT_PARENT_INDEX) {
        // Add this Gtree node as child of parent (if not root)
        childvecIdx = this->treeNodes[parentIdx].addChild(treeIdx);
        this->treeNodes[treeIdx].setParentChildIdx(childvecIdx);
        
        // Add parent's Gtree path to this nodes Gtree path then add parent
        this->treeNodes[treeIdx].addGtreePathNodes(this->treeNodes[parentIdx].gtreePath);
    }
    this->treeNodes[treeIdx].addGtreePathNode(treeIdx);
    
    // Find and add children of this node
    if (subgraph.size() > this->maxLeafSize) {
        this->addChildren(treeIdx, subgraph, originalGraph,partitioner);
    } else {
        // If this has less than tau graph vertex then we stop partitions
        this->treeNodes[treeIdx].setLeafNode();
        this->leafIdxs.push_back(treeIdx);
    }
    
    // Determine and add borders for this Gtree node
    if (parentIdx != ROOT_PARENT_INDEX) {
        // We mark the starting index of this tree nodes borders
        // in the parent's childBordersVec - we can do this as
        // we are about find all these borders and also add 
        // consecutively to the parent's childBordersVec
        this->treeNodes[parentIdx].addChildOffsetInChildBorderVec();
        int adjListStart, adjListSize, leafVerticesVecIdx, leafBordersVecIdx;
        bool isBorder;
        
        for (VertexID vertex: subgraph) {
            // If this node is a leaf node, we also need to map 
            // all it's graph vertices to they're Gtree leaf index
            if (this->treeNodes[treeIdx].isLeafNode()) {
                // Add each vertex as a leaf vertex (this
                // will include borders)
                leafVerticesVecIdx = this->treeNodes[treeIdx].addLeafVertex(vertex);
                this->setLeafVerticesVecIdx(vertex,leafVerticesVecIdx);
                this->vertexLeafIdxs[vertex] = treeIdx;
            }
            
            adjListStart = originalGraph.getEdgeListStartIndex(vertex);
            adjListSize = originalGraph.getEdgeListSize(vertex);
            isBorder = false;
            for (int i = adjListStart; i < adjListSize; ++i) {
                if (subgraph.find(originalGraph.edges[i].first) == subgraph.end()) {
                    // This mean we have found a neighbour of that outside this subgraph i.e. 
                    // this is vertex is a border so we add to list of borders of this tree node
                    // and the set of child borders for the parent
                    if (!isBorder) {
                        leafBordersVecIdx = this->treeNodes[treeIdx].addBorder(vertex);
                        if (this->treeNodes[treeIdx].isLeafNode()) {
                            this->setLeafBordersVecIdx(vertex,leafBordersVecIdx);
                        }
                        this->treeNodes[parentIdx].addChildBorder(vertex);
                        isBorder = true; // So that we don't add it again
                    }
                    if (this->treeNodes[treeIdx].isLeafNode()) {
                        // Since this edge leads to a vertex that is not in the subgraph
                        this->setEdgeNotInSubgraph(i);
                    }
                }
            }
        }

    }
    
}

void Gtree::addChildren(int parentTreeIdx, std::unordered_set<VertexID>& parentGraph, Graph& originalGraph, GraphPartitioner& partitioner) {
    // Assuming original graph is undirected graph (if not we must make undirected
    // in order for METIS to function)
    
    std::vector<VertexID> vertexOrder;
    std::vector<int> parts;
    partitioner.partition(originalGraph,parentGraph,this->fanout,vertexOrder,parts);

    // Create sets of subgraphs from partitioned parent graph
    // and add them to Gtree (and this will recurse)
    std::vector<std::unordered_set<VertexID>> childGraphs(this->fanout);

    for (std::size_t i = 0; i < vertexOrder.size(); ++i) {
        childGraphs[parts[i]].insert(vertexOrder[i]);
    }
    
    // Also add children of this node to Gtree
    for (std::size_t i = 0; i < childGraphs.size(); ++i) {
        this->addNode(parentTreeIdx,childGraphs[i],originalGraph,partitioner);
        // Note: We release child graph at this point as it is not need again
        // and this would improve overall memory usage (it won't affect vector size)
        utility::releaseSTLCollection(childGraphs[i]);
    }
}

void Gtree::computeDistanceMatrix(Graph& graph)
{
    std::vector<std::vector<int>> treeLevelIdxs = this->getTreeNodesByLevel();
    
    DijkstraSearch dijkstra;
    BinaryMinHeap<EdgeWeight,VertexID> *pqueue = new BinaryMinHeap<EdgeWeight,VertexID>();
    int currentIdx;
    DynamicGraph tempGraph(graph);
    std::unordered_set<VertexID>  *targets;
    std::vector<VertexID> adjVertices;
    std::vector<EdgeWeight> adjVertexWgts;
    std::vector<VertexID> *sourcesVec, *targetsVec;
    
    for (int i = treeLevelIdxs.size()-1; i >= 0; --i) {
        // Clear memory in unordered_map or it will continue to grow
        // Note: According to the paper, total number of borders at each level should be O(n)
        // so if we clear this map for each level then we should it's total size should be O(n)
            
        for (std::size_t j = 0; j < treeLevelIdxs[i].size(); ++j) {
            currentIdx = treeLevelIdxs[i][j];
            if (this->treeNodes[currentIdx].isLeafNode()) {
                // In a leaf we find distances from border to all leaf vertices
                sourcesVec = &this->treeNodes[currentIdx].getBorders();
                targets = &this->treeNodes[currentIdx].getLeafVerticesUset();
                targetsVec = &this->treeNodes[currentIdx].getLeafVertices();
            } else {
                // In a non-leaf we get distances from all child borders to all other child borders
                // Note: It's possible some of these distances have already been computed
                sourcesVec = &this->treeNodes[currentIdx].getChildBorders();
                targets = &this->treeNodes[currentIdx].getChildBordersUset();
                targetsVec = &this->treeNodes[currentIdx].getChildBorders();
            }
            
            std::unordered_map<VertexID,EdgeWeight> siblingBorderDistances;
            siblingBorderDistances.reserve(targetsVec->size());
            //Using single std::unordered_map
            int rowLength = targetsVec->size();
            this->treeNodes[currentIdx].matrixRowLength = rowLength;
            this->treeNodes[currentIdx].distanceMatrix.reserve(sourcesVec->size()*rowLength);
            if (this->treeNodes[currentIdx].isLeafNode()) {
                // If it is a leaf node we can search using original Graph
                // who's data structure is faster than DynamicGraph
                for (std::size_t i = 0; i < sourcesVec->size(); ++i) {
                    pqueue->clear();
                    dijkstra.findSSMTDistances(graph,(*sourcesVec)[i],(*targets),siblingBorderDistances,pqueue);
                    for (std::size_t j = 0; j < targetsVec->size(); ++j) {
                        this->treeNodes[currentIdx].distanceMatrix.push_back(siblingBorderDistances[(*targetsVec)[j]]);
                    }
                }
            } else {
                for (std::size_t i = 0; i < sourcesVec->size(); ++i) {
                    pqueue->clear();
                    dijkstra.findSSMTDistances(tempGraph,(*sourcesVec)[i],(*targets),siblingBorderDistances,pqueue);
                    for (std::size_t j = 0; j < targetsVec->size(); ++j) {
                        this->treeNodes[currentIdx].distanceMatrix.push_back(siblingBorderDistances[(*targetsVec)[j]]);
                    }
                }
            }
            
            // All future searches will be on this nodes border set (using closure property in paper)
            // because these will be the parent nodes child borders-> Therefore if we need only remove
            // unnecessary edges from these vertices (unimportant vertices will become disconnected)
            sourcesVec = &this->treeNodes[currentIdx].getBorders();

            VertexID border;
            int sourceIdx, targetIdx, distMatrixRowOffset;
            for (std::size_t i = 0; i < sourcesVec->size(); ++i) {
                border = (*sourcesVec)[i];
                if (this->treeNodes[currentIdx].isLeafNode()) {
                    sourceIdx = i;
                    distMatrixRowOffset = sourceIdx*this->treeNodes[currentIdx].matrixRowLength;
                } else {
                    sourceIdx = this->treeNodes[currentIdx].getBorderIdxInChildBorderVec(i);
                    distMatrixRowOffset = sourceIdx*this->treeNodes[currentIdx].matrixRowLength;
                }
                
                // Preserve edges to outside this gtree node (i.e. subgraph)
                adjVertices.clear();
                adjVertexWgts.clear();
                for (std::size_t i = 0; i < tempGraph.vertices[border].adjVertices.size(); ++i) {
                    if (targets->find(tempGraph.vertices[border].adjVertices[i]) == targets->end()) {
                        // This check whether the adj vertex is within the current gtree node
                        // if it is we do not need to preserve the edge
                        adjVertices.push_back(tempGraph.vertices[border].adjVertices[i]);
                        adjVertexWgts.push_back(tempGraph.vertices[border].adjVertexWgts[i]);
                    }
                }
                // Note: That this will make the graph disconnected but this doesn't
                // matter as removing disconnected vertex won't change search results and
                // we are only disconnecting vertices that are not borders of the subgraph
                // (i.e. they will not be needed again at higher levels)
                
                tempGraph.vertices[border].adjVertices = std::move(adjVertices);
                tempGraph.vertices[border].adjVertexWgts = std::move(adjVertexWgts);
                for (std::size_t j = 0; j < sourcesVec->size(); ++j) {
                    targetIdx = this->treeNodes[currentIdx].getBorderIdxInChildBorderVec(j); // We will return dist matrix idx whether leaf or not
                    if (border != (*sourcesVec)[j]) {
                        tempGraph.insertImaginaryNonInvertibleEdge(border,(*sourcesVec)[j],this->treeNodes[currentIdx].distanceMatrix[distMatrixRowOffset+targetIdx]);
                    }
                }
            }
        }
    }
    
    delete pqueue;
}

void Gtree::initialiseGtreeQueryStructure()
{
    // Allocate Memory for Gtree Query
    this->sourceToTreeNodeBorderDist.resize(this->treeNodes.size());
    for (size_t i = 0; i < this->treeNodes.size(); ++i) {
        this->sourceToTreeNodeBorderDist[i].resize(this->treeNodes[i].bordersVec.size());
    }
}

int Gtree::getNumLevels()
{
    return this->getTreeNodesByLevel().size();
}

int Gtree::getNumBorders()
{
    int numBorders = 0;
    for (std::size_t i = 0; i < this->treeNodes.size(); ++i) {
        numBorders += this->treeNodes[i].getNumBorders();
    }
    return numBorders;
}

int Gtree::getBorderToBorderRelationships()
{
    int numB2BRelationships = 0;
    for (std::size_t i = 0; i < this->treeNodes.size(); ++i) {
        numB2BRelationships += this->treeNodes[i].getNumBorders()*this->treeNodes[i].getNumBorders();
    }
    return numB2BRelationships;
}

int Gtree::getAvgPathCost(int treeIdx)
{
    int avgPathCost = 0;
    if (this->treeNodes[treeIdx].isLeafNode()) {
        avgPathCost = this->treeNodes[treeIdx].getNumBorders();
    } else {
        int totalChildPathCosts = 0, totalParentChildCost = 0;
        for (std::size_t i = 0; i < this->treeNodes[treeIdx].children.size(); ++i) {
            totalParentChildCost += this->treeNodes[treeIdx].getNumBorders()*this->treeNodes[this->treeNodes[treeIdx].children[i]].getNumBorders();
            totalChildPathCosts += this->getAvgPathCost(this->treeNodes[treeIdx].children[i]);
        }
        avgPathCost += totalParentChildCost/this->treeNodes[treeIdx].children.size();
        avgPathCost += totalChildPathCosts/this->treeNodes[treeIdx].children.size();
    }
    return avgPathCost;
}

std::vector< std::vector< int > > Gtree::getTreeNodesByLevel()
{
    std::vector<std::vector<int>> treeNodeLevel;
    std::vector<int> currentLevel, nextLevel;
    nextLevel.push_back(0);

    while(nextLevel.size() != 0){
            treeNodeLevel.push_back(nextLevel);
            currentLevel.swap(nextLevel);
            nextLevel.clear();
            for (int i: currentLevel){
                for (std::size_t j = 0; j < this->treeNodes[i].children.size(); j++ ){
                    nextLevel.push_back( this->treeNodes[i].children[j] );
                }
            }
    }

    return treeNodeLevel;
}

EdgeWeight Gtree::getShortestPathDistance(Graph& graph, VertexID u, VertexID v)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("computations_materialized",0);
    this->stats.initialiseStatistic("computations_executed",0);
    this->stats.initialiseStatistic("computations_total",0);
#endif    
    EdgeWeight spDist = 0;
    int uLeaf = this->getLeafIndex(u);
    int vLeaf = this->getLeafIndex(v);
    
    if (uLeaf == vLeaf) {
        // This mean they are both in the same node
        spDist = this->SPDistLeaf(u,v,uLeaf,graph);
    } else {
        int LCAIdx = 0;
        std::vector<int>& uPathFromRoot = this->treeNodes[uLeaf].getGtreePathFromRoot();
        std::vector<int>& vPathFromRoot = this->treeNodes[vLeaf].getGtreePathFromRoot();

        // Since the two vertices are in different leaf nodes we can guarantee that
        // there is at least one more node (the parent of both) in the G-tree path
        // which we call the the lowest common ancestor (LCA)
        
        // Search down Gtree from root until we find first ancestor that is different
        // this means the previous ancestor is the LCA and the indexes point its children
        unsigned int i, j;
        for (i = 0, j = 0; i < uPathFromRoot.size() && j < vPathFromRoot.size(); ++i, ++j) {
            if (uPathFromRoot[i] != vPathFromRoot[j]) {
                // When i = 0 and j = 0 it is referring to the root so in that
                // case uPathFromRoot[i] does equal vPathFromRoot[j]. This means
                // when they are equal i > 0, so i-1 is safe here
                LCAIdx = uPathFromRoot[i-1];
                break;
            }
            // Note: We can guarantee that LCAIdx will be set here. The only situation
            // it would not be set is both u and v were in the same leaf but we have
            // guarnateed that is not the case in the if/else statement
        }
        
        // Now search up G-tree (from source leaf) until we reach i and j, then we 
        // search down (to target leaf) computing shortest path distance
        
        // From source to source leaf
        this->SPDistToSourceLeafNode(u,uLeaf);
        
        // From source leaf to first child of LCA
        int x = i; // This is safe, depth of tree is O(log(n)) and n is at most 24 million in US dataset
        for (int k = uPathFromRoot.size()-1; k > x; --k) {
            // Since k > x and x is at worst 0, k-1 is safe here
            this->SPDistToParentNode(uPathFromRoot[k],uPathFromRoot[k-1],false);
        }
        
        // From first child of LCA to second child of LCA
        this->SPDistToSiblingNode(uPathFromRoot[i],vPathFromRoot[j],LCAIdx,false);
    
        // From second child of LCA to target leaf
        for (std::size_t k = j; k < vPathFromRoot.size()-1; ++k) {
            // Note the size()-1 in the above condition
            this->SPDistToChildNode(vPathFromRoot[k+1],vPathFromRoot[k],false);
        }
        
        // We assume target has not been visited
        spDist = this->SPDistToLeafTarget(v,vLeaf);
    }
    
    return spDist;
}

EdgeWeight Gtree::SPDistLeaf(VertexID u, VertexID v, int leafNode, Graph& graph)
{
    // Assume u and v are in the same leaf node
    EdgeWeight borderToBorderDist, spDist = 0;
    
    spDist = this->DijkstraDist(u,v,leafNode,graph);
    borderToBorderDist = this->BorderDist(u,v,leafNode);
    if (borderToBorderDist < spDist) {
        spDist = borderToBorderDist;
    }
    
    return spDist;
}

EdgeWeight Gtree::BorderDist(VertexID u, VertexID v, int leafNode)
{
    // Assume that u and v ARE in the same leaf node
    
    EdgeWeight borderDist, minDist = 0;
    
    int uIdx = this->getIdxInLeafVerticesVec(u);
    int vIdx = this->getIdxInLeafVerticesVec(v);
    int intermediateBorderIdx;
    if (this->treeNodes[leafNode].bordersVec.size() > 0) {
        intermediateBorderIdx = this->treeNodes[leafNode].getBorderIdxInChildBorderVec(0);
        minDist = this->treeNodes[leafNode].distanceMatrix[uIdx] 
            + this->treeNodes[leafNode].distanceMatrix[intermediateBorderIdx] 
            + this->treeNodes[leafNode].distanceMatrix[vIdx];
        for (std::size_t i = 0; i < this->treeNodes[leafNode].bordersVec.size(); ++i) {
            for (std::size_t j = 0; j < this->treeNodes[leafNode].bordersVec.size(); ++j) {
                // ASSUMING UNDIRECTED GRAPH (I.E. BORDER TO LEAF VERTEX DISTANCE = LEAF VERTEX TO BORDER DISTANCE
                // THIS IS BECAUSE WE HAVEN'T CALCULATED LEAF VERTEX TO BORDER DISTANCES
                if (i != j) {
                    // BorderDist return the distance between two vertices using outside vertices
                    // if i == j then we are entering and leaving through the same border
                    // which won't be smaller than the Dijkstra dist so we ignore it
                    intermediateBorderIdx = this->treeNodes[leafNode].getBorderIdxInChildBorderVec(j);
                    borderDist = this->treeNodes[leafNode].distanceMatrix[i*this->treeNodes[leafNode].matrixRowLength+uIdx] 
                        + this->treeNodes[leafNode].distanceMatrix[i*this->treeNodes[leafNode].matrixRowLength+intermediateBorderIdx] 
                        + this->treeNodes[leafNode].distanceMatrix[j*this->treeNodes[leafNode].matrixRowLength+vIdx];
                    if (borderDist < minDist) {
                        minDist = borderDist;
                    }
                }
            }
        }
    }

#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[leafNode].bordersVec.size()*this->treeNodes[leafNode].bordersVec.size());
#endif

    return minDist;
}

EdgeWeight Gtree::DijkstraDist(VertexID u, VertexID v, int leafNode, Graph& graph)
{
    // Assume that u and v ARE in the same leaf node
    
    DijkstraSearch dijk;
    
    return dijk.findShortestPathDistanceSubgraph(graph,u,v,this->edgeInLeafSubgraph);
}

// Note: The returned path does not include the LCA
std::vector<int> Gtree::getGtreePath(VertexID u, VertexID v, int& firstLCAChild)
{
    // Note: If u and v are in the same leaf node then path will be 
    // size 0 and LCA will contain leaf node
    std::vector<int> path;
    
    int uLeaf = this->getLeafIndex(u);
    int vLeaf = this->getLeafIndex(v);
    
    std::vector<int>& uPathFromRoot = this->treeNodes[uLeaf].getGtreePathFromRoot();
    std::vector<int>& vPathFromRoot = this->treeNodes[vLeaf].getGtreePathFromRoot();

    // The only case where firstLCAChild is not set is when both vertices are in the same leaf
    // Note: This will be overwritten later if it needs to be
    firstLCAChild = uPathFromRoot.back();
    // Note: Using back() is safe as uPathFromRoot will have at least root in it
    
    // Search down Gtree from root until we find first ancestor that is different
    // this means the previous ancestor is the lowest common ancestor
    unsigned int i, j;
    for (i = 0, j = 0; i < uPathFromRoot.size() && j < vPathFromRoot.size(); ++i, ++j) {
        if (uPathFromRoot[i] != vPathFromRoot[j]) {
            firstLCAChild = uPathFromRoot[i];
            break;
        }
    }
    
    // Note: That we do not include LCA in return path
    int x = static_cast<int>(i); // This is OK, i is greater than 0 as above
    for (int k = uPathFromRoot.size()-1; k >= x; --k) {
        // i represents the first non-common ancestor
        // so this will be the part of the path
        path.push_back(uPathFromRoot[k]);
    }

    // The path to v is in reverse order (we descend down)
    for (std::size_t k = j; k < vPathFromRoot.size(); ++k) {
        // j represents the first non-common ancestor
        // so this will be the part of the path
        path.push_back(vPathFromRoot[k]);
    }

    return path;

}

void Gtree::getKNNs(OccurenceList& occList, unsigned int k, VertexID queryVertex, std::vector<VertexID>& kNNs, 
                    std::vector<EdgeWeight>& kNNDistances, Graph& graph)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("computations_materialized",0);
    this->stats.initialiseStatistic("computations_executed",0);
    this->stats.initialiseStatistic("computations_total",0);
#endif
    // Note: We store both G-tree node indexes and road network VertexIDs in the
    // same priority queue. We  map the G-tree node indexes to a VertexID that 
    // is larger than the last vertex in the graph (and we check whenever 
    // a dequeue happens to determine what type min element is).
    BinaryMinHeap<EdgeWeight,VertexID> pqueue;
    
    int sourceLeaf = this->getLeafIndex(queryVertex);
    int Tn, prevTn, treeIdx;
    EdgeWeight Tmin, dist, minKey;
    VertexID minNode;
    
    Tn = sourceLeaf;
    Tmin = this->SPDistToSourceLeafNode(queryVertex,Tn);
    
    // Note: In order to use a priority queue with primitive types we
    // map the tree indexes as VertexID greater than the largest VertexID
    // for the current graph. Whenever we deque we can tell if the element
    // is a tree node or graph vertex by check if it's larger than the
    // number of vertices in the graph (since they are number 0 to n-1)
    VertexID firstMappedTreeNodeID = static_cast<VertexID>(this->numVertices);
    
    // If the source leaf node contains any object we must find them first
    if (occList.leafOccurenceList[sourceLeaf].size() > 0) {
        if (this->getSourceLeafkNNsByINE(queryVertex,k,occList.leafOccurenceSet[sourceLeaf],sourceLeaf,graph,kNNs,kNNDistances,pqueue)) {
            // If the in-leaf INE search was able to find definite k objects then we exist now
            return;
        }     
    }
    
    while (pqueue.size() > 0 || Tn != 0) {
        // UpdateT
        if (pqueue.size() == 0) {
            prevTn = Tn;
            Tn = this->getParentIndex(Tn);
            
            if (Tn != 0) {
                // We don't set Tmin if Tn is the root because this represents
                // infinite instead - it's faster to check Tn != 0 than Tmin < max int
                Tmin = this->SPDistToParentNode(prevTn,Tn);
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,prevTn));
#endif
            }
            
            // Tn is a leaf only at start, but pqueue is size 1 at this point
            for(int childIdx: occList.nonLeafOccurenceList[Tn]) {
                if (childIdx != prevTn) {
                    // We only add children that we haven't already visited
                    dist = this->SPDistToSiblingNode(prevTn,childIdx,Tn);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,prevTn));
#endif
                    pqueue.insert(childIdx+this->numVertices,dist);
                }
            }
            
        }
        
        if (pqueue.size() > 0) {
            // This check is different to paper, but is necessary
            // as above UpdateT will not guarantee pqueue will have
            // size greater than 0 - in the case where the previous
            // Tn parent doesn't have any children with objects
            minKey = pqueue.getMinKey();
            minNode = pqueue.extractMinElement();
            if (Tn != 0 && minKey > Tmin) {
                // Tn == 0 is equivalent to Tmin being infinity
                // So using Tn != 0 we can avoid setting Tmin to infinity (i.e. max int value)
                // which was seen to slow down comparisons

                // UpdateT
                prevTn = Tn;
                Tn = this->getParentIndex(Tn);
                if (Tn != 0) {
                    // We don't set Tmin if Tn is the root because this represents
                    // infinite instead - it's faster to check Tn != 0 than Tmin < max int
                    Tmin = this->SPDistToParentNode(prevTn,Tn);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,prevTn));
#endif
                }
                
                // Tn is guaranteed to be a non-leaf here because we move to parent
                for(int childIdx: occList.nonLeafOccurenceList[Tn]) {
                    if (childIdx != prevTn) {
                        dist = this->SPDistToSiblingNode(prevTn,childIdx,Tn);
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,prevTn));
#endif
                        pqueue.insert(childIdx+this->numVertices,dist);
                    }
                }
                
                pqueue.insert(minNode,minKey);
            } else if (minNode < firstMappedTreeNodeID) {
                // Min element is a graph vertex (i.e. object)
                kNNs.push_back(minNode);
                kNNDistances.push_back(minKey);
                if (kNNs.size() == k) {
                    // If this is the kth nearest neighbour object
                    // then we no longer need to iterate
                    break;
                }
            } else {
                treeIdx = minNode-this->numVertices;
                if (this->treeNodes[treeIdx].isLeafNode()) {
                    for(VertexID leafObject: occList.leafOccurenceList[treeIdx]) {
                        dist = this->SPDistToLeafTarget(leafObject,treeIdx);
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,treeIdx));
#endif
                        pqueue.insert(leafObject,dist);
                    }
                } else {
                    for(int childIdx: occList.nonLeafOccurenceList[treeIdx]) {
                        dist = this->SPDistToChildNode(childIdx,treeIdx);
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,treeIdx));
#endif
                        pqueue.insert(childIdx+this->numVertices,dist);
                    }
                }
            }
        }
    }
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_total",this->stats.getStatistic("computations_materialized")+this->stats.getStatistic("computations_executed"));
#endif
}

bool Gtree::getSourceLeafkNNsByINE(VertexID queryVertex, unsigned int k, std::unordered_set<VertexID>& targets, int leafNode, Graph& graph, std::vector<VertexID>& kNNs, 
                                   std::vector<EdgeWeight>& kNNDistances, BinaryMinHeap<EdgeWeight,VertexID>& pqueue)
{
    // Note: We assumpe vectors passed to this method are empty

    BinaryMinHeap<EdgeWeight,VertexID> localQueue;
    EdgeWeight minDist;
    VertexID minDistVertexID, adjVertex;
    // We use unordered_set because we are only searching subgraph (no benefit 
    // from allocating std::vector<bool> for all vertices, just overhead)
    std::unordered_set<VertexID> settledVertexIDs;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with query vertex ID
    localQueue.insert(queryVertex,0);
    
    bool borderEncountered = false;
    unsigned int targetsFound = 0;
    
    // We also exist if all target vertices have been found
    while (localQueue.size() > 0 && targetsFound < targets.size()) {
        // Extract and remove vertex with smallest distance from query point
        // and mark it as "settled" so we do not inspect again
        minDist = localQueue.getMinKey();
        minDistVertexID = localQueue.extractMinElement();
        if (settledVertexIDs.find(minDistVertexID) == settledVertexIDs.end()) {
            settledVertexIDs.insert(minDistVertexID);
            
            if (targets.find(minDistVertexID) != targets.end()) {
                ++targetsFound;
                // If the minimum is an object we have found a kNN
                if (!borderEncountered) {
                    // If we have not encounted a border so far then
                    // all objects we encounter in the source leaf
                    // cannot be bettered by a object outside the source leaf
                    // so we can count these as kNN results
                    kNNs.push_back(minDistVertexID);
                    kNNDistances.push_back(minDist);
                } else {
                    pqueue.insert(minDistVertexID,minDist);
                }
                if (kNNs.size() == k) {
                    return true;
                } else if (targetsFound == k) {
                    // We do not need to add more than k candidates within leaf
                    // because at best these are the actual kNN neighbours and at
                    // worst some of them will be replaced by closer kNNs outside this leaf
                    break;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);
            // Note: We have to make sure we don't exceed size of graph.edges vector
            
            bool isBorder = false;
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                if (this->isEdgeInLeafSubgraph(i)) {
                    if (settledVertexIDs.find(adjVertex) == settledVertexIDs.end()) {
                        // This is the first time we have seen this vertex (i.e. distance infinity)
                        localQueue.insert(adjVertex,minDist+graph.edges[i].second);
                    }
                } else {
                    if (!isBorder) {
                        isBorder = true;
                    }
                    if (!borderEncountered) {
                        borderEncountered = true;
                    }
                }
            }
            
            // If it is a border we add each of the other borders for this leaf
            // node and the distance from this border (assuming they have not 
            // already been settled - if they are settled we have already found
            // shortest path to them so the path through this border cannot
            // be an improvement)
            if (isBorder) {
                int borderIdx = this->getIdxInLeafBordersVec(minDistVertexID);
                for (std::size_t i = 0; i < this->treeNodes[leafNode].bordersVec.size(); ++i) {
                    // Note since this is leaf the above functions actual gets idx in leafVerticesVec
                    // Recall that the dimension of the distance matrix is bordersVec.size()*leafVerticesVec.size()
                    // and stores the distances from borders to leaf vertices so we can use these values
                    // compute the distances between two borders in the leaf vertice
                    if (settledVertexIDs.find(this->treeNodes[leafNode].bordersVec[i]) == settledVertexIDs.end()) {
                        int targetBorderIdx = this->treeNodes[leafNode].getBorderIdxInChildBorderVec(i);
                        localQueue.insert(this->treeNodes[leafNode].bordersVec[i],minDist+this->treeNodes[leafNode].distanceMatrix[borderIdx*this->treeNodes[leafNode].matrixRowLength+targetBorderIdx]);
                    }
                }
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("computations_executed",this->treeNodes[leafNode].bordersVec.size());
#endif
            }
        }
    }
    return false;
}

EdgeWeight Gtree::SPDistToSourceLeafNode(VertexID u, int sourceLeafIdx)
{
    EdgeWeight spDist = 0, sourceToNextBorderDist;
    int uIdx = this->getIdxInLeafVerticesVec(u);
    if (this->treeNodes[sourceLeafIdx].bordersVec.size() > 0) {
        spDist = this->treeNodes[sourceLeafIdx].distanceMatrix[uIdx];
        this->sourceToTreeNodeBorderDist[sourceLeafIdx][0] = spDist;
        for (std::size_t i = 1; i < this->treeNodes[sourceLeafIdx].bordersVec.size(); ++i) {
            sourceToNextBorderDist = this->treeNodes[sourceLeafIdx].distanceMatrix[i*this->treeNodes[sourceLeafIdx].matrixRowLength+uIdx];
            // ASSUMING DIRECTED GRAPH (I.E. BORDER TO LEAF VERTEX DISTANCE = LEAF VERTEX TO BORDER DISTANCE
            // THIS IS BECAUSE WE HAVEN'T CALCULATED LEAF VERTEX TO BORDER DISTANCES
            this->sourceToTreeNodeBorderDist[sourceLeafIdx][i] = sourceToNextBorderDist;
            if (sourceToNextBorderDist < spDist) {
                spDist = sourceToNextBorderDist;
            }
        }
    }
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[sourceLeafIdx].bordersVec.size());
#endif
    return spDist;
}

EdgeWeight Gtree::SPDistToParentNode(int childTreeIdx, int parentTreeIdx, bool computeSPDist)
{
    EdgeWeight spDist = 0, sourceToNextBorderDist, sourceToChildBorderDist, sourceToParentBorderDist;
    
    // Note: We assume we have already calculated distance to childTreeIdx from source
    
    int childPos, parentBorderIdx, childBorderOffset, childBorderIdx;

    childPos = this->treeNodes[childTreeIdx].getParentChildIdx();
    childBorderOffset = this->treeNodes[parentTreeIdx].getChildOffsetInChildBorderVec(childPos);
    
    if (this->treeNodes[childTreeIdx].bordersVec.size() > 0) {
        childBorderIdx = childBorderOffset;
        sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[childTreeIdx][0];
        for (std::size_t k = 0; k < this->treeNodes[parentTreeIdx].bordersVec.size(); ++k) {
            parentBorderIdx = this->treeNodes[parentTreeIdx].getBorderIdxInChildBorderVec(k);
            sourceToParentBorderDist = sourceToChildBorderDist + this->treeNodes[parentTreeIdx].distanceMatrix[childBorderIdx*this->treeNodes[parentTreeIdx].matrixRowLength+parentBorderIdx];
            // I.e. this is the first time we are computing distances to parent's border set so we initialise
            this->sourceToTreeNodeBorderDist[parentTreeIdx][k] = sourceToParentBorderDist;
        }
        for (std::size_t j = 1; j < this->treeNodes[childTreeIdx].bordersVec.size(); ++j) {
            childBorderIdx = childBorderOffset + j;
            sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[childTreeIdx][j];
            for (std::size_t k = 0; k < this->treeNodes[parentTreeIdx].bordersVec.size(); ++k) {
                parentBorderIdx = this->treeNodes[parentTreeIdx].getBorderIdxInChildBorderVec(k);
                sourceToParentBorderDist = sourceToChildBorderDist + this->treeNodes[parentTreeIdx].distanceMatrix[childBorderIdx*this->treeNodes[parentTreeIdx].matrixRowLength+parentBorderIdx];
                if (sourceToParentBorderDist < this->sourceToTreeNodeBorderDist[parentTreeIdx][k]) {
                    this->sourceToTreeNodeBorderDist[parentTreeIdx][k] = sourceToParentBorderDist;
                }
            }
        }
    }
    
    if (this->sourceToTreeNodeBorderDist[parentTreeIdx].size() > 0 && computeSPDist) {
        spDist = this->sourceToTreeNodeBorderDist[parentTreeIdx][0];
        for (std::size_t i = 1; i < this->sourceToTreeNodeBorderDist[parentTreeIdx].size(); ++i) {
            sourceToNextBorderDist = this->sourceToTreeNodeBorderDist[parentTreeIdx][i];
            if (sourceToNextBorderDist < spDist) {
                spDist = sourceToNextBorderDist;
            }
        }
    }

#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[childTreeIdx].bordersVec.size()*this->treeNodes[parentTreeIdx].bordersVec.size());
#endif
    return spDist;
}

EdgeWeight Gtree::SPDistToSiblingNode(int firstLCAChildIdx, int targetLCAChildIdx, int LCAIdx, bool computeSPDist)
{
    EdgeWeight spDist = 0, sourceToNextBorderDist, sourceToChildBorderDist, sourceToBorderDist;

    int childPos, childBorderOffset, childBorderIdx;
    int targetChildPos, targetBorderOffset, targetBorderIdx;
    
    childPos = this->treeNodes[firstLCAChildIdx].getParentChildIdx();
    childBorderOffset = this->treeNodes[LCAIdx].getChildOffsetInChildBorderVec(childPos);
    targetChildPos = this->treeNodes[targetLCAChildIdx].getParentChildIdx();
    targetBorderOffset = this->treeNodes[LCAIdx].getChildOffsetInChildBorderVec(targetChildPos);
    
    if (this->treeNodes[firstLCAChildIdx].bordersVec.size() > 0) {
        childBorderIdx = childBorderOffset;
        sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[firstLCAChildIdx][0];
        for (std::size_t k = 0; k < this->treeNodes[targetLCAChildIdx].bordersVec.size(); ++k) {
            targetBorderIdx = targetBorderOffset + k;
            sourceToBorderDist = sourceToChildBorderDist + this->treeNodes[LCAIdx].distanceMatrix[childBorderIdx*this->treeNodes[LCAIdx].matrixRowLength+targetBorderIdx];
            // I.e. this is the first time we are computing distances to tthe target set
            this->sourceToTreeNodeBorderDist[targetLCAChildIdx][k] = sourceToBorderDist;
        }
        for (std::size_t j = 1; j < this->treeNodes[firstLCAChildIdx].bordersVec.size(); ++j) {
            childBorderIdx = childBorderOffset + j;
            sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[firstLCAChildIdx][j];
            for (std::size_t k = 0; k < this->treeNodes[targetLCAChildIdx].bordersVec.size(); ++k) {
                targetBorderIdx = targetBorderOffset + k;
                sourceToBorderDist = sourceToChildBorderDist + this->treeNodes[LCAIdx].distanceMatrix[childBorderIdx*this->treeNodes[LCAIdx].matrixRowLength+targetBorderIdx];
                if (sourceToBorderDist < this->sourceToTreeNodeBorderDist[targetLCAChildIdx][k]) {
                    this->sourceToTreeNodeBorderDist[targetLCAChildIdx][k] = sourceToBorderDist;
                }
            }
        }
    }

    if (this->sourceToTreeNodeBorderDist[targetLCAChildIdx].size() > 0 && computeSPDist) {
        spDist = this->sourceToTreeNodeBorderDist[targetLCAChildIdx][0];
        for (std::size_t i = 1; i < this->sourceToTreeNodeBorderDist[targetLCAChildIdx].size(); ++i) {
            sourceToNextBorderDist = this->sourceToTreeNodeBorderDist[targetLCAChildIdx][i];
            if (sourceToNextBorderDist < spDist) {
                spDist = sourceToNextBorderDist;
            }
        }
    }
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[firstLCAChildIdx].bordersVec.size()*this->treeNodes[targetLCAChildIdx].bordersVec.size());
#endif

    return spDist;
}

EdgeWeight Gtree::SPDistToChildNode(int childTreeIdx, int parentTreeIdx, bool computeSPDist)
{
    EdgeWeight spDist = 0, sourceToNextBorderDist, sourceToChildBorderDist, sourceToBorderDist;

    int childPos, parentBorderIdx, childBorderOffset, childBorderIdx;
    childPos = this->treeNodes[childTreeIdx].getParentChildIdx();
    childBorderOffset = this->treeNodes[parentTreeIdx].getChildOffsetInChildBorderVec(childPos);

    if (this->treeNodes[parentTreeIdx].bordersVec.size() > 0) {
        parentBorderIdx = this->treeNodes[parentTreeIdx].getBorderIdxInChildBorderVec(0);
        sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[parentTreeIdx][0];
        for (std::size_t k = 0; k < this->treeNodes[childTreeIdx].bordersVec.size(); ++k) {
            childBorderIdx = childBorderOffset + k;
            sourceToBorderDist = sourceToChildBorderDist + this->treeNodes[parentTreeIdx].distanceMatrix[parentBorderIdx*this->treeNodes[parentTreeIdx].matrixRowLength+childBorderIdx];
            // I.e. this is the first time we are computing distances to tthe target set
            this->sourceToTreeNodeBorderDist[childTreeIdx][k] = sourceToBorderDist;
        }
        for (std::size_t j = 1; j < this->treeNodes[parentTreeIdx].bordersVec.size(); ++j) {
            parentBorderIdx = this->treeNodes[parentTreeIdx].getBorderIdxInChildBorderVec(j);
            sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[parentTreeIdx][j];
            for (std::size_t k = 0; k < this->treeNodes[childTreeIdx].bordersVec.size(); ++k) {
                childBorderIdx = childBorderOffset + k;
                sourceToBorderDist = sourceToChildBorderDist + this->treeNodes[parentTreeIdx].distanceMatrix[parentBorderIdx*this->treeNodes[parentTreeIdx].matrixRowLength+childBorderIdx];
                if (sourceToBorderDist < this->sourceToTreeNodeBorderDist[childTreeIdx][k]) {
                    this->sourceToTreeNodeBorderDist[childTreeIdx][k] = sourceToBorderDist;
                }
            }
        }
    }

    if (this->sourceToTreeNodeBorderDist[childTreeIdx].size() > 0 && computeSPDist) {
        spDist = this->sourceToTreeNodeBorderDist[childTreeIdx][0];
        for (std::size_t i = 1; i < this->sourceToTreeNodeBorderDist[childTreeIdx].size(); ++i) {
            sourceToNextBorderDist = this->sourceToTreeNodeBorderDist[childTreeIdx][i];
            if (sourceToNextBorderDist < spDist) {
                spDist = sourceToNextBorderDist;
            }
        }
    }
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[parentTreeIdx].bordersVec.size()*this->treeNodes[childTreeIdx].bordersVec.size());
#endif
    
    return spDist;
    
}

EdgeWeight Gtree::SPDistToLeafTarget(VertexID target, int leafIdx)
{
    EdgeWeight spDist = 0, candidateSourceToTargetDist;
    
    // Now find the border which has the shortest distance to the target in the leaf node
    int vIdx = this->getIdxInLeafVerticesVec(target);

    if (this->treeNodes[leafIdx].bordersVec.size() > 0) {
        spDist = this->sourceToTreeNodeBorderDist[leafIdx][0] + this->treeNodes[leafIdx].distanceMatrix[vIdx];
        for (std::size_t i = 1; i < this->treeNodes[leafIdx].bordersVec.size(); ++i) {
            candidateSourceToTargetDist = this->sourceToTreeNodeBorderDist[leafIdx][i] + this->treeNodes[leafIdx].distanceMatrix[i*this->treeNodes[leafIdx].matrixRowLength+vIdx];
            if (candidateSourceToTargetDist < spDist) {
                spDist = candidateSourceToTargetDist;
            }
        }
    }
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[leafIdx].bordersVec.size());
#endif
    
    return spDist;
}

int Gtree::getIdxInLeafVerticesVec(VertexID u)
{
    return this->vertexIDLeafVerticesVecIdx[u];
}

void Gtree::setLeafVerticesVecIdx(VertexID u, int leafVerticesVecIdx)
{
    this->vertexIDLeafVerticesVecIdx[u] = leafVerticesVecIdx;
}

int Gtree::getIdxInLeafBordersVec(VertexID u)
{
    return this->vertexIDLeafBordersVecIdx[u];
}

void Gtree::setLeafBordersVecIdx(VertexID u, int leafBordersVecIdx)
{
    this->vertexIDLeafBordersVecIdx[u] = leafBordersVecIdx;
}

bool Gtree::isEdgeInLeafSubgraph(VertexID edge)
{
    return this->edgeInLeafSubgraph[edge];
}

void Gtree::setEdgeNotInSubgraph(VertexID edge)
{
    this->edgeInLeafSubgraph[edge] = false;
}

double Gtree::computeIndexSize()
{
    double memoryUsage = 0;
    for (std::size_t i = 0; i < this->treeNodes.size(); ++i) {
        memoryUsage += this->treeNodes[i].computeIndexSizeBytes();
    }
    memoryUsage += sizeof(int)*this->leafIdxs.size();
    memoryUsage += sizeof(std::vector<EdgeWeight>)*this->sourceToTreeNodeBorderDist.size();
    for (std::size_t i = 0; i < sourceToTreeNodeBorderDist.size(); ++i) {
        memoryUsage += sizeof(EdgeWeight)*this->sourceToTreeNodeBorderDist[i].size();
    }
    memoryUsage += sizeof(int)*this->vertexIDLeafVerticesVecIdx.size();
    memoryUsage += sizeof(int)*this->vertexIDLeafBordersVecIdx.size();
    memoryUsage += (this->edgeInLeafSubgraph.size()+7)/8; // std::vector<bool> uses 1 bit per element
    memoryUsage += sizeof(int)*this->vertexLeafIdxs.size();
    return memoryUsage/(1024*1024);
}

EdgeWeight Gtree::getRepeatedShortestPathDistance(Graph& graph, VertexID u, VertexID v, std::vector<bool>& visited, std::vector<EdgeWeight>& leafVertexDistances, 
                                                             BinaryMinHeap<EdgeWeight,VertexID>& pqueue, std::unordered_set<VertexID>& leafVertexVisited)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("computations_materialized",0);
    this->stats.initialiseStatistic("computations_executed",0);
    this->stats.initialiseStatistic("computations_total",0);
#endif
    EdgeWeight spDist = 0;
    int uLeaf = this->getLeafIndex(u);
    int vLeaf = this->getLeafIndex(v);
    
    if (uLeaf == vLeaf) {
        // This mean they are both in the same node
        if (leafVertexDistances.size() == 0) {
            // This means repeated leaf search has not be initialised yet
            leafVertexDistances.resize(this->treeNodes[uLeaf].leafVerticesVec.size(),0);
            pqueue.insert(u,0); // Insert source vertex into queue for leaf search
        }
        spDist = this->getRepeatedSourceLeafShortestPathDistance(graph,uLeaf,v,leafVertexDistances,pqueue,leafVertexVisited);
    } else {
        int LCAIdx = 0;
        std::vector<int>& uPathFromRoot = this->treeNodes[uLeaf].getGtreePathFromRoot();
        std::vector<int>& vPathFromRoot = this->treeNodes[vLeaf].getGtreePathFromRoot();

        // Since the two vertices are in different leaf nodes we can guarantee that
        // there is at least one more node (the parent of both) in the G-tree path
        // which we call the the lowest common ancestor (LCA)
        
        // Search down Gtree from root until we find first ancestor that is different
        // this means the previous ancestor is the LCA and the indexes point its children
        unsigned int i, j;
        for (i = 0, j = 0; i < uPathFromRoot.size() && j < vPathFromRoot.size(); ++i, ++j) {
            if (uPathFromRoot[i] != vPathFromRoot[j]) {
                // When i = 0 and j = 0 it is referring to the root so in that
                // case uPathFromRoot[i] does equal vPathFromRoot[j]. This means
                // when they are equal i > 0, so i-1 is safe here
                LCAIdx = uPathFromRoot[i-1];
                break;
            }
            // Note: We can guarantee that LCAIdx will be set here. The only situation
            // it would not be set is both u and v were in the same leaf but we have
            // guarnateed that is not the case in the if/else statement
        }
        
        // Now search up G-tree (from source leaf) until we reach i and j, then we 
        // search down (to target leaf) computing shortest path distance
        
        // From source to source leaf
        if (!visited[uLeaf]) {
            this->SPDistToSourceLeafNode(u,uLeaf);
            visited[uLeaf] = true;
        }
        
        // From source leaf to first child of LCA
        int x = i; // This is safe, depth of tree is O(log(n)) and n is at most 24 million in US dataset
        for (int k = uPathFromRoot.size()-1; k > x; --k) {
            // Since k > x and x is at worst 0, k-1 is safe here
            if (!visited[uPathFromRoot[k-1]]) {
                this->SPDistToParentNode(uPathFromRoot[k],uPathFromRoot[k-1],false);
                visited[uPathFromRoot[k-1]] = true;
            }
        }
        
        // From first child of LCA to second child of LCA
        if (!visited[vPathFromRoot[j]]) {
            this->SPDistToSiblingNode(uPathFromRoot[i],vPathFromRoot[j],LCAIdx,false);
            visited[vPathFromRoot[j]] = true;
        }

        // From second child of LCA to target leaf
        for (std::size_t k = j; k < vPathFromRoot.size()-1; ++k) {
            // Note the size()-1 in the above condition
            if (!visited[vPathFromRoot[k+1]]) {
                this->SPDistToChildNode(vPathFromRoot[k+1],vPathFromRoot[k],false);
                visited[vPathFromRoot[k+1]] = true;
            }
        }
        
        // We assume target has not been visited
        spDist = this->SPDistToLeafTarget(v,vLeaf);
    }
    
    return spDist;
}

EdgeWeight Gtree::getRepeatedSourceLeafShortestPathDistance(Graph& graph, int leafNode, VertexID v, std::vector<EdgeWeight>& leafVertexDistances, 
                                                            BinaryMinHeap<EdgeWeight,VertexID>& pqueue, std::unordered_set<VertexID>& leafVertexVisited)
{
    // Assume vectors have been resized for number of vertices in source leaf
    // and query vertex has been inserted into pqueue
    
    if (leafVertexVisited.find(v) != leafVertexVisited.end()) {
        // If it's been already visited in the Dijkstra's search
        // we can return the distance immediately
        return leafVertexDistances[this->vertexIDLeafVerticesVecIdx[v]];
    }
    
    EdgeWeight minDist;
    VertexID minDistVertexID, adjVertex;
    int adjListStart, nextAdjListStart;
    
    // We also exist if all target vertices have been found
    while (pqueue.size() > 0) {
        // Extract and remove vertex with smallest distance from query point
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistVertexID = pqueue.extractMinElement();
        if (leafVertexVisited.find(minDistVertexID) == leafVertexVisited.end()) {
            if (minDistVertexID == v) {
                // Re-insert it back into queue so that we can dequeue
                // and add neighbours in repeated search
                pqueue.insert(minDistVertexID,minDist);
                return minDist;
            }

            leafVertexVisited.insert(minDistVertexID);
            leafVertexDistances[this->vertexIDLeafVerticesVecIdx[minDistVertexID]] = minDist;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistVertexID);
            nextAdjListStart = graph.getEdgeListSize(minDistVertexID);
            // Note: We have to make sure we don't exceed size of graph.edges vector
            
            bool isBorder = false;
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjVertex = graph.edges[i].first;
                if (this->isEdgeInLeafSubgraph(i)) {
                    if (leafVertexVisited.find(adjVertex) == leafVertexVisited.end()) {
                        // This is the first time we have seen this vertex (i.e. distance infinity)
                        pqueue.insert(adjVertex,minDist+graph.edges[i].second);
                    }
                } else {
                    if (!isBorder) {
                        isBorder = true;
                    }
                }
            }
            
            // If it is a border we add each of the other borders for this leaf
            // node and the distance from this border (assuming they have not 
            // already been settled - if they are settled we have already found
            // shortest path to them so the path through this border cannot
            // be an improvement)
            if (isBorder) {
                int borderIdx = this->getIdxInLeafBordersVec(minDistVertexID);
                for (std::size_t i = 0; i < this->treeNodes[leafNode].bordersVec.size(); ++i) {
                    // Note since this is leaf the above functions actual gets idx in leafVerticesVec
                    // Recall that the dimension of the distance matrix is bordersVec.size()*leafVerticesVec.size()
                    // and stores the distances from borders to leaf vertices so we can use these values
                    // compute the distances between two borders in the leaf vertice
                    if (leafVertexVisited.find(this->treeNodes[leafNode].bordersVec[i]) == leafVertexVisited.end()) {
                        int targetBorderIdx = this->treeNodes[leafNode].getBorderIdxInChildBorderVec(i);
                        pqueue.insert(this->treeNodes[leafNode].bordersVec[i],minDist+this->treeNodes[leafNode].distanceMatrix[borderIdx*this->treeNodes[leafNode].matrixRowLength+targetBorderIdx]);
                    }
                }
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("computations_executed",this->treeNodes[leafNode].bordersVec.size());
#endif
            }
            
        }
    }
    std::cerr << "INTERNAL ERROR: vertex " << v << " was not found in its G-tree leaf " << leafNode << std::endl;
    std::exit(1);
}

int Gtree::getComputations(int leafIdx, int targetIdx)
{
    int firstLCAChild;
    std::vector<int> gtreePath = this->getGtreePath(leafIdx,targetIdx,firstLCAChild);
    
    int numComputations = 0;
    std::size_t x = 0;
    
    if (gtreePath.size() > 0) {
        RN_DEBUG_ASSERT(gtreePath.size() >= 2 && "Invalid path - must contain at least two nodes e.g. two leaves or leaf and parent");
        numComputations += this->treeNodes[gtreePath[0]].getNumBorders();
        
        // Add computation up the path until we reach firstLCAChild
        for (x = 1; x < gtreePath.size() && gtreePath[x-1] != firstLCAChild; ++x) {
            numComputations += this->treeNodes[gtreePath[x-1]].bordersVec.size()*this->treeNodes[gtreePath[x]].bordersVec.size();
        }
        
        // Add computation between the two children of the LCA
        // Note: If the parent of leafIdx is the firstLCAChild then this loop does nothing
        // but x is now smaller than the gtreePath (because targetIdx is in the same path
        // as the leafIdx) then the following code does nothing
        if (x < gtreePath.size()) {
            numComputations += this->treeNodes[gtreePath[x-1]].bordersVec.size()*this->treeNodes[gtreePath[x]].bordersVec.size();
            ++x;
        }
        
        // Add computation down the path until we reach the end
        for (; x < gtreePath.size(); ++x) {
            numComputations += this->treeNodes[gtreePath[x-1]].bordersVec.size()*this->treeNodes[gtreePath[x]].bordersVec.size();
        }
    }

    return numComputations;
}

/*
 * GtreeNode
 */

GtreeNode::GtreeNode(int parentIdx, int numVertices): 
    parentIdx(parentIdx), numVertices(numVertices), isLeaf(false) {}

GtreeNode::GtreeNode() {}
    
bool GtreeNode::isLeafNode()
{
    return this->isLeaf;
}

void GtreeNode::setLeafNode()
{
    this->isLeaf = true;
}

int GtreeNode::getParentIdx()
{
    return this->parentIdx;
}

int GtreeNode::getNumVertices()
{
    return this->numVertices;
}

int GtreeNode::getParentChildIdx()
{
    return this->parentChildIdx;
}

void GtreeNode::setParentChildIdx(int idx)
{
    this->parentChildIdx = idx;
}

int GtreeNode::addChild(int treeIdx)
{
    this->children.push_back(treeIdx);
    return this->children.size()-1;
}

int GtreeNode::addBorder(VertexID vertex)
{
    // Borders are only ever added to G-tree node once
    this->bordersVec.push_back(vertex);
    int bordersVecIdx = this->bordersVec.size()-1;
    if (!this->isLeaf) {
        // If this is not a leaf, it's childBordersVec will be full
        // so we want to remember the array index for this border in
        // childBordersVec
        this->borderOffsetsInChildBorderVec.push_back(childBorderToChildBorderVecIdx[vertex]);
    } else {
        this->borderOffsetsInChildBorderVec.push_back(leafVerticeToLeafVerticeVecIdx[vertex]);
    }
    return bordersVecIdx;
}

std::vector<VertexID>& GtreeNode::getBorders()
{
    return this->bordersVec;
}

int GtreeNode::getNumBorders()
{
    return this->bordersVec.size();
}

int GtreeNode::addLeafVertex(VertexID vertex)
{
    this->leafVerticesUset.insert(vertex);
    this->leafVerticesVec.push_back(vertex);
    int leafVerticesVecIdx = this->leafVerticesVec.size()-1;
    this->leafVerticeToLeafVerticeVecIdx[vertex] = leafVerticesVecIdx;
    return leafVerticesVecIdx;
}

std::vector<VertexID>& GtreeNode::getLeafVertices()
{
    return this->leafVerticesVec;
}

std::unordered_set<VertexID >& GtreeNode::getLeafVerticesUset()
{
    return this->leafVerticesUset;
}

void GtreeNode::addChildBorder(VertexID vertex)
{
    // Child borders are only ever added once as they only belong to one child
    this->childBordersUset.insert(vertex);
    this->childBordersVec.push_back(vertex);
    this->childBorderToChildBorderVecIdx[vertex] = this->childBordersVec.size()-1;
}

std::vector<VertexID>& GtreeNode::getChildBorders()
{
    return this->childBordersVec;
}

std::unordered_set<VertexID>& GtreeNode::getChildBordersUset()
{
    return this->childBordersUset;
}

void GtreeNode::addGtreePathNode(int nodeIdx)
{
    this->gtreePath.push_back(nodeIdx);
}

void GtreeNode::addGtreePathNodes(const std::vector<int>& parentPath)
{
    this->gtreePath.insert(this->gtreePath.end(),parentPath.begin(),parentPath.end());
}

std::vector<int>& GtreeNode::getGtreePathFromRoot()
{
    return this->gtreePath;
}

int GtreeNode::getBorderIdxInChildBorderVec(int borderIdx)
{
    return this->borderOffsetsInChildBorderVec[borderIdx];
}

int GtreeNode::getChildOffsetInChildBorderVec(int childIdx)
{
    return this->childOffsetsInChildBorderVec[childIdx];
}

void GtreeNode::addChildOffsetInChildBorderVec()
{
    this->childOffsetsInChildBorderVec.push_back(this->childBordersVec.size());
}

double GtreeNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(int)*4;
    memoryUsage += sizeof(bool);
    memoryUsage += sizeof(this->children);
    memoryUsage += sizeof(this->bordersVec);
    memoryUsage += sizeof(this->leafVerticesVec);
    memoryUsage += sizeof(this->gtreePath);
    memoryUsage += sizeof(this->distanceMatrix);
    memoryUsage += sizeof(this->borderOffsetsInChildBorderVec);
    memoryUsage += sizeof(this->childOffsetsInChildBorderVec);
    memoryUsage += sizeof(int)*this->children.size();
    memoryUsage += sizeof(VertexID)*this->bordersVec.size();
    memoryUsage += sizeof(VertexID)*this->leafVerticesVec.size();
    memoryUsage += sizeof(int)*this->gtreePath.size();
    memoryUsage += sizeof(EdgeWeight)*this->distanceMatrix.size();
    memoryUsage += sizeof(int)*this->borderOffsetsInChildBorderVec.size();
    memoryUsage += sizeof(int)*this->childOffsetsInChildBorderVec.size();
    return memoryUsage;
}

/*
 * OccurenceList
 */

OccurenceList::OccurenceList(std::string setType, double setDensity, int setVariable, int setSize):
                             objSetType(setType), objSetDensity(setDensity), objSetVariable(setVariable), objSetSize(setSize)
{}

OccurenceList::OccurenceList() {}

void OccurenceList::addLeafOccurence(int leafIdx, VertexID object)
{
    if (this->leafOccurenceSet.find(leafIdx) == this->leafOccurenceSet.end()) {
        std::unordered_set<VertexID> occurenceSet;
        occurenceSet.insert(object);
        std::vector<VertexID> occurenceList;
        occurenceList.push_back(object);
        this->leafOccurenceSet[leafIdx] = occurenceSet;
        this->leafOccurenceList[leafIdx] = occurenceList;
    } else {
        if (this->leafOccurenceSet[leafIdx].find(object) == this->leafOccurenceSet[leafIdx].end()) {
            this->leafOccurenceSet[leafIdx].insert(object);
            this->leafOccurenceList[leafIdx].push_back(object);
        }
    }
}

void OccurenceList::addParentOccurence(int parentIdx, int childIdx)
{
    if (this->nonLeafOccurenceSet.find(parentIdx) == this->nonLeafOccurenceSet.end()) {
        std::unordered_set<int> occurenceSet;
        occurenceSet.insert(childIdx);
        this->nonLeafOccurenceSet[parentIdx] = occurenceSet;
        std::vector<int> occurenceList;
        occurenceList.push_back(childIdx);
        this->nonLeafOccurenceList[parentIdx] = occurenceList;
    } else {
        if (this->nonLeafOccurenceSet[parentIdx].find(childIdx) == this->nonLeafOccurenceSet[parentIdx].end()) {
            this->nonLeafOccurenceSet[parentIdx].insert(childIdx);
            this->nonLeafOccurenceList[parentIdx].push_back(childIdx);
        }
    }
}

std::string OccurenceList::getObjSetType()
{
    return this->objSetType;
}

double OccurenceList::getObjSetDensity()
{
    return this->objSetDensity;
}

int OccurenceList::getObjSetSize()
{
    return this->objSetSize;
}

int OccurenceList::getObjSetVariable()
{
    return this->objSetVariable;
}

double OccurenceList::computeIndexSize()
{
    double memoryUsage = 0;
    for (auto it = this->leafOccurenceList.begin(); it != this->leafOccurenceList.end(); ++it) {
        memoryUsage += sizeof(it->first)+sizeof(it->second)+sizeof(VertexID)*it->second.size();
    }
    for (auto it = this->nonLeafOccurenceList.begin(); it != this->nonLeafOccurenceList.end(); ++it) {
        memoryUsage += sizeof(it->first)+sizeof(it->second)+sizeof(int)*it->second.size();
    }
    return memoryUsage/(1024*1024);
    // Note: We don't include leafOccurenceSet even though it is serialized, because
    // it contains the same data as leafOccurenceList and can easily reconstructed
    // from that instead.
}
