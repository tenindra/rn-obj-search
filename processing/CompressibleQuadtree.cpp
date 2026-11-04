// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "CompressibleQuadtree.h"

#include "../utility/geometry.h"
#include "../utility/utility.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <deque>

CompressibleQuadtree::CompressibleQuadtree(int maxLeafColours): maxLeafColours(maxLeafColours) {}

CompressibleQuadtree::~CompressibleQuadtree()
{
    // Recursively delete all quadtree nodes
    delete rootNode;
}

void CompressibleQuadtree::buildFromGraph(Graph& graph, std::vector<MortonCode>& mortonCodes, std::vector<VertexID>& colourMap)
{
    this->rootNode = new CompressibleQuadtreeNode(true,true);
    this->firstLeaf = this->rootNode;
    this->lastLeaf = this->rootNode;
    
    VertexID numVertices = graph.getNumVertices();
    for (VertexID vertex = 0; vertex < numVertices; ++vertex) {
        this->insert(vertex,mortonCodes,colourMap);
    }
}

void CompressibleQuadtree::insert(VertexID vertex, std::vector<MortonCode>& mortonCodes, std::vector<VertexID>& colourMap)
{
    CompressibleQuadtreeNode* currentNode = this->rootNode;
    
    int index = vertex;

    // Find leaf node in which this vertex is contained
    while(!currentNode->isLeafNode()) {
        // Determine which of current (non-leaf) nodes children contains vertex
        if (currentNode->sw->contains(mortonCodes[index])) {
            currentNode = currentNode->sw;
        } else if (currentNode->se->contains(mortonCodes[index])) {
            currentNode = currentNode->se;
        } else if (currentNode->nw->contains(mortonCodes[index])) {
            currentNode = currentNode->nw;
        } else if (currentNode->ne->contains(mortonCodes[index])) {
            currentNode = currentNode->ne;
        } else {
            std::cerr << "No child quadrant of current node contains point" << std::endl;
            exit(1);
        }
    }
    
    currentNode->addVertex(vertex,colourMap[vertex]);
    if (currentNode->getNumColours() > this->maxLeafColours) {
        this->split(currentNode,mortonCodes,colourMap);
    }
}

void CompressibleQuadtree::split(CompressibleQuadtreeNode* treeNode, std::vector<MortonCode>& mortonCodes, std::vector<VertexID>& colourMap)
{
    // Mark node is being no longer a leaf node
    treeNode->setLeaf(false);
    
    // Create and update tree node pointers for split
    treeNode->sw = new CompressibleQuadtreeNode(treeNode->mortonCode,SW,treeNode->prevLeaf);
    treeNode->se = new CompressibleQuadtreeNode(treeNode->mortonCode,SE,treeNode->sw);
    treeNode->nw = new CompressibleQuadtreeNode(treeNode->mortonCode,NW,treeNode->se);
    treeNode->ne = new CompressibleQuadtreeNode(treeNode->mortonCode,NE,treeNode->nw,treeNode->nextLeaf);
    
    // Note: We don't know the next leaf nodes until all of them have created
    treeNode->sw->setNextLeaf(treeNode->se);
    treeNode->se->setNextLeaf(treeNode->nw);
    treeNode->nw->setNextLeaf(treeNode->ne);
    
    if (treeNode->prevLeaf != NULL) {
        treeNode->prevLeaf->nextLeaf = treeNode->sw;
        treeNode->prevLeaf = NULL;
    }
    if (treeNode->nextLeaf != NULL) {
        treeNode->nextLeaf->prevLeaf = treeNode->ne;
        treeNode->nextLeaf = NULL;
    }
    
    treeNode->pushDownVertices(mortonCodes,colourMap);

    // If this node was the first or last leaf 
    // we need update these with the child nodes
    if (treeNode == this->firstLeaf) {
        this->firstLeaf = treeNode->sw;
    }
    if (treeNode == this->lastLeaf) {
        this->lastLeaf = treeNode->ne;
    }
    
    // It is possible that all of items were placed in a single quadrant
    // so we must check and recursively split those as well
    if (treeNode->sw->getNumColours() > this->maxLeafColours) {
        this->split(treeNode->sw,mortonCodes,colourMap);
    }
    if (treeNode->se->getNumColours() > this->maxLeafColours) {
        this->split(treeNode->se,mortonCodes,colourMap);
    }
    if (treeNode->nw->getNumColours() > this->maxLeafColours) {
        this->split(treeNode->nw,mortonCodes,colourMap);
    }
    if (treeNode->ne->getNumColours() > this->maxLeafColours) {
        this->split(treeNode->ne,mortonCodes,colourMap);
    }
}

int CompressibleQuadtree::countLeafNodes(int& totalVertices, int& nonEmptyNodes) {
    totalVertices = 0;
    nonEmptyNodes = 0;
    int count = 0;
    CompressibleQuadtreeNode* currentNode = this->firstLeaf;
    while (currentNode != NULL) {
        totalVertices += currentNode->getNumVertices();
        if (currentNode->getNumVertices() > 0) {
            ++nonEmptyNodes;
        }
        currentNode = currentNode->nextLeaf;
        ++count;
    }
    totalVertices = 0;
    nonEmptyNodes = 0;
    count = 0;
    currentNode = this->lastLeaf;
    while (currentNode != NULL) {
        totalVertices += currentNode->getNumVertices();
        if (currentNode->getNumVertices() > 0) {
            ++nonEmptyNodes;
        }
        currentNode = currentNode->prevLeaf;
        ++count;
    }
    return count;
}

double CompressibleQuadtree::computeIndexSize()
{
    double memoryUsage = 0;
    std::deque<CompressibleQuadtreeNode*> stack;
    stack.push_back(this->rootNode);
    while (stack.size() > 0) {
        CompressibleQuadtreeNode* currNode = stack.front();
        stack.pop_front();
        memoryUsage += currNode->computeIndexSizeBytes();
        if (!currNode->isLeafNode()) {
            stack.push_back(currNode->sw);
            stack.push_back(currNode->se);
            stack.push_back(currNode->nw);
            stack.push_back(currNode->ne);
        }
    }
    return memoryUsage/(1024*1024);
}

/*
 * CompressibleQuadtreeNode
 */

CompressibleQuadtreeNode::CompressibleQuadtreeNode(MortonCode& parentCode, Direction dir, CompressibleQuadtreeNode* prevLeaf /* = NULL */, CompressibleQuadtreeNode* nextLeaf /* = NULL */)
{
    this->mortonCode.setCode(parentCode,dir);
    this->isRoot = false; // Cannot be root because it has a parent
    this->isLeaf = true;
    this->sw = this->se = this->nw = this->ne = NULL;
    this->prevLeaf = prevLeaf;
    this->nextLeaf = nextLeaf;
    this->numVertices = 0;
    this->numColours = 0;
}

CompressibleQuadtreeNode::CompressibleQuadtreeNode(bool isRoot, bool isLeaf)
{
    this->isRoot = isRoot;
    this->isLeaf = isLeaf;
    this->sw = this->se = this->nw = this->ne = NULL;
    this->prevLeaf = NULL;
    this->nextLeaf = NULL;
    this->numVertices = 0;
    this->numColours = 0;
}

CompressibleQuadtreeNode::~CompressibleQuadtreeNode() {
    // Recursively delete all quadtree nodes
    delete this->sw;
    delete this->se;
    delete this->nw;
    delete this->ne;
}

void CompressibleQuadtreeNode::setNextLeaf(CompressibleQuadtreeNode *nextLeaf) {
    this->nextLeaf = nextLeaf;
}

bool CompressibleQuadtreeNode::isRootNode() {
    return this->isRoot;
}

int CompressibleQuadtreeNode::getNumVertices()
{
    return this->numVertices;
}

int CompressibleQuadtreeNode::getNumColours()
{
    return this->numColours;
}

bool CompressibleQuadtreeNode::isLeafNode() {
    return this->isLeaf;
}

void CompressibleQuadtreeNode::setLeaf(bool isLeaf) {
    this->isLeaf = isLeaf;
}

void CompressibleQuadtreeNode::addVertex(VertexID vertex, VertexID colour) {
    ++this->numVertices;
    this->vertices.push_back(vertex);
    if (this->colourIDs.find(colour) == this->colourIDs.end()) {
        ++this->numColours;
        this->colourIDs.insert(colour);
    }
}

void CompressibleQuadtreeNode::pushDownVertices(std::vector<MortonCode>& allMortonCodes, std::vector<VertexID>& colourMap) {
    for (std::size_t i = 0; i < this->vertices.size(); ++i) {
        this->pushDownVertex(this->vertices[i],allMortonCodes[this->vertices[i]],colourMap[this->vertices[i]]);
    }
    this->vertices.clear();
}

void CompressibleQuadtreeNode::pushDownVertex(VertexID vertex, MortonCode& mortonCode, VertexID colour) {
    if (this->sw->contains(mortonCode)) {
        this->sw->addVertex(vertex,colour);
    } else if (this->se->contains(mortonCode)) {
        this->se->addVertex(vertex,colour);
    } else if (this->nw->contains(mortonCode)) {
        this->nw->addVertex(vertex,colour);
    } else if (this->ne->contains(mortonCode)) {
        this->ne->addVertex(vertex,colour);
    } else {
        std::cerr << "Cannot push down, no child quadrant contains point" << std::endl;
        exit(1);
    }
}

bool CompressibleQuadtreeNode::contains(const MortonCode& point) {
    return this->mortonCode.contains(point);
}

double CompressibleQuadtreeNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(CompressibleQuadtreeNode*)*6;
    memoryUsage += sizeof(bool)*2;
    memoryUsage += sizeof(int)*2;
    memoryUsage += sizeof(this->vertices);
    memoryUsage += this->mortonCode.computeIndexSizeBytes();
    memoryUsage += sizeof(VertexID)*this->vertices.size();
    // The colour set counts as its object and its colours, like vertices; buckets and allocator overhead are left out, as in the other index sizes
    memoryUsage += sizeof(this->colourIDs);
    memoryUsage += sizeof(VertexID)*this->colourIDs.size();
    return memoryUsage;
}