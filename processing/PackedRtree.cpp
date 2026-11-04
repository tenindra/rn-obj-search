// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "PackedRtree.h"
#include "../utility/geometry.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <deque>

PackedRtree::PackedRtree(int branchFactor, std::string setType, double setDensity, int setVariable, int setSize): 
    branchFactor(branchFactor), objSetType(setType), objSetDensity(setDensity), objSetVariable(setVariable), objSetSize(setSize)
{}

void PackedRtree::bulkLoad(std::vector<VertexID>& objectIDs, std::vector<CoordinatePair>& objectCoords)
{
    // The algorithm used to bulkload object into the R-tree is based on:
    // Leutenegger, Scott T., Mario Lopez, and Jeffrey Edgington
    // STR: A simple and efficient algorithm for R-tree packing
    // Proceedings of the 13th IEEE International Conference on Data
    // Engineering (ICDE), pp. 497-506, 1997

    // With fewer than 2 entries per node, each level has as many nodes as the level below
    // and the tree would never reach a single root node
    if (this->branchFactor < 2) {
        std::cerr << "R-tree branch factor must be at least 2" << std::endl;
        std::exit(1);
    }
    if (objectIDs.empty()) {
        std::cerr << "Cannot build R-tree with no objects" << std::endl;
        std::exit(1);
    }

    // Create a (leaf) entry for each input object
    std::size_t numEntries = objectIDs.size();
    std::vector<PtrRtreeEntry> leafLevelEntries(numEntries);
    for (std::size_t i = 0; i < numEntries; ++i) {
        leafLevelEntries[i].setEntryData(objectCoords[i],NULL,objectIDs[i]);
    }
    
    // Recursively build all levels (starting at leaf)
    PtrRtreeNode* root = this->buildLevel(0, leafLevelEntries);
    vectorise(root);
    
    // Delete all pointers as there not needed again
    std::vector<PtrRtreeNode*> allNodePtrs = this->getAllNodePtrs(root);
    for (std::size_t i = 0; i < allNodePtrs.size(); ++i) {
        // Note: This will try to delete child pointers in leaf nodes
        // (which are null) but this is fine as they are simply no-ops
        delete allNodePtrs[i];
    }
}

PtrRtreeNode* PackedRtree::buildLevel(int level, std::vector<PtrRtreeEntry>& levelEntries)
{
    int numNodes = std::ceil(static_cast<double>(levelEntries.size())/this->branchFactor);
    int numSlices = std::ceil(std::sqrt(numNodes));

    // Order entries by their X coordinate and divide into numSlices slices
    // Note: Each slice contain numSlices*this->branchFactor entries
    std::sort(levelEntries.begin(),levelEntries.end(),PtrRtreeEntry::compareX);
    
    // Order each slice by the Y coordinate and create node out of every 
    // this->branchFactor group of entries (except the last may have fewer)
    std::vector<PtrRtreeNode*> nodes;
    std::vector<Rectangle> rectangles;
    std::vector<PtrRtreeEntry>::iterator sliceStartIt, sliceEndIt, nodeStartIt, nodeEndIt;
    std::size_t nextSliceStartIdx = 0, nextGroupStartIdx = 0;
    for (std::size_t i = 0; i < levelEntries.size(); ) {
        // Sort the next group of this->branchFactor entries by their Y value
        // This is essentially considering one of the numSlices partitions
        // and we sort it according to the y value
        sliceStartIt = levelEntries.begin();
        std::advance(sliceStartIt,i);
        int maxEntriesInSlice = numSlices*this->branchFactor;
        if (i + maxEntriesInSlice < levelEntries.size()) {
            sliceEndIt = levelEntries.begin();
            std::advance(sliceEndIt,i + maxEntriesInSlice); // Note: We advanced by i first (for next slice)
            nextSliceStartIdx = i + maxEntriesInSlice;
        } else {
            // This mean we have reached the last slices and there
            //  were not be enough to fill the slice completely
            sliceEndIt = levelEntries.end();
            nextSliceStartIdx = levelEntries.size();
        }
        std::sort(sliceStartIt,sliceEndIt,PtrRtreeEntry::compareY);
        
        // Create node for each group of entries in the current slice
        for (std::size_t j = i; j < nextSliceStartIdx;) {
            PtrRtreeNode* rtreeNode = new PtrRtreeNode(level);
            Rectangle directoryRect;
            int maxEntriesInGroup = this->branchFactor;
            nodeStartIt = levelEntries.begin();
            std::advance(nodeStartIt,j);
            if (j + maxEntriesInGroup < nextSliceStartIdx) {
                nodeEndIt = levelEntries.begin();
                std::advance(nodeEndIt,j + maxEntriesInGroup); // Note: We advanced by j first (for next group in slice)
                nextGroupStartIdx = j + maxEntriesInGroup;
            } else {
                // This means we reached the last group of entries in this slice and 
                // there were not enough of them for to fill the node completely
                nodeEndIt = sliceEndIt;
                nextGroupStartIdx = nextSliceStartIdx;
            }
            // Every group has an entry: j is below the slice's end, and build checks that the branch factor is at least 2
            RN_DEBUG_ASSERT(nodeStartIt != nodeEndIt);
            // Initialise directory rectangle to cover first entry
            directoryRect.setDimensions(nodeStartIt->mbr);
            rtreeNode->addEntry(std::move(*nodeStartIt));
            ++nodeStartIt;
            for (; nodeStartIt != nodeEndIt; ++nodeStartIt) {
                directoryRect.expand(nodeStartIt->mbr);
                rtreeNode->addEntry(std::move(*nodeStartIt));
            }
            nodes.push_back(rtreeNode);
            rectangles.push_back(directoryRect);
            j = nextGroupStartIdx;
        }
        i = nextSliceStartIdx;
    }
    
    PtrRtreeNode* rootNode = NULL;
    std::size_t numEntries = nodes.size();
    if (numEntries == 1) {
        // This mean all entries were added to the root
        // and we can stop building the tree
        rootNode = nodes[0];
    } else {
        // Create the entries for each node and recursively build next level
        std::vector<PtrRtreeEntry> nextLevelEntries(numEntries);
        for (std::size_t i = 0; i < numEntries; ++i) {
            nextLevelEntries[i].setEntryData(std::move(rectangles[i]),nodes[i]);
        }
        rootNode = this->buildLevel(level+1,nextLevelEntries);
    }
    return rootNode;
}

BinaryMinHeap<EuclideanDistanceSqrd,PackedRtreeDataTuple> PackedRtree::getKNNs(unsigned int k, Coordinate queryPointX, Coordinate queryPointY, 
                                                                         std::vector<VertexID>& kNNs, std::vector<double>& kNNDistancesSqrd)
{
    BinaryMinHeap<EuclideanDistanceSqrd,PackedRtreeDataTuple> heap;
    EuclideanDistanceSqrd minDist;
    
    PackedRtreeDataTuple rootElement;
    rootElement.setNodeIdx(0);
    heap.insert(rootElement,0);
    while (heap.size() > 0) {
        minDist = heap.getMinKey();
        PackedRtreeDataTuple minElement = heap.extractMinElement();
        if (minElement.isNode) {
            if(this->nodes[minElement.nodeIdx].isLeaf()) {
                // Entries in leaf node are objects, compute Euclidean dist to it
                for (std::size_t i = 0; i < this->nodes[minElement.nodeIdx].entries.size(); ++i) {
                    minDist = this->nodes[minElement.nodeIdx].entries[i].mbr.getObjectDistSqrd(queryPointX,queryPointY);
                    PackedRtreeDataTuple element;
                    element.setObjectID(this->nodes[minElement.nodeIdx].entries[i].id);
                    heap.insert(element,minDist);
                }
            } else {
                // Entries in non-leaf nodes are other tree nodes, compute mindist
                for (std::size_t i = 0; i < this->nodes[minElement.nodeIdx].entries.size(); ++i) {
                    minDist = this->nodes[minElement.nodeIdx].entries[i].mbr.getMinDistSqrd(queryPointX,queryPointY);
                    PackedRtreeDataTuple element;
                    element.setNodeIdx(this->nodes[minElement.nodeIdx].entries[i].nodeIdx);
                    heap.insert(element,minDist);
                }
            }
        } else {
            // Then we have found object and the minimum distance to it is the actual distance
            // and since there is no closer node or object, it is the next nearest neighbour
            kNNs.push_back(minElement.objectID);
            kNNDistancesSqrd.push_back(minDist);
            if (kNNs.size() == k) {
                break;
            }
        }
    }
    
    return heap;
}

bool PackedRtree::getNextNearestNeighbour(BinaryMinHeap<EuclideanDistanceSqrd,PackedRtreeDataTuple>& heap, Coordinate queryPointX, 
                                          Coordinate queryPointY, VertexID& nn, EuclideanDistanceSqrd& nnDistSqrd)
{
    EuclideanDistanceSqrd minDist;
    while (heap.size() > 0) {
        minDist = heap.getMinKey();
        PackedRtreeDataTuple minElement = heap.extractMinElement();
        if (minElement.isNode) {
            if(this->nodes[minElement.nodeIdx].isLeaf()) {
                // Entries in leaf node are objects, compute Euclidean dist to it
                for (std::size_t i = 0; i < this->nodes[minElement.nodeIdx].entries.size(); ++i) {
                    minDist = this->nodes[minElement.nodeIdx].entries[i].mbr.getObjectDistSqrd(queryPointX,queryPointY);
                    PackedRtreeDataTuple element;
                    element.setObjectID(this->nodes[minElement.nodeIdx].entries[i].id);
                    heap.insert(element,minDist);
                }
            } else {
                // Entries in non-leaf nodes are other tree nodes, compute mindist
                for (std::size_t i = 0; i < this->nodes[minElement.nodeIdx].entries.size(); ++i) {
                    minDist = this->nodes[minElement.nodeIdx].entries[i].mbr.getMinDistSqrd(queryPointX,queryPointY);
                    PackedRtreeDataTuple element;
                    element.setNodeIdx(this->nodes[minElement.nodeIdx].entries[i].nodeIdx);
                    heap.insert(element,minDist);
                }
            }
        } else {
            // Then we have found object and the minimum distance to it is the actual distance
            // and since there is no closer node or object, it is the next nearest neighbour
            nn = minElement.objectID;
            nnDistSqrd = minDist;
            return true;
        }
    }
    // If we reach it means there are no objects to be found left
    return false;
}

BinaryMinHeap<EuclideanDistance,PackedRtreeDataTuple> PackedRtree::getAkNNs (unsigned int k, std::vector<Coordinate>& queryPointsX, std::vector<Coordinate>& queryPointsY, 
                                                                     std::vector<VertexID>& AkNNs, std::vector<EuclideanDistance>& AkNNDistances, 
                                                                     AggregateFunction< EuclideanDistance >& agg, std::vector<EuclideanDistance>& queryDistances)
{
    BinaryMinHeap<EuclideanDistance,PackedRtreeDataTuple> heap;
    EuclideanDistance aggDist;
    
    PackedRtreeDataTuple rootElement;
    rootElement.setNodeIdx(0);
    heap.insert(rootElement,0);
    while (heap.size() > 0) {
        aggDist = heap.getMinKey();
        PackedRtreeDataTuple minElement = heap.extractMinElement();
        if (minElement.isNode) {
            if(this->nodes[minElement.nodeIdx].isLeaf()) {
                // Entries in leaf node are objects, compute Euclidean dist to it
                for (std::size_t i = 0; i < this->nodes[minElement.nodeIdx].entries.size(); ++i) {
                    for (std::size_t j = 0; j < queryPointsX.size(); ++j) {
                        queryDistances[j] = std::sqrt(this->nodes[minElement.nodeIdx].entries[i].mbr.getObjectDistSqrd(queryPointsX[j],queryPointsY[j]));
                    }
                    aggDist = agg.aggregate(queryDistances);
                    PackedRtreeDataTuple element;
                    element.setObjectID(this->nodes[minElement.nodeIdx].entries[i].id);
                    heap.insert(element,aggDist);
                }
            } else {
                // Entries in non-leaf nodes are other tree nodes, compute mindist
                for (std::size_t i = 0; i < this->nodes[minElement.nodeIdx].entries.size(); ++i) {
                    for (std::size_t j = 0; j < queryPointsX.size(); ++j) {
                        queryDistances[j] = std::sqrt(this->nodes[minElement.nodeIdx].entries[i].mbr.getMinDistSqrd(queryPointsX[j],queryPointsY[j]));
                    }
                    aggDist = agg.aggregate(queryDistances);
                    PackedRtreeDataTuple element;
                    element.setNodeIdx(this->nodes[minElement.nodeIdx].entries[i].nodeIdx);
                    heap.insert(element,aggDist);
                }
            }
        } else {
            // Then we have found object and the minimum distance to it is the actual distance
            // and since there is no closer node or object, it is the next nearest neighbour
            AkNNs.push_back(minElement.objectID);
            AkNNDistances.push_back(aggDist);
            if (AkNNs.size() == k) {
                break;
            }
        }
    }
    
    return heap;
}

bool PackedRtree::getNextAggregateNearestNeighbour (BinaryMinHeap<EuclideanDistance,PackedRtreeDataTuple>& heap, std::vector<Coordinate>& queryPointsX, 
                                                std::vector< Coordinate >& queryPointsY, VertexID& ann, EuclideanDistance& annDist, 
                                                AggregateFunction<EuclideanDistance>& agg, std::vector<EuclideanDistance>& queryDistances)
{
    EuclideanDistance aggDist;
    while (heap.size() > 0) {
        aggDist = heap.getMinKey();
        PackedRtreeDataTuple minElement = heap.extractMinElement();
        if (minElement.isNode) {
            if(this->nodes[minElement.nodeIdx].isLeaf()) {
                // Entries in leaf node are objects, compute Euclidean dist to it
                for (std::size_t i = 0; i < this->nodes[minElement.nodeIdx].entries.size(); ++i) {
                    for (std::size_t j = 0; j < queryPointsX.size(); ++j) {
                        queryDistances[j] = std::sqrt(this->nodes[minElement.nodeIdx].entries[i].mbr.getObjectDistSqrd(queryPointsX[j],queryPointsY[j]));
                    }
                    aggDist = agg.aggregate(queryDistances);
                    PackedRtreeDataTuple element;
                    element.setObjectID(this->nodes[minElement.nodeIdx].entries[i].id);
                    heap.insert(element,aggDist);
                }
            } else {
                // Entries in non-leaf nodes are other tree nodes, compute mindist
                for (std::size_t i = 0; i < this->nodes[minElement.nodeIdx].entries.size(); ++i) {
                    for (std::size_t j = 0; j < queryPointsX.size(); ++j) {
                        queryDistances[j] = std::sqrt(this->nodes[minElement.nodeIdx].entries[i].mbr.getMinDistSqrd(queryPointsX[j],queryPointsY[j]));
                    }
                    aggDist = agg.aggregate(queryDistances);
                    PackedRtreeDataTuple element;
                    element.setNodeIdx(this->nodes[minElement.nodeIdx].entries[i].nodeIdx);
                    heap.insert(element,aggDist);
                }
            }
        } else {
            // Then we have found object and the minimum distance to it is the actual distance
            // and since there is no closer node or object, it is the next nearest neighbour
            ann = minElement.objectID;
            annDist = aggDist;
            return true;
        }
    }
    // If we reach it means there are no objects to be found left
    return false;

}

std::vector<VertexID> PackedRtree::windowQuery(Rectangle& rect)
{
    std::vector<VertexID> windowsPoints;
    std::vector<TreeNodeID> queue;
    
    queue.push_back(0);
    while (queue.size() > 0) {
        TreeNodeID currentNodeIdx = queue.back();
        queue.pop_back();
        if(this->nodes[currentNodeIdx].isLeaf()) {
            for (std::size_t i = 0; i < this->nodes[currentNodeIdx].entries.size(); ++i) {
                if (this->nodes[currentNodeIdx].entries[i].mbr.intersects(rect)) {
                    windowsPoints.push_back(this->nodes[currentNodeIdx].entries[i].id);
                }
            }
        } else {
            for (std::size_t i = 0; i < this->nodes[currentNodeIdx].entries.size(); ++i) {
                if (this->nodes[currentNodeIdx].entries[i].mbr.intersects(rect)) {
                    queue.push_back(this->nodes[currentNodeIdx].entries[i].nodeIdx);
                }
            }
        }
    }
    return windowsPoints;
}

std::vector<PtrRtreeNode*> PackedRtree::getAllNodePtrs(PtrRtreeNode* root)
{
    std::vector<PtrRtreeNode*> allNodePtrs;
    std::vector<PtrRtreeNode*> currentLevel, nextLevel;
    
    nextLevel.push_back(root);
    while(nextLevel.size() != 0){
        allNodePtrs.insert(allNodePtrs.end(), nextLevel.begin(), nextLevel.end());
        currentLevel.swap(nextLevel);
        nextLevel.clear();
        for (std::size_t i = 0; i < currentLevel.size(); ++i) {
            if (currentLevel[i] != NULL) {
                for (std::size_t j = 0; j < currentLevel[i]->entries.size(); ++j) {
                    nextLevel.push_back(currentLevel[i]->entries[j].rtreeNode);
                }
            }
        }
    }
    return allNodePtrs;
}

std::string PackedRtree::getObjSetType()
{
    return this->objSetType;
}

double PackedRtree::getObjSetDensity()
{
    return this->objSetDensity;
}

int PackedRtree::getObjSetSize()
{
    return this->objSetSize;
}

int PackedRtree::getObjSetVariable()
{
    return this->objSetVariable;
}

int PackedRtree::getBranchFactor()
{
    return this->branchFactor;
}

double PackedRtree::computeIndexSize()
{
    double memoryUsageBytes = 0;
    for (std::size_t i = 0; i < this->nodes.size(); ++i) {
        memoryUsageBytes += this->nodes[i].computeIndexSizeBytes();
    }
    return memoryUsageBytes/(1024*1024);
}

void PackedRtree::vectorise(PtrRtreeNode* root)
{
    typedef std::pair<PtrRtreeNode*,int> PtrIntPair;
    
    std::deque<PtrIntPair> queue;
    PackedRtreeNode rootNode(root->level);
    this->nodes.push_back(rootNode);
    queue.push_back(PtrIntPair(root,0));

    while (queue.size() > 0) {
        PtrIntPair currElement = queue.front();
        queue.pop_front();
        PtrRtreeNode* currNode = currElement.first;
        int currIdx = currElement.second;
        for (std::size_t i = 0; i < currNode->entries.size(); ++i) {
            PackedRtreeEntry entry;
            if (!currNode->isLeaf()) {
                int newNodeIdx = this->nodes.size();
                entry.setEntryData(currNode->entries[i].mbr,newNodeIdx,currNode->entries[i].id);
                PackedRtreeNode vecNode(currNode->entries[i].rtreeNode->level);
                this->nodes.push_back(vecNode);
                queue.push_back(PtrIntPair(currNode->entries[i].rtreeNode,newNodeIdx));
            } else {
                entry.setEntryData(currNode->entries[i].mbr,-1,currNode->entries[i].id);
            }
            this->nodes[currIdx].addEntry(entry);
        }
    }
}

/*
 * PackedRtreeNode
 */

PackedRtreeNode::PackedRtreeNode(int level): level(level)
{

}

void PackedRtreeNode::addEntry(PackedRtreeEntry entry)
{
    this->entries.push_back(entry);

}

bool PackedRtreeNode::isLeaf()
{
    return this->level == 0;
}


double PackedRtreeNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    for (std::size_t i = 0; i < this->entries.size(); ++i) {
        memoryUsage += this->entries[i].computeIndexSizeBytes();
    }
    memoryUsage += sizeof(int);
    return memoryUsage;

}

/*
 * PackedRtreeEntry
 */

void PackedRtreeEntry::setEntryData(Rectangle rect, int nodeIdx, VertexID id)
{
    this->mbr = rect;
    this->mbr.setCentroid(); 
    // Note: We assume centroid  has not been set (which allows us to 
    // delay computing it until all rectangle expansions are complete)
    this->nodeIdx = nodeIdx;
    this->id = id;
}

double PackedRtreeEntry::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(int);
    memoryUsage += sizeof(VertexID);
    memoryUsage += this->mbr.computeIndexSizeBytes();
    return memoryUsage;

}

/*
 * PtrRtreeNode
 */

PtrRtreeNode::PtrRtreeNode(int level): level(level)
{

}

void PtrRtreeNode::addEntry(PtrRtreeEntry entry)
{
    this->entries.push_back(entry);

}

bool PtrRtreeNode::isLeaf()
{
    return this->level == 0;
}

/*
 * PtrRtreeEntry
 */

void PtrRtreeEntry::setEntryData(CoordinatePair& ptCoords, PtrRtreeNode* rtreeNode, VertexID id)
{
    this->mbr.setDimensions(ptCoords.first,ptCoords.second,ptCoords.first,ptCoords.second);
    this->mbr.setCentroid();
    this->rtreeNode = rtreeNode;
    this->id = id;
}

void PtrRtreeEntry::setEntryData(Rectangle rect, PtrRtreeNode* rtreeNode, VertexID id)
{
    this->mbr = rect;
    this->mbr.setCentroid(); 
    // Note: We assume centroid  has not been set (which allows us to 
    // delay computing it until all rectangle expansions are complete)
    this->rtreeNode = rtreeNode;
    this->id = id;
}

bool PtrRtreeEntry::compareX(const PtrRtreeEntry& entry1, const PtrRtreeEntry& entry2)
{
    return entry1.mbr.centreX < entry2.mbr.centreX;
}

bool PtrRtreeEntry::compareY(const PtrRtreeEntry& entry1, const PtrRtreeEntry& entry2)
{
    return entry1.mbr.centreY < entry2.mbr.centreY;
}

/* 
 * Rectangle
 */

void Rectangle::setDimensions(Coordinate left, Coordinate bottom, Coordinate right, Coordinate top)
{
    this->left = left;
    this->bottom = bottom;
    this->right = right;
    this->top = top;
}

void Rectangle::setDimensions(Rectangle& rect)
{
    this->left = rect.left;
    this->bottom = rect.bottom;
    this->right = rect.right;
    this->top = rect.top;
}

void Rectangle::setCentroid()
{
    this->centreX = (this->left+this->right)/2;
    this->centreY = (this->bottom+this->top)/2;
}

void Rectangle::expand(Rectangle& rect)
{
    // Expand the current rectangle to all cover given rectangle
    this->left = this->left > rect.left ? rect.left : this->left;
    this->bottom = this->bottom > rect.bottom ? rect.bottom : this->bottom;
    this->right = this->right < rect.right ? rect.right : this->right;
    this->top = this->top < rect.top ? rect.top : this->top;
}

EuclideanDistanceSqrd Rectangle::getMinDistSqrd(Coordinate x, Coordinate y)
{
    Coordinate nearestX;
    if (x < this->left) {
        nearestX = this->left;
    } else if (x > this->right) {
        nearestX = this->right;
    } else {
        nearestX = x;
    }
    
    Coordinate nearestY;
    if (y < this->bottom) {
        nearestY = this->bottom;
    } else if (y > this->top) {
        nearestY = this->top;
    } else {
        nearestY = y;
    }
    
    return std::pow(nearestX-x,2) + std::pow(nearestY-y,2);
}

EuclideanDistanceSqrd Rectangle::getMinMaxDistSqrd(Coordinate x, Coordinate y)
{
    Coordinate furthestX, furthestY, nearestX, nearestY;

    // Get the corner of the rectangle that is further from the pointers
    // Note: This must necessarily be one of the corners
    EuclideanDistance midX = left + (right-left)/2, midY = bottom + (top-bottom)/2;

    if (x >= midX) {
        furthestX = left;
    } else {
        furthestX = right;
    }
    
    if (y >= midY) {
        furthestY = bottom;
    } else {
        furthestY = top;
    }
    
    // Compute the distance to this further point
    EuclideanDistanceSqrd maxDistInXSqrd, maxDistInYSqrd, maxDistSqrd, minMaxDistSqrd, currMinMaxDistSqrd;
    maxDistInXSqrd = std::pow(furthestX-x,2);
    maxDistInYSqrd = std::pow(furthestY-y,2);
    maxDistSqrd = maxDistInXSqrd + maxDistInYSqrd;
    
    // Remove the x (and then y) component of the distance to the furthest point
    // and replace it with the x (and then y component) of the nearest point
    // in the x (and then y) dimension) to calculate minmaxdist
    if (x <= midX) {
        nearestX = left;
    } else {
        nearestX = right;
    }
    minMaxDistSqrd = maxDistSqrd - maxDistInXSqrd + std::pow(nearestX-x,2);
    
    if (y <= midY) {
        nearestY = bottom;
    } else {
        nearestY = top;
    }
    currMinMaxDistSqrd = maxDistSqrd - maxDistInYSqrd + std::pow(nearestY-y,2);
    
    if (currMinMaxDistSqrd < minMaxDistSqrd) {
        minMaxDistSqrd = currMinMaxDistSqrd;
    }
    
    return minMaxDistSqrd;
}

EuclideanDistanceSqrd Rectangle::getMaxDistSqrd(Coordinate x, Coordinate y)
{
    Coordinate furthestX, furthestY;

    // Get the corner of the rectangle that is further from the point
    // Note: This must necessarily be one of the corners
    EuclideanDistance midX = left + (right-left)/2, midY = bottom + (top-bottom)/2;

    if (x >= midX) {
        furthestX = left;
    } else {
        furthestX = right;
    }
    
    if (y >= midY) {
        furthestY = bottom;
    } else {
        furthestY = top;
    }

    return std::pow(furthestX-x,2) + std::pow(furthestY-y,2);
}

// This assumes that the rectage is a point (therefore both corners are the same point)
EuclideanDistanceSqrd Rectangle::getObjectDistSqrd(Coordinate x, Coordinate y)
{
    return std::pow(left-x,2) + std::pow(bottom-y,2);
}

bool Rectangle::intersects(Rectangle& rect)
{
    return left <= rect.right && right >= rect.left &&
           bottom <= rect.top && top >= rect.bottom;
}

double Rectangle::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(Coordinate)*6;
    return memoryUsage;
}
