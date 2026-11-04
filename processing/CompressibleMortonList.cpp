// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "CompressibleMortonList.h"

#include "DijkstraSearch.h"
#include "../utility/geometry.h"
#include "../utility/utility.h"

#include <utility>
#include <bitset>
#include <iomanip>
#include <iostream>

/*
 * CompressibleMortonList
 */

void CompressibleMortonList::buildCompressibleMortonList(CompressibleQuadtree& qt) {
    std::vector<CompressibleMortonBlock> mortonList;
    CompressibleQuadtreeNode* currNode = qt.firstLeaf;
    
    // Note: We assume quadtree leaves contain at most one colour
    // and here we are simply converting to Morton List
    
    while (currNode != NULL) {
        // Note: A quadtree leaf without vertices has no colour, and is left out of the Morton list
        if (currNode->colourIDs.size() == 1) {
            CompressibleMortonBlock block(currNode->mortonCode,*(currNode->colourIDs.begin()));
            mortonList.push_back(block);
        } else if (currNode->colourIDs.size() > 1) {
            std::cout << "Quadtree leaf node has more than one colour, unable to build CompressibleMortonList!\n";
            std::cout << "Num Colours = " << currNode->colourIDs.size() << "\n";
            std::cout << "Colours are: ";
            for (auto it = currNode->colourIDs.begin(); it != currNode->colourIDs.end(); ++it) {
                std::cout << *it << ", ";
            }
            std::cout << "\n";
            std::exit(1);
        } /* else {
            // If it has no colour the it is implicitly dropped from MortonList so as not to waste space
        } */
        currNode = currNode->nextLeaf;
    }
    
    // We know that binary search tree will contain all the elements of list
    this->binarySearchTreeCompressibleMortonList.reserve(mortonList.size());
    // Build binary tree from linear quadtree Morton list
    // Note: We exploit the fact that Morton list is sorted
    // Each level holds the (left,right) ranges of its nodes. The nodes are added in level order, which for this
    // complete tree puts the children of position p at 2p+1 and 2p+2 (getLeftChild, getRightChild)
    std::vector<std::pair<int,int>> currentLevel, nextLevel;
    unsigned int left = 0, right = mortonList.size()-1, mid, size, largestPowerOf2;
    nextLevel.push_back(std::make_pair(left,right));
    while (nextLevel.size() > 0) {
        currentLevel = nextLevel;
        nextLevel.clear();
        for (auto it = currentLevel.begin(); it != currentLevel.end(); ++it) {
            left = it->first;
            right = it->second;
            size = right - left + 1; // right and left inclusive
            if (left == right) {
                this->binarySearchTreeCompressibleMortonList.push_back(std::move(mortonList[left]));
            } else if (left < right) {
                // We are building a complete binary search tree for space efficiency
                // and ease of inserts. This means that all levels are full up until
                // the final level where the elements are as far left as possible.
                
                // If n is the size of a binary search then the number elements on the left
                // hand side of the root is 2^x - 1 when the left handside is full (where 2^x
                // is the largest power of 2 that is less than n). However if the left-hand
                // side is not full the number of elements is n - 1 - (2^(x-1) - 1). This
                // is the total number of elements minus the root and minus the total number 
                // of elements in the right hand. Note that the right hand is always full, 
                // except the bottom level (which is why we use (x-1). Recall that we
                // always which fill from the left so since in thise case the left-hand side
                // is not full, the right hand side must be empty. n - 1 - (2^(x-1) - 1 can be
                // simplified to n - 2^(x-1), if n is full this expression will always greater
                // than or equal to 2^x - 1, because we are not considering the right-hand side
                // elements that may be present when the left-hand side is full. So we can
                // simply choose the minimum of the two and obtain the correct number of elments
                // on the left hand side.
                
                largestPowerOf2 = geometry::maxPower2(size); // Less than size
                mid = left + std::min(largestPowerOf2-1,size-(largestPowerOf2 >> 1));

                this->binarySearchTreeCompressibleMortonList.push_back(std::move(mortonList[mid]));

                nextLevel.push_back(std::make_pair(left,mid-1));
                nextLevel.push_back(std::make_pair(mid+1,right));
            }
        }
    }
    
    if (this->binarySearchTreeCompressibleMortonList.size() != mortonList.size()) {
        std::cerr << "INTERNAL ERROR: " << this->binarySearchTreeCompressibleMortonList.size() << " of the " << mortonList.size()
                  << " NVD blocks reached the search array" << std::endl;
        std::exit(1);
    }
}

VertexID CompressibleMortonList::getColourID(const MortonCode& targetCode) {
    // Binary search morton list to find block containing target
    int listIndex = this->getBlockIdxContainingTargetByBST(targetCode);
    return this->binarySearchTreeCompressibleMortonList[listIndex].getColourID();
}

std::size_t CompressibleMortonList::getBlockIdxContainingTargetByBST(const MortonCode& targetCode) {
    bool found = false;
    std::size_t parent = 0; // Start at root element
    int status;
    while(parent < this->binarySearchTreeCompressibleMortonList.size()) {
        status = this->binarySearchTreeCompressibleMortonList[parent].compareMortonCode(targetCode);
        if (status == 0) {
            // The target is contained within this morton block
            found = true;
            break;
        } else if (status == 1) {
            // Target is larger so it is in the right branch
            parent = this->getRightChild(parent);
        } else /* if (status == = -1)*/ {
            // Target is smaller is it is in the left branch
            parent = this->getLeftChild(parent);
        }
    }
    if (found) {
        return parent;
    } else {
        std::cerr << "getBlockIdxContainingTargetByBST() could not find containing block of target morton code in morton list" << std::endl;
        exit(1);
    }
}

std::size_t CompressibleMortonList::getLeftChild(std::size_t parent)
{
    return parent*2 + 1;
}

std::size_t CompressibleMortonList::getRightChild(std::size_t parent)
{
    return parent*2 + 2;
}

bool CompressibleMortonList::isTreeSorted(int pos)
{
    bool sorted = true;
    std::size_t leftChild = this->getLeftChild(pos);
    MortonCode parentCode, leftCode, rightCode;
    if (leftChild < this->binarySearchTreeCompressibleMortonList.size()) {
        parentCode = this->binarySearchTreeCompressibleMortonList[pos].getMortonCode();
        leftCode = this->binarySearchTreeCompressibleMortonList[leftChild].getMortonCode();
        if (leftCode > parentCode) {
            std::cout << "Left child at " << leftChild << " greater than parent at " << pos << std::endl;
            sorted = false;
        } else {
            sorted = this->isTreeSorted(leftChild);
        }
        // We don't have to check right child if left child is bigger than vector
        // because left child always appears before right child
        std::size_t rightChild = this->getRightChild(pos);
        if (rightChild < this->binarySearchTreeCompressibleMortonList.size()) {
            rightCode = this->binarySearchTreeCompressibleMortonList[rightChild].getMortonCode();
            if (rightCode < parentCode) {
                std::cout << "Right child at " << rightChild << " smaller than parent at " << pos << std::endl;
                sorted = false;
            } else {
                sorted = this->isTreeSorted(rightChild);
            }
        } 
    }
    return sorted;
}

// Check left and right subtree confirm to binary search tree property
// i.e. all nodes on left are small and all nodes on right are larger
bool CompressibleMortonList::areSubtreesValid(int pos)
{
    bool isValid = true;
    MortonCode parentCode, leftCode, rightCode;
    std::size_t leftChild = this->getLeftChild(pos);
    if (leftChild < this->binarySearchTreeCompressibleMortonList.size()) {
        parentCode = this->binarySearchTreeCompressibleMortonList[pos].getMortonCode();
        leftCode = this->binarySearchTreeCompressibleMortonList[leftChild].getMortonCode();
        if (leftCode > parentCode) {
            std::cout << "Left child at " << leftChild << " greater than parent at " << pos << std::endl;
            isValid = false;
        } else {
            isValid = this->isLeftSubtreeValid(leftChild,parentCode);
        }
        // We don't have to check right child if left child is bigger than vector
        // because left child always appears before right child
        std::size_t rightChild = this->getRightChild(pos);
        if (rightChild < this->binarySearchTreeCompressibleMortonList.size()) {
            rightCode = this->binarySearchTreeCompressibleMortonList[rightChild].getMortonCode();
            if (rightCode < parentCode) {
                std::cout << "Right child at " << rightChild << " smaller than parent at " << pos << std::endl;
                isValid = false;
            } else {
                isValid = this->isRightSubtreeValid(rightChild,parentCode);
            }
        } 
    }
    return isValid;
}

// All nodes on the left subtree must be smaller than parent
bool CompressibleMortonList::isLeftSubtreeValid(int subtreeRoot, MortonCode& parentCode)
{
    bool isValid = true;
    MortonCode leftCode, rightCode;
    std::size_t leftChild = this->getLeftChild(subtreeRoot);
    if (leftChild < this->binarySearchTreeCompressibleMortonList.size()) {
        leftCode = this->binarySearchTreeCompressibleMortonList[leftChild].getMortonCode();
        if (leftCode > parentCode) {
            std::cout << "Tree node in left subtree of root parent at " << leftChild << " greater than root parent" << std::endl;
            parentCode.printCode();
            leftCode.printCode();
            isValid = false;
        } else {
            isValid = this->isLeftSubtreeValid(leftChild,parentCode);
        }
        // We don't have to check right child if left child is bigger than vector
        // because left child always appears before right child
        std::size_t rightChild = this->getRightChild(subtreeRoot);
        if (rightChild < this->binarySearchTreeCompressibleMortonList.size()) {
            rightCode = this->binarySearchTreeCompressibleMortonList[rightChild].getMortonCode();
            if (rightCode > parentCode) {
                std::cout << "Tree node in left subtree of root parent at " << rightChild << " greater than root parent" << std::endl;
                parentCode.printCode();
                leftCode.printCode();
                isValid = false;
            } else {
                isValid = this->isLeftSubtreeValid(rightChild,parentCode);
            }
        } 
    }
    return isValid;
}

// All nodes on the right subtree must be larger than parent
bool CompressibleMortonList::isRightSubtreeValid(int subtreeRoot, MortonCode& parentCode)
{
    bool isValid = true;
    MortonCode leftCode, rightCode;
    std::size_t leftChild = this->getLeftChild(subtreeRoot);
    if (leftChild < this->binarySearchTreeCompressibleMortonList.size()) {
        leftCode = this->binarySearchTreeCompressibleMortonList[leftChild].getMortonCode();
        if (leftCode < parentCode) {
            std::cout << "Tree node in right subtree of root parent at " << leftChild << " smaller than root parent" << std::endl;
            parentCode.printCode();
            leftCode.printCode();
            isValid = false;
        } else {
            isValid = this->isRightSubtreeValid(leftChild,parentCode);
        }
        // We don't have to check right child if left child is bigger than vector
        // because left child always appears before right child
        std::size_t rightChild = this->getRightChild(subtreeRoot);
        if (rightChild < this->binarySearchTreeCompressibleMortonList.size()) {
            rightCode = this->binarySearchTreeCompressibleMortonList[rightChild].getMortonCode();
            if (rightCode < parentCode) {
                std::cout << "Tree node in right subtree of root parent at " << rightChild << " smaller than root parent" << std::endl;
                parentCode.printCode();
                leftCode.printCode();
                isValid = false;
            } else {
                isValid = this->isRightSubtreeValid(rightChild,parentCode);
            }
        } 
    }
    return isValid;
}

int CompressibleMortonList::getNumBlocks() {
    return this->binarySearchTreeCompressibleMortonList.size();
}

void CompressibleMortonList::printBST()
{
    for (std::size_t i = 0; i < this->binarySearchTreeCompressibleMortonList.size(); ++i) {
        MortonCode code = this->binarySearchTreeCompressibleMortonList[i].getMortonCode();
        std::cout << std::setfill('0') << std::setw(3)  << i << ": ";
        code.printCode();
        std::cout << "Object ID: " << static_cast<VertexID>(this->binarySearchTreeCompressibleMortonList[i].getColourID()) << std::endl;        
    }
}

double CompressibleMortonList::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(this->binarySearchTreeCompressibleMortonList);
    for (std::size_t i = 0; i < this->binarySearchTreeCompressibleMortonList.size(); ++i) {
        memoryUsage += this->binarySearchTreeCompressibleMortonList[i].computeIndexSizeBytes();
    }
    return memoryUsage;
}

void CompressibleMortonList::clear()
{
    this->binarySearchTreeCompressibleMortonList.clear();
}

/*
 * CompressibleMortonBlock
 */

CompressibleMortonBlock::CompressibleMortonBlock(MortonCode code, VertexID colourID): zNumber(code.getZNumber()), len(code.getLength()), colourID(colourID) {}

VertexID CompressibleMortonBlock::getColourID() {
    return this->colourID;
}

void CompressibleMortonBlock::printBlock() {
    std::cout << "Morton Code = " << std::bitset<64>(this->zNumber) << " (Len=" << static_cast<int>(this->len) << ")" << std::endl;
    std::cout << "Object ID: " << this->colourID << std::endl;
}

MortonNumber CompressibleMortonBlock::getMortonNumber()
{
    return this->zNumber;
}

MortonLength CompressibleMortonBlock::getMortonLength()
{
    return this->len;
}

int CompressibleMortonBlock::compareMortonCode(const MortonCode& code) const {
    // We assume morton code is shorter than input
    RN_DEBUG_ASSERT(code.getLength() >= this->len && "Input code is shorter than this code");
    MortonNumber testNum = 0;
    if (this->len > 0) {
        // We must shift by number less than MAX_MORTON_LENGTH otherwise it's undefined in C++
        // (but in our case shifting by length is equivalent to setting to 0)
        testNum = code.getZNumber();
        // Clear least significant bits that are not used by this code
        testNum = testNum >> (quadtree_constants::MAX_MORTON_LENGTH - this->len);
        testNum = testNum << (quadtree_constants::MAX_MORTON_LENGTH - this->len);
    }
    
    if (testNum < this->zNumber) {
        return -1;
    } else if (testNum > this->zNumber) {
        return 1;
    } else /* if (testNum == this->zNumber */ {
        return 0;
    }    
}

MortonCode CompressibleMortonBlock::getMortonCode() {
    return MortonCode(this->zNumber,this->len);
}

void CompressibleMortonBlock::decompose(int maxRangeExp, Coordinate& x, Coordinate& y, int& width)
{
    MortonNumber unshiftedZ = this->zNumber >> (quadtree_constants::MAX_MORTON_LENGTH - (maxRangeExp << 1)); 
    // Note: We can't use this->length here because we don't know if it is point or block
    
    MortonNumber mx = unshiftedZ & quadtree_constants::B[0], my = (unshiftedZ >> 1) & quadtree_constants::B[0];

    mx = (mx | (mx >> quadtree_constants::S[0])) & quadtree_constants::B[1];
    mx = (mx | (mx >> quadtree_constants::S[1])) & quadtree_constants::B[2];
    mx = (mx | (mx >> quadtree_constants::S[2])) & quadtree_constants::B[3];
    mx = (mx | (mx >> quadtree_constants::S[3])) & quadtree_constants::B[4];
    mx = (mx | (mx >> quadtree_constants::S[4])) & quadtree_constants::B[5];

    my = (my | (my >> quadtree_constants::S[0])) & quadtree_constants::B[1];
    my = (my | (my >> quadtree_constants::S[1])) & quadtree_constants::B[2];
    my = (my | (my >> quadtree_constants::S[2])) & quadtree_constants::B[3];
    my = (my | (my >> quadtree_constants::S[3])) & quadtree_constants::B[4];
    my = (my | (my >> quadtree_constants::S[4])) & quadtree_constants::B[5];
    
    x = mx;
    y = my;
    
    unsigned int widthPower = maxRangeExp - (this->len >> 1); // No precision will be lost because this->len is always multiple of two
    width = 1;
    width = width << widthPower;
}

bool CompressibleMortonBlock::operator==(const MortonCode& code) const  {
    if (this->len == code.getLength() && this->zNumber == code.getZNumber()) {
        return true;
    } else {
        return false;
    }
}

double CompressibleMortonBlock::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(this->colourID);
    memoryUsage += sizeof(this->zNumber);
    memoryUsage += sizeof(this->len);
    return memoryUsage;
}
