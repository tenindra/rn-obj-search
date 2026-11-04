// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "MortonCode.h"

#include <bitset>
#include <cstdlib>
#include <iostream>

// For points (width is 1, so len represents whole maxRangeExp)
MortonCode::MortonCode(Coordinate x, Coordinate y, int maxRangeExp)
{
    // Need to check that max length can fit in this MortonCode
    RN_DEBUG_ASSERT(maxRangeExp <= quadtree_constants::MAX_MORTON_LENGTH/2 && "Morton code will not fit in 64-bit zNumber");
    if (x < 0 || y < 0) {
        // The coordinates are translated so that the smallest is 0, and a negative one would wrap
        std::cerr << "INTERNAL ERROR: Morton code requested for a point with a negative coordinate (" << x << "," << y << ")" << std::endl;
        exit(1);
    }
    
    MortonNumber mx = x, my = y;

    mx = (mx | (mx << quadtree_constants::S[4])) & quadtree_constants::B[4];
    mx = (mx | (mx << quadtree_constants::S[3])) & quadtree_constants::B[3];
    mx = (mx | (mx << quadtree_constants::S[2])) & quadtree_constants::B[2];
    mx = (mx | (mx << quadtree_constants::S[1])) & quadtree_constants::B[1];
    mx = (mx | (mx << quadtree_constants::S[0])) & quadtree_constants::B[0];

    my = (my | (my << quadtree_constants::S[4])) & quadtree_constants::B[4];
    my = (my | (my << quadtree_constants::S[3])) & quadtree_constants::B[3];
    my = (my | (my << quadtree_constants::S[2])) & quadtree_constants::B[2];
    my = (my | (my << quadtree_constants::S[1])) & quadtree_constants::B[1];
    my = (my | (my << quadtree_constants::S[0])) & quadtree_constants::B[0];
    
    this->len = maxRangeExp << 1; // * 2
    this->zNumber = mx | (my << 1);
    this->zNumber <<= (quadtree_constants::MAX_MORTON_LENGTH - this->len);
    // Note: We can use this->len here instead of maxRangeExp << 1
    // because this is a point, the length will always be maxRangeExp
}

MortonCode::MortonCode(MortonNumber zNumber, MortonLength length): zNumber(zNumber), len(length) {}

void MortonCode::setCode(const MortonCode& parentBlock, Direction dir)
{
    MortonNumber parentZNumber = parentBlock.getZNumber();
    MortonLength parentLen = parentBlock.getLength();
    if (parentLen+2 > quadtree_constants::MAX_MORTON_LENGTH) {
        // The bit masks below have one entry per length a block can be split at
        std::cerr << "INTERNAL ERROR: cannot split a quadtree block of Morton code length " << static_cast<int>(parentLen) << std::endl;
        exit(1);
    }
    if (dir == SW) {
        this->zNumber = parentZNumber | quadtree_constants::ARRAY_SW_BIT_MASK[parentLen];
        // Note: This is equivalent to SW_MC << 64 - 2 - parentLen
        // I.e. we are the morton code for this quadrant to the end
        // of the parent block's morton code
    } else if (dir == SE) {
        this->zNumber = parentZNumber | quadtree_constants::ARRAY_SE_BIT_MASK[parentLen];
    } else if (dir == NW) {
        this->zNumber = parentZNumber | quadtree_constants::ARRAY_NW_BIT_MASK[parentLen];
    } else if (dir == NE) {
        this->zNumber = parentZNumber | quadtree_constants::ARRAY_NE_BIT_MASK[parentLen];
    } else {
        std::cerr << "Invalid direction provided" << std::endl;
        exit(1);
    }
    this->len = parentLen + 2;
}

MortonNumber MortonCode::getZNumber() const
{
    return this->zNumber;
}

MortonLength MortonCode::getLength() const
{
    return this->len;
}

bool MortonCode::contains(const MortonCode& point) const
{
    RN_DEBUG_ASSERT(point.getLength() >= this->len && "Input code is shorter than this code");
    MortonNumber testNum = point.getZNumber();
    // Clear least significant bits that are not used by this code
    testNum = testNum >> (quadtree_constants::MAX_MORTON_LENGTH - this->len);
    testNum = testNum << (quadtree_constants::MAX_MORTON_LENGTH - this->len);
    // If these two numbers are the same we return true because
    // the test code is prefixed by this code
    return (this->zNumber == testNum);
}

void MortonCode::decompose(int maxRangeExp, Coordinate& x, Coordinate& y, int& width)
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

bool MortonCode::operator==(const MortonCode& code) const  {
    if (this->len == code.getLength() && this->zNumber == code.getZNumber()) {
        return true;
    } else {
        return false;
    }
}

bool MortonCode::operator<(const MortonCode &code) const {
    if (this->zNumber < code.getZNumber()) {
        return true;
    } else {
        return false;
    }
}

bool MortonCode::operator>(const MortonCode &code) const  {
    if (this->zNumber > code.getZNumber()) {
        return true;
    } else {
        return false;
    }
}

void MortonCode::printCode() {
     std::cout << "Morton Code = " << std::bitset<64>(this->zNumber) << " (Len=" << static_cast<int>(this->len) << ")" << std::endl;
}

double MortonCode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(this->zNumber);
    memoryUsage += sizeof(this->len);
    return memoryUsage;
}
