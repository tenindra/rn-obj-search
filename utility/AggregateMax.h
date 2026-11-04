// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _AGGREGATEMAX_H
#define _AGGREGATEMAX_H

#include "AggregateFunction.h"

#include <algorithm>

template <typename T>
class AggregateMax: public AggregateFunction<T> {

    public:
        T aggregate(std::vector<T>& distances) {
            return *std::max_element(std::begin(distances),std::end(distances));
        }
        
        T getThreshold(T distance, int numQueryVertices) {
            return distance;
        }
        
        T getLowerBound(T distance, int numQueryVertices) {
            return distance;
        }

        std::string getPrintString() {
            return "max";
        }

        // Assumes both vectors are sorted on edge-weights
        int getMinObjListIndex(std::vector<VertexEdgeWeightPair>& objectList, std::vector<EdgeWeight>& queryVertexLmkDists) {
            // Each query vertex represents the break points of a sum of absolute value functions
            // The minimum value is given by the median break point
            std::vector<EdgeWeight> queryVertexLmkDistsCopy = queryVertexLmkDists; // Make a copy so we leave original intact
            std::sort(queryVertexLmkDistsCopy.begin(),queryVertexLmkDistsCopy.end()); // Sort so we can find min and max dists (we only do this once per OL)
            EdgeWeight minDist = (queryVertexLmkDistsCopy[0]+queryVertexLmkDistsCopy[queryVertexLmkDistsCopy.size()-1]) >> 1;
            
            // Note: If the sum of min and max distances is odd, then minDist will have decimal
            // component of 0.5 dropped. E.g. if the minDist = 12.5 and there object distances of 11 and 13
            // bother distances appear to the floored value of 12, even though 13 is the real minimum.
            // To avoid this problem, we need to always choose the larger of object distance, when two
            // object distances appear as close to minDist. Note that if minDist really is 12, then both
            // 11 and 13 are valid, and always choosing is still correct.
            
            // After this left equal right, and points to element greater than or equal to median dist
            std::size_t left = 0, right = objectList.size(), mid;
            while(left != right) {
                mid = (right+left)/2;
                if (objectList[mid].second == minDist) {
                    return mid;
                } else if (objectList[mid].second > minDist) {
                    right = mid;
                } else /* (objectList[mid].second < minDist)*/ {
                    left = mid+1;
                }
            }
            
            if (right == objectList.size()) {
                return objectList.size()-1;
            } else {
                if (right > 0) {
                    left = right-1; // First value smaller than minDist (minDist not present if it reaches here)
                    EdgeWeight leftDiff, rightDiff;
                    leftDiff = (objectList[left].second > minDist) ? objectList[left].second - minDist : minDist - objectList[left].second;
                    rightDiff = (objectList[right].second > minDist) ? objectList[right].second - minDist : minDist - objectList[right].second;
                    if (leftDiff < rightDiff) {
                        // Only return left, if difference is strictly less than difference for right value
                        // to avoid error described above
                        return left;
                    } else {
                        return right;
                    }                        
                } else {
                    return right;
                }
            }
        }

};

#endif // _AGGREGATEMAX_H