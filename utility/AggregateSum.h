// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _AGGREGATESUM_H
#define _AGGREGATESUM_H

#include "AggregateFunction.h"

#include <algorithm>
#include <numeric>

template <typename T>
class AggregateSum: public AggregateFunction<T> {

    public:

        T aggregate(std::vector<T>& distances) {
            return std::accumulate(distances.begin(), distances.end(), 0.0);
        }
        
        T getThreshold(T distance, int numQueryVertices) {
            return distance/numQueryVertices;
        }
        
        T getLowerBound(T distance, int numQueryVertices) {
            return distance*numQueryVertices;
        }

        std::string getPrintString() {
            return "sum";
        }

        // Assumes both vectors are sorted on edge-weights and neither are empty
        int getMinObjListIndex(std::vector<VertexEdgeWeightPair>& objectList, std::vector<EdgeWeight>& queryVertexLmkDists) {
            // Each query vertex represents the break points of a sum of absolute value functions
            // The minimum value is given by the median break point
            std::vector<EdgeWeight> queryVertexLmkDistsCopy = queryVertexLmkDists; // Make a copy so we leave original intact
            std::sort(queryVertexLmkDistsCopy.begin(),queryVertexLmkDistsCopy.end()); // Sort so we can find median dist (we only do this once per OL)
            size_t n = queryVertexLmkDistsCopy.size() >> 1;
            EdgeWeight medianDist = queryVertexLmkDistsCopy[n];
            
            // Note: When there is an even number of elements in the object list, we choose the 
            // smaller query landmark distance as the median. This means, when we search the object list
            // we need to find the object list element with distance (second value in pair) equal to 
            // or greater than the median distance. When the number of elements is odd, then we find the 
            // element closest to the median distance. This ensure correctness in the object with 
            // the smallest aggregate lower-bound
            
            // After this left equal right, and points to element greater than or equal to median dist
            std::size_t left = 0, right = objectList.size(), mid;
            while(left != right) {
                mid = (right+left)/2;
                if (objectList[mid].second == medianDist) {
                    return mid;
                } else if (objectList[mid].second > medianDist) {
                    right = mid;
                } else /* (objectList[mid].second < medianDist)*/ {
                    left = mid+1;
                }
            }
            
            if (right == objectList.size()) {
                return objectList.size()-1;
            } else {
                if (objectList.size()%2 == 0) {
                    return right;
                } else {
                    if (right > 0) {
                        left = right-1; // First value smaller than medianDist (medianDist not present if it reaches here)
                        EdgeWeight leftDiff, rightDiff;
                        leftDiff = (objectList[left].second > medianDist) ? objectList[left].second - medianDist : medianDist - objectList[left].second;
                        rightDiff = (objectList[right].second > medianDist) ? objectList[right].second - medianDist : medianDist - objectList[right].second;
                        if (leftDiff < rightDiff) {
                            return left;
                        } else {
                            return right;
                        }                        
                    } else {
                        return right;
                    }
                }
            }
        }

};

#endif // _AGGREGATESUM_H
