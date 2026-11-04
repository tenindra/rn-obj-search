// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _GROUNDTRUTH_H
#define _GROUNDTRUTH_H

#include "../common.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

// The queries of one experiment list (e.g. "varyingk") of one query type
template <typename Query>
struct QueryList {
    std::string name;
    int aggFunc; // constants::AGG_MAX or constants::AGG_SUM for AkNN lists, -1 otherwise
    std::vector<Query> queries;
};

// The ground truth of one query list: the queries (with every value resolved) and the correct result of each.
// Results are (object, distance) pairs sorted by (distance, object) for kNN/AkNN/range and by (distance desc, object) for FkN.
// Range results have distance 0 because INE does not report range distances.
template <typename Query>
struct GroundTruthSet {
    std::string networkName;
    int numVertices = 0;
    int numEdges = 0;
    std::string queryType; // knn, aknn, fkn or range
    std::string listName;
    int aggFunc = -1;      // constants::AGG_MAX or constants::AGG_SUM for AkNN, -1 otherwise
    std::vector<Query> queries;
    std::vector<std::vector<VertexDistancePair>> results; // results[i] is the result of queries[i]

    // Boost Serialization
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & this->networkName;
        ar & this->numVertices;
        ar & this->numEdges;
        ar & this->queryType;
        ar & this->listName;
        ar & this->aggFunc;
        ar & this->queries;
        ar & this->results;
    }
};

template <typename Query>
void saveGroundTruth(GroundTruthSet<Query>& groundTruth, const std::string& filePath)
{
    std::ofstream fs(filePath, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!fs.is_open()) {
        std::cerr << "Cannot open ground truth file " << filePath << " for writing" << std::endl;
        exit(1);
    }
    boost::archive::binary_oarchive oa(fs);
    oa << groundTruth;
}

// Returns false if the file cannot be opened
template <typename Query>
bool loadGroundTruth(const std::string& filePath, GroundTruthSet<Query>& groundTruth)
{
    std::ifstream fs(filePath, std::ios::in | std::ios::binary);
    if (!fs.is_open()) {
        return false;
    }
    boost::archive::binary_iarchive ia(fs);
    ia >> groundTruth;
    return true;
}

#endif // _GROUNDTRUTH_H
