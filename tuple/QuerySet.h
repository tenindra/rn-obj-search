// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _QUERYSET_H
#define _QUERYSET_H

#include "../common.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

struct KnnQuery {
    VertexID queryVertex;
    unsigned int k;
    double density;
    std::string type;
    int var;
    unsigned int objSetNum;
    std::string objSetFilenamePrefix;
    int objSetSize; // Note: In some cases, this takes the value of the AVERAGE object set size (not the one associated with the objSetNum)
    KnnQuery() : queryVertex(0), k(0), density(0), var(0), objSetNum(0), objSetSize(0) {} // Required to deserialize
    KnnQuery(VertexID queryVertex, unsigned int k, double density, std::string type, int var, unsigned int objSetNum, std::string objSetFilenamePrefix, int objSetSize) :
        queryVertex(queryVertex), k(k), density(density), type(type), var(var), objSetNum(objSetNum), objSetFilenamePrefix(objSetFilenamePrefix), objSetSize(objSetSize) {}
    void printQuery() {
        std::stringstream ss;
        ss << "QVID: " << queryVertex << "\t";
        ss << "k: " << k << "\t";
        ss.precision(5);
        ss << std::fixed;
        ss << "Density: " << density << "\t";
        ss << "Type: " << type << "\n";
        ss << "Variable: " << var << "\n";
        ss << "Object Set #: " << objSetNum << "\n";
        ss << "Object Set Index Filename Prefix: " << objSetFilenamePrefix << "\n";
        ss << "Object Set Size: " << objSetSize << "\n";
        std::cout << ss.str();
    }
    // Boost Serialization (used to store ground truth queries)
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & this->queryVertex;
        ar & this->k;
        ar & this->density;
        ar & this->type;
        ar & this->var;
        ar & this->objSetNum;
        ar & this->objSetFilenamePrefix;
        ar & this->objSetSize;
    }
};

struct AknnQuery {
    std::vector<VertexID> queryVertices;
    unsigned int k;
    double density;
    int querySubnetSize;
    std::string type;
    int var;
    unsigned int objSetNum;
    std::string objSetFilenamePrefix;
    int objSetSize;
    AknnQuery() : k(0), density(0), querySubnetSize(0), var(0), objSetNum(0), objSetSize(0) {} // Required to deserialize
    AknnQuery(std::vector<VertexID> queryVertices, unsigned int k, double density, int querySubnetSize, std::string type, int var, unsigned int objSetNum, std::string objSetFilenamePrefix, int objSetSize) :
        queryVertices(queryVertices), k(k), density(density), querySubnetSize(querySubnetSize), type(type), var(var), objSetNum(objSetNum), objSetFilenamePrefix(objSetFilenamePrefix), objSetSize(objSetSize) {}
    void printQuery() {
        std::stringstream ss;
        ss << "QVIDs: ";
        for (std::size_t i = 0; i < queryVertices.size(); ++i) {
            ss << queryVertices[i] << ",";
        }
        ss << "\t";
        ss << "k: " << k << "\t";
        ss << "Query Subnet Size (%): " << querySubnetSize << "\t";
        ss.precision(5);
        ss << std::fixed;
        ss << "Density: " << density << "\t";
        ss << "Type: " << type << "\n";
        ss << "Variable: " << var << "\n";
        ss << "Object Set #: " << objSetNum << "\n";
        ss << "Object Set Index Filename Prefix: " << objSetFilenamePrefix << "\n";
        ss << "Object Set Size: " << objSetSize << "\n";
        std::cout << ss.str();
    }
    // Boost Serialization (used to store ground truth queries)
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & this->queryVertices;
        ar & this->k;
        ar & this->density;
        ar & this->querySubnetSize;
        ar & this->type;
        ar & this->var;
        ar & this->objSetNum;
        ar & this->objSetFilenamePrefix;
        ar & this->objSetSize;
    }
};



struct RangeQuery {
    VertexID queryVertex;
    unsigned int range;
    double rangePerc;
    double density;
    std::string type;
    int var;
    unsigned int objSetNum;
    std::string objSetFilenamePrefix;
    int objSetSize;
    RangeQuery() : queryVertex(0), range(0), rangePerc(0), density(0), var(0), objSetNum(0), objSetSize(0) {} // Required to deserialize
    RangeQuery(VertexID queryVertex, unsigned int range, double rangePerc, double density, std::string type, int var, unsigned int objSetNum, std::string objSetFilenamePrefix, int objSetSize) :
        queryVertex(queryVertex), range(range), rangePerc(rangePerc), density(density), type(type), var(var), objSetNum(objSetNum), objSetFilenamePrefix(objSetFilenamePrefix), objSetSize(objSetSize) {}
    void printQuery() {
        std::stringstream ss;
        ss << "QVID: " << queryVertex << "\t";
        ss << "Range: " << range << "\t";
        ss.precision(5);
        ss << std::fixed;
        ss << "Range (%): " << rangePerc << "\t";
        ss << "Density: " << density << "\t";
        ss << "Type: " << type << "\n";
        ss << "Variable: " << var << "\n";
        ss << "Object Set #: " << objSetNum << "\n";
        ss << "Object Set Index Filename Prefix: " << objSetFilenamePrefix << "\n";
        ss << "Object Set Size: " << objSetSize << "\n";
        std::cout << ss.str();
    }
    // Boost Serialization (used to store ground truth queries)
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & this->queryVertex;
        ar & this->range;
        ar & this->rangePerc;
        ar & this->density;
        ar & this->type;
        ar & this->var;
        ar & this->objSetNum;
        ar & this->objSetFilenamePrefix;
        ar & this->objSetSize;
    }
};

#endif // _QUERYSET_H


