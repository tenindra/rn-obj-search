// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _UTILITY_H
#define _UTILITY_H

#include "../common.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace utility
{
    std::string getFormattedUsageString(std::string commandName, std::string description);
    std::string implodeVerticesAndDistances(std::vector<VertexID> kNNs,std::vector<VertexID> kNNDistances);
    std::vector<std::string> splitByDelim(std::string input, char delimeter);
    std::vector<int> getIngetersFromStringList(std::string input, char delimeter, int minValue);
    std::size_t getNumWarmUpQueries(std::size_t numQueries);
    bool isQueryListOn(const std::unordered_map<std::string,std::string>& parameters, std::string listSwitch);
    template <typename STLCollection>
    void releaseSTLCollection(STLCollection& collection) {
        STLCollection tmp;
        collection.swap(tmp);
    }

    // Query results
    std::vector<VertexDistancePair> combineVerticesAndDistanceVectors(std::vector<VertexID>& vertices, std::vector<EdgeWeight>& distances);

    // Reading/Writing Point Sets
    std::vector<VertexID> getPointSetFromFile(std::string queryPointsFilePath);
    std::vector<VertexID> getPointSetFromFile(std::string queryPointsFilePath, std::string& setType, double& setDensity, int& setSize, int& setVariable);
    void writeSampleSet(std::string outputFile, std::string networkName, std::string setType, 
                        double setDensity, int setSize, int setVariable, std::vector<VertexID>& sampleSet);
    void writeVectorSetToFile(std::string outputFile, std::string networkName, std::string vectorType, 
                              double numVectors, int numElementsPerVector, std::vector<std::vector<VertexID>>& vectorSet);
    std::vector<std::vector<VertexID>> getVectorSetFromFile(std::string inputFile, std::string& vectorType, 
                                                          unsigned int& numElementsPerVector);
    
    // Index Name Related
    std::string constructIndexFileName(std::string networkName, std::string indexType, const std::vector<std::string>& parameterNames, const std::vector<std::string>& parameterValues);
    // Object index file name: the object set prefix, the parameters the index was built with (see getObjectIndexParameters) and the index type as extension
    std::string constructObjectIndexFileName(std::string objSetPrefix, std::string indexType, const std::vector<std::string>& parameterNames, const std::vector<std::string>& parameterValues);
    std::string constructObjectIndexFileName(std::string objSetPrefix, std::string indexType, const std::unordered_map<std::string,std::string>& parameterMap);
    // The parameters recorded in an object index file name, taken from the experiment parameters: the road network index it was
    // built from (gtree_*, alt_*, sltree_*), the R-tree branch factor or the COL-Tree list size (lambda); none for nvd and txt
    void getObjectIndexParameters(std::string indexType, const std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string>& parameterNames, std::vector<std::string>& parameterValues);
    std::string constructObjsectSetFileName(std::string networkName, std::string setType, double setDensity, int setVariable, int setIdx);
    std::string constructQuerySetFileName(std::string networkName, std::string type, int numVectors, int vectorLength, int distance);
    std::string constructObjIdxFilenamePrefix(std::string networkName, std::string setType, double setDensity, int setVariable, int setIdx);
    // The number of object sets of a density that the object set and index builders write: a set of density one holds every vertex,
    // so there is only one (set index 000 in the file names above). The query lists still run every set index at every density (at
    // density one, all read that set), so each list times the same number of queries
    std::size_t getNumObjectSets(double setDensity, std::size_t numSets);

    // Experimental Results Output Related
    std::string removeTrailingZeroes(std::string str);
}

#endif // _UTILITY_H
