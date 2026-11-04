// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "utility.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <unordered_set>

namespace utility
{
    std::vector<VertexDistancePair> combineVerticesAndDistanceVectors(std::vector<VertexID>& vertices, std::vector<EdgeWeight>& distances)
    {
        std::vector<VertexDistancePair> verticesAndDistances;
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            verticesAndDistances.push_back(VertexDistancePair(vertices[i],distances[i]));
        }
        return verticesAndDistances;
    }
    
    std::string getFormattedUsageString(std::string commandName, std::string description) {

        // <--leading white space--><--name--><--mid white space--><--description-->
        // <--2 chars--><--15 chars max--><--mid white space-->

        RN_DEBUG_ASSERT(commandName.length()+1 < constants::MAX_COMMAND_NAME_SIZE && "Command name is too long to format");
        
        int leadingSpaceSize = 2;
        int middleSpaceSize = leadingSpaceSize + constants::MAX_COMMAND_NAME_SIZE - commandName.size();
        
        std::string leadingSpace(leadingSpaceSize,' '); // "  "
        std::string middleSpace(middleSpaceSize,' '); // "  "
        
        return leadingSpace + commandName + middleSpace + description;
    }
    
    std::string implodeVerticesAndDistances(std::vector<VertexID> kNNs, std::vector<VertexID> kNNDistances)
    {
        std::string implodedString = "", pair = "";
        for (std::size_t i = 0; i < kNNs.size(); ++i) {
            if (i != 0) {
                implodedString += ",";
            }
            pair = "[" + std::to_string(kNNs[i]) + "," + std::to_string(kNNDistances[i]) + "]";
            implodedString += pair;
        }
        return implodedString;
    }
    
    std::vector<VertexID> getPointSetFromFile(std::string queryPointsFilePath)
    {
        std::vector<VertexID> queryVertexIDs;
        VertexID queryVertexID;
        std::string networkName, setType; 
        double setDensity; 
        int setSize, setVariable;

        std::ifstream queryPointFS(queryPointsFilePath, std::ios::in);
        if (queryPointFS.is_open()) {
            // First line is query point set information
            queryPointFS >> networkName >> setType >> setDensity >> setVariable >> setSize;

            while (queryPointFS >> queryVertexID)  {
                queryVertexIDs.push_back(queryVertexID);
            }
        } else {
            std::cerr << "Cannot open point set file " << queryPointsFilePath << std::endl;
            exit(1);
        }
        queryPointFS.close();
        return queryVertexIDs;
    }
    
    std::vector<VertexID> getPointSetFromFile(std::string queryPointsFilePath, std::string& setType, double& setDensity, int& setSize, int& setVariable)
    {
        std::vector<VertexID> queryVertexIDs;
        VertexID queryVertexID;
        std::string networkName; 

        std::ifstream queryPointFS(queryPointsFilePath, std::ios::in);
        if (queryPointFS.is_open()) {
            // First line is query point set information
            queryPointFS >> networkName >> setType >> setDensity >> setVariable >> setSize;

            while (queryPointFS >> queryVertexID)  {
                queryVertexIDs.push_back(queryVertexID);
            }
            if (setSize < 0 || queryVertexIDs.size() != static_cast<std::size_t>(setSize)) {
                std::cerr << "Point set file " << queryPointsFilePath << " has " << queryVertexIDs.size() << " vertices, but its header says "
                          << setSize << std::endl;
                exit(1);
            }
        } else {
            std::cerr << "Cannot open point set file " << queryPointsFilePath << std::endl;
            exit(1);
        }
        queryPointFS.close();
        return queryVertexIDs;
    }
        
    std::string constructIndexFileName(std::string networkName, std::string indexType, const std::vector<std::string>& parameterNames, const std::vector<std::string>& parameterValues)
    {
        std::stringstream fileName;
        fileName << networkName;
        for (std::size_t i = 0; i < parameterNames.size(); ++i) {
            fileName << "_" << parameterNames[i] << "=" << parameterValues[i];
        }
        fileName << "." << indexType; // Index type is the extension
        return fileName.str();
    }

    std::string constructObjectIndexFileName(std::string objSetPrefix, std::string indexType, const std::vector<std::string>& parameterNames, const std::vector<std::string>& parameterValues)
    {
        return constructIndexFileName(objSetPrefix,indexType,parameterNames,parameterValues);
    }

    std::string constructObjectIndexFileName(std::string objSetPrefix, std::string indexType, const std::unordered_map<std::string,std::string>& parameterMap)
    {
        std::vector<std::string> parameterNames, parameterValues;
        getObjectIndexParameters(indexType,parameterMap,parameterNames,parameterValues);
        return constructIndexFileName(objSetPrefix,indexType,parameterNames,parameterValues);
    }

    // Value of an experiment parameter needed to name an object index file (the first value of a colon-separated list, as the builders use)
    static std::string getObjectIndexParameter(const std::unordered_map<std::string,std::string>& parameterMap, const std::string& key, const std::string& indexType)
    {
        auto it = parameterMap.find(key);
        if (it == parameterMap.end() || it->second == "") {
            std::cerr << "Parameter " << key << " needed to name " << indexType << " object index files was not provided!" << std::endl;
            exit(1);
        }
        return splitByDelim(it->second,':')[0];
    }

    void getObjectIndexParameters(std::string indexType, const std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string>& parameterNames, std::vector<std::string>& parameterValues)
    {
        parameterNames.clear();
        parameterValues.clear();
        std::vector<std::pair<std::string,std::string>> keys; // File name key and the experiment parameter it comes from
        if (indexType == constants::OBJ_IDX_RTREE) {
            keys = {{"branchfactor","rtree_branchfactor"}};
        } else if (indexType == constants::OBJ_IDX_LNDMRK) {
            keys = {{"landmarks","alt_numlandmarks"},{"type","alt_landmarktype"}};
        } else if (indexType == constants::OBJ_IDX_GTREE) {
            keys = {{"fanout","gtree_fanout"},{"maxleafsize","gtree_maxleafsize"}};
        } else if (indexType == constants::OBJ_IDX_COLTREE) {
            keys = {{"fanout","sltree_fanout"},{"maxleafsize","sltree_maxleafsize"},{"numlmks","sltree_numlmks"},{"lmktype","sltree_lmktype"},{"lambda","colt_listsize"}};
        }
        for (std::size_t i = 0; i < keys.size(); ++i) {
            parameterNames.push_back(keys[i].first);
            parameterValues.push_back(getObjectIndexParameter(parameterMap,keys[i].second,indexType));
        }
    }

    std::string constructObjIdxFilenamePrefix(std::string networkName, std::string setType, double setDensity, int setVariable, int setIdx)
    {
        std::stringstream fileName;
        fileName << networkName << "_" ;
        fileName << setType << "_" ;
        fileName << setDensity << "_" ;
        fileName << setVariable << "_" ;
        if (setDensity == 1) {
            // If density is one we only need to create one object index
            fileName << std::setfill('0') << std::setw(3) << 0;
        } else {
            fileName << std::setfill('0') << std::setw(3) << setIdx;
        }
        return fileName.str();
    }
    
    std::string constructObjsectSetFileName(std::string networkName, std::string setType, double setDensity, int setVariable, int setIdx)
    {
        std::stringstream fileName;
        fileName << networkName << "_" ;
        fileName << setType << "_" ;
        fileName << setDensity << "_" ;
        fileName << setVariable << "_" ;
        if (setDensity == 1) {
            // If density is one we only need to create one object set
            fileName << std::setfill('0') << std::setw(3) << 0;
        } else {
            fileName << std::setfill('0') << std::setw(3) << setIdx;
        }
        fileName << ".txt";
        return fileName.str();
    }

    std::size_t getNumObjectSets(double setDensity, std::size_t numSets)
    {
        return setDensity == 1 ? std::min(numSets,std::size_t{1}) : numSets;
    }

    std::string constructQuerySetFileName(std::string networkName, std::string type, int numVectors, int vectorLength, int distance)
    {
        std::stringstream fileName;
        fileName << networkName << "_" ;
        fileName << type << "_" ;
        fileName << numVectors << "_" ;
        fileName << vectorLength << "_" ;
        fileName << distance;
        fileName << ".txt";
        return fileName.str();
    }

    std::vector<std::string> splitByDelim(std::string input, char delimeter)
    {
        std::vector<std::string> columns;
        std::stringstream ss(input);
        std::string column;
        while (std::getline(ss, column, delimeter)) {
            columns.push_back(column);
        }
        return columns;
    }
    
    std::vector<int> getIngetersFromStringList(std::string input, char delimeter, int minValue)
    {
        std::vector<std::string> strValues = utility::splitByDelim(input,delimeter);
        std::vector<int> intValues;
        for(std::size_t i = 0; i < strValues.size(); ++i) {
            int value = std::stoi(strValues[i]);
            if (value > minValue) {
                intValues.push_back(value);
            } else {
                std::cerr << "Value in string is less than minimum value of " << minValue << std::endl;
                exit(1);    
            }
        }
        return intValues;
    }

    void writeSampleSet(std::string outputFile, std::string networkName, std::string setType, 
                                                double setDensity, int setSize, int setVariable, std::vector<VertexID>& sampleSet)
    {
        std::ofstream outputFS(outputFile, std::ios::out | std::ios::trunc);
        if (outputFS.is_open()) {
            // Write header line with object set details
            outputFS << networkName << " " << setType << " " 
                << setDensity << " " << setVariable << " " << setSize << std::endl;
            
            // Write one object VertexID per line
            for (auto it = sampleSet.begin(); it != sampleSet.end(); ++it) {
                outputFS << *it << std::endl;
            }
        } else {
            std::cerr << "Cannot open output file " << outputFile << std::endl;
        }
        outputFS.close();
    }


    // Number of untimed warm-up queries for a query list of numQueries (WARM_UP_QUERY_PERCENTAGE, at least one)
    std::size_t getNumWarmUpQueries(std::size_t numQueries)
    {
        std::size_t numWarmUpQueries = numQueries*constants::WARM_UP_QUERY_PERCENTAGE/100;
        return numWarmUpQueries > 0 ? numWarmUpQueries : 1;
    }

    // Whether an experiment runs the query list whose switch (e.g. aknn_varyingk) is in parameters: 1 runs it, any other value
    // skips it, and a missing switch runs it
    bool isQueryListOn(const std::unordered_map<std::string,std::string>& parameters, std::string listSwitch)
    {
        auto it = parameters.find(listSwitch);
        return it == parameters.end() || it->second == "1";
    }
    
    void writeVectorSetToFile(std::string outputFile, std::string networkName, std::string vectorType, 
                           double numVectors, int numElementsPerVector, std::vector<std::vector<VertexID>>& vectorSet)
    {
        std::ofstream outputFS(outputFile, std::ios::out | std::ios::trunc);
        if (outputFS.is_open()) {
            // Write header line with object set details
            outputFS << networkName << " " << vectorType << " " 
                     << numElementsPerVector << " " << numVectors << "\n";
            
            // Write one object VertexID per line
            for (std::size_t i = 0; i < vectorSet.size(); ++i) {
                for (std::size_t j = 0; j < vectorSet[i].size(); ++j) {
                    if (j != 0) {
                        outputFS << " ";
                    }
                    outputFS << vectorSet[i][j];
                }
                outputFS << "\n";
            }
        } else {
            std::cerr << "Cannot open output file " << outputFile << std::endl;
        }
        outputFS.close();
    }

    std::vector<std::vector<VertexID>> getVectorSetFromFile(std::string inputFile, std::string& vectorType, 
                                                          unsigned int& numElementsPerVector)
    {
        std::vector<std::vector<VertexID>> vectorSet;
        VertexID termID;
        std::string networkName, uselessField, line; 

        std::ifstream inputFS(inputFile, std::ios::in);
        if (inputFS.is_open()) {
            // First line is query point set information
            std::getline(inputFS, line);
            std::stringstream fss(line);
            fss >> networkName >> vectorType >> numElementsPerVector >> uselessField;

            while (std::getline(inputFS, line))  {
                std::stringstream ss(line);
                std::vector<VertexID> vector;
                while (ss >> termID) {
                    vector.push_back(termID);
                }
                if (vector.size() != numElementsPerVector) {
                    std::cout << "Vector with index " << vector.size()+1 << " has " << vector.size()
                              << " terms (should have " << numElementsPerVector << ")\n";
                }
                vectorSet.push_back(vector);
            }
        } else {
            std::cerr << "Cannot open input file " << inputFile << std::endl;
            exit(1);
        }
        inputFS.close();
        return vectorSet;
    }

    std::string removeTrailingZeroes(std::string str) {
        str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
        str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
        return str;
    }
    
}
