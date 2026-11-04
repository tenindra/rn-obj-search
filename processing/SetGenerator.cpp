// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "SetGenerator.h"

#include <random>
#include <unordered_set>
#include <iostream>
#include <limits>

// Status -1 indicates indicates ignore outdegrees, status 0 indicates exclude given degree, status 1 indicates include only given outdegree.
// Draws from the given engine, so a caller that seeds it gets the same sample every time
std::vector<VertexID> SetGenerator::generateRandomSampleSet(Graph& graph, unsigned int numObjects, std::mt19937& engine, int status, int outdegree)
{
    std::vector<bool> testedVertex(graph.getNumVertices(),false);
    std::vector<VertexID> sampleSet;

    int n = graph.getNumVertices()-1;

    std::uniform_int_distribution<int> uniform_dist(0, n);

    VertexID randomVertexID;

    // Randomly select vertices until object set is complete
    int numEdges;
    std::size_t numVerticesTested = 0;
    while (sampleSet.size() < numObjects) {
        randomVertexID = uniform_dist(engine);
        if (!testedVertex[randomVertexID]) {
            if (status == -1) {
                testedVertex[randomVertexID] = true;
                sampleSet.push_back(randomVertexID);
            } else if (status == 0) {
                numEdges = graph.getEdgeListSize(randomVertexID)-graph.getEdgeListStartIndex(randomVertexID);
                testedVertex[randomVertexID] = true; // Whether we pick this vertex or not we don't need to test it again
                if (numEdges != outdegree) {
                    // We exclude the given outdegree (assumes there are enough vertices
                    // with other outdegrees to satisfy numObjects)
                    sampleSet.push_back(randomVertexID);
                }
            } else if (status == 1) {
                numEdges = graph.getEdgeListSize(randomVertexID)-graph.getEdgeListStartIndex(randomVertexID);
                testedVertex[randomVertexID] = true; // Whether we pick this vertex or not we don't need to test it again
                if (numEdges == outdegree) {
                    // We only include the given outdegree (assumes there are enough vertices
                    // with outdegree to satisfy numObjects)
                    sampleSet.push_back(randomVertexID);
                }
            } else {
                std::cerr << "Invalid status provided" << std::endl;
                std::exit(1);
            }
            ++numVerticesTested;
            if (numVerticesTested == static_cast<std::size_t>(graph.getNumVertices()) && sampleSet.size() < numObjects) {
                std::cerr << "Only " << sampleSet.size() << " vertices satisfy the sample criteria but " << numObjects << " are required" << std::endl;
                std::exit(1);
            }
        }
    }
 
    return sampleSet;
}

// Draws from the given engine, so a caller that seeds it gets the same sample every time
std::vector<VertexID> SetGenerator::generateRandomSampleSet(int numVertices, unsigned int numObjects, std::mt19937& engine)
{
    std::unordered_set<VertexID> verticesInSet;
    std::vector<VertexID> sampleSet;
    
    if (numObjects > static_cast<unsigned int>(numVertices)) {
        std::cerr << "Cannot sample " << numObjects << " unique vertices from " << numVertices << " vertices" << std::endl;
        std::exit(1);
    }

    int n = numVertices-1;

    std::uniform_int_distribution<int> uniform_dist(0, n);    

    VertexID randomVertexID;

    // Randomly select vertices until object set is complete
    while (sampleSet.size() < numObjects) {
        randomVertexID = uniform_dist(engine);
        if (verticesInSet.find(randomVertexID) == verticesInSet.end()) {
            verticesInSet.insert(randomVertexID);
            sampleSet.push_back(randomVertexID);
        }
    }

    return sampleSet;
}

std::mt19937 SetGenerator::getEngine(unsigned int setSeed, unsigned int setKind, std::vector<unsigned int> setParameters)
{
    std::vector<unsigned int> seedValues = {setSeed,setKind};
    seedValues.insert(seedValues.end(),setParameters.begin(),setParameters.end());
    std::seed_seq seed(seedValues.begin(),seedValues.end());
    return std::mt19937(seed);
}

unsigned int SetGenerator::parseSeed(std::string seedValue, std::string seedOption)
{
    if (seedValue == "") {
        std::cerr << "No set seed given (" << seedOption << "): the sets are drawn from it, so that the same seed gives the same sets" << std::endl;
        std::exit(1);
    }
    // At most 10 digits, so that the conversion cannot overflow before the range check
    if (seedValue.find_first_not_of("0123456789") != std::string::npos || seedValue.size() > 10
        || std::stoull(seedValue) > std::numeric_limits<unsigned int>::max()) {
        std::cerr << "Invalid set seed " << seedValue << " (" << seedOption << "): it must be an integer from 0 to "
                  << std::numeric_limits<unsigned int>::max() << std::endl;
        std::exit(1);
    }
    return static_cast<unsigned int>(std::stoull(seedValue));
}
