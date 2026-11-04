// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "GenerateVertexVectorCommand.h"

#include "../processing/Graph.h"
#include "../processing/SetGenerator.h"

#include "../utility/StopWatch.h"
#include "../utility/utility.h"
#include "../utility/serialization.h"

#include <random>
#include <set>
#include <deque>

void GenerateVertexVectorCommand::execute(int argc, char* argv[])
{
    std::string vectorType = "";
    int numVectors = 0;
    int numElementsPerVector = 0;
    std::string subnetworkSizes = ""; // in percentages
    std::string outputFilePath = "";
    std::string outputFilePathPrefix = "";
    std::string bgrFilePath = "";
    std::string sourceBgrFilePath = "";
    std::string setSeed = "";


    /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "t:n:v:o:b:p:s:a:r:")) != -1) {
        switch (opt) {
            case 't':
                vectorType = optarg;
                break;
            case 'r':
                setSeed = optarg;
                break;
            case 'n':
                numVectors = std::stoi(optarg);
                break;
            case 'v':
                numElementsPerVector = std::stoi(optarg);
                break;
            case 's':
                subnetworkSizes = optarg;
                break;
            case 'o':
                outputFilePath = optarg;
                break;
            case 'b':
                bgrFilePath = optarg;
                break;
            case 'a':
                sourceBgrFilePath = optarg;
                break;
            case 'p':
                outputFilePathPrefix = optarg;
                break;
            default:
                std::cerr << "Unknown option(s) provided!\n\n";
                showCommandUsage(argv[0]);
                exit(1);
        }
    }

    // Validate Command Line Arguments
    if (argc == 5 && vectorType != "") {
        // This is 5 so that user can just enter method to find out what parameters are required for method
        this->showSetTypeUsage(vectorType,argv[0]);
        exit(1);
    }    
    
    if (argc < 11) {
        // Arguments: -c sample_set -t random -n <num vectors> -v <num terms per vector> -i <binary keyword index file>
        std::cerr << "Too few arguments!\n\n";
        this->showSetTypeUsage(vectorType,argv[0]);
        exit(1);
    }
    
    if (vectorType == "" || numVectors <= 0 || numElementsPerVector <= 0) {
        std::cerr << "Invalid argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }

    if (sourceBgrFilePath != "" && vectorType != SUBNET) {
        std::cerr << "Only " << SUBNET << " sets can be copied from another dataset (-a)!\n";
        exit(1);
    }

    /*
     * Serialize KeywordIndex from Binary File
     */
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);

    if (vectorType == RANDOM) {

        // Process Set Type Command-Line Arguments
        if (argc < 13) {
            // Arguments: -o <output file>
            std::cerr << "Too few arguments!\n\n";
            this->showSetTypeUsage(vectorType,argv[0]);
            exit(1);
        } else if (outputFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showSetTypeUsage(vectorType,argv[0]);
            exit(1);
        }
        
        /*
         * Create Randomly Selected Vertex Vectors
         */

        if (numElementsPerVector > graph.getNumVertices()) {
            std::cerr << "Cannot create vectors of " << numElementsPerVector << " unique vertices from " << graph.getNumVertices() << " vertices!\n";
            exit(1);
        }

        unsigned int seed = SetGenerator::parseSeed(setSeed,"-r");

        StopWatch sw;
        sw.start();
        int n = graph.getNumVertices()-1;
        std::mt19937 e1 = SetGenerator::getEngine(seed,SetGenerator::VERTEX_VECTOR_SET,{static_cast<unsigned int>(numVectors),static_cast<unsigned int>(numElementsPerVector)});
        std::uniform_int_distribution<int> uniform_dist(0,n);
        std::vector<std::vector<VertexID>> vectorSet;
        
        for (int i = 0; i < numVectors; ++i) {
            VertexID randomVertexID;
            std::vector<VertexID> vertexVector;
            std::set<VertexID> vertexIDsInVector;
            while ( vertexVector.size() < static_cast<std::size_t>(numElementsPerVector)) {
                randomVertexID = uniform_dist(e1);
                if (vertexIDsInVector.find(randomVertexID) == vertexIDsInVector.end()) {
                    // If term not already in vector, append term to it
                    vertexVector.push_back(randomVertexID);
                    vertexIDsInVector.insert(randomVertexID);
                }
            }
            vectorSet.push_back(vertexVector);
            
        }
        sw.stop();
        
        // Write term vector to file
        utility::writeVectorSetToFile(outputFilePath,graph.getNetworkName(),vectorType,numVectors,numElementsPerVector,vectorSet);

        std::cout << vectorSet.size() << " random query vertex vector sets of size " << numElementsPerVector << " successfully created in " << sw.getTimeMs() << "ms!" << std::endl;

    } else if (vectorType == SUBNET) {

        // Process Set Type Command-Line Arguments
        if (argc < 15) {
            // Arguments: -s <list of subnetwork size percentages as integers>
            std::cerr << "Too few arguments!\n\n";
            this->showSetTypeUsage(vectorType,argv[0]);
            exit(1);
        } else if (subnetworkSizes == "" || outputFilePathPrefix == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showSetTypeUsage(vectorType,argv[0]);
            exit(1);
        }

        if (sourceBgrFilePath != "") {
            // Copy the sets generated for another dataset with the same vertex IDs and neighbours (the subnetworks are grown by
            // BFS, so they do not depend on the edge weights); every set is read and checked before any is written
            Graph sourceGraph = serialization::getIndexFromBinaryFile<Graph>(sourceBgrFilePath);
            std::string difference;
            if (!graph.hasSameNeighbours(sourceGraph,difference)) {
                std::cerr << "Cannot copy the sets of " << sourceGraph.getNetworkName() << " to " << graph.getNetworkName() << ": " << difference << std::endl;
                exit(1);
            }
            std::string sourceNetworkName = sourceGraph.getNetworkName();
            std::vector<std::string> strSubnetworkSizes = utility::splitByDelim(subnetworkSizes,',');
            std::vector<std::vector<std::vector<VertexID>>> vectorSets;
            for (std::size_t i = 0; i < strSubnetworkSizes.size(); ++i) {
                std::string sourceFilePath = outputFilePathPrefix + "/query_sets/" + utility::constructQuerySetFileName(sourceNetworkName,vectorType,numVectors,numElementsPerVector,std::stoi(strSubnetworkSizes[i]));
                std::string sourceVectorType;
                unsigned int sourceNumElementsPerVector;
                vectorSets.push_back(utility::getVectorSetFromFile(sourceFilePath,sourceVectorType,sourceNumElementsPerVector));
                bool validVectorSet = sourceVectorType == vectorType && sourceNumElementsPerVector == static_cast<unsigned int>(numElementsPerVector)
                                      && vectorSets[i].size() == static_cast<std::size_t>(numVectors);
                for (std::size_t j = 0; j < vectorSets[i].size(); ++j) {
                    if (vectorSets[i][j].size() != static_cast<std::size_t>(numElementsPerVector)) {
                        validVectorSet = false;
                    }
                }
                if (!validVectorSet) {
                    std::cerr << sourceFilePath << " does not hold " << numVectors << " " << vectorType << " vertex vectors of size " << numElementsPerVector << std::endl;
                    exit(1);
                }
            }
            for (std::size_t i = 0; i < strSubnetworkSizes.size(); ++i) {
                outputFilePath = outputFilePathPrefix + "/query_sets/" + utility::constructQuerySetFileName(graph.getNetworkName(),vectorType,numVectors,numElementsPerVector,std::stoi(strSubnetworkSizes[i]));
                utility::writeVectorSetToFile(outputFilePath,graph.getNetworkName(),vectorType,numVectors,numElementsPerVector,vectorSets[i]);
            }
            std::cout << "Subnetwork vertex vector sets of size " << numElementsPerVector << " for " << strSubnetworkSizes.size() << " subnetwork sizes successfully copied from " << sourceNetworkName << "!" << std::endl;
            return;
        }

        unsigned int seed = SetGenerator::parseSeed(setSeed,"-r");

        StopWatch sw;
        sw.start();

        /*
         * Create subnetwork of the graph around random vertex for each subnetwork size
         */

        std::vector<std::string> strSubnetworkSizes = utility::splitByDelim(subnetworkSizes,',');

        int n = graph.getNumVertices()-1;
        std::uniform_int_distribution<int> uniform_dist1(0,n);

        VertexID centreVertexID;

        for(std::size_t i = 0; i < strSubnetworkSizes.size(); ++i) {

            int subnetworkSize = std::stoi(strSubnetworkSizes[i]);
            if (subnetworkSize <= 0 || subnetworkSize > 100) {
                std::cerr << "Invalid subnetwork size percentage in list provided, must be greater than 0 and at most 100!\n";
                exit(1);
            }

            // Each subnetwork size has its own engine, so its centre vertex and vectors do not depend on the other sizes in the list
            std::mt19937 e1 = SetGenerator::getEngine(seed,SetGenerator::SUBNET_VECTOR_SET,{static_cast<unsigned int>(numVectors),static_cast<unsigned int>(numElementsPerVector),static_cast<unsigned int>(subnetworkSize)});

            std::size_t requiredNumVertices = subnetworkSize/100.0 * graph.getNumVertices();

            // Find all vertex IDs in subnetwork
            std::vector<VertexID> subnetVertices;
            
            // Randomly select a centre vertex
            centreVertexID = uniform_dist1(e1);
            
            // Conduct a BFS until required subnetwork size is method
            std::deque<VertexID> queue;
            std::vector<bool> visited(graph.getNumVertices());
            
            queue.push_back(centreVertexID);
            
            VertexID currentVertex, adjVertex;
            int adjListStart, nextAdjListStart;
            
            while (queue.size() > 0) {
                currentVertex = queue.front();
                queue.pop_front();
                if (!visited[currentVertex]) {
                    visited[currentVertex] = true;
                    subnetVertices.push_back(currentVertex);
                    if (subnetVertices.size() == requiredNumVertices) {
                        break;
                    }
                    
                    // Inspect each neighbour and update pqueue using edge weights
                    adjListStart = graph.getEdgeListStartIndex(currentVertex);
                    nextAdjListStart = graph.getEdgeListSize(currentVertex);

                    for (int i = adjListStart; i < nextAdjListStart; ++i) {
                        adjVertex = graph.edges[i].first;
                        if (!visited[adjVertex]) {
                            queue.push_back(adjVertex);
                        }
                    }
                }
            }

            // Check enough unique vertex IDs available to fill at least one vertex vector
            if (subnetVertices.size() > static_cast<std::size_t>(numElementsPerVector)) {
                // Randomly select vertices from within the subnetwork to fill vector sets
                int maxIdx = subnetVertices.size()-1;
                std::uniform_int_distribution<int> uniform_dist2(0,maxIdx);
                std::vector<std::vector<VertexID>> vectorSet;
                
                for (int i = 0; i < numVectors; ++i) {
                    VertexID randomVertexIdx, randomVertexID;
                    std::vector<VertexID> vertexVec;
                    std::set<VertexID> vertexIDsInVec;
                    while (vertexVec.size() < static_cast<std::size_t>(numElementsPerVector)) {
                        randomVertexIdx = uniform_dist2(e1);
                        randomVertexID = subnetVertices[randomVertexIdx];
                        if (vertexIDsInVec.find(randomVertexID) == vertexIDsInVec.end()) {
                            // If ID not already in vector, append ID to it
                            vertexVec.push_back(randomVertexID);
                            vertexIDsInVec.insert(randomVertexID);
                        }
                    }
                    vectorSet.push_back(vertexVec);
                }
                
                // Write term vector to file
                outputFilePath = outputFilePathPrefix + "/query_sets/" + utility::constructQuerySetFileName(graph.getNetworkName(),vectorType,numVectors,numElementsPerVector,subnetworkSize);
                utility::writeVectorSetToFile(outputFilePath,graph.getNetworkName(),vectorType,numVectors,numElementsPerVector,vectorSet);
            } else {
                std::cerr << "Not enough vertices in bult subnetwork for even one query vertex vector!\n";
                std::exit(1);                
            }
        }
        
        sw.stop();
        
        std::cout << "Subnetwork vertex vector sets of size " << numElementsPerVector << " for " << strSubnetworkSizes.size() << " subnework sizes successfully created in " << sw.getTimeMs() << "ms!" << std::endl;
        
    } else {
        std::cerr << "Invalid set type!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }

}

void GenerateVertexVectorCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::VERTEX_VECTOR_CMD + " -t <set type>\n\n"
              << "Set Type Options:\n"
              << utility::getFormattedUsageString(RANDOM,"Select vertices randomly") << "\n"
              << utility::getFormattedUsageString(SUBNET,"Select vertices by subgraph of some percentage") << "\n";
}

void GenerateVertexVectorCommand::showSetTypeUsage(std::string vectorType, std::string programName)
{
    if (vectorType == RANDOM) {
        std::cerr << "Usage: " << programName << " -c " + constants::VERTEX_VECTOR_CMD + 
                     " -t " + RANDOM + " -n <num vectors>\n-v <num terms per vector> -b <binary graph index file> -o <output file> -r <set seed>\n\n"
                  << "The vertices are drawn from the set seed (an integer) and the numbers of vectors and terms, so the same seed and\n"
                  << "arguments give the same vectors of a graph\n";
    } else if (vectorType == SUBNET) {
        std::cerr << "Usage: " << programName << " -c " + constants::VERTEX_VECTOR_CMD +
                     " -t " + SUBNET + " -n <num vectors>\n-v <num terms per vector> -s <list of subnetwork size percentages as integers> -b <binary graph index file> -p <output file path prefix>\n"
                  << "-r <set seed> [-a <binary graph of another dataset>]\n\n"
                  << "Each subnetwork size's centre vertex and vectors are drawn from the set seed (an integer), the numbers of vectors\n"
                  << "and terms and the subnetwork size, so the same seed and arguments give the same vectors of a graph.\n"
                  << "With -a, the seed is not needed and the sets are copied from the other dataset's sets of the same sizes instead of\n"
                  << "generated. The other dataset must have the same vertex IDs and neighbours (edge weights may differ)\n";
    } else {
        this->showCommandUsage(programName);
    }
}
