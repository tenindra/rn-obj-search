// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "GenerateSampleSetCommand.h"

#include "../processing/Graph.h"
#include "../processing/SetGenerator.h"
#include "../utility/StopWatch.h"
#include "../utility/utility.h"
#include "../utility/serialization.h"

#include <cmath>
#include <iostream>

void GenerateSampleSetCommand::execute(int argc, char* argv[])
{
    std::string setType = "";
    double setDensity = 0;
    int setVariable = 0;
    int setSize = 0;
    std::string bgrFilePath = "";
    std::string sourceBgrFilePath = "";
    std::string outputFilePath = "";
    std::string status = "ignore";
    int outdegree = 0;
    std::string setSeed = "";

    /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "t:i:o:d:n:s:e:v:a:r:")) != -1) {
        switch (opt) {
            case 't':
                setType = optarg;
                break;
            case 'r':
                setSeed = optarg;
                break;
            case 'i':
                bgrFilePath = optarg;
                break;
            case 'a':
                sourceBgrFilePath = optarg;
                break;
            case 'o':
                outputFilePath = optarg;
                break;
            case 'd':
                setDensity = std::stod(optarg);
                break;
            case 'v':
                setVariable = std::stoi(optarg);
                break;
            case 'n':
                setSize = std::stoi(optarg);
                break;
            case 's':
                status = optarg;
                break;
            case 'e':
                outdegree = std::stoi(optarg);
                break;
            default:
                std::cerr << "Unknown option(s) provided!\n\n";
                showCommandUsage(argv[0]);
                exit(1);
        }
    }

    // Validate Command Line Arguments
    if (argc == 5 && setType != "") {
        // This is 5 so that user can just enter method to find out what parameters are required for method
        this->showSetTypeUsage(setType,argv[0]);
        exit(1);
    }    
    
    if (setType == "" || bgrFilePath == "" || outputFilePath == "") {
        std::cerr << "Invalid argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }

    if (sourceBgrFilePath != "" && setType != constants::RAND_OBJ_SET) {
        std::cerr << "Only " << constants::RAND_OBJ_SET << " sets can be copied from another dataset (-a)!\n";
        exit(1);
    }

    /*
     * Serialize Graph from Binary File
     */
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);

    /*
     * Sample Graph to Generate Object Set
     */          
    int numVertices = graph.getNumVertices();
    SetGenerator sg;
    StopWatch sw;

    if (setType == constants::RAND_OBJ_SET) {
        /*
         * Create Randomly Selected Object Set
         */
        
        // Process Set Type Command-Line Arguments
        if (argc < 11) {
            // Arguments: -c sample_set -t random -i <binary graph file> -o <output file> 
            // (at least one of -d <density> or -n <size>) [-s <ignore|include|exclude> [-e <outdegree of vertices to include or exclude>]] [-v <set variable - currently not used>]
            std::cerr << "Too few arguments!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        } else if ((setDensity <= 0 || setDensity >= 1) && setSize <= 0) {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        }
        
        // Determine how many objects must be generated if density was provided
        if (setSize == 0) {
            setSize = std::ceil(numVertices*setDensity);
        } else {
            setDensity = setSize/numVertices;
        }

        if (sourceBgrFilePath != "") {
            // Copy the set generated for another dataset with the same vertex IDs and neighbours: the file in the same directory
            // whose name starts with the other dataset's network name instead of this one's
            Graph sourceGraph = serialization::getIndexFromBinaryFile<Graph>(sourceBgrFilePath);
            std::string difference;
            if (!graph.hasSameNeighbours(sourceGraph,difference)) {
                std::cerr << "Cannot copy the sets of " << sourceGraph.getNetworkName() << " to " << graph.getNetworkName() << ": " << difference << std::endl;
                exit(1);
            }
            std::string sourceNetworkName = sourceGraph.getNetworkName();
            std::string networkNamePrefix = graph.getNetworkName() + "_";
            std::size_t fileNamePos = outputFilePath.find_last_of('/') + 1;
            if (outputFilePath.compare(fileNamePos,networkNamePrefix.size(),networkNamePrefix) != 0) {
                std::cerr << "Cannot copy a set to " << outputFilePath << ", as its file name does not start with " << networkNamePrefix << std::endl;
                exit(1);
            }
            std::string sourceFilePath = outputFilePath.substr(0,fileNamePos) + sourceNetworkName + "_" + outputFilePath.substr(fileNamePos+networkNamePrefix.size());
            std::string sourceSetType;
            double sourceSetDensity;
            int sourceSetSize, sourceSetVariable;
            std::vector<VertexID> sampleSet = utility::getPointSetFromFile(sourceFilePath,sourceSetType,sourceSetDensity,sourceSetSize,sourceSetVariable);
            if (sourceSetType != setType || sourceSetSize != setSize || sampleSet.size() != static_cast<std::size_t>(setSize)) {
                std::cerr << sourceFilePath << " is not a " << setType << " set of size " << setSize << std::endl;
                exit(1);
            }
            utility::writeSampleSet(outputFilePath,graph.getNetworkName(),sourceSetType,sourceSetDensity,sourceSetSize,sourceSetVariable,sampleSet);
            std::cout << "Random sample set of size " << sampleSet.size() << " successfully copied from " << sourceFilePath << "!" << std::endl;
            return;
        }

        unsigned int seed = SetGenerator::parseSeed(setSeed,"-r");

        // Generat set and write to file
        sw.start();
        std::vector<VertexID> sampleSet;
        std::string message = "";
        int statusValue;
        if (status == "include") {
            if (outdegree <= 0) {
                std::cerr << "Invalid out degree, must be greater than zero\n\n";
                this->showSetTypeUsage(setType,argv[0]);
                exit(1);
            }
            statusValue = 1;
            message = " including only vertices with degree " + std::to_string(outdegree);
        } else if (status == "exclude") {
            if (outdegree <= 0) {
                std::cerr << "Invalid out degree, must be greater than zero\n\n";
                this->showSetTypeUsage(setType,argv[0]);
                exit(1);
            }
            statusValue = 0;
            message = " excluding all vertices with degree " + std::to_string(outdegree);
        } else if (status == "ignore") {
            statusValue = -1;
            outdegree = 0; // Not a sampling criterion, so it does not select the engine either
            message = " including vertices with any degree";
        } else {
            std::cerr << "-s option argument does not have a valid value, it must be either \"include\", \"exclude\" or \"ignore\" if provided at all\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        }
        // The set's engine is seeded with its size and sampling criteria
        std::mt19937 engine = SetGenerator::getEngine(seed,SetGenerator::SAMPLE_SET,{static_cast<unsigned int>(setSize),static_cast<unsigned int>(statusValue+1),static_cast<unsigned int>(outdegree)});
        sampleSet = sg.generateRandomSampleSet(graph,setSize,engine,statusValue,outdegree);
        sw.stop();
        utility::writeSampleSet(outputFilePath,graph.getNetworkName(),setType,setDensity,setSize,setVariable,sampleSet);
        std::cout << "Random sample set of size " << sampleSet.size() << message << " successfully created in " << sw.getTimeMs() << "ms!" << std::endl;
        
    } else {
        std::cerr << "Invalid set type!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }

}

void GenerateSampleSetCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::SAMPLE_SET_CMD + " -t <set type>\n\n"
              << "Set Type Options:\n"
              << utility::getFormattedUsageString(constants::RAND_OBJ_SET,"Randomly sample graph vertices") << "\n";
}

void GenerateSampleSetCommand::showSetTypeUsage(std::string setType, std::string programName)
{
    if (setType == constants::RAND_OBJ_SET) {
        std::cerr << "Usage: " << programName << " -c " + constants::SAMPLE_SET_CMD + 
                     " -t " +constants::RAND_OBJ_SET + " -i <binary graph file>"
                  << " -o <set output file> (one option and argument from group 1)\n "
                  << "-r <set seed> [-s <ignore|include|exclude> [-e <outdegree of vertices to include or exclude>]]\n"
                  << "[-a <binary graph of another dataset>]\n\n"
                  << "Group 1:\n"
                  << utility::getFormattedUsageString("-d <fraction>","Density of sampled set") << "\n"
                  << utility::getFormattedUsageString("-n <size>","Size of sampled set") << "\n\n"
                  << "The vertices are drawn from the set seed (an integer), the set size and the sampling criteria, so the same\n"
                  << "seed and arguments give the same set of a graph. With -a, the seed is not needed and the set is copied from the other dataset instead of sampled: from the file in the same directory\n"
                  << "whose name starts with the other dataset's network name instead of this one's. The other dataset must have\n"
                  << "the same vertex IDs and neighbours (edge weights may differ)\n";
    } else {
        this->showCommandUsage(programName);
    }
}
