// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "LMKExperimentsCommand.h"

#include "../processing/INE.h"
#include "../processing/IER.h"
#include "../processing/SetGenerator.h"
#include "../processing/ILBR.h"
#include "../tuple/IndexTuple.h"
#include "../tuple/ObjectIndexTuple.h"
#include "../tuple/KnnQueryTuple.h"
#include "../tuple/QuerySet.h"
#include "../common.h"
#include "../utility/StopWatch.h"
#include "../utility/Statistics.h"
#include <fstream>
#include "../utility/utility.h"
#include "../utility/serialization.h"

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <climits>
#include <iomanip>
#include <sstream>
#include <omp.h>

void LMKExperimentsCommand::execute(int argc, char* argv[])
{
    std::string experiment = "";
    std::string bgrFilePath = "";
    std::string filePathPrefix = "";
    std::string statsOutputFile = "";
    std::string parameters = "";
    std::string objDensities = "";
    std::string objTypes = "";
    std::string objVariable = "";
    unsigned int numSets = 0;
    std::string queryVertexFile = "";
    std::string kValues = "";
    std::string rwObjSetNames = "";
    unsigned int minLambda = 0;
    std::string fanouts = "";
    std::string numLandmarks = "";
    std::string landmarkTypes = "";
    std::string lambdas = "";
    std::string sourceBgrFilePath = "";
    
    /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "e:g:p:f:s:n:d:t:q:k:v:r:u:w:x:y:z:a:")) != -1) {
        switch (opt) {
            case 'e':
                experiment = optarg;
                break;
            case 'g':
                bgrFilePath = optarg;
                break;
            case 'p':
                parameters = optarg;
                break;
            case 'f':
                filePathPrefix = optarg;
                break;
            case 's':
                statsOutputFile = optarg;
                break;
            case 'n':
                numSets = std::stoul(optarg);
                break;
            case 'd':
                objDensities = optarg;
                break;
            case 't':
                objTypes = optarg;
                break;
            case 'q':
                queryVertexFile = optarg;
                break;
            case 'k':
                kValues = optarg;
                break;
            case 'v':
                objVariable = optarg;
                break;
            case 'r':
                rwObjSetNames = optarg;
                break;
            case 'u':
                minLambda = std::stoul(optarg);
                break;
            case 'w':
                fanouts = optarg;
                break;
            case 'x':
                numLandmarks = optarg;
                break;
            case 'y':
                landmarkTypes = optarg;
                break;
            case 'z':
                lambdas = optarg;
                break;
            case 'a':
                sourceBgrFilePath = optarg;
                break;
            default:
                std::cerr << "Unknown option(s) provided!\n\n";
                showCommandUsage(argv[0]);
                exit(1);
        }
    }

    // Validate Command Line Arguments
    if (argc == 5) {
        // This is 5 so that user can just enter method to find out what parameters are required for method
        this->showPhaseUsage(experiment,argv[0]);
        exit(1);
    } 
    
    if (experiment == "") {
        std::cerr << "Invalid argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
    if (experiment == this->BUILD_RN_INDEXES) {
        if (argc < 13) {
            // Arguments: -g <binary graph file> -p <parameters> -f <index output file path prefix> -s <stats output file>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->buildIndexes(bgrFilePath,parameters,filePathPrefix,statsOutputFile);
        this->buildExternalIndexes(bgrFilePath,parameters,filePathPrefix,statsOutputFile);
        
    } else if (experiment == this->BUILD_PARTITIONS) {
        if (argc < 13) {
            // Arguments: -g <binary graph file> -p <parameters> -f <index output file path prefix> -s <stats output file>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->buildPartitions(bgrFilePath,parameters,filePathPrefix,statsOutputFile);
        
    } else if (experiment == this->RUN_SDL_EXPERIMENT) {
        if (argc < 13) {
            // Arguments: -g <binary graph file> -p <parameters> -f <index output file path prefix> -s <stats output file>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->runSDLExperiment(bgrFilePath,parameters,filePathPrefix,statsOutputFile);
        
    } else if (experiment == this->BUILD_OBJ_INDEXES) {
        if (argc < 21) {
            // Arguments: -g <binary graph file>  -p <parameters> -n <num sets> -d <list of object densities> 
            // -t <list of object types> -v <list of some object variable> -f <index output file path prefix> 
            // -s <stats output file>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == ""
            || numSets == 0 || objDensities == "" || objTypes == "" || objVariable == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->buildObjectIndexes(bgrFilePath,parameters,numSets,objDensities,objTypes,objVariable,filePathPrefix,statsOutputFile,sourceBgrFilePath);
        
    } else if (experiment == this->BUILD_RW_OBJ_INDEXES) {
        if (argc < 15) {
            // Arguments: -g <binary graph file>  -p <parameters> -n <num sets> -d <list of object densities> 
            // -t <list of object types> -v <list of some object variable> -f <index output file path prefix> 
            // -s <stats output file> -r <rw object set file names>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == "" || rwObjSetNames == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->buildRealWorldObjectIndexes(bgrFilePath,parameters,filePathPrefix,statsOutputFile,rwObjSetNames);
        
    } else if (experiment == this->RUN_KNN_QUERIES) {
        if (argc < 23) {
            // Arguments: -g <binary graph file> -q <query vertex file> -k <k values> -p <parameters> -n <num sets> -d <list of object densities> -t <list of object types> 
            // -v <list of some object variable> -f <index output file path prefix> [-r <rw object set file names>] 
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || numSets == 0 || objDensities == "" ||
            objTypes == "" || queryVertexFile == "" || kValues == "" || objVariable == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->runQueries2(bgrFilePath,queryVertexFile,kValues,parameters,numSets,objDensities,objTypes,objVariable,rwObjSetNames,filePathPrefix);

    } else if (experiment == this->RUN_KNN_SINGLE_QUERIES) {
        if (argc < 25) {
            // Arguments: -g <binary graph file> -q <query vertex file> -k <k values> -p <parameters> -n <num sets> -d <list of object densities> -t <list of object types> 
            // -v <list of some object variable> -f <index output file path prefix> -s <stats output file>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == "" || numSets == 0 || objDensities == "" || 
            objTypes == "" || queryVertexFile == "" || kValues == "" || objVariable == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->runSingleQueries2(bgrFilePath,queryVertexFile,kValues,parameters,numSets,objDensities,objTypes,objVariable,filePathPrefix,statsOutputFile);

    } else if (experiment == this->BUILD_PARAMETER_TEST_RN_INDEXES) {
        if (argc < 15) {
            // Arguments: -g <binary graph file> -f <index output file path prefix> 
            // -u <minimum lambda> -w <fanouts> -x <num landmarks> -y <landmark types>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || filePathPrefix == "" || minLambda == 0 || fanouts == "" ||
            numLandmarks == "" || landmarkTypes == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->buildSULTreeParameterIndexes(bgrFilePath,parameters,filePathPrefix,fanouts,numLandmarks,landmarkTypes,minLambda);

    } else if (experiment == this->BUILD_PARAMETER_TEST_OBJ_INDEXES) {
        if (argc < 27) {
            // Arguments: -g <binary graph file> -p <parameters> -f <index output file path prefix>
            // -u <minimum lambda> -w <fanouts> -x <num landmarks> -y <landmark types> -z <lambdas>
            // -n <num sets> -d <list of object densities> -t <list of object types> -v <list of some object variable> 
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == ""|| minLambda == 0 || 
            fanouts == "" || numLandmarks == "" || landmarkTypes == "" || lambdas == "" || 
            numSets == 0 || objDensities == "" || objTypes == "" || objVariable == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->buildCOLTParameterIndexes(bgrFilePath,parameters,filePathPrefix,numSets,objDensities,objTypes,objVariable,fanouts,numLandmarks,landmarkTypes,lambdas,minLambda);

    } else if (experiment == this->RUN_PARAMETER_TEST_KNN_QUERIES) {
        if (argc < 31) {
            // Arguments: -g <binary graph file> -p <parameters> -f <index output file path prefix>
            // -q <query vertex file> -u <minimum lambda> -w <fanouts> -x <num landmarks> -y <landmark types> -z <lambdas>
            // -n <num sets> -d <list of object densities> -t <list of object types> -v <list of some object variable> 
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || queryVertexFile == "" ||
            minLambda == 0 || fanouts == "" || numLandmarks == "" || landmarkTypes == "" || lambdas == "" || 
            numSets == 0 || objDensities == "" || objTypes == "" || objVariable == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->runCOLTParameterTestQueries(bgrFilePath,parameters,filePathPrefix,queryVertexFile,numSets,objDensities,objTypes,objVariable,fanouts,numLandmarks,landmarkTypes,lambdas,minLambda);

    } else {
        std::cout << "Input Experiment Phase: " << experiment << "\n\n";
        std::cerr << "Invalid experimental step!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
}

void LMKExperimentsCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::LMK_EXPERIMENTS_CMD + " -e <experimental step>\n\n"
              << "Steps:\n"
              << utility::getFormattedUsageString(this->BUILD_PARTITIONS,"Build the hierarchical graph partitions the SUL-Trees are built from") + "\n"
              << utility::getFormattedUsageString(this->BUILD_RN_INDEXES,"Build all indexes") + "\n"
              << utility::getFormattedUsageString(this->BUILD_OBJ_INDEXES,"Generate object sets and build object indexes") + "\n"
              << utility::getFormattedUsageString(this->BUILD_RW_OBJ_INDEXES,"Build object indexes for supplied real-world object sets") + "\n"
              << utility::getFormattedUsageString(this->RUN_KNN_QUERIES,"Run all landmark kNN experiments") + "\n"
              << utility::getFormattedUsageString(this->RUN_KNN_SINGLE_QUERIES,"Run single landmark kNN experiments") + "\n"
              << utility::getFormattedUsageString(this->RUN_SDL_EXPERIMENT,"Time the subgraph distance list (SDL) computation with each method") + "\n"
              << utility::getFormattedUsageString(this->BUILD_PARAMETER_TEST_RN_INDEXES,"Build road network indexes for varying parameters experiments") + "\n"
              << utility::getFormattedUsageString(this->BUILD_PARAMETER_TEST_OBJ_INDEXES,"Build object indexes for varying parameters experiments") + "\n"
              << utility::getFormattedUsageString(this->RUN_PARAMETER_TEST_KNN_QUERIES,"Run kNN queries for parameter experiments") + "\n";
}

void LMKExperimentsCommand::showPhaseUsage(std::string method, std::string programName)
{
    if (method == this->BUILD_RN_INDEXES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_EXPERIMENTS_CMD
                  << " -e " + this->BUILD_RN_INDEXES + " -g <binary graph file>\n"
                  << "-p <parameter key-value pairs> -f <index output file path prefix> -s <stats output file>\n";
    } else if (method == this->BUILD_PARTITIONS ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_EXPERIMENTS_CMD
                  << " -e " + this->BUILD_PARTITIONS + " -g <binary graph file>\n"
                  << "-p <parameter key-value pairs> -f <index output file path prefix> -s <stats output file>\n";
    } else if (method == this->RUN_SDL_EXPERIMENT ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_EXPERIMENTS_CMD
                  << " -e " + this->RUN_SDL_EXPERIMENT + " -g <binary graph file>\n"
                  << "-p <parameter key-value pairs> -f <index file path prefix> -s <stats output file>\n";
    } else if (method == this->BUILD_OBJ_INDEXES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_EXPERIMENTS_CMD
                  << " -e " + this->BUILD_OBJ_INDEXES + " -g <binary graph file>\n"
                  << "-p <parameter key-value pairs> -n <num sets> -d <list of object set densities>\n"
                  << "-t <list of object types>  -v <list of some object variable> -f <index output file path prefix> -s <stats output file>\n"
                  << "[-a <binary graph of another dataset>]\n\n"
                  << "With obj_sets=1, the random object sets are generated, each drawn from set_seed (an integer, in -p), its size,\n"
                  << "variable and index, so the same seed gives the same sets of a graph. With -a and obj_sets=1, the seed is not\n"
                  << "needed and the random object sets are copied from the other dataset's sets of the same names instead of\n"
                  << "generated. The other dataset must have the same vertex IDs and neighbours (edge weights may differ)\n";
    } else if (method == this->BUILD_RW_OBJ_INDEXES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_EXPERIMENTS_CMD
                  << " -e " + this->BUILD_RW_OBJ_INDEXES + " -g <binary graph file>\n"
                  << "-p <parameter key-value pairs> -n <num sets> -r <list of real-world object set names in pre-defined directory>\n"
                  << "-f <index output file path prefix> -s <stats output file>\n";
    } else if (method == this->RUN_KNN_QUERIES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_EXPERIMENTS_CMD
                  << " -e " + this->RUN_KNN_QUERIES + " -g <binary graph file>\n"
                  << "-q <query vertex file> -k <k values> -p <parameters> -n <num sets> -d <list of object set densities>\n"
                  << "-t <list of object types>  -v <list of some object variable> -f <index output file path prefix>\n"
                  << "[-r <rw object set file names>]\n";
    } else if (method == this->RUN_KNN_SINGLE_QUERIES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_EXPERIMENTS_CMD
                  << " -e " + this->RUN_KNN_SINGLE_QUERIES + " -g <binary graph file>\n"
                  << "-q <query vertex file> -k <k values> -p <parameters> -n <num sets> -d <list of object set densities>\n"
                  << "-t <list of object types>  -v <list of some object variable> -f <index output file path prefix>\n"
                  << "-s <stats output file>\n";
    } else if (method == this->BUILD_PARAMETER_TEST_RN_INDEXES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_EXPERIMENTS_CMD
                  << " -e " + this->BUILD_PARAMETER_TEST_RN_INDEXES + " -g <binary graph file> -f <index output file path prefix>\n"
                  << "-u <minimum lambda> -w <fanouts> -x <num landmarks> -y <landmark types>\n";
    } else if (method == this->BUILD_PARAMETER_TEST_OBJ_INDEXES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_EXPERIMENTS_CMD
                  << " -e " + this->BUILD_PARAMETER_TEST_OBJ_INDEXES + " -g <binary graph file> -p <parameters> -f <index output file path prefix>\n"
                  << "-u <minimum lambda> -w <fanouts> -x <num landmarks> -y <landmark types> -z <lambdas>\n"
                  << "-n <num sets> -d <object densities> -t <object types> -v <object variables>";
    } else if (method == this->RUN_PARAMETER_TEST_KNN_QUERIES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_EXPERIMENTS_CMD
                  << " -e " + this->RUN_PARAMETER_TEST_KNN_QUERIES + " -g <binary graph file> -p <parameters> -f <index output file path prefix>\n"
                  << "-q <query vertex file> -u <minimum lambda> -w <fanouts> -x <num landmarks> -y <landmark types> -z <lambdas>\n"
                  << "-n <num sets> -d <object densities> -t <object types> -v <object variables>";
    } else {
        std::cerr << "Invalid experiment phase!" << std::endl;
        this->showCommandUsage(programName);
    }
}

std::unordered_map<std::string, std::string> LMKExperimentsCommand::getParameters(std::string parameters)
{
    std::unordered_map<std::string,std::string> parameterMap;
    std::vector<std::string> pairs = utility::splitByDelim(parameters,',');
    
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        std::vector<std::string> pair = utility::splitByDelim(pairs[i],'=');
        if (pair.size() == 2) {
            parameterMap[pair[0]] = pair[1];
        } else {
            std::cerr << "Invalid key-value pair in parameter string" << std::endl;
            exit(1);
        }
    }
    return parameterMap;
}


void LMKExperimentsCommand::buildPartitions(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile)
{
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);

    if (parameterMap["sltree"] != "1") {
        std::cout << "SULTree index is not enabled, no partition to build" << std::endl;
        return;
    }

    // Find all additional fields we need to add to index stats tuples
    std::vector<std::string> specialFields;
    int field = 0;
    while (true) {
        std::string key = "special_field_" + std::to_string(field);
        if (parameterMap.find(key) != parameterMap.end()) {
            specialFields.push_back(parameterMap[key]);
            ++field;
        } else {
            break;
        }
    }

    int fanout = std::stoi(parameterMap["sltree_fanout"]);
    std::size_t maxLeafSize = std::stoi(parameterMap["sltree_maxleafsize"]);
    if (fanout < 2 || maxLeafSize < 32) {
        std::cerr << "Invalid SULTree partitioning parameters!\n";
        exit(1);
    }

    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);
    std::vector<std::string> parameterKeys, parameterValues;
    parameterKeys.push_back("fanout");
    parameterValues.push_back(parameterMap["sltree_fanout"]);
    parameterKeys.push_back("maxleafsize");
    parameterValues.push_back(parameterMap["sltree_maxleafsize"]);
    std::string idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PARTITION_CMD,parameterKeys,parameterValues);
    HierarchicalPartition partition;
    this->buildPartition(graph,fanout,maxLeafSize,idxOutputFile,statsOutputFile,specialFields,partition);
}

void LMKExperimentsCommand::buildPartition(Graph& graph, int fanout, std::size_t maxLeafSize, std::string idxOutputFile, std::string statsOutputFile,
                                           std::vector<std::string> specialFields, HierarchicalPartition& partition)
{
    // The partitioning method is chosen here (METIS is the only one so far)
    METISGraphPartitioner partitioner(graph.getNumVertices(),graph.getNumEdges());
    partition.build(graph,partitioner,fanout,maxLeafSize);

    IndexTuple stats(partition.getNetworkName(),partition.getNumVertices(),partition.getNumEdges(),
                     constants::IDX_PARTITION_CMD,partition.getBuildTimeMs(),partition.computeIndexSize());
    stats.addSupplementaryFields("fanout",std::to_string(fanout));
    stats.addSupplementaryFields("max_leaf_size",std::to_string(maxLeafSize));
    stats.addSupplementaryFields("method",partition.getMethodName());
    stats.addSupplementaryFields("levels",std::to_string(partition.getNumLevels()));
    stats.addSupplementaryFields("tree_nodes",std::to_string(static_cast<unsigned int>(partition.nodes.size())));
    stats.setAdditionalFields(specialFields);

    serialization::outputIndexToBinaryFile<HierarchicalPartition>(partition,idxOutputFile);

    this->outputCommandStats(statsOutputFile,stats.getTupleString());

    std::cout << "Hierarchical partition (f=" + std::to_string(fanout) + ",lambda=" + std::to_string(maxLeafSize) + ") successfully created!" << std::endl;
}

// Dijkstra search spillover (gamma) of the subgraph distance list computation: vertices settled vs vertices in the subgraph, as the
// average over all searches weighted by subgraph size (all vertices settled over all subgraph vertices) and the minimum and maximum
void LMKExperimentsCommand::computeSpillFactors(SULTree& sltree, double& avg, double& min, double& max)
{
    double sumSettled = 0, sumVertices = 0;
    double factor;
    max = 0;
    min = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < sltree.nodeSpilloverFactors.size(); ++i) {
        int numNodeVertices = sltree.treeNodes[i].getNumVertices();
        for (std::size_t j = 0; j < sltree.nodeSpilloverFactors[i].size(); ++j) {
            factor = sltree.nodeSpilloverFactors[i][j];
            sumSettled += factor*numNodeVertices; // Number of vertices that were settled by Dijkstra
            sumVertices += numNodeVertices; // Ideal number of vertices that should be settled
            if (factor > max) {
                max = factor;
            }
            if (factor < min) {
                min = factor;
            }
        }
    }
    avg = sumSettled/sumVertices;
}

// Statistics of the Dijkstra search spillover (gamma) of the subgraph distance list computation (see computeSpillFactors)
void LMKExperimentsCommand::addSpillFactorFields(SULTree& sltree, IndexTuple& stats)
{
    double avg, min, max;
    this->computeSpillFactors(sltree,avg,min,max);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(3) << max;
    stats.addSupplementaryFields("spill_factor_max",stream.str());
    stream.str("");
    stream << std::fixed << std::setprecision(3) << min;
    stats.addSupplementaryFields("spill_factor_min",stream.str());
    stream.str("");
    stream << std::fixed << std::setprecision(3) << avg;
    stats.addSupplementaryFields("spill_factor_avg",stream.str());
}

// SDL experiment: build the SUL-Tree (in memory only) with each subgraph distance list method and record the
// construction time, the phase times and the search spillover of each as one stats row per method
void LMKExperimentsCommand::runSDLExperiment(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile)
{
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);

    // Find all additional fields we need to add to index stats tuples
    std::vector<std::string> specialFields;
    int field = 0;
    while (true) {
        std::string key = "special_field_" + std::to_string(field);
        if (parameterMap.find(key) != parameterMap.end()) {
            specialFields.push_back(parameterMap[key]);
            ++field;
        } else {
            break;
        }
    }

    int fanout = std::stoi(parameterMap["sltree_fanout"]);
    std::size_t maxLeafSize = std::stoi(parameterMap["sltree_maxleafsize"]);
    int numNodeLandmarks = std::stoi(parameterMap["sltree_numlmks"]);
    std::string landmarkTypeStr = parameterMap["sltree_lmktype"];
    if (fanout < 2 || maxLeafSize < 32 || numNodeLandmarks < 1 || landmarkTypeStr == "") {
        std::cerr << "Invalid SULTree parameters!\n";
        exit(1);
    }

    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);

    std::vector<std::string> parameterKeys, parameterValues;
    parameterKeys.push_back("fanout");
    parameterValues.push_back(parameterMap["sltree_fanout"]);
    parameterKeys.push_back("maxleafsize");
    parameterValues.push_back(parameterMap["sltree_maxleafsize"]);
    std::string partitionFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PARTITION_CMD,parameterKeys,parameterValues);
    if (!std::ifstream(partitionFile).good()) {
        std::cerr << "Partition file " << partitionFile << " not found, build it first with -e " << this->BUILD_PARTITIONS << " (bash buildPartitions)" << std::endl;
        exit(1);
    }
    HierarchicalPartition partition = serialization::getIndexFromBinaryFile<HierarchicalPartition>(partitionFile);
    double vertexNodePositionsSize = partition.computeVertexNodePositionsSize(); // In MB
    std::cout << "Vertex node positions without subgraph vertex ordering would need " << vertexNodePositionsSize << "MB" << std::endl;

    // The methods in the order they were developed (each row improves on the previous one)
    const std::vector<SDL_METHOD> methods = {SDL_METHOD::SDL_HASH_TABLE, SDL_METHOD::SDL_DIJKSTRA, SDL_METHOD::SDL_BORDER_SET};
    const std::vector<std::string> methodNames = {"hash_table", "dijkstra", "border_set"}; // In the stats rows (sdl_<name>)
    const std::vector<std::string> methodDescriptions = {"SDL construction with Dijkstra using hash table",
                                                         "SDL construction with Dijkstra without hash table",
                                                         "SDL construction with border-vertex-restricted-Dijkstra without hash table"};
    for (std::size_t m = 0; m < methods.size(); ++m) {
        std::cout << "\nBuilding SULTree: " << methodDescriptions[m] << "..." << std::endl;
        SULTree sltree(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),numNodeLandmarks);
        bool validMethod = false;
        COLTREE_LANDMARK_TYPE landmarkType = sltree.getLandmarkTypeFromString(landmarkTypeStr,validMethod);
        if (!validMethod) {
            std::cerr << "Invalid landmark selection method!" << std::endl;
            std::exit(1);
        }
        sltree.buildSULTree(graph,partition,landmarkType,methods[m]);
        // The methods differ only in the distance lists of the non-root nodes, so the experiment reports their computation time
        // (the partition, the root's lists and the landmark selection are the same for every method)
        double sdlTimeMs = sltree.getDistanceListTimeMs();

        double indexSize = sltree.computeIndexSize(); // In MB

        IndexTuple stats(sltree.getNetworkName(),sltree.getNumVertices(),sltree.getNumEdges(),"sdl_" + methodNames[m],sdlTimeMs,indexSize);
        stats.addSupplementaryFields("fanout",std::to_string(fanout));
        stats.addSupplementaryFields("max_leaf_size",std::to_string(maxLeafSize));
        stats.addSupplementaryFields("num_lmks",std::to_string(numNodeLandmarks));
        stats.addSupplementaryFields("lmk_type",std::to_string(landmarkType));
        std::stringstream stream;
        stream << std::fixed << std::setprecision(3) << partition.getBuildTimeMs();
        stats.addSupplementaryFields("partition_time",stream.str());
        stream.str("");
        stream << std::fixed << std::setprecision(3) << sltree.getLandmarkTimeMs();
        stats.addSupplementaryFields("landmark_time",stream.str());
        stream.str("");
        stream << std::fixed << std::setprecision(3) << sltree.getDistanceListTimeMs();
        stats.addSupplementaryFields("sdl_time",stream.str());
        stream.str("");
        stream << std::fixed << std::setprecision(3) << sltree.getAllocationTimeMs();
        stats.addSupplementaryFields("allocation_time",stream.str());
        this->addSpillFactorFields(sltree,stats);
        stats.addSupplementaryFields("tree_nodes",std::to_string(static_cast<unsigned int>(sltree.treeNodes.size())));
        // Without the subgraph vertex ordering, the index would also need the vertex node positions (its size is in the memory column)
        stream.str("");
        stream << std::fixed << std::setprecision(3) << vertexNodePositionsSize;
        stats.addSupplementaryFields("vertex_node_positions_size",stream.str());
        stats.setAdditionalFields(specialFields);

        this->outputCommandStats(statsOutputFile,stats.getTupleString());

        double spillAvg, spillMin, spillMax;
        this->computeSpillFactors(sltree,spillAvg,spillMin,spillMax);
        std::stringstream spill; // Local stream, so the precision does not carry over to later output
        spill << std::fixed << std::setprecision(3) << "avg " << spillAvg << ", min " << spillMin << ", max " << spillMax;
        std::stringstream size;
        size << std::fixed << std::setprecision(3) << indexSize << "MB (" << indexSize+vertexNodePositionsSize << "MB with vertex node positions)";
        std::cout << methodDescriptions[m] << ": SDL computation time " << sdlTimeMs << "ms, spillover factor (gamma) " << spill.str()
                  << ", index size " << size.str() << std::endl;
    }
}

void LMKExperimentsCommand::buildIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile)
{
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    
    bool buildGtreeIdx = parameterMap["gtree"] == "1";
    bool buildALTIdx = parameterMap["alt"] == "1";
    bool buildSULTreeIdx = parameterMap["sltree"] == "1";
    
    // Find all additional fields we need to add to index stats tuples
    std::vector<std::string> specialFields;
    int field = 0;
    while (true) {
        std::string key = "special_field_" + std::to_string(field);
        if (parameterMap.find(key) != parameterMap.end()) {
            specialFields.push_back(parameterMap[key]);
            ++field;
        } else {
            break;
        }
    }    
    
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);

    // The SUL-Tree is built from a saved hierarchical partition (-e partitions), so check it exists before building any index
    std::string partitionFile = "";
    if (buildSULTreeIdx) {
        const std::vector<std::string> partitionKeys = {"fanout","maxleafsize"};
        const std::vector<std::string> partitionValues = {parameterMap["sltree_fanout"],parameterMap["sltree_maxleafsize"]};
        partitionFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PARTITION_CMD,partitionKeys,partitionValues);
        if (!std::ifstream(partitionFile).good()) {
            std::cerr << "Partition file " << partitionFile << " not found, build it first with -e " << this->BUILD_PARTITIONS << " (bash buildPartitions)" << std::endl;
            exit(1);
        }
    }
    
    // Put all of these into functions (except input validation) so indexes get released after serialization
    
    if (buildGtreeIdx) {
        int fanout = std::stoi(parameterMap["gtree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["gtree_maxleafsize"]);
        if (fanout < 2 || maxLeafSize < 32) {
            std::cerr << "Invalid Gtree parameters!\n";
            exit(1);    
        }
        const std::vector<std::string> gtreeKeys = {"fanout","maxleafsize"};
        const std::vector<std::string> gtreeValues = {parameterMap["gtree_fanout"],parameterMap["gtree_maxleafsize"]};
        std::string idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,gtreeKeys,gtreeValues);
        this->buildGtree(graph,fanout,maxLeafSize,idxOutputFile,statsOutputFile,specialFields);
    }
    
    if (buildALTIdx) {
        int numLandmarks = std::stoi(parameterMap["alt_numlandmarks"]);
        std::string selectionMethod = parameterMap["alt_landmarktype"];
        if (numLandmarks < 1 || selectionMethod == "") {
            std::cerr << "Invalid ALT parameters!\n";
            exit(1);    
        }
        const std::vector<std::string> altKeys = {"landmarks","type"};
        const std::vector<std::string> altValues = {parameterMap["alt_numlandmarks"],parameterMap["alt_landmarktype"]};
        std::string idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,altKeys,altValues);
        this->buildALT(graph,numLandmarks,selectionMethod,idxOutputFile,statsOutputFile,specialFields);
    }
    
    if (buildSULTreeIdx) {
        int fanout = std::stoi(parameterMap["sltree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["sltree_maxleafsize"]);
        std::size_t numNodeLandmarks = std::stoi(parameterMap["sltree_numlmks"]);
        std::string landmarkType = parameterMap["sltree_lmktype"];
        if (fanout < 2 || maxLeafSize < 32 || numNodeLandmarks < 1 || landmarkType == "") {
            std::cerr << "Invalid SULTree parameters!\n";
            exit(1);    
        }
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {parameterMap["sltree_fanout"],parameterMap["sltree_maxleafsize"],parameterMap["sltree_numlmks"],parameterMap["sltree_lmktype"]};
        std::string idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        HierarchicalPartition partition = serialization::getIndexFromBinaryFile<HierarchicalPartition>(partitionFile);
        this->buildSULTree(graph,partition,numNodeLandmarks,landmarkType,idxOutputFile,statsOutputFile,specialFields);
    }
    
    std::cout << "Index building complete!" << std::endl;
}

void LMKExperimentsCommand::buildExternalIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile)
{
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    
    bool buildPHLIdx = parameterMap["phl"] == "1";
    
    // Find all additional fields we need to add to index stats tuples
    std::vector<std::string> specialFields;
    int field = 0;
    while (true) {
        std::string key = "special_field_" + std::to_string(field);
        if (parameterMap.find(key) != parameterMap.end()) {
            specialFields.push_back(parameterMap[key]);
            ++field;
        } else {
            break;
        }
    }
    
    // Put all of these into functions (except input validation) so indexes get released after serialization
    std::string idxOutputFile = "";
    std::vector<std::string> parameterKeys, parameterValues;

    // Some additional file paths to create intermediary and additional files for some indexes
    std::string graphDataOutputFile = "";
    
    std::string networkName = "";
    int numVertices = 0, numEdges = 0;
    if (buildPHLIdx) {
        
        // Convert graph to input file PHL index construction expects then build PHL
        {
            // Use scope to destroy graph data structure to have extra memory for methods that don't need it
            Graph tempGraph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);
            networkName = tempGraph.getNetworkName();
            numVertices = tempGraph.getNumVertices();
            numEdges = tempGraph.getNumEdges();
            graphDataOutputFile = filePathPrefix + "/data/" + utility::constructIndexFileName(networkName,"tsv",parameterKeys,parameterValues);
            tempGraph.outputToTSVFile(graphDataOutputFile);
        }

        idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(networkName,constants::IDX_PHL_CMD,parameterKeys,parameterValues);
        this->buildPHL(networkName,numVertices,numEdges,graphDataOutputFile,idxOutputFile,statsOutputFile,specialFields);

        // Delete intermediary text file as it is not needed
        std::remove(graphDataOutputFile.c_str());
    }

    std::cout << "External index building complete!" << std::endl;
}

void LMKExperimentsCommand::buildGtree(Graph& graph, int fanout, std::size_t maxLeafSize, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields)
{
    StopWatch sw;
    
    sw.start();
    Gtree gtree(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),fanout,maxLeafSize);
    gtree.buildGtree(graph);
    sw.stop();

    double processingTimeMs = sw.getTimeMs();
    // Note: G-tree uses Graph to search for object within leaf (only edges and edge weights same parts as INE)
    double memoryUsage = gtree.computeIndexSize() + graph.computeINEIndexSize();
    
    IndexTuple stats(gtree.getNetworkName(),gtree.getNumVertices(),gtree.getNumEdges(),
                     constants::IDX_GTREE_CMD,processingTimeMs,memoryUsage);
    stats.addSupplementaryFields("max_leaf_size",std::to_string(maxLeafSize));
    stats.addSupplementaryFields("fanout",std::to_string(fanout));
    stats.addSupplementaryFields("levels",std::to_string(gtree.getNumLevels()));
    int numTreeNodes = gtree.getTreeSize(), numBorders = gtree.getNumBorders(), b2bRelationships = gtree.getBorderToBorderRelationships(), avgPathCost = gtree.getAvgPathCost();
    stats.addSupplementaryFields("num_tree_nodes",std::to_string(numTreeNodes));
    stats.addSupplementaryFields("num_borders",std::to_string(numBorders));
    stats.addSupplementaryFields("num_b2b_rel",std::to_string(b2bRelationships));
    stats.addSupplementaryFields("avg_path_cost",std::to_string(avgPathCost));
    stats.setAdditionalFields(specialFields);
    serialization::outputIndexToBinaryFile<Gtree>(gtree,idxOutputFile);

    this->outputCommandStats(statsOutputFile,stats.getTupleString());
    
    std::cout << "Gtree index successfully created!" << std::endl;
}

void LMKExperimentsCommand::buildPHL(std::string networkName, int numVertices, int numEdges, std::string dataOutputFile, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields)
{
    PrunedHighwayLabeling phl;
    phl.ConstructLabel(dataOutputFile.c_str());
    
    double processingTimeMs = phl.getConstructionTime()*1000; // convert to ms
    double memoryUsage = phl.computeIndexSize();

#if defined(PARALLEL_BUILDS)
    // rn_obj_search_par builds the labels with several threads, and they differ slightly from rn_obj_search's, so its
    // construction time is printed but not recorded
    std::cout << "PHL labels built with " << omp_get_max_threads() << " threads in " << processingTimeMs << "ms (" << memoryUsage
              << "MB); the time is not recorded in " << statsOutputFile << std::endl;
#else
    IndexTuple stats(networkName,numVertices,numEdges,constants::IDX_PHL_CMD,processingTimeMs,memoryUsage);
    stats.setAdditionalFields(specialFields);
    this->outputCommandStats(statsOutputFile,stats.getTupleString());
#endif

    phl.StoreLabel(idxOutputFile.c_str());
    
    std::cout << "PHL index successfully created!" << std::endl;
}

void LMKExperimentsCommand::buildALT(Graph& graph, int numLandmarks, std::string selectionMethod, std::string idxOutputFile, 
                                  std::string statsOutputFile, std::vector<std::string> specialFields)
{
    StopWatch sw;
    bool validMethod;
    
    sw.start();
    ALT alt(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges());
    LANDMARK_TYPE landmarkType = alt.getLandmarkType(selectionMethod,validMethod);
    if (!validMethod) {
        std::cerr << "Invalid landmark selection method!" << std::endl;
        std::exit(1);
    }
    alt.buildALT(graph,landmarkType,numLandmarks);
    sw.stop();
    
    double processingTimeMs = sw.getTimeMs();    
    double memoryUsage = alt.computeIndexSize();

    IndexTuple stats(alt.getNetworkName(),alt.getNumVertices(),alt.getNumEdges(),constants::IDX_ALT_CMD,processingTimeMs,memoryUsage);
    stats.setAdditionalFields(specialFields);
    stats.addSupplementaryFields("landmarks",std::to_string(numLandmarks));
    stats.addSupplementaryFields("type",selectionMethod);

    serialization::outputIndexToBinaryFile<ALT>(alt,idxOutputFile);
    
    this->outputCommandStats(statsOutputFile,stats.getTupleString());

    std::cout << "ALT index successfully created!" << std::endl;
}

void LMKExperimentsCommand::buildSULTree(Graph& graph, const HierarchicalPartition& partition, int numNodeLandmarks, std::string landmarkTypeStr,
                                        std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields)
{
    StopWatch sw;
    
    sw.start();
    SULTree sltree(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),numNodeLandmarks);
    bool validMethod = false;
    COLTREE_LANDMARK_TYPE landmarkType = sltree.getLandmarkTypeFromString(landmarkTypeStr,validMethod);
    if (!validMethod) {
        std::cerr << "Invalid landmark selection method!" << std::endl;
        std::exit(1);
    }
    sltree.buildSULTree(graph,partition,landmarkType);
    sw.stop();

    double processingTimeMs = sw.getTimeMs() + partition.getBuildTimeMs(); // The partitioning is part of the SUL-Tree construction
    double memoryUsage = sltree.computeIndexSize();
    
    IndexTuple stats(sltree.getNetworkName(),sltree.getNumVertices(),sltree.getNumEdges(),
                     constants::IDX_SULTREE_CMD,processingTimeMs,memoryUsage);
    stats.addSupplementaryFields("fanout",std::to_string(partition.getFanout()));
    stats.addSupplementaryFields("max_leaf_size",std::to_string(partition.getMaxLeafSize()));
    stats.addSupplementaryFields("num_lmks",std::to_string(numNodeLandmarks));
    stats.addSupplementaryFields("lmk_type",std::to_string(landmarkType));
    stats.addSupplementaryFields("levels",std::to_string(sltree.getNumLevels()));
    this->addSpillFactorFields(sltree,stats);
    stats.addSupplementaryFields("tree_nodes",std::to_string(static_cast<unsigned int>(sltree.treeNodes.size())));
    std::stringstream stream;
    stream << std::fixed << std::setprecision(3) << partition.getBuildTimeMs();
    stats.addSupplementaryFields("partition_time",stream.str());
    stats.setAdditionalFields(specialFields);

    serialization::outputIndexToBinaryFile<SULTree>(sltree,idxOutputFile);

    this->outputCommandStats(statsOutputFile,stats.getTupleString());
    
    std::cout << "SULTree index (f=" + std::to_string(partition.getFanout()) + ",lambda=" +std::to_string(partition.getMaxLeafSize()) 
        + ",m=" + std::to_string(numNodeLandmarks)+ ",type=" + landmarkTypeStr + ") successfully created!" << std::endl;
}
void LMKExperimentsCommand::buildObjectIndexes(std::string bgrFileName, std::string parameters, std::size_t numSets, std::string objDensities, 
                                               std::string objTypes, std::string objVariable, std::string filePathPrefix, std::string statsOutputFile, std::string sourceBgrFilePath)
{
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    
    // Find all additional fields we need to add to index stats tuples
    std::vector<std::string> specialFields;
    int field = 0;
    while (true) {
        std::string key = "special_field_" + std::to_string(field);
        if (parameterMap.find(key) != parameterMap.end()) {
            specialFields.push_back(parameterMap[key]);
            ++field;
        } else {
            break;
        }
    }

    bool generateObjSets = parameterMap["obj_sets"] == "1";
    bool buildOccList = parameterMap["occ_list"] == "1";
    bool buildRtree = parameterMap["rtree"] == "1";
    bool buildObjList = parameterMap["obj_list"] == "1";
    bool buildNVD = parameterMap["nvd"] == "1";
    bool buildCOLT = parameterMap["colt"] == "1";
    if (sourceBgrFilePath != "" && !generateObjSets) {
        std::cerr << "Object sets can only be copied from another dataset (-a) with obj_sets=1!\n";
        exit(1);
    }
    unsigned int setSeed = 0;
    if (generateObjSets && sourceBgrFilePath == "") {
        setSeed = SetGenerator::parseSeed(parameterMap["set_seed"],"set_seed in -p");
    }

    std::vector<std::string> strDensities = utility::splitByDelim(objDensities,',');
    std::vector<double> densities;
    for(std::size_t i = 0; i < strDensities.size(); ++i) {
        double density = std::stod(strDensities[i]);
        if (density > 0) {
            densities.push_back(density);
        } else {
            std::cerr << "Invalid density in list provided!\n";
            exit(1);
        }
    }
    std::vector<std::string> strVariables = utility::splitByDelim(objVariable,',');
    std::vector<int> variables;
    for(std::size_t i = 0; i < strVariables.size(); ++i) {
        int variable = std::stoi(strVariables[i]);
        if (variable > 0) {
            variables.push_back(variable);
        } else {
            std::cerr << "Invalid variable in list provided (must be greater than zero)!\n";
            exit(1);    
        }
    }
    std::vector<std::string> types = utility::splitByDelim(objTypes,',');
    if (densities.size() == 0 || types.size() == 0 || variables.size() == 0) {
        std::cerr << "Not enough densities or types provided!\n";
        exit(1);
    }

    if (generateObjSets) {
        // Generate and Output Object Sets
        this->generateObjectSets(graph,numSets,densities,types,variables,filePathPrefix,statsOutputFile,setSeed,sourceBgrFilePath);
    }
    
    if (buildOccList) {
        int fanout = std::stoi(parameterMap["gtree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["gtree_maxleafsize"]);
        if (fanout < 2 || maxLeafSize < 32) {
            std::cerr << "Invalid Gtree parameters!\n";
            exit(1);
        }
        const std::vector<std::string> gtreeKeys = {"fanout","maxleafsize"};
        const std::vector<std::string> gtreeValues = {parameterMap["gtree_fanout"],parameterMap["gtree_maxleafsize"]};
        std::string gtreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,gtreeKeys,gtreeValues);
        this->buildOccurenceLists(gtreeIdxFile,numSets,densities,types,variables,filePathPrefix,statsOutputFile,parameterMap,specialFields);
    }
    
    if (buildRtree) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() != 1) {
            std::cerr << "One and only one R-tree branch factor can be provided in LMK experiments!\n";
            exit(1);
        }
        this->buildRtrees(graph,branchFactors[0],numSets,densities,types,variables,filePathPrefix,statsOutputFile,parameterMap,specialFields);
    }

    if (buildObjList) {
        int numLandmarks = std::stoi(parameterMap["alt_numlandmarks"]);
        if (numLandmarks < 1) {
            std::cerr << "Invalid ALT parameters!\n";
            exit(1);
        }
        const std::vector<std::string> altKeys = {"landmarks","type"};
        const std::vector<std::string> altValues = {parameterMap["alt_numlandmarks"],parameterMap["alt_landmarktype"]};
        std::string altIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,altKeys,altValues);
        this->buildObjectLists(graph,altIdxFile,numSets,densities,types,variables,filePathPrefix,statsOutputFile,parameterMap,specialFields);
    }
    
    if (buildNVD) {
        int numLandmarks = std::stoi(parameterMap["alt_numlandmarks"]);
        if (numLandmarks < 1) {
            std::cerr << "Invalid ALT parameters!\n";
            exit(1);
        }
        this->buildNVD(graph,numSets,densities,types,variables,filePathPrefix,statsOutputFile,parameterMap,specialFields);
    }

    if (buildCOLT) {
        int fanout = std::stoi(parameterMap["sltree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["sltree_maxleafsize"]);
        std::size_t numNodeLandmarks = std::stoi(parameterMap["sltree_numlmks"]);
        std::string landmarkType = parameterMap["sltree_lmktype"];
        if (fanout < 2 || maxLeafSize < 32 || numNodeLandmarks < 1 || landmarkType == "") {
            std::cerr << "Invalid SULTree parameters!\n";
            exit(1);    
        }
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {parameterMap["sltree_fanout"],parameterMap["sltree_maxleafsize"],parameterMap["sltree_numlmks"],parameterMap["sltree_lmktype"]};
        std::string sltreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        
        std::vector<int> listSizes = utility::getIngetersFromStringList(parameterMap["colt_listsize"],':',0);
        if (listSizes.size() < 1) {
            std::cerr << "Invalid number of list size for COLTree!\n";
            exit(1);
        } else if (listSizes.size() > 1) {
            std::cerr << "Only first list size will be used for COLTree construction!\n";
        }
        
        this->buildCOLT(graph,sltreeIdxFile,listSizes[0],numSets,densities,types,variables,filePathPrefix,statsOutputFile,parameterMap,specialFields);
    }
}

void LMKExperimentsCommand::generateObjectSets(Graph& graph, std::size_t numSets, std::vector<double> objDensities, 
                                            std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile,
                                            unsigned int setSeed, std::string sourceBgrFilePath)
{
    SetGenerator sg;
    INE ine;
    int numVertices = graph.getNumVertices(), setSize;
    std::vector<VertexID> sampleSet;
    std::string objSetOutputFile = "";
    StopWatch sw;
    double totalINETime, totalINEMemory;
    
    // With another dataset's binary graph, its random object sets of the same names are copied instead of generated (as the sets
    // are drawn from the vertex IDs alone, they are valid for any dataset with the same vertex IDs and neighbours)
    std::string sourceNetworkName = "", sourceSetFile, sourceSetType;
    double sourceSetDensity;
    int sourceSetSize, sourceSetVariable;
    if (sourceBgrFilePath != "") {
        Graph sourceGraph = serialization::getIndexFromBinaryFile<Graph>(sourceBgrFilePath);
        std::string difference;
        if (!graph.hasSameNeighbours(sourceGraph,difference)) {
            std::cerr << "Cannot copy the sets of " << sourceGraph.getNetworkName() << " to " << graph.getNetworkName() << ": " << difference << std::endl;
            exit(1);
        }
        sourceNetworkName = sourceGraph.getNetworkName();
        // Every set to copy must exist before any is written (the set of density one holds every vertex, so it is not copied)
        std::size_t numMissing = 0;
        for (std::size_t i = 0; i < objTypes.size(); ++i) {
            if (objTypes[i] != constants::RAND_OBJ_SET) {
                std::cerr << "Only " << constants::RAND_OBJ_SET << " object sets can be copied from another dataset (-a)!\n";
                exit(1);
            }
            for (std::size_t j = 0; j < objDensities.size(); ++j) {
                for (std::size_t l = 0; l < objVariable.size() && objDensities[j] != 1; ++l) {
                    for (std::size_t k = 0; k < numSets; ++k) {
                        sourceSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(sourceNetworkName,objTypes[i],objDensities[j],objVariable[l],k);
                        if (!std::ifstream(sourceSetFile).good()) {
                            if (numMissing < 5) {
                                std::cerr << "Missing " << sourceSetFile << std::endl;
                            }
                            ++numMissing;
                        }
                    }
                }
            }
        }
        if (numMissing > 0) {
            std::cerr << numMissing << " object set file(s) of " << sourceNetworkName << " are missing, so none was copied" << std::endl;
            exit(1);
        }
    }

    // Generate and Output Object Sets
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        if (objTypes[i] == constants::RAND_OBJ_SET) {
            for (std::size_t j = 0; j < objDensities.size(); ++j) {
                const std::size_t totalSets = utility::getNumObjectSets(objDensities[j],numSets);
                for (std::size_t l = 0; l < objVariable.size(); ++l) {
                    setSize = std::ceil(numVertices*objDensities[j]);
                    totalINETime = 0;
                    totalINEMemory = 0;
                    for (std::size_t k = 0; k < totalSets; ++k) {
                        if (objDensities[j] != 1 && sourceNetworkName != "") {
                            sourceSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(sourceNetworkName,objTypes[i],objDensities[j],objVariable[l],k);
                            sampleSet = utility::getPointSetFromFile(sourceSetFile,sourceSetType,sourceSetDensity,sourceSetSize,sourceSetVariable);
                            if (sourceSetType != objTypes[i] || sourceSetDensity != objDensities[j] || sourceSetVariable != objVariable[l]
                                || sourceSetSize != setSize || sampleSet.size() != static_cast<std::size_t>(setSize)) {
                                std::cerr << sourceSetFile << " is not a " << objTypes[i] << " object set of density " << objDensities[j] << " and size " << setSize << std::endl;
                                exit(1);
                            }
                        } else if (objDensities[j] != 1) {
                            // Each object set has its own engine, seeded with its size (for the density), variable and index
                            std::mt19937 engine = SetGenerator::getEngine(setSeed,SetGenerator::OBJECT_SET,{static_cast<unsigned int>(setSize),static_cast<unsigned int>(objVariable[l]),static_cast<unsigned int>(k)});
                            sampleSet = sg.generateRandomSampleSet(numVertices,setSize,engine);
                        } else {
                            sampleSet = graph.getVertexIDsVector();
                        }

                        objSetOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                        utility::writeSampleSet(objSetOutputFile,graph.getNetworkName(),objTypes[i],objDensities[j],setSize,objVariable[l],sampleSet);
                        // We only include the size of the input to the other object indexes and the time taken
                        // to load load this information into the graph as the baselines for INE
                        sw.reset();
                        sw.start();
                        std::vector<bool> vertexObjectStatus = ine.buildVertexObjectStatusVector(graph,sampleSet);
                        sw.stop();
                        totalINETime += sw.getTimeUs();
                        totalINEMemory += static_cast<double>(sizeof(bool)*vertexObjectStatus.size())/(1024*1024);
                    }
                    
                    // Average stats over all object sets and then write to stats file
                    double processingTimeMs = totalINETime/totalSets;
                    double memoryUsage = totalINEMemory/totalSets*1024; // in KB
                    
                    ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),
                                            constants::INE_KNN_QUERY,processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[l],setSize);

                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        } else {
            std::cerr << "Invalid object set type provided!\n";
            exit(1);
        }
    }
    if (sourceNetworkName != "") {
        std::cout << "Object sets successfully copied from " << sourceNetworkName << " for " << graph.getNetworkName() << std::endl;
    } else {
        std::cout << "Object sets successfully generated for " << graph.getNetworkName() << std::endl;
    }
}

// The object index builders read each object set file by the name made from its type and density, so its header must have the
// same ones. The name and the header both write the density with <<, so the densities are compared as written
void LMKExperimentsCommand::checkObjectSetHeader(std::string objSetFile, std::string objType, double objDensity, std::string setType, double setDensity)
{
    std::stringstream nameDensity, headerDensity;
    nameDensity << objDensity;
    headerDensity << setDensity;
    if (setType != objType || headerDensity.str() != nameDensity.str()) {
        std::cerr << "Object set file " << objSetFile << " has type " << setType << " and density " << headerDensity.str()
                  << " in its header, but its name has type " << objType << " and density " << nameDensity.str() << std::endl;
        std::exit(1);
    }
}

void LMKExperimentsCommand::buildOccurenceLists(std::string gtreeIdxFile, std::size_t numSets, std::vector<double> objDensities, 
                                             std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                                             std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string> specialFields)
{
    Gtree gtree = serialization::getIndexFromBinaryFile<Gtree>(gtreeIdxFile);
    
    int setSize = 0, setVariable;
    std::string setType;
    double setDensity;
    std::string objSetFile, objIdxOutputFile, objIdxOutputFilePrefix;
    
    StopWatch sw;
    double totalTime = 0, totalMemory = 0;
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            const std::size_t totalSets = utility::getNumObjectSets(objDensities[j],numSets);
            for (std::size_t l = 0; l < objVariable.size(); ++l) {
                totalTime = 0;
                totalMemory = 0;
                
                for (std::size_t k = 0; k < totalSets; ++k) {
                    objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(gtree.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                    std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                    this->checkObjectSetHeader(objSetFile,objTypes[i],objDensities[j],setType,setDensity);

                    sw.start();
                    OccurenceList occList(setType,setDensity,setVariable,setSize);
                    for (auto objIt = objectVertices.begin(); objIt != objectVertices.end(); ++objIt) {
                        // Find the Gtree leaf index for this object and add it to occurence list
                        // for that leaf (create list if it doesn't exist)
                        int leafIdx = gtree.getLeafIndex(*objIt);
                        occList.addLeafOccurence(leafIdx,*objIt);
                        
                        // Propagate this to occurence lists of parents of leaf node
                        int childIdx = leafIdx;
                        int parentIdx = gtree.getParentIndex(leafIdx);
                        while (parentIdx != -1) {
                            occList.addParentOccurence(parentIdx, childIdx);
                            
                            // Go up a level (until we reach root)
                            childIdx = parentIdx;
                            parentIdx = gtree.getParentIndex(childIdx);
                        }
                    }
                    sw.stop();
                    totalTime += sw.getTimeUs();
                    totalMemory += occList.computeIndexSize()*1024; // in KB

                    objIdxOutputFilePrefix = utility::constructObjIdxFilenamePrefix(gtree.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                    objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxOutputFilePrefix,constants::OBJ_IDX_GTREE,parameterMap);
                    serialization::outputIndexToBinaryFile<OccurenceList>(occList,objIdxOutputFile);
                    
                }
                
                // Average stats over the object sets and then write to stats file
                double processingTimeMs = totalTime/totalSets;
                double memoryUsage = totalMemory/totalSets;
                
                ObjectIndexTuple stats(gtree.getNetworkName(),gtree.getNumVertices(),gtree.getNumEdges(),
                                    constants::OBJ_IDX_GTREE,processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[l],setSize);
                stats.setAdditionalFields(specialFields);

                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
    }
    std::cout << "Occurrence lists successfully created for " << gtree.getNetworkName() << std::endl;
}

void LMKExperimentsCommand::buildRtrees(Graph& graph, int branchFactor, std::size_t numSets, std::vector<double> objDensities, 
                                     std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string> specialFields)
{
    int setSize = 0, setVariable;
    std::string setType;
    double setDensity;
    std::string objSetFile, objIdxOutputFile, objIdxOutputFilePrefix;

    StopWatch sw;
    double totalTime, totalMemory;
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            const std::size_t totalSets = utility::getNumObjectSets(objDensities[j],numSets);
            for (std::size_t l = 0; l < objVariable.size(); ++l) {
                totalTime = 0;
                totalMemory = 0;

                for (std::size_t k = 0; k < totalSets; ++k) {
                    objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                    std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                    this->checkObjectSetHeader(objSetFile,objTypes[i],objDensities[j],setType,setDensity);

                    sw.start();
                    std::vector<CoordinatePair> objectCoords;
                    for (std::size_t i = 0; i < objectVertices.size(); ++i) {
                        CoordinatePair objectCoordPair;
                        graph.getCoordinates(objectVertices[i],objectCoordPair.first,objectCoordPair.second);
                        objectCoords.push_back(objectCoordPair);
                    }
                    PackedRtree rtree(branchFactor,setType,setDensity,setVariable,setSize);
                    rtree.bulkLoad(objectVertices,objectCoords);
                    sw.stop();
                    totalTime += sw.getTimeUs();
                    totalMemory += rtree.computeIndexSize()*1024; // in KB

                    objIdxOutputFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                    objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxOutputFilePrefix,constants::OBJ_IDX_RTREE,parameterMap);
                    serialization::outputIndexToBinaryFile<PackedRtree>(rtree,objIdxOutputFile);
                }

                // Average stats over the object sets and then write to stats file
                double processingTimeMs = totalTime/totalSets;
                double memoryUsage = totalMemory/totalSets;
                
                ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),
                                    constants::OBJ_IDX_RTREE,processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[l],setSize);
                stats.setAdditionalFields(specialFields);
                stats.addSupplementaryFields("rtree_branch_factor",std::to_string(branchFactor));

                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
    }
    std::cout << "Rtrees successfully created for " << graph.getNetworkName() << std::endl;
}

void LMKExperimentsCommand::buildObjectLists(Graph& graph, std::string altIdxFile, std::size_t numSets, std::vector<double> objDensities, 
                                          std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, 
                                          std::string statsOutputFile, std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string> specialFields)
{
    ALT alt = serialization::getIndexFromBinaryFile<ALT>(altIdxFile);
    
    int setSize = 0, setVariable;
    std::string setType;
    double setDensity;
    std::string objSetFile, objIdxOutputFile, objIdxOutputFilePrefix;
    
    StopWatch sw;
    double totalTime = 0, totalMemory = 0;
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            const std::size_t totalSets = utility::getNumObjectSets(objDensities[j],numSets);
            for (std::size_t l = 0; l < objVariable.size(); ++l) {
                totalTime = 0;
                totalMemory = 0;
                
                for (std::size_t k = 0; k < totalSets; ++k) {
                    objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(alt.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                    std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                    this->checkObjectSetHeader(objSetFile,objTypes[i],objDensities[j],setType,setDensity);

                    sw.start();
                    ObjectList ol(setType,setDensity,setVariable,setSize);
                    ol.buildObjectList(alt,objectVertices);
                    sw.stop();
                    totalTime += sw.getTimeUs();
                    totalMemory += ol.computeIndexSize()*1024; // in KB

                    objIdxOutputFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                    objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxOutputFilePrefix,constants::OBJ_IDX_LNDMRK,parameterMap);
                    serialization::outputIndexToBinaryFile<ObjectList>(ol,objIdxOutputFile);
                    
                }
                
                // Average stats over the object sets and then write to stats file
                double processingTimeMs = totalTime/totalSets;
                double memoryUsage = totalMemory/totalSets;
                
                ObjectIndexTuple stats(alt.getNetworkName(),alt.getNumVertices(),alt.getNumEdges(),
                                       constants::OBJ_IDX_LNDMRK,processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[l],setSize);
                stats.setAdditionalFields(specialFields);

                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
    }
    std::cout << "Object Lists successfully created for " << alt.getNetworkName() << std::endl;
}

void LMKExperimentsCommand::buildNVD(Graph& graph, std::size_t numSets, std::vector<double> objDensities, std::vector<std::string> objTypes, 
                                  std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string> specialFields)
{
    int setSize = 0, setVariable;
    std::string setType;
    double setDensity;
    std::string objSetFile, objIdxOutputFile;
    
    StopWatch sw;
    double totalTime = 0, totalMemory = 0;
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            const std::size_t totalSets = utility::getNumObjectSets(objDensities[j],numSets);
            for (std::size_t l = 0; l < objVariable.size(); ++l) {
                totalTime = 0;
                totalMemory = 0;
                
                for (std::size_t k = 0; k < totalSets; ++k) {
                    objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                    std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                    this->checkObjectSetHeader(objSetFile,objTypes[i],objDensities[j],setType,setDensity);

                    sw.start();
                    NVDQuadtree nvd(setType,setDensity,setVariable,setSize);
                    nvd.buildNVD(graph,objectVertices,graph.mortonCodes,1);
                    sw.stop();
                    totalTime += sw.getTimeUs();
                    totalMemory += nvd.computeIndexSize()*1024; // in KB

                    objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k),constants::OBJ_IDX_NVD,parameterMap);
                    serialization::outputIndexToBinaryFile<NVDQuadtree>(nvd,objIdxOutputFile);
                    
                }
                
                // Average stats over the object sets and then write to stats file
                double processingTimeMs = totalTime/totalSets;
                double memoryUsage = totalMemory/totalSets;
                
                ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),
                                       constants::OBJ_IDX_NVD,processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[l],setSize);
                stats.setAdditionalFields(specialFields);

                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
    }
    std::cout << "Network Voronoi Diagrams successfully created for " << graph.getNetworkName() << std::endl;
}

void LMKExperimentsCommand::buildCOLT(Graph& graph, std::string sltreeIdxFile, int objListSize, std::size_t numSets, std::vector<double> objDensities, 
                                      std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                                      std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string> specialFields)
{
    SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(sltreeIdxFile);
    
    int setSize = 0, setVariable;
    std::string setType;
    double setDensity;
    std::string objSetFile, objIdxOutputFile, objIdxOutputFilePrefix;
    
    StopWatch sw;
    double totalTime = 0, totalMemory = 0;
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            const std::size_t totalSets = utility::getNumObjectSets(objDensities[j],numSets);
            for (std::size_t l = 0; l < objVariable.size(); ++l) {
                totalTime = 0;
                totalMemory = 0;
                
                for (std::size_t k = 0; k < totalSets; ++k) {
                    objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(sltree.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                    std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                    this->checkObjectSetHeader(objSetFile,objTypes[i],objDensities[j],setType,setDensity);

                    sw.start();
                    COLTree colt(setType,setDensity,setVariable,setSize);
                    colt.buildCOLTTopDown(sltree,objectVertices,objListSize);
                    sw.stop();
                    totalTime += sw.getTimeUs();
                    totalMemory += colt.computeIndexSize()*1024; // in KB

                    objIdxOutputFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                    objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxOutputFilePrefix,constants::OBJ_IDX_COLTREE,parameterMap);
                    serialization::outputIndexToBinaryFile<COLTree>(colt,objIdxOutputFile);
                    
                }
                
                // Average stats over the object sets and then write to stats file
                double processingTimeMs = totalTime/totalSets;
                double memoryUsage = totalMemory/totalSets;
                
                ObjectIndexTuple stats(sltree.getNetworkName(),sltree.getNumVertices(),sltree.getNumEdges(),
                                       constants::OBJ_IDX_COLTREE,processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[l],setSize);
                stats.setAdditionalFields(specialFields);
                stats.addSupplementaryFields("listsize",std::to_string(objListSize));
                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
    }
    std::cout << "COLTree indexes successfully created for " << sltree.getNetworkName() << std::endl;
}

void LMKExperimentsCommand::buildRealWorldObjectIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile, std::string rwObjSetNames)
{
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);

    bool buildOccList = parameterMap["occ_list"] == "1";
    bool buildRtree = parameterMap["rtree"] == "1";
    bool buildObjList = parameterMap["obj_list"] == "1";
    bool buildNVD = parameterMap["nvd"] == "1";
    bool buildCOLT = parameterMap["colt"] == "1";

    std::vector<std::string> rwObjSetNamesVec = utility::splitByDelim(rwObjSetNames,',');
    if (rwObjSetNamesVec.size() == 0) {
        std::cerr << "No real-world object set names provided!\n";
        std::exit(1);
    }

    // The POI sets are those of the network whose text files the graph was built from (rw_network in the parameters), which is the
    // graph's own name unless the graph was given another one (e.g. a second dataset with other edge weights)
    std::string rwNetworkName = parameterMap.find("rw_network") != parameterMap.end() ? parameterMap["rw_network"] : graph.getNetworkName();
    std::string rwPOIDirectory = filePathPrefix + "/real_world_pois/" + rwNetworkName + "/";
    
    StopWatch sw;
    double totalTime = 0, totalMemory = 0;
    int setSize, setVariable = 0, totalSetSize;
    std::string setType;
    double setDensity = 0;
    
    // Copy text files to object index folder with expected format 
    for (std::size_t j = 0; j < rwObjSetNamesVec.size(); ++j) {
        std::string objSetFile = rwPOIDirectory + rwObjSetNamesVec[j] + ".txt";
        std::string objIdxOutputFile = filePathPrefix + "/obj_indexes/" + graph.getNetworkName() + "_" + rwObjSetNamesVec[j] + ".txt";
        std::ifstream  input(objSetFile, std::ios::binary);
        std::ofstream  output(objIdxOutputFile, std::ios::binary);
        output << input.rdbuf();
    }
    
    if (buildOccList) {
        int fanout = std::stoi(parameterMap["gtree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["gtree_maxleafsize"]);
        if (fanout < 2 || maxLeafSize < 32) {
            std::cerr << "Invalid Gtree parameters!\n";
            exit(1);
        }
        const std::vector<std::string> gtreeKeys = {"fanout","maxleafsize"};
        const std::vector<std::string> gtreeValues = {parameterMap["gtree_fanout"],parameterMap["gtree_maxleafsize"]};
        std::string gtreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,gtreeKeys,gtreeValues);
        Gtree gtree = serialization::getIndexFromBinaryFile<Gtree>(gtreeIdxFile);
        totalTime = 0;
        totalMemory = 0;
        totalSetSize = 0;
        for (std::size_t j = 0; j < rwObjSetNamesVec.size(); ++j) {
            std::string objSetFile = rwPOIDirectory + rwObjSetNamesVec[j] + ".txt";
            std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
            
            sw.start();
            OccurenceList occList(setType,setDensity,setVariable,setSize);
            for (auto objIt = objectVertices.begin(); objIt != objectVertices.end(); ++objIt) {
                // Find the Gtree leaf index for this object and add it to occurence list
                // for that leaf (create list if it doesn't exist)
                int leafIdx = gtree.getLeafIndex(*objIt);
                occList.addLeafOccurence(leafIdx,*objIt);
                
                // Propagate this to occurence lists of parents of leaf node
                int childIdx = leafIdx;
                int parentIdx = gtree.getParentIndex(leafIdx);
                while (parentIdx != -1) {
                    occList.addParentOccurence(parentIdx, childIdx);
                    
                    // Go up a level (until we reach root)
                    childIdx = parentIdx;
                    parentIdx = gtree.getParentIndex(childIdx);
                }
            }
            sw.stop();
            totalTime += sw.getTimeUs();
            totalMemory += occList.computeIndexSize()*1024; // in KB
            totalSetSize += occList.getObjSetSize();
            
            // Serialize object index to output file
            std::string objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(graph.getNetworkName() + "_" + rwObjSetNamesVec[j],constants::OBJ_IDX_GTREE,parameterMap);
            serialization::outputIndexToBinaryFile<OccurenceList>(occList,objIdxOutputFile);
        }
        // Collect stats and return to output to file
        ObjectIndexTuple stats(gtree.getNetworkName(),gtree.getNumVertices(),gtree.getNumEdges(),
                            constants::OBJ_IDX_GTREE,totalTime,totalMemory,"all_rw_sets",setDensity,setVariable,totalSetSize);

        this->outputCommandStats(statsOutputFile,stats.getTupleString());
    }
    
    if (buildRtree) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() != 1) {
            std::cerr << "One and only one R-tree branch factor can be provided in LMK experiments!\n";
            exit(1);
        }
        totalTime = 0;
        totalMemory = 0;
        totalSetSize = 0;
        for (std::size_t j = 0; j < rwObjSetNamesVec.size(); ++j) {
            // Note: We don't bother setting branch factor in file name for real-world POIs
            std::string objSetFile = rwPOIDirectory + rwObjSetNamesVec[j] + ".txt";
            std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
            
            sw.start();
            std::vector<CoordinatePair> objectCoords;
            for (std::size_t i = 0; i < objectVertices.size(); ++i) {
                CoordinatePair objectCoordPair;
                graph.getCoordinates(objectVertices[i],objectCoordPair.first,objectCoordPair.second);
                objectCoords.push_back(objectCoordPair);
            }
            PackedRtree rtree(branchFactors[0],setType,setDensity,setVariable,setSize);
            rtree.bulkLoad(objectVertices,objectCoords);
            sw.stop();
            totalTime += sw.getTimeUs();
            totalMemory += rtree.computeIndexSize()*1024; // in KB
            totalSetSize += rtree.getObjSetSize();

            // Serialize object index to output file
            std::string objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(graph.getNetworkName() + "_" + rwObjSetNamesVec[j],constants::OBJ_IDX_RTREE,parameterMap);
            serialization::outputIndexToBinaryFile<PackedRtree>(rtree,objIdxOutputFile);
        }
        // Collect stats and return to output to file
        ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),
                            constants::OBJ_IDX_RTREE,totalTime,totalMemory,"all_rw_sets",setDensity,setVariable,totalSetSize);
        stats.addSupplementaryFields("rtree_branch_factor",std::to_string(branchFactors[0]));

        this->outputCommandStats(statsOutputFile,stats.getTupleString());
    }
    
    if (buildObjList) {
        int numLandmarks = std::stoi(parameterMap["alt_numlandmarks"]);
        if (numLandmarks < 1) {
            std::cerr << "Invalid ALT parameters!\n";
            exit(1);
        }
        const std::vector<std::string> altKeys = {"landmarks","type"};
        const std::vector<std::string> altValues = {parameterMap["alt_numlandmarks"],parameterMap["alt_landmarktype"]};
        std::string altIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,altKeys,altValues);
        ALT alt = serialization::getIndexFromBinaryFile<ALT>(altIdxFile);
        totalTime = 0;
        totalMemory = 0;
        totalSetSize = 0;
        for (std::size_t j = 0; j < rwObjSetNamesVec.size(); ++j) {
            std::string objSetFile = rwPOIDirectory + rwObjSetNamesVec[j] + ".txt";
            std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
            
            sw.start();
            ObjectList ol(setType,setDensity,setVariable,setSize);
            ol.buildObjectList(alt,objectVertices);
            sw.stop();
            totalTime += sw.getTimeUs();
            totalMemory += ol.computeIndexSize()*1024; // in KB
            totalSetSize += ol.getObjSetSize();
            
            // Serialize object index to output file
            std::string objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(graph.getNetworkName() + "_" + rwObjSetNamesVec[j],constants::OBJ_IDX_LNDMRK,parameterMap);
            serialization::outputIndexToBinaryFile<ObjectList>(ol,objIdxOutputFile);
        }
        // Collect stats and return to output to file
        ObjectIndexTuple stats(alt.getNetworkName(),alt.getNumVertices(),alt.getNumEdges(),
                                constants::OBJ_IDX_LNDMRK,totalTime,totalMemory,"all_rw_sets",setDensity,setVariable,totalSetSize);

        this->outputCommandStats(statsOutputFile,stats.getTupleString());
    }
    
    if (buildNVD) {
        totalTime = 0;
        totalMemory = 0;
        totalSetSize = 0;
        for (std::size_t j = 0; j < rwObjSetNamesVec.size(); ++j) {
            std::string objSetFile = rwPOIDirectory + rwObjSetNamesVec[j] + ".txt";
            std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
            
            sw.start();
            NVDQuadtree nvd(setType,setDensity,setVariable,setSize);
            nvd.buildNVD(graph,objectVertices,graph.mortonCodes,1);
            sw.stop();
            totalTime += sw.getTimeUs();
            totalMemory += nvd.computeIndexSize()*1024; // in KB
            totalSetSize += nvd.getObjSetSize();
            
            // Serialize object index to output file
            std::string objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(graph.getNetworkName() + "_" + rwObjSetNamesVec[j],constants::OBJ_IDX_NVD,parameterMap);
            serialization::outputIndexToBinaryFile<NVDQuadtree>(nvd,objIdxOutputFile);
        }
        // Collect stats and return to output to file
        ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),
                                constants::OBJ_IDX_NVD,totalTime,totalMemory,"all_rw_sets",setDensity,setVariable,totalSetSize);

        this->outputCommandStats(statsOutputFile,stats.getTupleString());
    }
    
    if (buildCOLT) {
        int fanout = std::stoi(parameterMap["sltree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["sltree_maxleafsize"]);
        std::size_t numNodeLandmarks = std::stoi(parameterMap["sltree_numlmks"]);
        std::string landmarkType = parameterMap["sltree_lmktype"];
        if (fanout < 2 || maxLeafSize < 32 || numNodeLandmarks < 1 || landmarkType == "") {
            std::cerr << "Invalid SULTree parameters!\n";
            exit(1);    
        }
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {parameterMap["sltree_fanout"],parameterMap["sltree_maxleafsize"],parameterMap["sltree_numlmks"],parameterMap["sltree_lmktype"]};
        std::string sltreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(sltreeIdxFile);
        
        std::vector<int> listSizes = utility::getIngetersFromStringList(parameterMap["colt_listsize"],':',0);
        if (listSizes.size() < 1) {
            std::cerr << "Invalid number of list size for COLTree!\n";
            exit(1);
        } else if (listSizes.size() > 1) {
            std::cerr << "Only first list size will be used for COLTree construction!\n";
        }
        
        totalTime = 0;
        totalMemory = 0;
        totalSetSize = 0;
        for (std::size_t j = 0; j < rwObjSetNamesVec.size(); ++j) {
            std::string objSetFile = rwPOIDirectory + rwObjSetNamesVec[j] + ".txt";
            std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
            
            sw.start();
            COLTree colt(setType,setDensity,setVariable,setSize);
            colt.buildCOLTTopDown(sltree,objectVertices,listSizes[0]);
            sw.stop();
            totalTime += sw.getTimeUs();
            totalMemory += colt.computeIndexSize()*1024; // in KB         
            totalSetSize += colt.getObjSetSize();   
            sw.reset();
            
            // Serialize object index to output file
            std::string objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(graph.getNetworkName() + "_" + rwObjSetNamesVec[j],constants::OBJ_IDX_COLTREE,parameterMap);
            serialization::outputIndexToBinaryFile<COLTree>(colt,objIdxOutputFile);
        }
        // Collect stats and return to output to file
        ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),
                               constants::OBJ_IDX_COLTREE,totalTime,totalMemory,"all_rw_sets",setDensity,setVariable,totalSetSize);

        this->outputCommandStats(statsOutputFile,stats.getTupleString());
    }
    
    std::cout << "Object indexes for real-world object sets successfully created for " << graph.getNetworkName() << std::endl;

}

void LMKExperimentsCommand::writeBlankLines(std::string fileName)
{
    if (fileName == "") {
        return; // A switched-off query list has no stats file
    }
    std::ofstream fs;
    fs.open(fileName, std::ios_base::app | std::ios_base::out);
    fs << "\n\n"; //  GNU Plot requires double line separation between data points for different plot lines
    fs.close();
}

void LMKExperimentsCommand::runQueries2(std::string bgrFileName, std::string queryVertexFile, std::string kValues, std::string parameters, std::size_t numSets, 
                                        std::string objDensities, std::string objTypes, std::string objVariable, std::string rwObjSetNames, std::string filePathPrefix)
{
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);

    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    
    std::vector<VertexID> queryVertices = utility::getPointSetFromFile(queryVertexFile);
    if (queryVertices.size() == 0) {
        std::cerr << "No query points were provided!\n";
        exit(1);    
    }
    if (numSets == 0) {
        std::cerr << "No object sets were provided!\n";
        exit(1);    
    }
    
    // Get Default Values for Variables
    unsigned int defaultK = std::stoul(parameterMap["knn_def_k"]);
    double defaultDensity = std::stod(parameterMap["knn_def_density"]);
    std::string defaultObjType = parameterMap["knn_def_type"];
    int defaultObjVar = std::stoi(parameterMap["knn_def_var"]);
    if (defaultObjType == "") {
        std::cerr << "Invalid default object set type of " << defaultObjType << "!\n";
        std::exit(1);    
    }
    
    /*
     * Create Query Sets for Each Experiment
     */
    
    // Each query list runs unless its switch in the parameters is off (utility::isQueryListOn). A switched-off list is not built, and
    // its stats file name is left empty, so that neither its runs nor the blank separator lines write to its stats file
    bool runVaryingK = utility::isQueryListOn(parameterMap,"knn_varyingk");
    bool runVaryingDensity = utility::isQueryListOn(parameterMap,"knn_varyingdensity");
    
    // Varying K Experiment
    std::vector<std::vector<KnnQuery>> varyingKQuerySets;

    std::vector<std::string> strKValuesVec = utility::splitByDelim(kValues,',');
    std::vector<unsigned int> kValuesVec;
    for(std::size_t i = 0; i < strKValuesVec.size(); ++i) {
        unsigned int k = std::stoul(strKValuesVec[i]);
        if (k > 0) {
            kValuesVec.push_back(k);
        } else {
            std::cerr << "Invalid k value " << k << " in list provided!\n";
            std::exit(1);    
        }
    }
    
    if (runVaryingK) {
        for (std::size_t i = 0; i < kValuesVec.size(); ++i) {
            std::vector<KnnQuery> kQuerySet;
            for (std::size_t j = 0; j < numSets; ++j) {
                std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
                for (std::size_t k = 0; k < queryVertices.size(); ++k) {
                    KnnQuery query(queryVertices[k],kValuesVec[i],defaultDensity,defaultObjType,defaultObjVar,j,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*defaultDensity));
                    kQuerySet.push_back(query);
                }
            }
            varyingKQuerySets.push_back(kQuerySet);
        }
    }
    
    if (parameterMap.find("knn_varyingk_stats_file") == parameterMap.end() || parameterMap["knn_varyingk_stats_file"] == "") {
        std::cerr << "Varying k stats file not provided in parameters";
        std::exit(1);
    }
    std::string varyingKStatsFile = runVaryingK ? filePathPrefix + "/stats/" + parameterMap["knn_varyingk_stats_file"] : "";
    
    // Varying Density Experiment
    std::vector<std::vector<KnnQuery>> varyingDensityQuerySets;

    std::vector<std::string> strObjDensitiesVec = utility::splitByDelim(objDensities,',');
    std::vector<double> objDensitiesVec;
    for(std::size_t i = 0; i < strObjDensitiesVec.size(); ++i) {
        double density = std::stod(strObjDensitiesVec[i]);
        if (density > 0) {
            objDensitiesVec.push_back(density);
        } else {
            std::cerr << "Invalid density in list provided!\n";
            exit(1);    
        }
    }
    
    if (runVaryingDensity) {
        for (std::size_t i = 0; i < objDensitiesVec.size(); ++i) {
            std::vector<KnnQuery> densityQuerySet;
            for (std::size_t j = 0; j < numSets; ++j) {
                std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,objDensitiesVec[i],defaultObjVar,j);
                for (std::size_t k = 0; k < queryVertices.size(); ++k) {
                    KnnQuery query(queryVertices[k],defaultK,objDensitiesVec[i],defaultObjType,defaultObjVar,j,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*objDensitiesVec[i]));
                    densityQuerySet.push_back(query);
                }
            }
            varyingDensityQuerySets.push_back(densityQuerySet);
        }
    }
    
    if (parameterMap.find("knn_varyingdensity_stats_file") == parameterMap.end() || parameterMap["knn_varyingdensity_stats_file"] == "") {
        std::cerr << "Varying density stats file not provided in parameters";
        std::exit(1);
    }
    std::string varyingDensityStatsFile = runVaryingDensity ? filePathPrefix + "/stats/" + parameterMap["knn_varyingdensity_stats_file"] : "";
    
    // Varying Real-World Object Set Experiment
    std::vector<std::vector<KnnQuery>> varyingRWObjSetQuerySets;

    std::vector<std::string> rwObjSetNamesVec = utility::splitByDelim(rwObjSetNames,',');
    bool runRealWorld = parameterMap["knn_rw"] == "1";
    if (runRealWorld) {
        if (rwObjSetNamesVec.size() == 0) {
            std::cerr << "No real-world object set names provided!\n";
            std::exit(1);
        }
        // Each real-world object set is stored under its name, so a name may not repeat or be that of a synthetic set
        std::unordered_set<std::string> objSetFilenamePrefixes;
        for (std::size_t i = 0; i < objDensitiesVec.size(); ++i) {
            for (std::size_t j = 0; j < numSets; ++j) {
                objSetFilenamePrefixes.insert(utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,objDensitiesVec[i],defaultObjVar,j));
            }
        }
        for (std::size_t i = 0; i < rwObjSetNamesVec.size(); ++i) {
            if (!objSetFilenamePrefixes.insert(graph.getNetworkName() + "_" + rwObjSetNamesVec[i]).second) {
                std::cerr << "Real-world object set " << rwObjSetNamesVec[i] << " already exists, either duplicated or name-clash with synthetic set!\n";
                std::exit(1);
            }
        }
        // Each real-world POI type has only one object set, so we use a separate query vertex file with
        // (number of object sets x number of query vertices) query vertices to run the same number of queries
        // as each synthetic object set experiment
        if (parameterMap.find("knn_rw_qv_file") == parameterMap.end() || parameterMap["knn_rw_qv_file"] == "") {
            std::cerr << "Real-world query vertex file not provided in parameters!\n";
            std::exit(1);
        }
        std::vector<VertexID> rwQueryVertices = utility::getPointSetFromFile(filePathPrefix + "/query_sets/" + parameterMap["knn_rw_qv_file"]);
        if (rwQueryVertices.size() != numSets*queryVertices.size()) {
            std::cerr << "Real-world query vertex file has " << rwQueryVertices.size() << " query vertices but " << numSets*queryVertices.size()
                      << " (number of object sets x number of query vertices) are required!\n";
            std::exit(1);
        }

        for (std::size_t i = 0; i < rwObjSetNamesVec.size(); ++i) {
            std::vector<KnnQuery> rwPOITypeQuerySet;
            std::string objIdxFilenamePrefix = graph.getNetworkName() + "_" + rwObjSetNamesVec[i];
            for (std::size_t k = 0; k < rwQueryVertices.size(); ++k) {
                KnnQuery query(rwQueryVertices[k],defaultK,0,rwObjSetNamesVec[i],0,0,objIdxFilenamePrefix,0);
                rwPOITypeQuerySet.push_back(query);
            }
            varyingRWObjSetQuerySets.push_back(rwPOITypeQuerySet);
        }
    }

    if (parameterMap.find("knn_varyingrwset_stats_file") == parameterMap.end() || parameterMap["knn_varyingrwset_stats_file"] == "") {
        std::cerr << "Varying real-world stats file not provided in parameters";
        std::exit(1);
    }
    std::string varyingRWObjSetStatsFile = filePathPrefix + "/stats/" + parameterMap["knn_varyingrwset_stats_file"];
    
    bool queryINE = parameterMap["ine"] == "1";
    bool queryGtree = parameterMap["gtree"] == "1";
    bool queryIERPHL = parameterMap["ier_phl"] == "1";
    bool queryOLPHL = parameterMap["ol_phl"] == "1";
    bool queryNVDPHL = parameterMap["nvd_phl"] == "1";
    bool queryLLBPHL = parameterMap["llb_phl"] == "1";
    bool queryCOLTPHL = parameterMap["colt_phl"] == "1";
    
    std::string setType;
    double setDensity;
    int setSize, setVariable;

    if (queryINE) {
        INE ine;

        // Load all object sets into memory
        std::unordered_map<std::string,std::vector<bool>> fileNameToObjIdx;
        for (std::size_t i = 0; i < objDensitiesVec.size(); ++i) {
            for (std::size_t j = 0; j < numSets; ++j) {
                std::string objSetFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,objDensitiesVec[i],defaultObjVar,j);
                std::string objSetFilePath = filePathPrefix + "/obj_indexes/" + objSetFilenamePrefix + ".txt";
                std::vector<VertexID> objSetVertices = utility::getPointSetFromFile(objSetFilePath,setType,setDensity,setSize,setVariable);
                std::vector<bool> vertexObjectStatus = ine.buildVertexObjectStatusVector(graph,objSetVertices);
                // Do we need to check it hasn't already been inserted?
                fileNameToObjIdx[objSetFilenamePrefix] = vertexObjectStatus;
            }
        }
        
        if (runRealWorld) {
            for (std::size_t i = 0; i < rwObjSetNamesVec.size(); ++i) {
                std::string objSetFilenamePrefix = graph.getNetworkName() + "_" + rwObjSetNamesVec[i];
                std::string objSetFilePath = filePathPrefix + "/obj_indexes/" + objSetFilenamePrefix + ".txt";
                std::vector<VertexID> objSetVertices = utility::getPointSetFromFile(objSetFilePath,setType,setDensity,setSize,setVariable);
                std::vector<bool> vertexObjectStatus = ine.buildVertexObjectStatusVector(graph,objSetVertices);
                fileNameToObjIdx[objSetFilenamePrefix] = vertexObjectStatus;
            }
        }

        this->runINEQueries2(graph,varyingKQuerySets,fileNameToObjIdx,varyingKStatsFile,"varyingk");
        this->runINEQueries2(graph,varyingDensityQuerySets,fileNameToObjIdx,varyingDensityStatsFile,"varyingdensity");
        if (runRealWorld) {
            this->runINEQueries2(graph,varyingRWObjSetQuerySets,fileNameToObjIdx,varyingRWObjSetStatsFile,"varyingrwset");
        }

        this->writeBlankLines(varyingKStatsFile);
        this->writeBlankLines(varyingDensityStatsFile);
    }
    
    if (queryGtree) {
        int fanout = std::stoi(parameterMap["gtree_fanout"]);
        unsigned int maxLeafSize = std::stoul(parameterMap["gtree_maxleafsize"]);
        if (fanout < 2 || maxLeafSize < 32) {
            std::cerr << "Invalid Gtree parameters!\n";
            exit(1);    
        }
        const std::vector<std::string> gtreeKeys = {"fanout","maxleafsize"};
        const std::vector<std::string> gtreeValues = {parameterMap["gtree_fanout"],parameterMap["gtree_maxleafsize"]};
        std::string gtreeIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,gtreeKeys,gtreeValues);
        
        Gtree gtree = serialization::getIndexFromBinaryFile<Gtree>(gtreeIdxFilePath);

        // Load all object indexes into memory
        std::unordered_map<std::string,OccurenceList> fileNameToObjIdx;
        for (std::size_t i = 0; i < objDensitiesVec.size(); ++i) {
            for (std::size_t j = 0; j < numSets; ++j) {
                std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,objDensitiesVec[i],defaultObjVar,j);
                std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_GTREE,parameterMap);
                OccurenceList objIdx = serialization::getIndexFromBinaryFile<OccurenceList>(objIdxFilePath);
                fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
            }
        }
        
        if (runRealWorld) {
            for (std::size_t i = 0; i < rwObjSetNamesVec.size(); ++i) {
                std::string objIdxFilenamePrefix = graph.getNetworkName() + "_" + rwObjSetNamesVec[i];
                std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_GTREE,parameterMap);
                OccurenceList objIdx = serialization::getIndexFromBinaryFile<OccurenceList>(objIdxFilePath);
                fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
            }
        }
        
        this->runGtreeQueries2(graph,gtree,varyingKQuerySets,fileNameToObjIdx,varyingKStatsFile,"varyingk");
        this->runGtreeQueries2(graph,gtree,varyingDensityQuerySets,fileNameToObjIdx,varyingDensityStatsFile,"varyingdensity");
        if (runRealWorld) {
            this->runGtreeQueries2(graph,gtree,varyingRWObjSetQuerySets,fileNameToObjIdx,varyingRWObjSetStatsFile,"varyingrwset");
        }
    
        this->writeBlankLines(varyingKStatsFile);
        this->writeBlankLines(varyingDensityStatsFile);
    }    
    
    if (queryIERPHL || queryOLPHL || queryNVDPHL || queryLLBPHL || queryCOLTPHL) {
        std::string phlIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PHL_CMD,{},{});
        PrunedHighwayLabeling phl;
        phl.LoadLabel(phlIdxFilePath.c_str());
        
        if (queryIERPHL) {

            // Load all object indexes into memory
            std::unordered_map<std::string,PackedRtree> fileNameToObjIdx;
            for (std::size_t i = 0; i < objDensitiesVec.size(); ++i) {
                for (std::size_t j = 0; j < numSets; ++j) {
                    std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,objDensitiesVec[i],defaultObjVar,j);
                    std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_RTREE,parameterMap);
                    PackedRtree objIdx = serialization::getIndexFromBinaryFile<PackedRtree>(objIdxFilePath);
                    fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                }
            }
        
            if (runRealWorld) {
                for (std::size_t i = 0; i < rwObjSetNamesVec.size(); ++i) {
                    std::string objIdxFilenamePrefix = graph.getNetworkName() + "_" + rwObjSetNamesVec[i];
                    std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_RTREE,parameterMap);
                    PackedRtree objIdx = serialization::getIndexFromBinaryFile<PackedRtree>(objIdxFilePath);
                    fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                }
            }

            this->runIERPHL(graph,phl,varyingKQuerySets,fileNameToObjIdx,varyingKStatsFile,"varyingk");
            this->runIERPHL(graph,phl,varyingDensityQuerySets,fileNameToObjIdx,varyingDensityStatsFile,"varyingdensity");
            if (runRealWorld) {
                this->runIERPHL(graph,phl,varyingRWObjSetQuerySets,fileNameToObjIdx,varyingRWObjSetStatsFile,"varyingrwset");
            }

            this->writeBlankLines(varyingKStatsFile);
            this->writeBlankLines(varyingDensityStatsFile);
        }
        
        if (queryOLPHL || queryNVDPHL || queryLLBPHL || queryCOLTPHL) {
            const std::vector<std::string> altKeys = {"landmarks","type"};
            const std::vector<std::string> altValues = {parameterMap["alt_numlandmarks"],parameterMap["alt_landmarktype"]};
            std::string altIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,altKeys,altValues);
            ALT alt = serialization::getIndexFromBinaryFile<ALT>(altIdxFilePath);
            
            if (queryOLPHL) {
                // Load all object indexes into memory
                std::unordered_map<std::string,ObjectList> fileNameToObjIdx;
                for (std::size_t i = 0; i < objDensitiesVec.size(); ++i) {
                    for (std::size_t j = 0; j < numSets; ++j) {
                        std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,objDensitiesVec[i],defaultObjVar,j);
                        std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_LNDMRK,parameterMap);
                        ObjectList objIdx = serialization::getIndexFromBinaryFile<ObjectList>(objIdxFilePath);
                        fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                    }
                }

                if (runRealWorld) {
                    for (std::size_t i = 0; i < rwObjSetNamesVec.size(); ++i) {
                        std::string objIdxFilenamePrefix = graph.getNetworkName() + "_" + rwObjSetNamesVec[i];
                        std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_LNDMRK,parameterMap);
                        ObjectList objIdx = serialization::getIndexFromBinaryFile<ObjectList>(objIdxFilePath);
                        fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                    }
                }

                this->runOLPHL(graph,phl,alt,varyingKQuerySets,fileNameToObjIdx,varyingKStatsFile,"varyingk");
                this->runOLPHL(graph,phl,alt,varyingDensityQuerySets,fileNameToObjIdx,varyingDensityStatsFile,"varyingdensity");
                if (runRealWorld) {
                    this->runOLPHL(graph,phl,alt,varyingRWObjSetQuerySets,fileNameToObjIdx,varyingRWObjSetStatsFile,"varyingrwset");
                }

                this->writeBlankLines(varyingKStatsFile);
                this->writeBlankLines(varyingDensityStatsFile);
            }
            
            if (queryNVDPHL) {
                // Load all object indexes into memory
                std::unordered_map<std::string,NVDQuadtree> fileNameToObjIdx;
                for (std::size_t i = 0; i < objDensitiesVec.size(); ++i) {
                    for (std::size_t j = 0; j < numSets; ++j) {
                        std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,objDensitiesVec[i],defaultObjVar,j);
                        std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_NVD,parameterMap);
                        NVDQuadtree objIdx = serialization::getIndexFromBinaryFile<NVDQuadtree>(objIdxFilePath);
                        fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                    }
                }

                if (runRealWorld) {
                    for (std::size_t i = 0; i < rwObjSetNamesVec.size(); ++i) {
                        std::string objIdxFilenamePrefix = graph.getNetworkName() + "_" + rwObjSetNamesVec[i];
                        std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_NVD,parameterMap);
                        NVDQuadtree objIdx = serialization::getIndexFromBinaryFile<NVDQuadtree>(objIdxFilePath);
                        fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                    }
                }

                this->runNVDPHL(graph,phl,alt,varyingKQuerySets,fileNameToObjIdx,varyingKStatsFile,"varyingk");
                this->runNVDPHL(graph,phl,alt,varyingDensityQuerySets,fileNameToObjIdx,varyingDensityStatsFile,"varyingdensity");
                if (runRealWorld) {
                    this->runNVDPHL(graph,phl,alt,varyingRWObjSetQuerySets,fileNameToObjIdx,varyingRWObjSetStatsFile,"varyingrwset");
                }

                this->writeBlankLines(varyingKStatsFile);
                this->writeBlankLines(varyingDensityStatsFile);
            }
            
            if (queryLLBPHL) {
                // Load all object indexes into memory
                std::unordered_map<std::string,std::vector<VertexID>> fileNameToObjIdx;
                for (std::size_t i = 0; i < objDensitiesVec.size(); ++i) {
                    for (std::size_t j = 0; j < numSets; ++j) {
                        std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,objDensitiesVec[i],defaultObjVar,j);
                        std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + objIdxFilenamePrefix + ".txt";
                        std::vector<VertexID> objIdx = utility::getPointSetFromFile(objIdxFilePath,setType,setDensity,setSize,setVariable);
                        fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                    }
                }
                if (runRealWorld) {
                    for (std::size_t i = 0; i < rwObjSetNamesVec.size(); ++i) {
                        std::string objIdxFilenamePrefix = graph.getNetworkName() + "_" + rwObjSetNamesVec[i];
                        std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + objIdxFilenamePrefix + ".txt";
                        std::vector<VertexID> objIdx = utility::getPointSetFromFile(objIdxFilePath,setType,setDensity,setSize,setVariable);
                        fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                    }
                }

                this->runLLBPHL(graph,phl,alt,varyingKQuerySets,fileNameToObjIdx,varyingKStatsFile,"varyingk");
                this->runLLBPHL(graph,phl,alt,varyingDensityQuerySets,fileNameToObjIdx,varyingDensityStatsFile,"varyingdensity");
                if (runRealWorld) {
                    this->runLLBPHL(graph,phl,alt,varyingRWObjSetQuerySets,fileNameToObjIdx,varyingRWObjSetStatsFile,"varyingrwset");
                }

                this->writeBlankLines(varyingKStatsFile);
                this->writeBlankLines(varyingDensityStatsFile);
            }
            
            if (queryCOLTPHL) {
                const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
                const std::vector<std::string> sltreeValues = {parameterMap["sltree_fanout"],parameterMap["sltree_maxleafsize"],parameterMap["sltree_numlmks"],parameterMap["sltree_lmktype"]};
                std::string sltreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
                StopWatch sw;
                sw.start();
                SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(sltreeIdxFile);
                sw.stop();
                std::cout << "SULTree load time = " << sw.getTimeMs() << "\n";
        
                // Load all object indexes into memory
                std::unordered_map<std::string,COLTree> fileNameToObjIdx;
                for (std::size_t i = 0; i < objDensitiesVec.size(); ++i) {
                    for (std::size_t j = 0; j < numSets; ++j) {
                        std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,objDensitiesVec[i],defaultObjVar,j);
                        std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_COLTREE,parameterMap);
                        COLTree objIdx = serialization::getIndexFromBinaryFile<COLTree>(objIdxFilePath);
                        fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                    }
                }

                if (runRealWorld) {
                    for (std::size_t i = 0; i < rwObjSetNamesVec.size(); ++i) {
                        std::string objIdxFilenamePrefix = graph.getNetworkName() + "_" + rwObjSetNamesVec[i];
                        std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_COLTREE,parameterMap);
                        COLTree objIdx = serialization::getIndexFromBinaryFile<COLTree>(objIdxFilePath);
                        fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                    }
                }

                this->runCOLTPHL(graph,phl,alt,sltree,varyingKQuerySets,fileNameToObjIdx,varyingKStatsFile,"varyingk");
                this->runCOLTPHL(graph,phl,alt,sltree,varyingDensityQuerySets,fileNameToObjIdx,varyingDensityStatsFile,"varyingdensity");
                if (runRealWorld) {
                    this->runCOLTPHL(graph,phl,alt,sltree,varyingRWObjSetQuerySets,fileNameToObjIdx,varyingRWObjSetStatsFile,"varyingrwset");
                }

                this->writeBlankLines(varyingKStatsFile);
                this->writeBlankLines(varyingDensityStatsFile);
            }
        }
    }
    
}

void LMKExperimentsCommand::runSingleQueries2(std::string bgrFileName, std::string queryVertexFile, std::string kValues, std::string parameters, 
                                              std::size_t numSets, std::string objDensities, std::string objTypes, std::string objVariable, 
                                              std::string filePathPrefix, std::string statsOutputFile)
{
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);

    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    
    std::vector<VertexID> queryVertices = utility::getPointSetFromFile(queryVertexFile);
    if (queryVertices.size() == 0) {
        std::cerr << "No query points were provided!\n";
        exit(1);    
    }
    if (numSets == 0) {
        std::cerr << "No object sets were provided!\n";
        exit(1);    
    }
    
    // Get Default Values for Variables
    unsigned int defaultK = std::stoul(parameterMap["knn_def_k"]);
    double defaultDensity = std::stod(parameterMap["knn_def_density"]);
    std::string defaultObjType = parameterMap["knn_def_type"];
    int defaultObjVar = std::stoi(parameterMap["knn_def_var"]);
    if (defaultObjType == "") {
        std::cerr << "Invalid default object set type of " << defaultObjType << "!\n";
        std::exit(1);    
    }

    /*
     * Create set of queries to run on each method
     */
    std::vector<std::vector<KnnQuery>> singleQuerySet;
    std::vector<KnnQuery> querySet;
    for (std::size_t i = 0; i < numSets; ++i) {
        std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,i);
        for (std::size_t j = 0; j < queryVertices.size(); ++j) {
            KnnQuery query(queryVertices[j],defaultK,defaultDensity,defaultObjType,defaultObjVar,i,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*defaultDensity));
            querySet.push_back(query);
        }
    }
    singleQuerySet.push_back(querySet);
    
    bool queryINE = parameterMap["ine"] == "1";
    bool queryGtree = parameterMap["gtree"] == "1";
    bool queryIERPHL = parameterMap["ier_phl"] == "1";
    bool queryOLPHL = parameterMap["ol_phl"] == "1";
    bool queryNVDPHL = parameterMap["nvd_phl"] == "1";
    bool queryLLBPHL = parameterMap["llb_phl"] == "1";
    bool queryCOLTPHL = parameterMap["colt_phl"] == "1";
    
    std::string setType;
    double setDensity;
    int setSize, setVariable;

    if (queryINE) {
        INE ine;

        // Load all object sets into memory
        std::unordered_map<std::string,std::vector<bool>> fileNameToObjIdx;
        for (std::size_t j = 0; j < numSets; ++j) {
            std::string objSetFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
            std::string objSetFile = filePathPrefix + "/obj_indexes/" + objSetFilenamePrefix + ".txt";
            std::vector<VertexID> objSetVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
            std::vector<bool> vertexObjectStatus = ine.buildVertexObjectStatusVector(graph,objSetVertices);
            // Do we need to check it hasn't already been inserted?
            fileNameToObjIdx[objSetFilenamePrefix] = vertexObjectStatus;
        }

        this->runINEQueries2(graph,singleQuerySet,fileNameToObjIdx,statsOutputFile,"single_query_set");
    }
    
    if (queryGtree) {
        int fanout = std::stoi(parameterMap["gtree_fanout"]);
        unsigned int maxLeafSize = std::stoul(parameterMap["gtree_maxleafsize"]);
        if (fanout < 2 || maxLeafSize < 32) {
            std::cerr << "Invalid Gtree parameters!\n";
            exit(1);    
        }
        const std::vector<std::string> gtreeKeys = {"fanout","maxleafsize"};
        const std::vector<std::string> gtreeValues = {parameterMap["gtree_fanout"],parameterMap["gtree_maxleafsize"]};
        std::string gtreeIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,gtreeKeys,gtreeValues);
        
        Gtree gtree = serialization::getIndexFromBinaryFile<Gtree>(gtreeIdxFilePath);

        // Load all object indexes into memory
        std::unordered_map<std::string,OccurenceList> fileNameToObjIdx;
        for (std::size_t j = 0; j < numSets; ++j) {
            std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
            std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_GTREE,parameterMap);
            OccurenceList objIdx = serialization::getIndexFromBinaryFile<OccurenceList>(objIdxFilePath);
            fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
        }
        
        this->runGtreeQueries2(graph,gtree,singleQuerySet,fileNameToObjIdx,statsOutputFile,"single_query_set");
    }
    
    if (queryIERPHL || queryOLPHL || queryNVDPHL || queryLLBPHL || queryCOLTPHL) {
        std::string phlIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PHL_CMD,{},{});
        PrunedHighwayLabeling phl;
        phl.LoadLabel(phlIdxFilePath.c_str());
        
        if (queryIERPHL) {

            // Load all object indexes into memory
            std::unordered_map<std::string,PackedRtree> fileNameToObjIdx;
            for (std::size_t j = 0; j < numSets; ++j) {
                std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
                std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_RTREE,parameterMap);
                PackedRtree objIdx = serialization::getIndexFromBinaryFile<PackedRtree>(objIdxFilePath);
                fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
            }
        
            this->runIERPHL(graph,phl,singleQuerySet,fileNameToObjIdx,statsOutputFile,"single_query_set");
        }
        
        if (queryOLPHL || queryNVDPHL || queryLLBPHL || queryCOLTPHL) {
            const std::vector<std::string> altKeys = {"landmarks","type"};
            const std::vector<std::string> altValues = {parameterMap["alt_numlandmarks"],parameterMap["alt_landmarktype"]};
            std::string altIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,altKeys,altValues);
            ALT alt = serialization::getIndexFromBinaryFile<ALT>(altIdxFilePath);
            
            if (queryOLPHL) {
                // Load all object indexes into memory
                std::unordered_map<std::string,ObjectList> fileNameToObjIdx;
                for (std::size_t j = 0; j < numSets; ++j) {
                    std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
                    std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_LNDMRK,parameterMap);
                    ObjectList objIdx = serialization::getIndexFromBinaryFile<ObjectList>(objIdxFilePath);
                    fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                }
            
                this->runOLPHL(graph,phl,alt,singleQuerySet,fileNameToObjIdx,statsOutputFile,"single_query_set");
            }
            
            if (queryNVDPHL) {
                // Load all object indexes into memory
                std::unordered_map<std::string,NVDQuadtree> fileNameToObjIdx;
                for (std::size_t j = 0; j < numSets; ++j) {
                    std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
                    std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_NVD,parameterMap);
                    NVDQuadtree objIdx = serialization::getIndexFromBinaryFile<NVDQuadtree>(objIdxFilePath);
                    fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                }
            
                this->runNVDPHL(graph,phl,alt,singleQuerySet,fileNameToObjIdx,statsOutputFile,"single_query_set");
            }
            
            if (queryLLBPHL) {
                // Load all object indexes into memory
                std::unordered_map<std::string,std::vector<VertexID>> fileNameToObjIdx;
                for (std::size_t j = 0; j < numSets; ++j) {
                    std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
                    std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + objIdxFilenamePrefix + ".txt";
                    std::vector<VertexID> objIdx = utility::getPointSetFromFile(objIdxFilePath,setType,setDensity,setSize,setVariable);
                    fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                }
            
                this->runLLBPHL(graph,phl,alt,singleQuerySet,fileNameToObjIdx,statsOutputFile,"single_query_set");
            }
            
            if (queryCOLTPHL) {
                const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
                const std::vector<std::string> sltreeValues = {parameterMap["sltree_fanout"],parameterMap["sltree_maxleafsize"],parameterMap["sltree_numlmks"],parameterMap["sltree_lmktype"]};
                std::string sltreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
                SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(sltreeIdxFile);

                // Load all object indexes into memory
                std::unordered_map<std::string,COLTree> fileNameToObjIdx;
                for (std::size_t j = 0; j < numSets; ++j) {
                    std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
                    std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_COLTREE,parameterMap);
                    COLTree objIdx = serialization::getIndexFromBinaryFile<COLTree>(objIdxFilePath);
                    fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                }
            
                this->runCOLTPHL(graph,phl,alt,sltree,singleQuerySet,fileNameToObjIdx,statsOutputFile,"single_query_set");
            }
        }
    }
}


void LMKExperimentsCommand::runINEQueries2(Graph& graph, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,std::vector<bool>>& objectIndexes, 
                                           std::string statsOutputFile, std::string desc, std::vector<std::string> specialFields) 
{
    INE ine;

    StopWatch sw;
    double totalQueryTime;
    int totalQueries;
    std::vector<VertexID> kNNs;
    std::vector<EdgeWeight> kNNDistances;
    
    // Warm-up: the first WARM_UP_QUERY_PERCENTAGE percent of the first list (at least one query) runs untimed, so that the cold start of a fresh
    // process or of freshly loaded indexes stays out of the timings (skipped in rn_obj_search_stats, whose statistics do not depend on it)
#if !defined(COLLECT_STATISTICS)
    if (querySets.size() > 0) {
        std::size_t numWarmUpQueries = utility::getNumWarmUpQueries(querySets[0].size());
        for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
            ine.getKNNs(graph,objectIndexes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertex,kNNs,kNNDistances);
            kNNs.clear();
            kNNDistances.clear();
        }
    }
#endif

    for (std::size_t i = 0; i < querySets.size(); ++i) {
        totalQueries = querySets[i].size();
        sw.reset();
        sw.start();
        for (std::size_t j = 0; j < querySets[i].size(); ++j) {
            ine.getKNNs(graph,objectIndexes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertex,kNNs,kNNDistances);
            kNNs.clear();
            kNNDistances.clear();
        }
        sw.stop();
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;
        
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::INE_KNN_QUERY,
                              querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].var,querySets[i][0].objSetSize,kNNs,kNNDistances);
        results.setAdditionalFields(specialFields);
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "INE kNN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}


void LMKExperimentsCommand::runGtreeQueries2(Graph& graph, Gtree& gtree, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,OccurenceList>& objectIndexes, 
                                             std::string statsOutputFile, std::string desc, std::vector<std::string> specialFields) 
{
    StopWatch sw;
    double totalQueryTime;
    int totalQueries;
    std::vector<VertexID> kNNs;
    std::vector<EdgeWeight> kNNDistances;
    
    Statistics knnStats;
    
    // Warm-up: the first WARM_UP_QUERY_PERCENTAGE percent of the first list (at least one query) runs untimed, so that the cold start of a fresh
    // process or of freshly loaded indexes stays out of the timings (skipped in rn_obj_search_stats, whose statistics do not depend on it)
#if !defined(COLLECT_STATISTICS)
    if (querySets.size() > 0) {
        std::size_t numWarmUpQueries = utility::getNumWarmUpQueries(querySets[0].size());
        for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
            gtree.getKNNs(objectIndexes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertex,kNNs,kNNDistances,graph);
            kNNs.clear();
            kNNDistances.clear();
        }
    }
#endif

    for (std::size_t i = 0; i < querySets.size(); ++i) {
        totalQueries = querySets[i].size();
#if defined(COLLECT_STATISTICS)
        knnStats.clear();
#endif
        sw.reset();
        sw.start();
        for (std::size_t j = 0; j < querySets[i].size(); ++j) {
            gtree.getKNNs(objectIndexes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertex,kNNs,kNNDistances,graph);
            kNNs.clear();
            kNNDistances.clear();
#if defined(COLLECT_STATISTICS)
            knnStats.mergeStatistics(gtree.stats);
#endif
        }
        sw.stop();
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;  
#if defined(COLLECT_STATISTICS)
        knnStats.normalizeStatistics(totalQueries);
#endif
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::GTREE_KNN_QUERY,
                            querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].var,querySets[i][0].objSetSize,kNNs,kNNDistances);
        results.setAdditionalFields(specialFields);
#if defined(COLLECT_STATISTICS)
        knnStats.populateTupleFields(results,0);
#endif
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "G-tree kNN queries successfully executed for " << gtree.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKExperimentsCommand::runIERPHL(Graph& graph, PrunedHighwayLabeling& phl, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,PackedRtree>& objectIndexes, 
                                      std::string statsOutputFile, std::string desc, std::vector<std::string> specialFields) 
{
    IER ier;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries;
    std::vector<VertexID> kNNs;
    std::vector<EdgeWeight> kNNDistances;

    Statistics knnStats;
                    
    // Warm-up: the first WARM_UP_QUERY_PERCENTAGE percent of the first list (at least one query) runs untimed, so that the cold start of a fresh
    // process or of freshly loaded indexes stays out of the timings (skipped in rn_obj_search_stats, whose statistics do not depend on it)
#if !defined(COLLECT_STATISTICS)
    if (querySets.size() > 0) {
        std::size_t numWarmUpQueries = utility::getNumWarmUpQueries(querySets[0].size());
        if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
            for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
                ier.getKNNsByPHL(objectIndexes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertex,kNNs,kNNDistances,graph,phl);
                kNNs.clear();
                kNNDistances.clear();
            }
        } else {
            for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
                ier.getKNNsByPHLTravelTimes(objectIndexes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertex,kNNs,kNNDistances,graph,phl);
                kNNs.clear();
                kNNDistances.clear();
            }
        }
    }
#endif

    for (std::size_t i = 0; i < querySets.size(); ++i) {
        totalQueries = querySets[i].size();
#if defined(COLLECT_STATISTICS)
        knnStats.clear();
#endif
        if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
            sw.reset();
            sw.start();
            for (std::size_t j = 0; j < querySets[i].size(); ++j) {
                ier.getKNNsByPHL(objectIndexes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertex,kNNs,kNNDistances,graph,phl);
                kNNs.clear();
                kNNDistances.clear();
#if defined(COLLECT_STATISTICS)
                knnStats.mergeStatistics(ier.stats);
#endif
            }
            sw.stop();
        } else {
            sw.reset();
            sw.start();
            for (std::size_t j = 0; j < querySets[i].size(); ++j) {
                ier.getKNNsByPHLTravelTimes(objectIndexes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertex,kNNs,kNNDistances,graph,phl);
                kNNs.clear();
                kNNDistances.clear();
#if defined(COLLECT_STATISTICS)
                knnStats.mergeStatistics(ier.stats);
#endif
            }
            sw.stop();
        }
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;  
#if defined(COLLECT_STATISTICS)
        knnStats.normalizeStatistics(totalQueries);
#endif
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::IER_PHL_KNN_QUERY,
                            querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].var,querySets[i][0].objSetSize,kNNs,kNNDistances);
        results.setAdditionalFields(specialFields);
#if defined(COLLECT_STATISTICS)
        knnStats.populateTupleFields(results,0);
#endif
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "IER-PHL kNN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKExperimentsCommand::runOLPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,ObjectList>& objectIndexes, 
                                     std::string statsOutputFile, std::string desc, std::vector<std::string> specialFields) 
{
    ILBR ilbr;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries;
    std::vector<VertexID> kNNs;
    std::vector<EdgeWeight> kNNDistances;

    Statistics knnStats;
    
    // Warm-up: the first WARM_UP_QUERY_PERCENTAGE percent of the first list (at least one query) runs untimed, so that the cold start of a fresh
    // process or of freshly loaded indexes stays out of the timings (skipped in rn_obj_search_stats, whose statistics do not depend on it)
#if !defined(COLLECT_STATISTICS)
    if (querySets.size() > 0) {
        std::size_t numWarmUpQueries = utility::getNumWarmUpQueries(querySets[0].size());
        for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
            ilbr.getKNNsByPHLImprovedByQueryLandmark(alt,objectIndexes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertex,kNNs,kNNDistances,phl);
            kNNs.clear();
            kNNDistances.clear();
        }
    }
#endif

    for (std::size_t i = 0; i < querySets.size(); ++i) {
        totalQueries = querySets[i].size();
#if defined(COLLECT_STATISTICS)
        knnStats.clear();
#endif
        sw.reset();
        sw.start();
        for (std::size_t j = 0; j < querySets[i].size(); ++j) {
            ilbr.getKNNsByPHLImprovedByQueryLandmark(alt,objectIndexes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertex,kNNs,kNNDistances,phl);
            kNNs.clear();
            kNNDistances.clear();
#if defined(COLLECT_STATISTICS)
            knnStats.mergeStatistics(ilbr.stats);
#endif
        }
        sw.stop();
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;  
#if defined(COLLECT_STATISTICS)
        knnStats.normalizeStatistics(totalQueries);
#endif
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::OL_ILBR_QLM_PHL_KNN_QUERY,
                            querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].var,querySets[i][0].objSetSize,kNNs,kNNDistances);
        results.setAdditionalFields(specialFields);
#if defined(COLLECT_STATISTICS)
        knnStats.populateTupleFields(results,0);
#endif
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "OL-PHL kNN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKExperimentsCommand::runNVDPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,NVDQuadtree>& objectIndexes, 
                                      std::string statsOutputFile, std::string desc, std::vector<std::string> specialFields) 
{
    ILBR ilbr;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries;
    std::vector<VertexID> kNNs;
    std::vector<EdgeWeight> kNNDistances;

    Statistics knnStats;
    
    // Warm-up: the first WARM_UP_QUERY_PERCENTAGE percent of the first list (at least one query) runs untimed, so that the cold start of a fresh
    // process or of freshly loaded indexes stays out of the timings (skipped in rn_obj_search_stats, whose statistics do not depend on it)
#if !defined(COLLECT_STATISTICS)
    if (querySets.size() > 0) {
        std::size_t numWarmUpQueries = utility::getNumWarmUpQueries(querySets[0].size());
        for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
            ilbr.getKNNByNVDQuadtreePHL(alt,objectIndexes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertex,kNNs,kNNDistances,phl,graph);
            kNNs.clear();
            kNNDistances.clear();
        }
    }
#endif

    for (std::size_t i = 0; i < querySets.size(); ++i) {
        totalQueries = querySets[i].size();
#if defined(COLLECT_STATISTICS)
        knnStats.clear();
#endif
        sw.reset();
        sw.start();
        for (std::size_t j = 0; j < querySets[i].size(); ++j) {
            ilbr.getKNNByNVDQuadtreePHL(alt,objectIndexes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertex,kNNs,kNNDistances,phl,graph);
            kNNs.clear();
            kNNDistances.clear();
#if defined(COLLECT_STATISTICS)
            knnStats.mergeStatistics(ilbr.stats);
#endif
        }
        sw.stop();
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;  
#if defined(COLLECT_STATISTICS)
        knnStats.normalizeStatistics(totalQueries);
#endif
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::NVD_ILBR_PHL_KNN_QUERY,
                            querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].var,querySets[i][0].objSetSize,kNNs,kNNDistances);
        results.setAdditionalFields(specialFields);
#if defined(COLLECT_STATISTICS)
        knnStats.populateTupleFields(results,0);
#endif
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "NVD-PHL kNN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKExperimentsCommand::runLLBPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,std::vector<VertexID>>& objectIndexes, 
                                      std::string statsOutputFile, std::string desc, std::vector<std::string> specialFields) 
{
    ILBR ilbr;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries;
    std::vector<VertexID> kNNs;
    std::vector<EdgeWeight> kNNDistances;

    Statistics knnStats;
    
    // Warm-up: the first WARM_UP_QUERY_PERCENTAGE percent of the first list (at least one query) runs untimed, so that the cold start of a fresh
    // process or of freshly loaded indexes stays out of the timings (skipped in rn_obj_search_stats, whose statistics do not depend on it)
#if !defined(COLLECT_STATISTICS)
    if (querySets.size() > 0) {
        std::size_t numWarmUpQueries = utility::getNumWarmUpQueries(querySets[0].size());
        for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
            ilbr.getKNNsByPHLNaive(alt,objectIndexes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertex,kNNs,kNNDistances,phl);
            kNNs.clear();
            kNNDistances.clear();
        }
    }
#endif

    for (std::size_t i = 0; i < querySets.size(); ++i) {
        totalQueries = querySets[i].size();
#if defined(COLLECT_STATISTICS)
        knnStats.clear();
#endif
        sw.reset();
        sw.start();
        for (std::size_t j = 0; j < querySets[i].size(); ++j) {
            ilbr.getKNNsByPHLNaive(alt,objectIndexes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertex,kNNs,kNNDistances,phl);
            kNNs.clear();
            kNNDistances.clear();
#if defined(COLLECT_STATISTICS)
            knnStats.mergeStatistics(ilbr.stats);
#endif
        }
        sw.stop();
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;  
#if defined(COLLECT_STATISTICS)
        knnStats.normalizeStatistics(totalQueries);
#endif
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::LLB_ILBR_PHL_KNN_QUERY,
                            querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].var,querySets[i][0].objSetSize,kNNs,kNNDistances);
        results.setAdditionalFields(specialFields);
#if defined(COLLECT_STATISTICS)
        knnStats.populateTupleFields(results,0);
#endif
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "LLB-PHL kNN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKExperimentsCommand::runCOLTPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, SULTree& sltree, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,COLTree>& objectIndexes, 
                                       std::string statsOutputFile, std::string desc, std::vector<std::string> specialFields, std::vector<std::pair<std::string,std::string>> supplementaryFields)
{
    ILBR ilbr;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries;
    std::vector<VertexID> kNNs;
    std::vector<EdgeWeight> kNNDistances;

    Statistics knnStats;
    
    // Warm-up: the first WARM_UP_QUERY_PERCENTAGE percent of the first list (at least one query) runs untimed, so that the cold start of a fresh
    // process or of freshly loaded indexes stays out of the timings (skipped in rn_obj_search_stats, whose statistics do not depend on it)
#if !defined(COLLECT_STATISTICS)
    if (querySets.size() > 0) {
        std::size_t numWarmUpQueries = utility::getNumWarmUpQueries(querySets[0].size());
        for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
            ilbr.getKNNByCOLTPHL(alt,objectIndexes[querySets[0][j].objSetFilenamePrefix],sltree,querySets[0][j].k,querySets[0][j].queryVertex,kNNs,kNNDistances,phl);
            kNNs.clear();
            kNNDistances.clear();
        }
    }
#endif

    for (std::size_t i = 0; i < querySets.size(); ++i) {
        totalQueries = querySets[i].size();
#if defined(COLLECT_STATISTICS)
        knnStats.clear();
#endif
        sw.reset();
        sw.start();
        for (std::size_t j = 0; j < querySets[i].size(); ++j) {
            ilbr.getKNNByCOLTPHL(alt,objectIndexes[querySets[i][j].objSetFilenamePrefix],sltree,querySets[i][j].k,querySets[i][j].queryVertex,kNNs,kNNDistances,phl);
            kNNs.clear();
            kNNDistances.clear();
#if defined(COLLECT_STATISTICS)
            knnStats.mergeStatistics(ilbr.stats);
#endif
        }
        sw.stop();
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;  
#if defined(COLLECT_STATISTICS)
        knnStats.normalizeStatistics(totalQueries);
#endif
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::COLTREE_ILBR_PHL_KNN_QUERY,
                            querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].var,querySets[i][0].objSetSize,kNNs,kNNDistances);
        results.setAdditionalFields(specialFields);
        for (std::size_t f = 0; f < supplementaryFields.size(); ++f) {
            results.addSupplementaryFields(supplementaryFields[f].first,supplementaryFields[f].second);
        }
#if defined(COLLECT_STATISTICS)
        knnStats.populateTupleFields(results,0);
#endif
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "COLT-PHL kNN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKExperimentsCommand::buildSULTreeParameterIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, 
                                                        std::string strFanouts, std::string strNumLandmarks, std::string strLandmarkTypes, 
                                                        int minLambda)
{
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    // Parameter testing directory for the SUL-Tree and COL-Tree variants (PHL, ALT and the object sets are read from filePathPrefix)
    std::string parFilePathPrefix = parameterMap["par_path"];
    if (parFilePathPrefix == "") {
        std::cerr << "Parameter testing path (par_path) was not provided!\n";
        exit(1);
    }

    // Find all additional fields we need to add to index stats tuples
    std::vector<std::string> specialFields;
    int field = 0;
    while (true) {
        std::string key = "special_field_" + std::to_string(field);
        if (parameterMap.find(key) != parameterMap.end()) {
            specialFields.push_back(parameterMap[key]);
            ++field;
        } else {
            break;
        }
    }
    
    // Process input variables
    std::vector<std::string> strFanoutsVec = utility::splitByDelim(strFanouts,',');
    std::vector<int> fanouts;
    for(std::size_t i = 0; i < strFanoutsVec.size(); ++i) {
        int fanout = std::stoi(strFanoutsVec[i]);
        if (fanout > 0) {
            fanouts.push_back(fanout);
        } else {
            std::cerr << "Invalid fanout " << strFanoutsVec[i] << " in list provided!\n";
            exit(1);
        }
    }

    std::vector<std::string> strNumLandmarksVec = utility::splitByDelim(strNumLandmarks,',');
    std::vector<int> numLandmarks;
    for(std::size_t i = 0; i < strNumLandmarksVec.size(); ++i) {
        int landmarks = std::stoi(strNumLandmarksVec[i]);
        if (landmarks > 0) {
            numLandmarks.push_back(landmarks);
        } else {
            std::cerr << "Invalid number of landmarks " << strNumLandmarksVec[i] << " in list provided!\n";
            exit(1);
        }
    }

    std::vector<std::string> landmarkTypes = utility::splitByDelim(strLandmarkTypes,',');

    if (fanouts.size() == 0 || numLandmarks.size() == 0 || landmarkTypes.size() == 0) {
        std::cerr << "Not enough fanouots, number of landmarks or landmark types provided!\n";
        exit(1);
    }

    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);

#if defined(PARALLEL_BUILDS)
    // Threads of the SUL-Tree builds (rn_obj_search_par only, see PARALLEL_BUILDS_PRAGMA in common.h; rn_obj_search builds serially),
    // split below Phase 1
    int numThreads = omp_get_max_threads(); // The available cores, or OMP_NUM_THREADS if set
#endif

    // Phase 1: all SUL-Trees with the same fanout (and maximum leaf size) are built from one partition, loaded if it was built before.
    // The partitions are built one after another, also in rn_obj_search_par: METIS draws from libc's process-wide rand(), so partitions
    // built in parallel differ from run to run, while these are the partitions lmk_exp -e partitions builds for the same fanout and
    // maximum leaf size (SUL-Tree landmarks do not depend on the order either, see SULTree::getLandmarkEngine)
    std::vector<HierarchicalPartition> partitions(fanouts.size());
    for (std::size_t i = 0; i < fanouts.size(); ++i) {
        std::vector<std::string> partitionKeys, partitionValues;
        partitionKeys.push_back("fanout");
        partitionValues.push_back(std::to_string(fanouts[i]));
        partitionKeys.push_back("maxleafsize");
        partitionValues.push_back(std::to_string(minLambda));
        std::string partitionFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PARTITION_CMD,partitionKeys,partitionValues);
        if (std::ifstream(partitionFile).good()) {
            std::cout << "Loading partition " << partitionFile << std::endl;
            partitions[i] = serialization::getIndexFromBinaryFile<HierarchicalPartition>(partitionFile);
        } else {
            std::cout << "Building partition (f=" + std::to_string(fanouts[i]) + ",lambda=" + std::to_string(minLambda) + ")..." << std::endl;
            std::string partitionStatsFile = filePathPrefix + "/stats/" + graph.getNetworkName() + "_fanout=" + std::to_string(fanouts[i]) + "_lambda="
                + std::to_string(minLambda) + "_parameter_partition_stats.txt";
            this->buildPartition(graph,fanouts[i],minLambda,partitionFile,partitionStatsFile,specialFields,partitions[i]);
        }
    }

#if defined(PARALLEL_BUILDS)
    // Each SUL-Tree processes its tree nodes with 4 threads (the nested region in SULTree::buildSULTree), up to which the node loop scales
    // well, and as many SUL-Trees are built at once as the other threads allow (2 x 4 on 8 cores). Each SUL-Tree being built holds its
    // distance lists in memory, so if the largest would not fit that many times, fewer are built at once and their threads go to the node
    // loops instead (on NW, 1 SUL-Tree with 8 threads was 9% slower than 2 with 4, while 4 SUL-Trees with 1 thread each took twice as long).
    // A SUL-Tree's memory is estimated from its partition (the distance lists, plus the vertex orderings and border lists it copies) with
    // 10% added for the searches (on USA it used 8% more than its index size), against 90% of the memory available after Phase 1
    int maxNumLandmarks = *std::max_element(numLandmarks.begin(),numLandmarks.end());
    long maxBuildBytes = 0;
    for (std::size_t i = 0; i < partitions.size(); ++i) {
        long listEntries = 0, borders = 0;
        for (std::size_t t = 0; t < partitions[i].nodes.size(); ++t) {
            listEntries += partitions[i].nodes[t].numVertices;
            borders += partitions[i].nodes[t].borderVertexIDs.size();
        }
        long copiedIDs = 2*static_cast<long>(partitions[i].reorderedToOriginal.size()) + borders;
        long buildBytes = (listEntries*maxNumLandmarks*static_cast<long>(sizeof(EdgeWeight)) + copiedIDs*static_cast<long>(sizeof(VertexID)))*11/10;
        maxBuildBytes = std::max(maxBuildBytes,buildBytes);
    }
    long availableBytes = 0;
    std::ifstream meminfo("/proc/meminfo");
    std::string meminfoLine;
    while (std::getline(meminfo,meminfoLine)) {
        if (meminfoLine.compare(0,13,"MemAvailable:") == 0) {
            availableBytes = std::stol(meminfoLine.substr(13))*1024; // The value is in kB
        }
    }
    int numBuilds = static_cast<int>(fanouts.size()*numLandmarks.size()*landmarkTypes.size());
    int memoryLimit = numThreads;
    if (availableBytes == 0) {
        std::cerr << "WARNING: MemAvailable not found in /proc/meminfo, the number of SUL-Trees built at once is not limited by memory" << std::endl;
    } else if (maxBuildBytes > 0) {
        memoryLimit = static_cast<int>(std::min<long>(numThreads,availableBytes*9/10/maxBuildBytes));
        if (memoryLimit == 0) {
            std::cerr << "WARNING: The largest SUL-Tree may not fit in the available memory" << std::endl;
        }
    }
    int outerThreads = std::max(1,std::min(numBuilds,std::min(numThreads/4,memoryLimit)));
    int innerThreads = std::max(1,numThreads/outerThreads);
    omp_set_max_active_levels(2); // Nested regions: the SUL-Tree loop below and the node loop in SULTree::buildSULTree
    std::ostringstream memoryInfo; // Formatted separately so that std::cout keeps its format for the times printed later
    memoryInfo << std::fixed << std::setprecision(1) << "largest SUL-Tree about " << maxBuildBytes/1e9 << " GB, " << availableBytes/1e9 << " GB of memory available";
    std::cout << "Building " << outerThreads << " SUL-Tree(s) at a time with " << innerThreads << " thread(s) each (" << memoryInfo.str() << ")" << std::endl;
#endif

    // Phase 2: every SUL-Tree in one loop over all fanouts, so that a parallel region is not split (and its threads left idle) by fanout.
    // Note: schedule(dynamic) because the build time grows with the number of landmarks, which a static schedule gives to the same threads
    PARALLEL_BUILDS_PRAGMA(omp parallel for collapse(3) schedule(dynamic) num_threads(outerThreads))
    for (std::size_t i = 0; i < fanouts.size(); ++i) {
        for (std::size_t j = 0; j < numLandmarks.size(); ++j) {
            for (std::size_t k = 0; k < landmarkTypes.size(); ++k) {
                std::cout << "Building SULTree index (f=" + std::to_string(fanouts[i]) + ",lambda=" +std::to_string(minLambda) 
                    + ",m=" + std::to_string(numLandmarks[j])+ ",t=" + landmarkTypes[k] + ")..." << std::endl;
                std::vector<std::string> parameterKeys, parameterValues;
                parameterKeys.push_back("fanout");
                parameterValues.push_back(std::to_string(fanouts[i]));
                parameterKeys.push_back("maxleafsize");
                parameterValues.push_back(std::to_string(minLambda));
                parameterKeys.push_back("numlmks");
                parameterValues.push_back(std::to_string(numLandmarks[j]));
                parameterKeys.push_back("lmktype");
                parameterValues.push_back(landmarkTypes[k]);
                std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,parameterKeys,parameterValues);
                std::string statsOutputFile = filePathPrefix + "/stats/" + graph.getNetworkName() + "_fanout=" + std::to_string(fanouts[i])  + "_lambda="
                   + std::to_string(minLambda) + "_numlmks=" + std::to_string(numLandmarks[j]) + "_lmktype=" + landmarkTypes[k] + "_parameter_index_stats.txt";
#if defined(PARALLEL_BUILDS)
                omp_set_num_threads(innerThreads); // Threads of this SUL-Tree's node loop (a nested region in SULTree::buildSULTree)
#endif
                this->buildSULTree(graph,partitions[i],numLandmarks[j],landmarkTypes[k],idxOutputFile,statsOutputFile,specialFields);
            }
        }
    }

}

void LMKExperimentsCommand::buildCOLTParameterIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::size_t numSets, 
                                                      std::string strObjectDensities, std::string objTypes, std::string objVariable, std::string strFanouts, 
                                                      std::string strNumLandmarks, std::string strLandmarkTypes, std::string strLambdas, int minLambda)
{
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    // Parameter testing directory for the SUL-Tree and COL-Tree variants (PHL, ALT and the object sets are read from filePathPrefix)
    std::string parFilePathPrefix = parameterMap["par_path"];
    if (parFilePathPrefix == "") {
        std::cerr << "Parameter testing path (par_path) was not provided!\n";
        exit(1);
    }
    int setSize, setVariable;
    std::string setType;
    double setDensity;
    // Find all additional fields we need to add to index stats tuples
    std::vector<std::string> specialFields;
    int field = 0;
    while (true) {
        std::string key = "special_field_" + std::to_string(field);
        if (parameterMap.find(key) != parameterMap.end()) {
            specialFields.push_back(parameterMap[key]);
            ++field;
        } else {
            break;
        }
    }

    // Process object set input variables
    std::vector<std::string> strObjectDensitiesVec = utility::splitByDelim(strObjectDensities,',');
    std::vector<double> densities;
    for(std::size_t i = 0; i < strObjectDensitiesVec.size(); ++i) {
        double density = std::stod(strObjectDensitiesVec[i]);
        if (density > 0) {
            densities.push_back(density);
        } else {
            std::cerr << "Invalid density " << strObjectDensitiesVec[i] << " in list provided!\n";
            exit(1);
        }
    }
    std::vector<std::string> types = utility::splitByDelim(objTypes,',');
    if (densities.size() == 0 || types.size() == 0) {
        std::cerr << "Not enough densities or types provided!\n";
        exit(1);
    }
    double defaultDensity = std::stod(parameterMap["par_def_density"]);
    std::string defaultObjType = parameterMap["par_def_type"];
    int defaultObjVariable = std::stoi(parameterMap["par_def_variable"]);
    if (defaultDensity == 0 || defaultObjType == "") {
        std::cerr << "Invalid detault deensity or object type!\n";
        exit(1);
    }

     // Process parametere testing input variables
    std::vector<std::string> strFanoutsVec = utility::splitByDelim(strFanouts,',');
    std::vector<int> fanouts;
    for(std::size_t i = 0; i < strFanoutsVec.size(); ++i) {
        int fanout = std::stoi(strFanoutsVec[i]);
        if (fanout > 0) {
            fanouts.push_back(fanout);
        } else {
            std::cerr << "Invalid fanout " << strFanoutsVec[i] << " in list provided!\n";
            exit(1);
        }
    }
    std::vector<std::string> strNumLandmarksVec = utility::splitByDelim(strNumLandmarks,',');
    std::vector<int> numLandmarks;
    for(std::size_t i = 0; i < strNumLandmarksVec.size(); ++i) {
        int landmarks = std::stoi(strNumLandmarksVec[i]);
        if (landmarks > 0) {
            numLandmarks.push_back(landmarks);
        } else {
            std::cerr << "Invalid number of landmarks " << strNumLandmarksVec[i] << " in list provided!\n";
            exit(1);
        }
    }
    std::vector<std::string> strLambdasVec = utility::splitByDelim(strLambdas,',');
    std::vector<int> objectListSizes; // A.k.a. lambda
    for(std::size_t i = 0; i < strLambdasVec.size(); ++i) {
        int lambda = std::stoi(strLambdasVec[i]);
        if (lambda > 0) {
            objectListSizes.push_back(lambda);
        } else {
            std::cerr << "Invalid object list size (lambda) " << strLambdasVec[i] << " in list provided!\n";
            exit(1);
        }
    }
    std::vector<std::string> landmarkTypes = utility::splitByDelim(strLandmarkTypes,',');
    if (fanouts.size() == 0 || numLandmarks.size() == 0 || landmarkTypes.size() == 0) {
        std::cerr << "Not enough fanouots, number of landmarks or landmark types provided!\n";
        exit(1);
    }
    int defFanout = std::stoi(parameterMap["par_def_fanout"]);
    int defNumLandmarks = std::stoi(parameterMap["par_def_numlmks"]);
    int defLambda = std::stoi(parameterMap["par_def_lambda"]);
    std::string defLandmarkType = parameterMap["par_def_lmktype"];
    if (defFanout == 0 || defNumLandmarks == 0 || defLambda == 0 || defLandmarkType == "") {
        std::cerr << "Invalid detault fanout, num landmarks, lambda, or landmark type!\n";
        exit(1);
    }

    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);

    // Note: Lambda variable should be inner most parameter becausee you only need one SULTree

    // Load each required default object set once for reuse across parameter sweeps
    const std::size_t defaultTotalSets = utility::getNumObjectSets(defaultDensity,numSets);
    std::vector<std::vector<VertexID>> defaultObjectSets;
    for (std::size_t i = 0; i < defaultTotalSets; ++i) {
        std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVariable,i);
        std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
        defaultObjectSets.push_back(objectVertices);
    }
    
    StopWatch sw;
    double totalTime = 0, totalMemory = 0;

    // Build object indexes for fanout vs lambda experiment
    if (parameterMap.find("fanout_vs_lambda_stats_file_path") == parameterMap.end() || parameterMap["fanout_vs_lambda_stats_file_path"] == "") {
        std::cerr << "Fanout vs lamda experiment stats file path not provided in parameters";
        std::exit(1);
    }
    std::string statsOutputFile = parameterMap["fanout_vs_lambda_stats_file_path"];
    for (std::size_t i = 0; i < fanouts.size(); ++i) {
        // Load SULTree for given fanout, min lambda, and defaults for other variables
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(fanouts[i]),std::to_string(minLambda),std::to_string(defNumLandmarks),defLandmarkType};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t j = 0; j < objectListSizes.size(); ++j) {
            totalTime = 0;
            totalMemory = 0;
            for (std::size_t k = 0; k < defaultTotalSets; ++k) {
                sw.start();
                COLTree colt(defaultObjType,defaultDensity,defaultObjVariable,defaultObjectSets[k].size());
                colt.buildCOLTTopDown(sltree,defaultObjectSets[k],objectListSizes[j]);
                sw.stop();
                totalTime += sw.getTimeUs();
                totalMemory += colt.computeIndexSize()*1024; // in KB
                sw.reset();

                std::string objIdxOutputFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVariable,k);
                std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
                coltKeys.push_back("lambda");
                coltValues.push_back(std::to_string(objectListSizes[j]));
                std::string objIdxOutputFile = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxOutputFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                serialization::outputIndexToBinaryFile<COLTree>(colt,objIdxOutputFile);          
            }
            // Average stats over all object sets and then write to stats file
            double processingTimeMs = totalTime/defaultTotalSets;
            double memoryUsage = totalMemory/defaultTotalSets;
            
            ObjectIndexTuple stats(sltree.getNetworkName(),sltree.getNumVertices(),sltree.getNumEdges(),constants::OBJ_IDX_COLTREE,
                                   processingTimeMs,memoryUsage,defaultObjType,defaultDensity,defaultObjVariable,defaultObjectSets[0].size());
            stats.setAdditionalFields(specialFields);
            stats.addSupplementaryFields("fanout",std::to_string(fanouts[i]));
            stats.addSupplementaryFields("numlmks",std::to_string(defNumLandmarks));
            stats.addSupplementaryFields("lmktype",defLandmarkType);
            stats.addSupplementaryFields("lambda",std::to_string(objectListSizes[j]));
            this->outputCommandStats(statsOutputFile,stats.getTupleString());
        }
    }
    std::cout << "COLTree indexes successfully created for " << graph.getNetworkName() << " fanout vs. lambda parameter experiment" << std::endl;

    // Build object indexes for num landmarks vs lambda experiment
    if (parameterMap.find("numlmks_vs_lambda_stats_file_path") == parameterMap.end() || parameterMap["numlmks_vs_lambda_stats_file_path"] == "") {
        std::cerr << "Num landmarks vs lamda experiment stats file path not provided in parameters";
        std::exit(1);
    }
    statsOutputFile = parameterMap["numlmks_vs_lambda_stats_file_path"];
    for (std::size_t i = 0; i < numLandmarks.size(); ++i) {
        // Load SULTree for given landmark number, min lambda, and defaults for other variables
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(defFanout),std::to_string(minLambda),std::to_string(numLandmarks[i]),defLandmarkType};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t j = 0; j < objectListSizes.size(); ++j) {
            totalTime = 0;
            totalMemory = 0;
            for (std::size_t k = 0; k < defaultTotalSets; ++k) {
                sw.start();
                COLTree colt(defaultObjType,defaultDensity,defaultObjVariable,defaultObjectSets[k].size());
                colt.buildCOLTTopDown(sltree,defaultObjectSets[k],objectListSizes[j]);
                sw.stop();
                totalTime += sw.getTimeUs();
                totalMemory += colt.computeIndexSize()*1024; // in KB
                sw.reset();

                std::string objIdxOutputFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVariable,k);
                std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
                coltKeys.push_back("lambda");
                coltValues.push_back(std::to_string(objectListSizes[j]));
                std::string objIdxOutputFile = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxOutputFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                serialization::outputIndexToBinaryFile<COLTree>(colt,objIdxOutputFile);
            }
            // Average stats over all object sets and then write to stats file
            double processingTimeMs = totalTime/defaultTotalSets;
            double memoryUsage = totalMemory/defaultTotalSets;
            
            ObjectIndexTuple stats(sltree.getNetworkName(),sltree.getNumVertices(),sltree.getNumEdges(),constants::OBJ_IDX_COLTREE,
                                   processingTimeMs,memoryUsage,defaultObjType,defaultDensity,defaultObjVariable,defaultObjectSets[0].size());
            stats.setAdditionalFields(specialFields);
            stats.addSupplementaryFields("fanout",std::to_string(defFanout));
            stats.addSupplementaryFields("numlmks",std::to_string(numLandmarks[i]));
            stats.addSupplementaryFields("lmktype",defLandmarkType);
            stats.addSupplementaryFields("lambda",std::to_string(objectListSizes[j]));
            this->outputCommandStats(statsOutputFile,stats.getTupleString());
        }
    }
    std::cout << "COLTree indexes successfully created for " << graph.getNetworkName() << " number of landmarks vs. lambda parameter experiment" << std::endl;

    // Build object indexes for num landmarks vs landmark type experiment
    if (parameterMap.find("numlmks_vs_lmktype_stats_file_path") == parameterMap.end() || parameterMap["numlmks_vs_lmktype_stats_file_path"] == "") {
        std::cerr << "Num landmarks vs landmark type experiment stats file path not provided in parameters";
        std::exit(1);
    }
    statsOutputFile = parameterMap["numlmks_vs_lmktype_stats_file_path"];
    for (std::size_t i = 0; i < numLandmarks.size(); ++i) {
        for (std::size_t j = 0; j < landmarkTypes.size(); ++j) {
            // Load SULTree for given landmark number, type, and defaults for other variables
            const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
            const std::vector<std::string> sltreeValues = {std::to_string(defFanout),std::to_string(minLambda),std::to_string(numLandmarks[i]),landmarkTypes[j]};
            std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
            SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

            totalTime = 0;
            totalMemory = 0;
            for (std::size_t k = 0; k < defaultTotalSets; ++k) {
                sw.start();
                COLTree colt(defaultObjType,defaultDensity,defaultObjVariable,defaultObjectSets[k].size());
                colt.buildCOLTTopDown(sltree,defaultObjectSets[k],defLambda);
                sw.stop();
                totalTime += sw.getTimeUs();
                totalMemory += colt.computeIndexSize()*1024; // in KB
                sw.reset();

                std::string objIdxOutputFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVariable,k);
                std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
                coltKeys.push_back("lambda");
                coltValues.push_back(std::to_string(defLambda));
                std::string objIdxOutputFile = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxOutputFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                serialization::outputIndexToBinaryFile<COLTree>(colt,objIdxOutputFile);
            }
            // Average stats over all object sets and then write to stats file
            double processingTimeMs = totalTime/defaultTotalSets;
            double memoryUsage = totalMemory/defaultTotalSets;
            
            ObjectIndexTuple stats(sltree.getNetworkName(),sltree.getNumVertices(),sltree.getNumEdges(),constants::OBJ_IDX_COLTREE,
                                   processingTimeMs,memoryUsage,defaultObjType,defaultDensity,defaultObjVariable,defaultObjectSets[0].size());
            stats.setAdditionalFields(specialFields);
            stats.addSupplementaryFields("fanout",std::to_string(defFanout));
            stats.addSupplementaryFields("numlmks",std::to_string(numLandmarks[i]));
            stats.addSupplementaryFields("lmktype",landmarkTypes[j]);
            stats.addSupplementaryFields("lambda",std::to_string(defLambda));
            this->outputCommandStats(statsOutputFile,stats.getTupleString());
        }
    }
    std::cout << "COLTree indexes successfully created for " << graph.getNetworkName() << " number of landmarks vs. landmark type parameter experiment" << std::endl;

    // Build object indexes for density vs num landmarks experiment
    if (parameterMap.find("density_vs_numlmks_stats_file_path") == parameterMap.end() || parameterMap["density_vs_numlmks_stats_file_path"] == "") {
        std::cerr << "Density vs num landmarks experiment stats file path not provided in parameters";
        std::exit(1);
    }
    statsOutputFile = parameterMap["density_vs_numlmks_stats_file_path"];
    for (std::size_t i = 0; i < numLandmarks.size(); ++i) {
        // Load SULTree for given fanout, min lambda, and defaults for other variables
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(defFanout),std::to_string(minLambda),std::to_string(numLandmarks[i]),defLandmarkType};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t j = 0; j < densities.size(); ++j) {
            const std::size_t totalSets = utility::getNumObjectSets(densities[j],numSets);
            totalTime = 0;
            totalMemory = 0;
            int objectSizeSet = 0;
            for (std::size_t k = 0; k < totalSets; ++k) {
                // Need load object set files (since previous load for oone density)
                std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                objectSizeSet = objectVertices.size();
                
                sw.start();
                COLTree colt(defaultObjType,densities[j],defaultObjVariable,objectVertices.size());
                colt.buildCOLTTopDown(sltree,objectVertices,defLambda);
                sw.stop();
                totalTime += sw.getTimeUs();
                totalMemory += colt.computeIndexSize()*1024; // in KB
                sw.reset();

                std::string objIdxOutputFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
                coltKeys.push_back("lambda");
                coltValues.push_back(std::to_string(defLambda));
                std::string objIdxOutputFile = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxOutputFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                serialization::outputIndexToBinaryFile<COLTree>(colt,objIdxOutputFile);
            }
            // Average stats over all object sets and then write to stats file
            double processingTimeMs = totalTime/totalSets;
            double memoryUsage = totalMemory/totalSets;
            
            ObjectIndexTuple stats(sltree.getNetworkName(),sltree.getNumVertices(),sltree.getNumEdges(),constants::OBJ_IDX_COLTREE,
                                   processingTimeMs,memoryUsage,defaultObjType,densities[j],defaultObjVariable,objectSizeSet);
            stats.setAdditionalFields(specialFields);
            stats.addSupplementaryFields("fanout",std::to_string(defFanout));
            stats.addSupplementaryFields("numlmks",std::to_string(numLandmarks[i]));
            stats.addSupplementaryFields("lmktype",defLandmarkType);
            stats.addSupplementaryFields("lambda",std::to_string(defLambda));
            stats.addSupplementaryFields("density",std::to_string(densities[j]));
            this->outputCommandStats(statsOutputFile,stats.getTupleString());
        }
    }
    std::cout << "COLTree indexes successfully created for " << graph.getNetworkName() << " number of landmarks vs. density parameter experiment" << std::endl;

    // Build object indexes for density vs landmark types experiment
    if (parameterMap.find("density_vs_lmktype_stats_file_path") == parameterMap.end() || parameterMap["density_vs_lmktype_stats_file_path"] == "") {
        std::cerr << "Density vs landmark type experiment stats file path not provided in parameters";
        std::exit(1);
    }
    statsOutputFile = parameterMap["density_vs_lmktype_stats_file_path"];
    for (std::size_t i = 0; i < landmarkTypes.size(); ++i) {
        // Load SULTree for given fanout, min lambda, and defaults for other variables
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(defFanout),std::to_string(minLambda),std::to_string(defNumLandmarks),landmarkTypes[i]};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t j = 0; j < densities.size(); ++j) {
            const std::size_t totalSets = utility::getNumObjectSets(densities[j],numSets);
            totalTime = 0;
            totalMemory = 0;
            int objectSizeSet = 0;
            for (std::size_t k = 0; k < totalSets; ++k) {
                // Need load object set files (since previous load for oone density)
                std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                objectSizeSet = objectVertices.size();
                
                sw.start();
                COLTree colt(defaultObjType,densities[j],defaultObjVariable,objectVertices.size());
                colt.buildCOLTTopDown(sltree,objectVertices,defLambda);
                sw.stop();
                totalTime += sw.getTimeUs();
                totalMemory += colt.computeIndexSize()*1024; // in KB
                sw.reset();

                std::string objIdxOutputFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
                coltKeys.push_back("lambda");
                coltValues.push_back(std::to_string(defLambda));
                std::string objIdxOutputFile = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxOutputFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                serialization::outputIndexToBinaryFile<COLTree>(colt,objIdxOutputFile);
            }
            // Average stats over all object sets and then write to stats file
            double processingTimeMs = totalTime/totalSets;
            double memoryUsage = totalMemory/totalSets;
            
            ObjectIndexTuple stats(sltree.getNetworkName(),sltree.getNumVertices(),sltree.getNumEdges(),constants::OBJ_IDX_COLTREE,
                                   processingTimeMs,memoryUsage,defaultObjType,densities[j],defaultObjVariable,objectSizeSet);
            stats.setAdditionalFields(specialFields);
            stats.addSupplementaryFields("fanout",std::to_string(defFanout));
            stats.addSupplementaryFields("numlmks",std::to_string(defNumLandmarks));
            stats.addSupplementaryFields("lmktype",landmarkTypes[i]);
            stats.addSupplementaryFields("lambda",std::to_string(defLambda));
            stats.addSupplementaryFields("density",std::to_string(densities[j]));
            this->outputCommandStats(statsOutputFile,stats.getTupleString());
        }
    }
    std::cout << "COLTree indexes successfully created for " << graph.getNetworkName() << " landmark type vs. density parameter experiment" << std::endl;

    // Build object indexes for density vs lambda experiment
    {
        // Use scope to release the default SULTree before the density vs fanout experiment loads its own
        if (parameterMap.find("density_vs_lambda_stats_file_path") == parameterMap.end() || parameterMap["density_vs_lambda_stats_file_path"] == "") {
            std::cerr << "Density vs lambda type experiment stats file path not provided in parameters";
            std::exit(1);
        }
        statsOutputFile = parameterMap["density_vs_lambda_stats_file_path"];

        // Load relevant SULTree (it does not vary for lambda)
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(defFanout),std::to_string(minLambda),std::to_string(defNumLandmarks),defLandmarkType};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t i = 0; i < objectListSizes.size(); ++i) {
            for (std::size_t j = 0; j < densities.size(); ++j) {
                const std::size_t totalSets = utility::getNumObjectSets(densities[j],numSets);
                totalTime = 0;
                totalMemory = 0;
                int objectSizeSet = 0;
                for (std::size_t k = 0; k < totalSets; ++k) {
                    // Need load object set files (since previous load for oone density)
                    std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                    std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                    objectSizeSet = objectVertices.size();

                    sw.start();
                    COLTree colt(defaultObjType,densities[j],defaultObjVariable,objectVertices.size());
                    colt.buildCOLTTopDown(sltree,objectVertices,objectListSizes[i]);
                    sw.stop();
                    totalTime += sw.getTimeUs();
                    totalMemory += colt.computeIndexSize()*1024; // in KB
                    sw.reset();

                    std::string objIdxOutputFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                    std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
                    coltKeys.push_back("lambda");
                    coltValues.push_back(std::to_string(objectListSizes[i]));
                    std::string objIdxOutputFile = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxOutputFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                    serialization::outputIndexToBinaryFile<COLTree>(colt,objIdxOutputFile);
                }
                // Average stats over all object sets and then write to stats file
                double processingTimeMs = totalTime/totalSets;
                double memoryUsage = totalMemory/totalSets;

                ObjectIndexTuple stats(sltree.getNetworkName(),sltree.getNumVertices(),sltree.getNumEdges(),constants::OBJ_IDX_COLTREE,
                                       processingTimeMs,memoryUsage,defaultObjType,densities[j],defaultObjVariable,objectSizeSet);
                stats.setAdditionalFields(specialFields);
                stats.addSupplementaryFields("fanout",std::to_string(defFanout));
                stats.addSupplementaryFields("numlmks",std::to_string(defNumLandmarks));
                stats.addSupplementaryFields("lmktype",defLandmarkType);
                stats.addSupplementaryFields("lambda",std::to_string(objectListSizes[i]));
                stats.addSupplementaryFields("density",std::to_string(densities[j]));
                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
        std::cout << "COLTree indexes successfully created for " << graph.getNetworkName() << " lambda vs. density parameter experiment" << std::endl;
    }

    // Build object indexes for density vs fanout experiment
    if (parameterMap.find("density_vs_fanout_stats_file_path") == parameterMap.end() || parameterMap["density_vs_fanout_stats_file_path"] == "") {
        std::cerr << "Density vs fanout type experiment stats file path not provided in parameters";
        std::exit(1);
    }
    statsOutputFile = parameterMap["density_vs_fanout_stats_file_path"];
    for (std::size_t i = 0; i < fanouts.size(); ++i) {
        // Load relevant SULTree
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(fanouts[i]),std::to_string(minLambda),std::to_string(defNumLandmarks),defLandmarkType};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t j = 0; j < densities.size(); ++j) {
            const std::size_t totalSets = utility::getNumObjectSets(densities[j],numSets);
            totalTime = 0;
            totalMemory = 0;
            int objectSizeSet = 0;
            for (std::size_t k = 0; k < totalSets; ++k) {
                // Need load object set files (since previous load for oone density)
                std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                std::vector<VertexID> objectVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                objectSizeSet = objectVertices.size();

                sw.start();
                COLTree colt(defaultObjType,densities[j],defaultObjVariable,objectVertices.size());
                colt.buildCOLTTopDown(sltree,objectVertices,defLambda);
                sw.stop();
                totalTime += sw.getTimeUs();
                totalMemory += colt.computeIndexSize()*1024; // in KB
                sw.reset();

                std::string objIdxOutputFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
                coltKeys.push_back("lambda");
                coltValues.push_back(std::to_string(defLambda));
                std::string objIdxOutputFile = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxOutputFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                serialization::outputIndexToBinaryFile<COLTree>(colt,objIdxOutputFile);
            }
            // Average stats over all object sets and then write to stats file
            double processingTimeMs = totalTime/totalSets;
            double memoryUsage = totalMemory/totalSets;

            ObjectIndexTuple stats(sltree.getNetworkName(),sltree.getNumVertices(),sltree.getNumEdges(),constants::OBJ_IDX_COLTREE,
                                   processingTimeMs,memoryUsage,defaultObjType,densities[j],defaultObjVariable,objectSizeSet);
            stats.setAdditionalFields(specialFields);
            stats.addSupplementaryFields("fanout",std::to_string(fanouts[i]));
            stats.addSupplementaryFields("numlmks",std::to_string(defNumLandmarks));
            stats.addSupplementaryFields("lmktype",defLandmarkType);
            stats.addSupplementaryFields("lambda",std::to_string(defLambda));
            stats.addSupplementaryFields("density",std::to_string(densities[j]));
            this->outputCommandStats(statsOutputFile,stats.getTupleString());
        }
    }
    std::cout << "COLTree indexes successfully created for " << graph.getNetworkName() << " fanout vs. density parameter experiment" << std::endl;
}

void LMKExperimentsCommand::runCOLTParameterTestQueries(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string queryVertexFile, 
                                                        std::size_t numSets, std::string strObjectDensities, std::string objTypes, std::string objVariable, 
                                                        std::string strFanouts, std::string strNumLandmarks, std::string strLandmarkTypes, std::string strLambdas, 
                                                        int minLambda)
{
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    // Parameter testing directory for the SUL-Tree and COL-Tree variants (PHL, ALT and the object sets are read from filePathPrefix)
    std::string parFilePathPrefix = parameterMap["par_path"];
    if (parFilePathPrefix == "") {
        std::cerr << "Parameter testing path (par_path) was not provided!\n";
        exit(1);
    }

    // Find all additional fields we need to add to index stats tuples
    std::vector<std::string> specialFields;
    int field = 0;
    while (true) {
        std::string key = "special_field_" + std::to_string(field);
        if (parameterMap.find(key) != parameterMap.end()) {
            specialFields.push_back(parameterMap[key]);
            ++field;
        } else {
            break;
        }
    }

    // Process object set input variables
    std::vector<std::string> strObjectDensitiesVec = utility::splitByDelim(strObjectDensities,',');
    std::vector<double> densities;
    for(std::size_t i = 0; i < strObjectDensitiesVec.size(); ++i) {
        double density = std::stod(strObjectDensitiesVec[i]);
        if (density > 0) {
            densities.push_back(density);
        } else {
            std::cerr << "Invalid density " << strObjectDensitiesVec[i] << " in list provided!\n";
            exit(1);
        }
    }
    std::vector<std::string> types = utility::splitByDelim(objTypes,',');
    if (densities.size() == 0 || types.size() == 0) {
        std::cerr << "Not enough densities or types provided!\n";
        exit(1);
    }
    double defaultDensity = std::stod(parameterMap["par_def_density"]);
    std::string defaultObjType = parameterMap["par_def_type"];
    int defaultObjVariable = std::stoi(parameterMap["par_def_variable"]);
    int defaultK = std::stoi(parameterMap["par_def_knn_def_k"]);
    if (defaultDensity == 0 || defaultObjType == "" || defaultK == 0) {
        std::cerr << "Invalid detault density or object type!\n";
        exit(1);
    }

     // Process parametere testing input variables
    std::vector<std::string> strFanoutsVec = utility::splitByDelim(strFanouts,',');
    std::vector<int> fanouts;
    for(std::size_t i = 0; i < strFanoutsVec.size(); ++i) {
        int fanout = std::stoi(strFanoutsVec[i]);
        if (fanout > 0) {
            fanouts.push_back(fanout);
        } else {
            std::cerr << "Invalid fanout " << strFanoutsVec[i] << " in list provided!\n";
            exit(1);
        }
    }
    std::vector<std::string> strNumLandmarksVec = utility::splitByDelim(strNumLandmarks,',');
    std::vector<int> numLandmarks;
    for(std::size_t i = 0; i < strNumLandmarksVec.size(); ++i) {
        int landmarks = std::stoi(strNumLandmarksVec[i]);
        if (landmarks > 0) {
            numLandmarks.push_back(landmarks);
        } else {
            std::cerr << "Invalid number of landmarks " << strNumLandmarksVec[i] << " in list provided!\n";
            exit(1);
        }
    }
    std::vector<std::string> strLambdasVec = utility::splitByDelim(strLambdas,',');
    std::vector<int> objectListSizes; // A.k.a. lambda
    for(std::size_t i = 0; i < strLambdasVec.size(); ++i) {
        int lambda = std::stoi(strLambdasVec[i]);
        if (lambda > 0) {
            objectListSizes.push_back(lambda);
        } else {
            std::cerr << "Invalid object list size (lambda) " << strLambdasVec[i] << " in list provided!\n";
            exit(1);
        }
    }
    std::vector<std::string> landmarkTypes = utility::splitByDelim(strLandmarkTypes,',');
    if (fanouts.size() == 0 || numLandmarks.size() == 0 || landmarkTypes.size() == 0) {
        std::cerr << "Not enough fanouts, number of landmarks or landmark types provided!\n";
        exit(1);
    }
    int defFanout = std::stoi(parameterMap["par_def_fanout"]);
    int defNumLandmarks = std::stoi(parameterMap["par_def_numlmks"]);
    int defLambda = std::stoi(parameterMap["par_def_lambda"]);
    std::string defLandmarkType = parameterMap["par_def_lmktype"];
    if (defFanout == 0 || defNumLandmarks == 0 || defLambda == 0 || defLandmarkType == "") {
        std::cerr << "Invalid detault fanout, num landmarks, lambda, or landmark type!\n";
        exit(1);
    }

    std::vector<VertexID> queryVertices = utility::getPointSetFromFile(queryVertexFile);

    // Load supporting road network index files
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);

    std::string phlIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PHL_CMD,{},{});
    PrunedHighwayLabeling phl;
    phl.LoadLabel(phlIdxFilePath.c_str());

    const std::vector<std::string> altKeys = {"landmarks","type"};
    const std::vector<std::string> altValues = {parameterMap["alt_numlandmarks"],parameterMap["alt_landmarktype"]};
    std::string altIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,altKeys,altValues);
    ALT alt = serialization::getIndexFromBinaryFile<ALT>(altIdxFilePath);

    // Query lists as runQueries2 builds them: the default queries (default k and density) and one list per density for the density experiments
    std::vector<std::vector<KnnQuery>> defaultQuerySets(1), varyingDensityQuerySets;
    for (std::size_t j = 0; j < numSets; ++j) {
        std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVariable,j);
        for (std::size_t k = 0; k < queryVertices.size(); ++k) {
            defaultQuerySets[0].push_back(KnnQuery(queryVertices[k],defaultK,defaultDensity,defaultObjType,defaultObjVariable,j,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*defaultDensity)));
        }
    }
    for (std::size_t i = 0; i < densities.size(); ++i) {
        std::vector<KnnQuery> densityQuerySet;
        for (std::size_t j = 0; j < numSets; ++j) {
            std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,densities[i],defaultObjVariable,j);
            for (std::size_t k = 0; k < queryVertices.size(); ++k) {
                densityQuerySet.push_back(KnnQuery(queryVertices[k],defaultK,densities[i],defaultObjType,defaultObjVariable,j,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*densities[i])));
            }
        }
        varyingDensityQuerySets.push_back(densityQuerySet);
    }

    // For each tested parameter combination: load the SUL-Tree and the COL-Trees of all object sets built with it, time the query
    // lists on them with the same runner as the query experiments and record the combination in the rows as supplementary fields

    // Run queries for fanout vs lambda experiment
    if (parameterMap.find("fanout_vs_lambda_stats_file_path") == parameterMap.end() || parameterMap["fanout_vs_lambda_stats_file_path"] == "") {
        std::cerr << "Fanout vs lambda experiment stats file path not provided in parameters";
        std::exit(1);
    }
    std::string statsOutputFile = parameterMap["fanout_vs_lambda_stats_file_path"];
    for (std::size_t i = 0; i < fanouts.size(); ++i) {
        // Load SULTree for given fanout, min lambda, and defaults for other variables
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(fanouts[i]),std::to_string(minLambda),std::to_string(defNumLandmarks),defLandmarkType};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t j = 0; j < objectListSizes.size(); ++j) {
            std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
            coltKeys.push_back("lambda");
            coltValues.push_back(std::to_string(objectListSizes[j]));
            std::unordered_map<std::string,COLTree> fileNameToObjIdx;
            for (std::size_t k = 0; k < numSets; ++k) {
                std::string objIdxFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVariable,k);
                std::string objIdxFilePath = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                COLTree objIdx = serialization::getIndexFromBinaryFile<COLTree>(objIdxFilePath);
                fileNameToObjIdx[objIdxFilePrefix] = objIdx;
            }
            std::vector<std::pair<std::string,std::string>> testedParameters = {{"fanout",std::to_string(fanouts[i])},{"numlmks",std::to_string(defNumLandmarks)},{"lmktype",defLandmarkType},{"lambda",std::to_string(objectListSizes[j])}};
            this->runCOLTPHL(graph,phl,alt,sltree,defaultQuerySets,fileNameToObjIdx,statsOutputFile,"fanout vs lambda",specialFields,testedParameters);
        }
    }

    // Run queries for num landmarks vs lambda experiment
    if (parameterMap.find("numlmks_vs_lambda_stats_file_path") == parameterMap.end() || parameterMap["numlmks_vs_lambda_stats_file_path"] == "") {
        std::cerr << "Num landmarks vs lambda experiment stats file path not provided in parameters";
        std::exit(1);
    }
    statsOutputFile = parameterMap["numlmks_vs_lambda_stats_file_path"];
    for (std::size_t i = 0; i < numLandmarks.size(); ++i) {
        // Load SULTree for given landmark number, min lambda, and defaults for other variables
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(defFanout),std::to_string(minLambda),std::to_string(numLandmarks[i]),defLandmarkType};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t j = 0; j < objectListSizes.size(); ++j) {
            std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
            coltKeys.push_back("lambda");
            coltValues.push_back(std::to_string(objectListSizes[j]));
            std::unordered_map<std::string,COLTree> fileNameToObjIdx;
            for (std::size_t k = 0; k < numSets; ++k) {
                std::string objIdxFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVariable,k);
                std::string objIdxFilePath = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                COLTree objIdx = serialization::getIndexFromBinaryFile<COLTree>(objIdxFilePath);
                fileNameToObjIdx[objIdxFilePrefix] = objIdx;
            }
            std::vector<std::pair<std::string,std::string>> testedParameters = {{"fanout",std::to_string(defFanout)},{"numlmks",std::to_string(numLandmarks[i])},{"lmktype",defLandmarkType},{"lambda",std::to_string(objectListSizes[j])}};
            this->runCOLTPHL(graph,phl,alt,sltree,defaultQuerySets,fileNameToObjIdx,statsOutputFile,"num landmarks vs lambda",specialFields,testedParameters);
        }
    }

    // Run queries for num landmarks vs landmark type experiment
    if (parameterMap.find("numlmks_vs_lmktype_stats_file_path") == parameterMap.end() || parameterMap["numlmks_vs_lmktype_stats_file_path"] == "") {
        std::cerr << "Num landmarks vs landmark type experiment stats file path not provided in parameters";
        std::exit(1);
    }
    statsOutputFile = parameterMap["numlmks_vs_lmktype_stats_file_path"];
    for (std::size_t i = 0; i < numLandmarks.size(); ++i) {
        for (std::size_t j = 0; j < landmarkTypes.size(); ++j) {
            // Load SULTree for given landmark number, type, and defaults for other variables
            const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
            const std::vector<std::string> sltreeValues = {std::to_string(defFanout),std::to_string(minLambda),std::to_string(numLandmarks[i]),landmarkTypes[j]};
            std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
            SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

            std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
            coltKeys.push_back("lambda");
            coltValues.push_back(std::to_string(defLambda));
            std::unordered_map<std::string,COLTree> fileNameToObjIdx;
            for (std::size_t k = 0; k < numSets; ++k) {
                std::string objIdxFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVariable,k);
                std::string objIdxFilePath = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                COLTree objIdx = serialization::getIndexFromBinaryFile<COLTree>(objIdxFilePath);
                fileNameToObjIdx[objIdxFilePrefix] = objIdx;
            }
            std::vector<std::pair<std::string,std::string>> testedParameters = {{"fanout",std::to_string(defFanout)},{"numlmks",std::to_string(numLandmarks[i])},{"lmktype",landmarkTypes[j]},{"lambda",std::to_string(defLambda)}};
            this->runCOLTPHL(graph,phl,alt,sltree,defaultQuerySets,fileNameToObjIdx,statsOutputFile,"num landmarks vs landmark type",specialFields,testedParameters);
        }
    }

    // Run queries for density vs num landmarks experiment
    if (parameterMap.find("density_vs_numlmks_stats_file_path") == parameterMap.end() || parameterMap["density_vs_numlmks_stats_file_path"] == "") {
        std::cerr << "Density vs num landmarks experiment stats file path not provided in parameters";
        std::exit(1);
    }
    statsOutputFile = parameterMap["density_vs_numlmks_stats_file_path"];
    for (std::size_t i = 0; i < numLandmarks.size(); ++i) {
        // Load SULTree for given landmark number, min lambda, and defaults for other variables
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(defFanout),std::to_string(minLambda),std::to_string(numLandmarks[i]),defLandmarkType};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t j = 0; j < densities.size(); ++j) {
            std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
            coltKeys.push_back("lambda");
            coltValues.push_back(std::to_string(defLambda));
            std::unordered_map<std::string,COLTree> fileNameToObjIdx;
            for (std::size_t k = 0; k < numSets; ++k) {
                std::string objIdxFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                std::string objIdxFilePath = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                COLTree objIdx = serialization::getIndexFromBinaryFile<COLTree>(objIdxFilePath);
                fileNameToObjIdx[objIdxFilePrefix] = objIdx;
            }
            std::vector<std::pair<std::string,std::string>> testedParameters = {{"fanout",std::to_string(defFanout)},{"numlmks",std::to_string(numLandmarks[i])},{"lmktype",defLandmarkType},{"lambda",std::to_string(defLambda)},{"density",std::to_string(densities[j])}};
            std::vector<std::vector<KnnQuery>> densityQuerySets(1,varyingDensityQuerySets[j]);
            this->runCOLTPHL(graph,phl,alt,sltree,densityQuerySets,fileNameToObjIdx,statsOutputFile,"density vs num landmarks",specialFields,testedParameters);
        }
    }

    // Run queries for density vs landmark type experiment
    if (parameterMap.find("density_vs_lmktype_stats_file_path") == parameterMap.end() || parameterMap["density_vs_lmktype_stats_file_path"] == "") {
        std::cerr << "Density vs landmark type experiment stats file path not provided in parameters";
        std::exit(1);
    }
    statsOutputFile = parameterMap["density_vs_lmktype_stats_file_path"];
    for (std::size_t i = 0; i < landmarkTypes.size(); ++i) {
        // Load SULTree for given landmark type, min lambda, and defaults for other variables
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(defFanout),std::to_string(minLambda),std::to_string(defNumLandmarks),landmarkTypes[i]};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t j = 0; j < densities.size(); ++j) {
            std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
            coltKeys.push_back("lambda");
            coltValues.push_back(std::to_string(defLambda));
            std::unordered_map<std::string,COLTree> fileNameToObjIdx;
            for (std::size_t k = 0; k < numSets; ++k) {
                std::string objIdxFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                std::string objIdxFilePath = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                COLTree objIdx = serialization::getIndexFromBinaryFile<COLTree>(objIdxFilePath);
                fileNameToObjIdx[objIdxFilePrefix] = objIdx;
            }
            std::vector<std::pair<std::string,std::string>> testedParameters = {{"fanout",std::to_string(defFanout)},{"numlmks",std::to_string(defNumLandmarks)},{"lmktype",landmarkTypes[i]},{"lambda",std::to_string(defLambda)},{"density",std::to_string(densities[j])}};
            std::vector<std::vector<KnnQuery>> densityQuerySets(1,varyingDensityQuerySets[j]);
            this->runCOLTPHL(graph,phl,alt,sltree,densityQuerySets,fileNameToObjIdx,statsOutputFile,"density vs landmark type",specialFields,testedParameters);
        }
    }

    // Run queries for density vs lambda experiment
    {
        // Use scope to release the default SULTree before the density vs fanout experiment loads its own
        if (parameterMap.find("density_vs_lambda_stats_file_path") == parameterMap.end() || parameterMap["density_vs_lambda_stats_file_path"] == "") {
            std::cerr << "Density vs lambda experiment stats file path not provided in parameters";
            std::exit(1);
        }
        statsOutputFile = parameterMap["density_vs_lambda_stats_file_path"];
        // Load the default SULTree (it does not vary with lambda)
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(defFanout),std::to_string(minLambda),std::to_string(defNumLandmarks),defLandmarkType};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t i = 0; i < objectListSizes.size(); ++i) {
            for (std::size_t j = 0; j < densities.size(); ++j) {
                std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
                coltKeys.push_back("lambda");
                coltValues.push_back(std::to_string(objectListSizes[i]));
                std::unordered_map<std::string,COLTree> fileNameToObjIdx;
                for (std::size_t k = 0; k < numSets; ++k) {
                    std::string objIdxFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                    std::string objIdxFilePath = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                    COLTree objIdx = serialization::getIndexFromBinaryFile<COLTree>(objIdxFilePath);
                    fileNameToObjIdx[objIdxFilePrefix] = objIdx;
                }
                std::vector<std::pair<std::string,std::string>> testedParameters = {{"fanout",std::to_string(defFanout)},{"numlmks",std::to_string(defNumLandmarks)},{"lmktype",defLandmarkType},{"lambda",std::to_string(objectListSizes[i])},{"density",std::to_string(densities[j])}};
                std::vector<std::vector<KnnQuery>> densityQuerySets(1,varyingDensityQuerySets[j]);
                this->runCOLTPHL(graph,phl,alt,sltree,densityQuerySets,fileNameToObjIdx,statsOutputFile,"density vs lambda",specialFields,testedParameters);
            }
        }
    }

    // Run queries for density vs fanout experiment
    if (parameterMap.find("density_vs_fanout_stats_file_path") == parameterMap.end() || parameterMap["density_vs_fanout_stats_file_path"] == "") {
        std::cerr << "Density vs fanout experiment stats file path not provided in parameters";
        std::exit(1);
    }
    statsOutputFile = parameterMap["density_vs_fanout_stats_file_path"];
    for (std::size_t i = 0; i < fanouts.size(); ++i) {
        // Load SULTree for given fanout, min lambda, and defaults for other variables
        const std::vector<std::string> sltreeKeys = {"fanout","maxleafsize","numlmks","lmktype"};
        const std::vector<std::string> sltreeValues = {std::to_string(fanouts[i]),std::to_string(minLambda),std::to_string(defNumLandmarks),defLandmarkType};
        std::string idxOutputFile = parFilePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SULTREE_CMD,sltreeKeys,sltreeValues);
        SULTree sltree = serialization::getIndexFromBinaryFile<SULTree>(idxOutputFile);

        for (std::size_t j = 0; j < densities.size(); ++j) {
            std::vector<std::string> coltKeys = sltreeKeys, coltValues = sltreeValues;
            coltKeys.push_back("lambda");
            coltValues.push_back(std::to_string(defLambda));
            std::unordered_map<std::string,COLTree> fileNameToObjIdx;
            for (std::size_t k = 0; k < numSets; ++k) {
                std::string objIdxFilePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,densities[j],defaultObjVariable,k);
                std::string objIdxFilePath = parFilePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilePrefix,constants::OBJ_IDX_COLTREE,coltKeys,coltValues);
                COLTree objIdx = serialization::getIndexFromBinaryFile<COLTree>(objIdxFilePath);
                fileNameToObjIdx[objIdxFilePrefix] = objIdx;
            }
            std::vector<std::pair<std::string,std::string>> testedParameters = {{"fanout",std::to_string(fanouts[i])},{"numlmks",std::to_string(defNumLandmarks)},{"lmktype",defLandmarkType},{"lambda",std::to_string(defLambda)},{"density",std::to_string(densities[j])}};
            std::vector<std::vector<KnnQuery>> densityQuerySets(1,varyingDensityQuerySets[j]);
            this->runCOLTPHL(graph,phl,alt,sltree,densityQuerySets,fileNameToObjIdx,statsOutputFile,"density vs fanout",specialFields,testedParameters);
        }
    }
}
