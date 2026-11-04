// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "LMKAggregateExperimentsCommand.h"

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
#include "../utility/utility.h"
#include "../utility/serialization.h"
#include "../utility/AggregateMax.h"
#include "../utility/AggregateSum.h"

#include <cstdio>

void LMKAggregateExperimentsCommand::execute(int argc, char* argv[])
{
    std::string experiment = "";
    std::string bgrFilePath = "";
    std::string filePathPrefix = "";
    std::string statsOutputFile = "";
    std::vector<std::string> statsOutputFiles;
    std::string parameters = "";
    std::string objDensities = "";
    std::string objTypes = "";
    std::string objVariable = "";
    unsigned int numSets = 0;
    std::string kValues = "";
    std::string queryVertexSetSizes = "";
    std::string queryVertexSetSubnetSizes = "";
    std::string rwObjSetNames;
    unsigned int minLambda = 0;
    std::string fanouts = "";
    std::string numLandmarks = "";
    std::string landmarkTypes = "";
    std::string lambdas = "";
    std::string defaultAggFunc = "";
    
    /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "e:g:p:f:s:n:d:t:k:m:v:r:y:u:w:x:h:i:j:")) != -1) {
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
                statsOutputFiles.push_back(statsOutputFile);
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
            case 'k':
                kValues = optarg;
                break;
            case 'm':
                queryVertexSetSizes = optarg;
                break;
            case 'y':
                queryVertexSetSubnetSizes = optarg;
                break;
            case 'v':
                objVariable = optarg;
                break;
            case 'r':
                rwObjSetNames = optarg;
                break;
            case 'u':
                minLambda = std::stoul(optarg);;
                break;
            case 'w':
                fanouts = optarg;
                break;
            case 'x':
                numLandmarks = optarg;
                break;
            case 'h':
                landmarkTypes = optarg;
                break;
            case 'i':
                lambdas = optarg;
                break;
            case 'j':
                defaultAggFunc = optarg;
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
    
    if (experiment == this->RUN_KNN_QUERIES) {
        if (argc < 25) {
            // Arguments: -g <binary graph file> -k <k values> -p <parameters> -n <num sets> -d <list of object densities> -t <list of object types> 
            // -v <list of some object variable> -f <index output file path prefix> [-r <rw object set file names>] -m <query vertex set sizes>
            // -y <query vertex set subnet sizes as perc>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || numSets == 0 || objDensities == "" ||
            objTypes == "" || kValues == "" || objVariable == "" || queryVertexSetSizes == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->runQueries2(bgrFilePath,kValues,parameters,numSets,objDensities,objTypes,objVariable,rwObjSetNames,queryVertexSetSizes,queryVertexSetSubnetSizes,filePathPrefix);

    } else if (experiment == this->RUN_KNN_SINGLE_QUERIES) {
        if (argc < 25) {
            // Arguments: -g <binary graph file> -k <k values> -p <parameters> -n <num sets> -d <list of object densities> -t <list of object types> 
            // -v <list of some object variable> -f <index output file path prefix> -s <max stats output file> -s <sum stats output file>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFiles.size() == 0 || numSets == 0 || objDensities == "" || 
            objTypes == "" || kValues == "" || objVariable == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->runSingleQueries2(bgrFilePath,kValues,parameters,numSets,objDensities,objTypes,objVariable,filePathPrefix,statsOutputFiles);

    } else if (experiment == this->RUN_PARAMETER_TEST_KNN_QUERIES) {
        if (argc < 29) {
            // Arguments: -g <binary graph file> -p <parameters> -f <index output file path prefix>
            // -u <minimum lambda> -w <fanouts> -x <num landmarks> -h <landmark types> -i <lambdas> 
            // -j <default aggregate function> -n <num sets> -d <list of object densities> -t <list of object types> -v <list of some object variable> 
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" ||
            minLambda == 0 || fanouts == "" || numLandmarks == "" || landmarkTypes == "" || lambdas == "" || 
            numSets == 0 || objDensities == "" || objTypes == "" || objVariable == "" || defaultAggFunc == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        if (defaultAggFunc == "sum") {
            AggregateSum<EdgeWeight> agg;
            this->runCOLTParameterTestQueries(bgrFilePath,parameters,filePathPrefix,numSets,objDensities,objTypes,objVariable,fanouts,numLandmarks,landmarkTypes,lambdas,minLambda,agg);
        } else if (defaultAggFunc == "max") {
            AggregateMax<EdgeWeight> agg;
            this->runCOLTParameterTestQueries(bgrFilePath,parameters,filePathPrefix,numSets,objDensities,objTypes,objVariable,fanouts,numLandmarks,landmarkTypes,lambdas,minLambda,agg);
        } else {
            std::cerr << "Invalid default aggregate function argument -j " << defaultAggFunc << ". Must be \"max\" or \"sum\"\n\n";
            exit(1);
        }

    } else {
        std::cerr << "Invalid experimental step!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
}

void LMKAggregateExperimentsCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::LMK_AGG_EXPERIMENTS_CMD + " -e <experimental step>\n\n"
              << "Steps:\n"
              << utility::getFormattedUsageString(this->RUN_KNN_QUERIES,"Run all landmark kNN experiments") + "\n"
              << utility::getFormattedUsageString(this->RUN_KNN_SINGLE_QUERIES,"Run single landmark kNN experiments") + "\n"
              << utility::getFormattedUsageString(this->RUN_PARAMETER_TEST_KNN_QUERIES,"Run kNN queries for parameter experiments") + "\n";
}

void LMKAggregateExperimentsCommand::showPhaseUsage(std::string method, std::string programName)
{
    if (method == this->RUN_KNN_QUERIES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_AGG_EXPERIMENTS_CMD
                  << " -e " + this->RUN_KNN_QUERIES + " -g <binary graph file>\n"
                  << "-k <k values> -p <parameters> -n <num sets> -d <list of object set densities>\n"
                  << "-t <list of object types>  -v <list of some object variable> [-r <rw object set file names>] -m <query vertex set sizes> -y <query vertex set subnet sizes as perc> -f <index output file path prefix>\n";
    } else if (method == this->RUN_KNN_SINGLE_QUERIES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_AGG_EXPERIMENTS_CMD
                  << " -e " + this->RUN_KNN_SINGLE_QUERIES + " -g <binary graph file>\n"
                  << "-k <k values> -p <parameters> -n <num sets> -d <list of object set densities>\n"
                  << "-t <list of object types>  -v <list of some object variable> -f <index output file path prefix>\n"
                  << "-s <max stats output file> -s <sum stats output file>\n";
    } else if (method == this->RUN_PARAMETER_TEST_KNN_QUERIES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_AGG_EXPERIMENTS_CMD
                  << " -e " + this->RUN_PARAMETER_TEST_KNN_QUERIES + " -g <binary graph file> -p <parameters> -f <index output file path prefix>\n"
                  << "-u <minimum lambda> -w <fanouts> -x <num landmarks> -h <landmark types> -i <lambdas> \n"
                  << "-j <default aggregate function> -n <num sets> -d <list of object densities> -t <list of object types> -v <list of some object variable> ";
    } else {
        std::cerr << "Invalid experiment phase!" << std::endl;
        this->showCommandUsage(programName);
    }
}

std::unordered_map<std::string, std::string> LMKAggregateExperimentsCommand::getParameters(std::string parameters)
{
    std::unordered_map<std::string,std::string> parameterMap;
    std::vector<std::string> pairs = utility::splitByDelim(parameters,',');
    
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        std::vector<std::string> pair = utility::splitByDelim(pairs[i],'=');
        if (pair.size() == 2) {
            parameterMap[pair[0]] = pair[1];
        } else {
            std::cerr << "Invalid key-value pair in parameter string" << std::endl;
            std::cerr << "Parameter string: " << pairs[i] << std::endl;
            exit(1);
        }
    }
    return parameterMap;
}

void LMKAggregateExperimentsCommand::writeBlankLines(std::string fileName)
{
    if (fileName == "") {
        return; // A switched-off query list has no stats file
    }
    std::ofstream fs;
    fs.open(fileName, std::ios_base::app | std::ios_base::out);
    fs << "\n\n"; //  GNU Plot requires double line separation between data points for different plot lines
    fs.close();
}


void LMKAggregateExperimentsCommand::runQueries2(std::string bgrFileName, std::string kValues, std::string parameters, std::size_t numSets, 
                                                 std::string objDensities, std::string objTypes, std::string objVariable, std::string rwObjSetNames, 
                                                 std::string queryVertexSetSizes, std::string queryVertexSetSubnetSizes, std::string filePathPrefix)
{
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);

    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    
    if (numSets == 0) {
        std::cerr << "No object sets were provided!\n";
        exit(1);    
    }
    
    // Get Default Values for Variables
    unsigned int defaultK = std::stoul(parameterMap["aknn_def_k"]);
    double defaultDensity = std::stod(parameterMap["aknn_def_density"]);
    std::string defaultObjType = parameterMap["aknn_def_type"];
    if (defaultObjType == "") {
        std::cerr << "Invalid default object set type of " << defaultObjType << "!\n";
        std::exit(1);    
    }
    int defaultObjVar = std::stoi(parameterMap["aknn_def_var"]);
    
    int defaultNumQueryVertices = std::stoi(parameterMap["aknn_qv_def_setsize"]);
    if (defaultNumQueryVertices <= 0) {
        std::cerr << "Invalid default query vertex set size " << defaultNumQueryVertices << "!\n";
        std::exit(1);    
    }
    
    int defaultQueryVertexSubnetSize = std::stoi(parameterMap["aknn_qv_def_setsubnetsize"]);
    int maxNumQueryVertices = std::stoi(parameterMap["aknn_qv_max_setsize"]);
    int numQueryVertexSets = std::stoi(parameterMap["aknn_num_qv_sets"]);
    if (defaultQueryVertexSubnetSize <= 0 || defaultQueryVertexSubnetSize > 100) {
        std::cerr << "Invalid default query vertex set subnet size " << defaultQueryVertexSubnetSize << "!\n";
        std::exit(1);
    }
    if (maxNumQueryVertices <= 1) {
        std::cerr << "Invalid max query vertex set size " << maxNumQueryVertices << "!\n";
        std::exit(1);
    }
    if (numQueryVertexSets <= 1) {
        std::cerr << "Invalid number of query vertex sets " << numQueryVertexSets << "!\n";
        std::exit(1);
    }

    std::string vectorType;
    unsigned int numElementsPerVec;
    
    // Load default subnetwork size query vertex sets
    std::string defSubnetQueryVertexFilePath = filePathPrefix + "/query_sets/" + utility::constructQuerySetFileName(graph.getNetworkName(),"subnet",numQueryVertexSets,maxNumQueryVertices,defaultQueryVertexSubnetSize);
    std::vector<std::vector<VertexID>> queryVertexSets = utility::getVectorSetFromFile(defSubnetQueryVertexFilePath,vectorType,numElementsPerVec);

    if (queryVertexSets.size() == 0) {
        std::cerr << "No query points were provided!\n";
        exit(1);    
    }
    
    // Create query vertex set vector for default set size    
    std::vector<std::vector<VertexID>> defaultSizeQueryVertexSetsVec;
    for (std::size_t i = 0; i < queryVertexSets.size(); ++i) {
        std::vector<VertexID> tempVec;
        for (std::size_t j = 0; j < static_cast<std::size_t>(defaultNumQueryVertices); ++j) {
            tempVec.push_back(queryVertexSets[i][j]);
        }
        defaultSizeQueryVertexSetsVec.push_back(tempVec);
    }
    
    /*
     * Create Query Sets for Each Experiment
     */
    
    // Each query list runs unless its switch in the parameters is off (utility::isQueryListOn). A switched-off list is not built, and
    // its stats file name is left empty, so that neither its runs nor the blank separator lines write to its stats file
    bool runVaryingK = utility::isQueryListOn(parameterMap,"aknn_varyingk");
    bool runVaryingDensity = utility::isQueryListOn(parameterMap,"aknn_varyingdensity");
    bool runVaryingQueryVertexSetSize = utility::isQueryListOn(parameterMap,"aknn_varyingqvsetsize");
    bool runVaryingQueryVertexSetSubnetSize = utility::isQueryListOn(parameterMap,"aknn_varyingqvsetsubnetsize");
    
    /*
     * Varying K Experiment
     */
    std::vector<std::vector<AknnQuery>> varyingKQuerySets;

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
            std::vector<AknnQuery> kQuerySet;
            for (std::size_t j = 0; j < numSets; ++j) {
                std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
                for (std::size_t k = 0; k < defaultSizeQueryVertexSetsVec.size(); ++k) {
                    AknnQuery query(defaultSizeQueryVertexSetsVec[k],kValuesVec[i],defaultDensity,defaultQueryVertexSubnetSize,defaultObjType,defaultObjVar,j,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*defaultDensity));
                    kQuerySet.push_back(query);
                }
            }
            varyingKQuerySets.push_back(kQuerySet);
        }
    }
    
    if (parameterMap.find("max_aknn_varyingk_stats_file") == parameterMap.end() || parameterMap["max_aknn_varyingk_stats_file"] == "") {
        std::cerr << "Max varying k stats file not provided in parameters";
        std::exit(1);
    }
    std::string maxVaryingKStatsFile = runVaryingK ? filePathPrefix + "/stats/" + parameterMap["max_aknn_varyingk_stats_file"] : "";
    
    if (parameterMap.find("sum_aknn_varyingk_stats_file") == parameterMap.end() || parameterMap["sum_aknn_varyingk_stats_file"] == "") {
        std::cerr << "Sum varying k stats file not provided in parameters";
        std::exit(1);
    }
    std::string sumVaryingKStatsFile = runVaryingK ? filePathPrefix + "/stats/" + parameterMap["sum_aknn_varyingk_stats_file"] : "";
    
    
    /* 
     * Varying Density Experiment
     */
    std::vector<std::vector<AknnQuery>> varyingDensityQuerySets;

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
            std::vector<AknnQuery> densityQuerySet;
            for (std::size_t j = 0; j < numSets; ++j) {
                std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,objDensitiesVec[i],defaultObjVar,j);
                for (std::size_t k = 0; k < defaultSizeQueryVertexSetsVec.size(); ++k) {
                    AknnQuery query(defaultSizeQueryVertexSetsVec[k],defaultK,objDensitiesVec[i],defaultQueryVertexSubnetSize,defaultObjType,defaultObjVar,j,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*objDensitiesVec[i]));
                    densityQuerySet.push_back(query);
                }
            }
            varyingDensityQuerySets.push_back(densityQuerySet);
        }
    }
    
    if (parameterMap.find("max_aknn_varyingdensity_stats_file") == parameterMap.end() || parameterMap["max_aknn_varyingdensity_stats_file"] == "") {
        std::cerr << "Max varying density stats file not provided in parameters";
        std::exit(1);
    }
    std::string maxVaryingDensityStatsFile = runVaryingDensity ? filePathPrefix + "/stats/" + parameterMap["max_aknn_varyingdensity_stats_file"] : "";

    if (parameterMap.find("sum_aknn_varyingdensity_stats_file") == parameterMap.end() || parameterMap["sum_aknn_varyingdensity_stats_file"] == "") {
        std::cerr << "Sum  varying density stats file not provided in parameters";
        std::exit(1);
    }
    std::string sumVaryingDensityStatsFile = runVaryingDensity ? filePathPrefix + "/stats/" + parameterMap["sum_aknn_varyingdensity_stats_file"] : "";

    
    /*
     * Varying Real-World Object Set Experiment
     */
    std::vector<std::vector<AknnQuery>> varyingRWObjSetQuerySets;

    std::vector<std::string> rwObjSetNamesVec = utility::splitByDelim(rwObjSetNames,',');
    bool runRealWorld = parameterMap["aknn_rw"] == "1";
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
        // Each real-world POI type has only one object set, so we use a separate query vertex set file with
        // (number of object sets x number of query vertex sets) query vertex sets to run the same number of
        // queries as each synthetic object set experiment
        // Note: The file must be a subnet query vertex set file for the default subnetwork size
        if (parameterMap.find("aknn_rw_qv_file") == parameterMap.end() || parameterMap["aknn_rw_qv_file"] == "") {
            std::cerr << "Real-world query vertex set file not provided in parameters!\n";
            std::exit(1);
        }
        unsigned int rwNumElementsPerVec;
        std::vector<std::vector<VertexID>> rwQueryVertexSets = utility::getVectorSetFromFile(filePathPrefix + "/query_sets/" + parameterMap["aknn_rw_qv_file"],vectorType,rwNumElementsPerVec);
        if (rwQueryVertexSets.size() != numSets*defaultSizeQueryVertexSetsVec.size()) {
            std::cerr << "Real-world query vertex set file has " << rwQueryVertexSets.size() << " query vertex sets but " << numSets*defaultSizeQueryVertexSetsVec.size()
                      << " (number of object sets x number of query vertex sets) are required!\n";
            std::exit(1);
        }
        if (rwNumElementsPerVec < static_cast<unsigned int>(defaultNumQueryVertices)) {
            std::cerr << "Real-world query vertex sets have " << rwNumElementsPerVec << " query vertices but at least " << defaultNumQueryVertices << " are required!\n";
            std::exit(1);
        }

        // Create query vertex set vector for default set size
        std::vector<std::vector<VertexID>> rwDefaultSizeQueryVertexSetsVec;
        for (std::size_t i = 0; i < rwQueryVertexSets.size(); ++i) {
            rwDefaultSizeQueryVertexSetsVec.push_back(std::vector<VertexID>(rwQueryVertexSets[i].begin(),rwQueryVertexSets[i].begin()+defaultNumQueryVertices));
        }

        for (std::size_t i = 0; i < rwObjSetNamesVec.size(); ++i) {
            std::vector<AknnQuery> rwPOITypeQuerySet;
            std::string objIdxFilenamePrefix = graph.getNetworkName() + "_" + rwObjSetNamesVec[i];
            for (std::size_t k = 0; k < rwDefaultSizeQueryVertexSetsVec.size(); ++k) {
                AknnQuery query(rwDefaultSizeQueryVertexSetsVec[k],defaultK,0,defaultQueryVertexSubnetSize,rwObjSetNamesVec[i],0,0,objIdxFilenamePrefix,0);
                rwPOITypeQuerySet.push_back(query);
            }
            varyingRWObjSetQuerySets.push_back(rwPOITypeQuerySet);
        }
    }
    
    if (parameterMap.find("max_aknn_varyingrwset_stats_file") == parameterMap.end() || parameterMap["max_aknn_varyingrwset_stats_file"] == "") {
        std::cerr << "Max varying real-world set stats file not provided in parameters";
        std::exit(1);
    }
    std::string maxVaryingRWObjSetStatsFile = filePathPrefix + "/stats/" + parameterMap["max_aknn_varyingrwset_stats_file"];
    
    if (parameterMap.find("sum_aknn_varyingrwset_stats_file") == parameterMap.end() || parameterMap["sum_aknn_varyingrwset_stats_file"] == "") {
        std::cerr << "Sum varying real-world set stats file not provided in parameters";
        std::exit(1);
    }
    std::string sumVaryingRWObjSetStatsFile = filePathPrefix + "/stats/" + parameterMap["sum_aknn_varyingrwset_stats_file"];

    
    /* 
     * Varying Query Vertex Set Experiment
     */
    std::vector<std::vector<AknnQuery>> varyingQueryVertexSetSizeQuerySets;
    std::vector<std::string> strSizeValues = utility::splitByDelim(queryVertexSetSizes,',');
    std::vector<unsigned int> sizeValuesVec;
    for(std::size_t i = 0; i < strSizeValues.size(); ++i) {
        unsigned int size = std::stoul(strSizeValues[i]);
        if (size > 0) {
            sizeValuesVec.push_back(size);
        } else {
            std::cerr << "Invalid query vertex set size '" << size << "' in list provided!\n";
            std::exit(1);    
        }
    }
    
    if (runVaryingQueryVertexSetSize) {
        for (std::size_t i = 0; i < sizeValuesVec.size(); ++i) {
            // Need first create query term vector set with only the 
            // the number of terms for the current query num terms value
            std::vector<std::vector<VertexID>> currentSizeQueryVertexSets;
            for (std::size_t j = 0; j < queryVertexSets.size(); ++j) {
                std::vector<VertexID> tempVec;
                for (std::size_t k = 0; k < sizeValuesVec[i]; ++k) {
                    tempVec.push_back(queryVertexSets[j][k]);
                }
                currentSizeQueryVertexSets.push_back(tempVec);
            }

            std::vector<AknnQuery> queryVertexSetSizeQuerySet;
            for (std::size_t j = 0; j < numSets; ++j) {
                std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
                for (std::size_t k = 0; k < currentSizeQueryVertexSets.size(); ++k) {
                    AknnQuery query(currentSizeQueryVertexSets[k],defaultK,defaultDensity,defaultQueryVertexSubnetSize,defaultObjType,defaultObjVar,j,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*defaultDensity));
                    queryVertexSetSizeQuerySet.push_back(query);
                }
            }
            varyingQueryVertexSetSizeQuerySets.push_back(queryVertexSetSizeQuerySet);
        }
    }

    if (parameterMap.find("max_aknn_varyingqvsetsize_stats_file") == parameterMap.end() || parameterMap["max_aknn_varyingqvsetsize_stats_file"] == "") {
        std::cerr << "Max varying query set size file not provided in parameters";
        std::exit(1);
    }
    std::string maxVaryingQueryVertexSetSizeStatsFile = runVaryingQueryVertexSetSize ? filePathPrefix + "/stats/" + parameterMap["max_aknn_varyingqvsetsize_stats_file"] : "";
    
    if (parameterMap.find("sum_aknn_varyingqvsetsize_stats_file") == parameterMap.end() || parameterMap["sum_aknn_varyingqvsetsize_stats_file"] == "") {
        std::cerr << "Sum varying query set size file not provided in parameters";
        std::exit(1);
    }
    std::string sumVaryingQueryVertexSetSizeStatsFile = runVaryingQueryVertexSetSize ? filePathPrefix + "/stats/" + parameterMap["sum_aknn_varyingqvsetsize_stats_file"] : "";
    
    /*
     * Varying Query Vertex Subnetwork Size Experiment
     */
    
    std::vector<std::vector<AknnQuery>> varyingQueryVertexSetSubnetSizeQuerySets;

    std::vector<std::string> strQueryVertexSetSubnetSizes = utility::splitByDelim(queryVertexSetSubnetSizes,',');
    std::vector<unsigned int> queryVertexSubnetSizes;
    for(std::size_t i = 0; i < strQueryVertexSetSubnetSizes.size(); ++i) {
        unsigned int subnetSize = std::stoul(strQueryVertexSetSubnetSizes[i]);
        if (subnetSize > 0 && subnetSize <= 100) {
            queryVertexSubnetSizes.push_back(subnetSize);
        } else {
            std::cerr << "Invalid default query vertex set subnet size " << subnetSize << "!\n";
            std::exit(1);    
        }
    }    
    
    if (runVaryingQueryVertexSetSubnetSize) {
        for (std::size_t i = 0; i < queryVertexSubnetSizes.size(); ++i) {
            std::vector<AknnQuery> distQuerySet;
        
            std::string distQueryVertexFilePath = filePathPrefix + "/query_sets/" + utility::constructQuerySetFileName(graph.getNetworkName(),"subnet",numQueryVertexSets,maxNumQueryVertices,queryVertexSubnetSizes[i]);
            std::vector<std::vector<VertexID>> distQueryVertexSets = utility::getVectorSetFromFile(distQueryVertexFilePath,vectorType,numElementsPerVec);

            // Create query vertex set vector for default set size
            std::vector<std::vector<VertexID>> defSizeSubnetSizeQueryVertexSet;
            for (std::size_t j = 0; j < distQueryVertexSets.size(); ++j) {
                std::vector<VertexID> tempVec;
                for (std::size_t k = 0; k < static_cast<std::size_t>(defaultNumQueryVertices); ++k) {
                    tempVec.push_back(distQueryVertexSets[j][k]);
                }
                defSizeSubnetSizeQueryVertexSet.push_back(tempVec);
            }
        
            for (std::size_t j = 0; j < numSets; ++j) {
                std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
                for (std::size_t k = 0; k < defSizeSubnetSizeQueryVertexSet.size(); ++k) {
                    AknnQuery query(defSizeSubnetSizeQueryVertexSet[k],defaultK,defaultDensity,queryVertexSubnetSizes[i],defaultObjType,defaultObjVar,j,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*defaultDensity));
                    distQuerySet.push_back(query);
                }
                if (distQuerySet.size() == 0) {
                    std::cerr << "No query sets for distance = " << queryVertexSubnetSizes[i] << "!\n";
                    std::exit(1);
                }
            }
            varyingQueryVertexSetSubnetSizeQuerySets.push_back(distQuerySet);
        }
    }

    if (parameterMap.find("max_aknn_varyingqvsetsubnetsize_stats_file") == parameterMap.end() || parameterMap["max_aknn_varyingqvsetsubnetsize_stats_file"] == "") {
        std::cerr << "Max varying query set subnet size file not provided in parameters";
        std::exit(1);
    }
    std::string maxVaryingQueryVertexSetSubnetSizeStatsFile = runVaryingQueryVertexSetSubnetSize ? filePathPrefix + "/stats/" + parameterMap["max_aknn_varyingqvsetsubnetsize_stats_file"] : "";
    
    if (parameterMap.find("sum_aknn_varyingqvsetsubnetsize_stats_file") == parameterMap.end() || parameterMap["sum_aknn_varyingqvsetsubnetsize_stats_file"] == "") {
        std::cerr << "Sum varying query set subnet size file not provided in parameters";
        std::exit(1);
    }
    std::string sumVaryingQueryVertexSetSubnetSizeStatsFile = runVaryingQueryVertexSetSubnetSize ? filePathPrefix + "/stats/" + parameterMap["sum_aknn_varyingqvsetsubnetsize_stats_file"] : "";
    
    bool queryIERPHL = parameterMap["ier_phl"] == "1";
    bool queryNVDPHL = parameterMap["nvd_phl"] == "1";
    bool queryCOLTPHL = parameterMap["colt_phl"] == "1";
    
    AggregateSum<EdgeWeight> aggSum;
    AggregateSum<EuclideanDistance> euclidAggSum;
    AggregateMax<EdgeWeight> aggMax;
    AggregateMax<EuclideanDistance> euclidAggMax;

    if (queryIERPHL || queryNVDPHL || queryCOLTPHL) {
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

            this->runIERPHL(graph,phl,varyingKQuerySets,fileNameToObjIdx,maxVaryingKStatsFile,"varyingk",aggMax,euclidAggMax);
            this->runIERPHL(graph,phl,varyingDensityQuerySets,fileNameToObjIdx,maxVaryingDensityStatsFile,"varyingdensity",aggMax,euclidAggMax);
            this->runIERPHL(graph,phl,varyingQueryVertexSetSizeQuerySets,fileNameToObjIdx,maxVaryingQueryVertexSetSizeStatsFile,"varyingqvsetsize",aggMax,euclidAggMax);
            this->runIERPHL(graph,phl,varyingQueryVertexSetSubnetSizeQuerySets,fileNameToObjIdx,maxVaryingQueryVertexSetSubnetSizeStatsFile,"varyingqvsetsubnetsize",aggMax,euclidAggMax);
            this->runIERPHL(graph,phl,varyingKQuerySets,fileNameToObjIdx,sumVaryingKStatsFile,"varyingk",aggSum,euclidAggSum);
            this->runIERPHL(graph,phl,varyingDensityQuerySets,fileNameToObjIdx,sumVaryingDensityStatsFile,"varyingdensity",aggSum,euclidAggSum);
            this->runIERPHL(graph,phl,varyingQueryVertexSetSizeQuerySets,fileNameToObjIdx,sumVaryingQueryVertexSetSizeStatsFile,"varyingqvsetsize",aggSum,euclidAggSum);
            this->runIERPHL(graph,phl,varyingQueryVertexSetSubnetSizeQuerySets,fileNameToObjIdx,sumVaryingQueryVertexSetSubnetSizeStatsFile,"varyingqvsetsubnetsize",aggSum,euclidAggSum);
            if (runRealWorld) {
                this->runIERPHL(graph,phl,varyingRWObjSetQuerySets,fileNameToObjIdx,maxVaryingRWObjSetStatsFile,"varyingrwset",aggMax,euclidAggMax);
                this->runIERPHL(graph,phl,varyingRWObjSetQuerySets,fileNameToObjIdx,sumVaryingRWObjSetStatsFile,"varyingrwset",aggSum,euclidAggSum);
            }

            this->writeBlankLines(maxVaryingKStatsFile);
            this->writeBlankLines(sumVaryingKStatsFile);
            this->writeBlankLines(maxVaryingDensityStatsFile);
            this->writeBlankLines(sumVaryingDensityStatsFile);
            this->writeBlankLines(maxVaryingQueryVertexSetSizeStatsFile);
            this->writeBlankLines(sumVaryingQueryVertexSetSizeStatsFile);
            this->writeBlankLines(maxVaryingQueryVertexSetSubnetSizeStatsFile);
            this->writeBlankLines(sumVaryingQueryVertexSetSubnetSizeStatsFile);
        }
        
        if (queryNVDPHL || queryCOLTPHL) {
            const std::vector<std::string> altKeys = {"landmarks","type"};
            const std::vector<std::string> altValues = {parameterMap["alt_numlandmarks"],parameterMap["alt_landmarktype"]};
            std::string altIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,altKeys,altValues);
            ALT alt = serialization::getIndexFromBinaryFile<ALT>(altIdxFilePath);
            
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

                this->runNVDPHL(graph,phl,alt,varyingKQuerySets,fileNameToObjIdx,maxVaryingKStatsFile,"varyingk",aggMax);
                this->runNVDPHL(graph,phl,alt,varyingDensityQuerySets,fileNameToObjIdx,maxVaryingDensityStatsFile,"varyingdensity",aggMax);
                this->runNVDPHL(graph,phl,alt,varyingQueryVertexSetSizeQuerySets,fileNameToObjIdx,maxVaryingQueryVertexSetSizeStatsFile,"varyingqvsetsize",aggMax);
                this->runNVDPHL(graph,phl,alt,varyingQueryVertexSetSubnetSizeQuerySets,fileNameToObjIdx,maxVaryingQueryVertexSetSubnetSizeStatsFile,"varyingqvsetsubnetsize",aggMax);
                this->runNVDPHL(graph,phl,alt,varyingKQuerySets,fileNameToObjIdx,sumVaryingKStatsFile,"varyingk",aggSum);
                this->runNVDPHL(graph,phl,alt,varyingDensityQuerySets,fileNameToObjIdx,sumVaryingDensityStatsFile,"varyingdensity",aggSum);
                this->runNVDPHL(graph,phl,alt,varyingQueryVertexSetSizeQuerySets,fileNameToObjIdx,sumVaryingQueryVertexSetSizeStatsFile,"varyingqvsetsize",aggSum);
                this->runNVDPHL(graph,phl,alt,varyingQueryVertexSetSubnetSizeQuerySets,fileNameToObjIdx,sumVaryingQueryVertexSetSubnetSizeStatsFile,"varyingqvsetsubnetsize",aggSum);
                if (runRealWorld) {
                    this->runNVDPHL(graph,phl,alt,varyingRWObjSetQuerySets,fileNameToObjIdx,maxVaryingRWObjSetStatsFile,"varyingrwset",aggMax);
                    this->runNVDPHL(graph,phl,alt,varyingRWObjSetQuerySets,fileNameToObjIdx,sumVaryingRWObjSetStatsFile,"varyingrwset",aggSum);
                }

                this->writeBlankLines(maxVaryingKStatsFile);
                this->writeBlankLines(sumVaryingKStatsFile);
                this->writeBlankLines(maxVaryingDensityStatsFile);
                this->writeBlankLines(sumVaryingDensityStatsFile);
                this->writeBlankLines(maxVaryingQueryVertexSetSizeStatsFile);
                this->writeBlankLines(sumVaryingQueryVertexSetSizeStatsFile);
                this->writeBlankLines(maxVaryingQueryVertexSetSubnetSizeStatsFile);
                this->writeBlankLines(sumVaryingQueryVertexSetSubnetSizeStatsFile);
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

                this->runCOLTPHL(graph,phl,alt,sltree,varyingKQuerySets,fileNameToObjIdx,maxVaryingKStatsFile,"varyingk",aggMax);
                this->runCOLTPHL(graph,phl,alt,sltree,varyingDensityQuerySets,fileNameToObjIdx,maxVaryingDensityStatsFile,"varyingdensity",aggMax);
                this->runCOLTPHL(graph,phl,alt,sltree,varyingQueryVertexSetSizeQuerySets,fileNameToObjIdx,maxVaryingQueryVertexSetSizeStatsFile,"varyingqvsetsize",aggMax);
                this->runCOLTPHL(graph,phl,alt,sltree,varyingQueryVertexSetSubnetSizeQuerySets,fileNameToObjIdx,maxVaryingQueryVertexSetSubnetSizeStatsFile,"varyingqvsetsubnetsize",aggMax);
                this->runCOLTPHL(graph,phl,alt,sltree,varyingKQuerySets,fileNameToObjIdx,sumVaryingKStatsFile,"varyingk",aggSum);
                this->runCOLTPHL(graph,phl,alt,sltree,varyingDensityQuerySets,fileNameToObjIdx,sumVaryingDensityStatsFile,"varyingdensity",aggSum);
                this->runCOLTPHL(graph,phl,alt,sltree,varyingQueryVertexSetSizeQuerySets,fileNameToObjIdx,sumVaryingQueryVertexSetSizeStatsFile,"varyingqvsetsize",aggSum);
                this->runCOLTPHL(graph,phl,alt,sltree,varyingQueryVertexSetSubnetSizeQuerySets,fileNameToObjIdx,sumVaryingQueryVertexSetSubnetSizeStatsFile,"varyingqvsetsubnetsize",aggSum);
                if (runRealWorld) {
                    this->runCOLTPHL(graph,phl,alt,sltree,varyingRWObjSetQuerySets,fileNameToObjIdx,maxVaryingRWObjSetStatsFile,"varyingrwset",aggMax);
                    this->runCOLTPHL(graph,phl,alt,sltree,varyingRWObjSetQuerySets,fileNameToObjIdx,sumVaryingRWObjSetStatsFile,"varyingrwset",aggSum);
                }

                this->writeBlankLines(maxVaryingKStatsFile);
                this->writeBlankLines(sumVaryingKStatsFile);
                this->writeBlankLines(maxVaryingDensityStatsFile);
                this->writeBlankLines(sumVaryingDensityStatsFile);
                this->writeBlankLines(maxVaryingQueryVertexSetSizeStatsFile);
                this->writeBlankLines(sumVaryingQueryVertexSetSizeStatsFile);
                this->writeBlankLines(maxVaryingQueryVertexSetSubnetSizeStatsFile);
                this->writeBlankLines(sumVaryingQueryVertexSetSubnetSizeStatsFile);
            }
        }
    }
    
}

void LMKAggregateExperimentsCommand::runSingleQueries2(std::string bgrFileName, std::string kValues, std::string parameters, 
                                              std::size_t numSets, std::string objDensities, std::string objTypes, std::string objVariable, 
                                              std::string filePathPrefix, std::vector<std::string>& statsOutputFiles)
{

    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);

    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    
    // Find all additional fields we need to add to stats tuples
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
    
    std::string vectorType;
    unsigned int numElementsPerVec;
    if (numSets == 0) {
        std::cerr << "No object sets were provided!\n";
        exit(1);    
    }
    
    if (statsOutputFiles.size() != 2) {
        std::cerr << "Not enough stats output files provide, need one for each aggregate function max and sum!\n";
        exit(1);    
    }
    std::string maxStatsOutputFile = statsOutputFiles[0];
    if (maxStatsOutputFile.find("max") == std::string::npos){
        std::cerr << "Max stats output file " << maxStatsOutputFile << " does not seem to have max in its name, is it correct?\n";
        exit(1);    
    }
    std::string sumStatsOutputFile = statsOutputFiles[1];
    if (sumStatsOutputFile.find("sum") == std::string::npos){
        std::cerr << "Sum stats output file " << sumStatsOutputFile << " does not seem to have sum in its name, is it correct?\n";
        exit(1);    
    }
    
    // Get Default Values for Variables
    unsigned int defaultK = std::stoul(parameterMap["aknn_def_k"]);
    double defaultDensity = std::stod(parameterMap["aknn_def_density"]);
    std::string defaultObjType = parameterMap["aknn_def_type"];
    int defaultObjVar = std::stoi(parameterMap["aknn_def_var"]);
    if (defaultObjType == "") {
        std::cerr << "Invalid default object set type of " << defaultObjType << "!\n";
        std::exit(1);    
    }
    int defaultNumQueryVertices = std::stoi(parameterMap["aknn_qv_def_setsize"]);
    if (defaultNumQueryVertices <= 0) {
        std::cerr << "Invalid default query vertex set size " << defaultNumQueryVertices << "!\n";
        std::exit(1);    
    }
    
    int defaultQueryVertexSubnetSize = std::stoi(parameterMap["aknn_qv_def_setsubnetsize"]);
    int maxNumQueryVertices = std::stoi(parameterMap["aknn_qv_max_setsize"]);
    int numQueryVertexSets = std::stoi(parameterMap["aknn_num_qv_sets"]);
    if (defaultQueryVertexSubnetSize <= 0 || defaultQueryVertexSubnetSize > 100) {
        std::cerr << "Invalid default query vertex set subnet size " << defaultQueryVertexSubnetSize << "!\n";
        std::exit(1);    
    }
    if (maxNumQueryVertices <= 1) {
        std::cerr << "Invalid max query vertex set size " << maxNumQueryVertices << "!\n";
        std::exit(1);    
    }
    if (numQueryVertexSets <= 1) {
        std::cerr << "Invalid number of query vertex sets " << numQueryVertexSets << "!\n";
        std::exit(1);    
    }

    // Load default subnetwork size query vertex sets (same sets as runQueries2)
    std::string defSubnetQueryVertexFilePath = filePathPrefix + "/query_sets/" + utility::constructQuerySetFileName(graph.getNetworkName(),"subnet",numQueryVertexSets,maxNumQueryVertices,defaultQueryVertexSubnetSize);
    std::vector<std::vector<VertexID>> queryVertexSets = utility::getVectorSetFromFile(defSubnetQueryVertexFilePath,vectorType,numElementsPerVec);
    if (queryVertexSets.size() == 0) {
        std::cerr << "No query points were provided!\n";
        exit(1);    
    }

    // Create query vertex set vector for default set size    
    std::vector<std::vector<VertexID>> defaultSizeQueryVertexSetsVec;
    for (std::size_t i = 0; i < queryVertexSets.size(); ++i) {
        std::vector<VertexID> tempVec;
        for (std::size_t j = 0; j < static_cast<std::size_t>(defaultNumQueryVertices); ++j) {
            tempVec.push_back(queryVertexSets[i][j]);
        }
        defaultSizeQueryVertexSetsVec.push_back(tempVec);
    }
    
    /*
     * Create set of queries to run on each method
     */
    std::vector<std::vector<AknnQuery>> singleQuerySet;
    std::vector<AknnQuery> querySet;
    for (std::size_t i = 0; i < numSets; ++i) {
        std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,i);
        for (std::size_t j = 0; j < defaultSizeQueryVertexSetsVec.size(); ++j) {
            AknnQuery query(defaultSizeQueryVertexSetsVec[j],defaultK,defaultDensity,defaultQueryVertexSubnetSize,defaultObjType,defaultObjVar,i,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*defaultDensity));
            querySet.push_back(query);
        }
    }
    singleQuerySet.push_back(querySet);
    
    bool queryIERPHL = parameterMap["ier_phl"] == "1";
    bool queryNVDPHL = parameterMap["nvd_phl"] == "1";
    bool queryCOLTPHL = parameterMap["colt_phl"] == "1";

    AggregateSum<EdgeWeight> aggSum;
    AggregateSum<EuclideanDistance> euclidAggSum;
    AggregateMax<EdgeWeight> aggMax;
    AggregateMax<EuclideanDistance> euclidAggMax;

    if (queryIERPHL || queryNVDPHL || queryCOLTPHL) {
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
        
            this->runIERPHL(graph,phl,singleQuerySet,fileNameToObjIdx,maxStatsOutputFile,"single_query_set",aggMax,euclidAggMax,specialFields);
            this->runIERPHL(graph,phl,singleQuerySet,fileNameToObjIdx,sumStatsOutputFile,"single_query_set",aggSum,euclidAggSum,specialFields);
        }
        
        if (queryNVDPHL || queryCOLTPHL) {
            const std::vector<std::string> altKeys = {"landmarks","type"};
            const std::vector<std::string> altValues = {parameterMap["alt_numlandmarks"],parameterMap["alt_landmarktype"]};
            std::string altIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,altKeys,altValues);
            ALT alt = serialization::getIndexFromBinaryFile<ALT>(altIdxFilePath);

            if (queryNVDPHL) {
                // Load all object indexes into memory
                std::unordered_map<std::string,NVDQuadtree> fileNameToObjIdx;
                for (std::size_t j = 0; j < numSets; ++j) {
                    std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
                    std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(objIdxFilenamePrefix,constants::OBJ_IDX_NVD,parameterMap);
                    NVDQuadtree objIdx = serialization::getIndexFromBinaryFile<NVDQuadtree>(objIdxFilePath);
                    fileNameToObjIdx[objIdxFilenamePrefix] = objIdx;
                }
                this->runNVDPHL(graph,phl,alt,singleQuerySet,fileNameToObjIdx,maxStatsOutputFile,"single_query_set",aggMax,specialFields);
                this->runNVDPHL(graph,phl,alt,singleQuerySet,fileNameToObjIdx,sumStatsOutputFile,"single_query_set",aggSum,specialFields);
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
                
                this->runCOLTPHL(graph,phl,alt,sltree,singleQuerySet,fileNameToObjIdx,maxStatsOutputFile,"single_query_set",aggMax,specialFields);
                this->runCOLTPHL(graph,phl,alt,sltree,singleQuerySet,fileNameToObjIdx,sumStatsOutputFile,"single_query_set",aggSum,specialFields);
            }
        }
    }
}

void LMKAggregateExperimentsCommand::runIERPHL(Graph& graph, PrunedHighwayLabeling& phl, std::vector<std::vector<AknnQuery>>& querySets, std::unordered_map<std::string,PackedRtree>& objectIndexes, 
                                      std::string statsOutputFile, std::string desc, AggregateFunction<EdgeWeight>& agg, AggregateFunction<EuclideanDistance>& euclidAgg, std::vector<std::string> specialFields) 
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
                ier.getAkNNsByPHL(objectIndexes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertices,kNNs,kNNDistances,graph,phl,agg,euclidAgg);
                kNNs.clear();
                kNNDistances.clear();
            }
        } else {
            for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
                ier.getAkNNsByPHLTravelTimes(objectIndexes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertices,kNNs,kNNDistances,graph,phl,agg,euclidAgg);
                kNNs.clear();
                kNNDistances.clear();
            }
        }
    }
#endif

    for (std::size_t i = 0; i < querySets.size(); ++i) {
        totalQueries = querySets[i].size();
#if defined(COLLECT_STATISTICS)
        // Statistics only (rn_obj_search runs the timed loops below): the queries of the list run in parallel, each thread with its
        // own IER and statistics, which are then summed (integer counts, so the totals do not depend on the order)
        knnStats.clear();
        sw.reset();
        sw.start();
        #pragma omp parallel
        {
            IER threadIER;
            Statistics threadStats;
            std::vector<VertexID> threadKNNs;
            std::vector<EdgeWeight> threadKNNDistances;
            #pragma omp for schedule(dynamic)
            for (std::size_t j = 0; j < querySets[i].size(); ++j) {
                PackedRtree& rtree = objectIndexes.at(querySets[i][j].objSetFilenamePrefix); // at() because operator[] is not thread-safe
                if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
                    threadIER.getAkNNsByPHL(rtree,querySets[i][j].k,querySets[i][j].queryVertices,threadKNNs,threadKNNDistances,graph,phl,agg,euclidAgg);
                } else {
                    threadIER.getAkNNsByPHLTravelTimes(rtree,querySets[i][j].k,querySets[i][j].queryVertices,threadKNNs,threadKNNDistances,graph,phl,agg,euclidAgg);
                }
                threadKNNs.clear();
                threadKNNDistances.clear();
                threadStats.mergeStatistics(threadIER.stats);
            }
            #pragma omp critical
            knnStats.mergeStatistics(threadStats);
        }
        sw.stop();
#else
        if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
            sw.reset();
            sw.start();
            for (std::size_t j = 0; j < querySets[i].size(); ++j) {
                ier.getAkNNsByPHL(objectIndexes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertices,kNNs,kNNDistances,graph,phl,agg,euclidAgg);
                kNNs.clear();
                kNNDistances.clear();
            }
            sw.stop();
        } else {
            sw.reset();
            sw.start();
            for (std::size_t j = 0; j < querySets[i].size(); ++j) {
                ier.getAkNNsByPHLTravelTimes(objectIndexes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertices,kNNs,kNNDistances,graph,phl,agg,euclidAgg);
                kNNs.clear();
                kNNDistances.clear();
            }
            sw.stop();
        }
#endif
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;  
#if defined(COLLECT_STATISTICS)
        knnStats.normalizeStatistics(totalQueries);
#endif
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::IER_PHL_KNN_QUERY,
                            querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].queryVertices.size(),querySets[i][0].objSetSize,kNNs,kNNDistances);
        std::vector<std::string> moreFields;
        moreFields.push_back(std::to_string(querySets[i][0].querySubnetSize));
        results.setAdditionalFields(moreFields);
        results.setAdditionalFields(specialFields);
        results.addSupplementaryFields("agg_func",agg.getPrintString());
#if defined(COLLECT_STATISTICS)
        knnStats.populateTupleFields(results,0);
#endif
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "IER-PHL AkNN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKAggregateExperimentsCommand::runNVDPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, std::vector<std::vector<AknnQuery>>& querySets, std::unordered_map<std::string,NVDQuadtree>& objectIndexes, 
                                      std::string statsOutputFile, std::string desc, AggregateFunction<EdgeWeight>& agg, std::vector<std::string> specialFields) 
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
            ilbr.getAkNNByNVDQuadtreePHLByConcExps(alt,objectIndexes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertices,kNNs,kNNDistances,phl,graph,agg);
            kNNs.clear();
            kNNDistances.clear();
        }
    }
#endif

    for (std::size_t i = 0; i < querySets.size(); ++i) {
        totalQueries = querySets[i].size();
#if defined(COLLECT_STATISTICS)
        // Statistics only (rn_obj_search runs the timed loop below): the queries of the list run in parallel, each thread with its
        // own ILBR and statistics, which are then summed (integer counts, so the totals do not depend on the order)
        knnStats.clear();
        sw.reset();
        sw.start();
        #pragma omp parallel
        {
            ILBR threadILBR;
            Statistics threadStats;
            std::vector<VertexID> threadKNNs;
            std::vector<EdgeWeight> threadKNNDistances;
            #pragma omp for schedule(dynamic)
            for (std::size_t j = 0; j < querySets[i].size(); ++j) {
                NVDQuadtree& nvd = objectIndexes.at(querySets[i][j].objSetFilenamePrefix); // at() because operator[] is not thread-safe
                threadILBR.getAkNNByNVDQuadtreePHLByConcExps(alt,nvd,querySets[i][j].k,querySets[i][j].queryVertices,threadKNNs,threadKNNDistances,phl,graph,agg);
                threadKNNs.clear();
                threadKNNDistances.clear();
                threadStats.mergeStatistics(threadILBR.stats);
            }
            #pragma omp critical
            knnStats.mergeStatistics(threadStats);
        }
        sw.stop();
#else
        sw.reset();
        sw.start();
        for (std::size_t j = 0; j < querySets[i].size(); ++j) {
            ilbr.getAkNNByNVDQuadtreePHLByConcExps(alt,objectIndexes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertices,kNNs,kNNDistances,phl,graph,agg);
            kNNs.clear();
            kNNDistances.clear();
        }
        sw.stop();
#endif
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;  
#if defined(COLLECT_STATISTICS)
        knnStats.normalizeStatistics(totalQueries);
#endif
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::NVD_ILBR_PHL_KNN_QUERY,
                            querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].queryVertices.size(),querySets[i][0].objSetSize,kNNs,kNNDistances);
        std::vector<std::string> moreFields;
        moreFields.push_back(std::to_string(querySets[i][0].querySubnetSize));
        results.setAdditionalFields(moreFields);
        results.setAdditionalFields(specialFields);
        results.addSupplementaryFields("agg_func",agg.getPrintString());
#if defined(COLLECT_STATISTICS)
        knnStats.populateTupleFields(results,0);
#endif
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "NVD-PHL AkNN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKAggregateExperimentsCommand::runCOLTPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, SULTree& sltree, std::vector<std::vector<AknnQuery>>& querySets, std::unordered_map<std::string,COLTree>& objectIndexes, 
                                                std::string statsOutputFile, std::string desc, AggregateFunction<EdgeWeight>& agg, std::vector<std::string> specialFields, std::vector<std::pair<std::string,std::string>> supplementaryFields)
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
            ilbr.getAkNNByCOLTPHL(alt,objectIndexes[querySets[0][j].objSetFilenamePrefix],sltree,querySets[0][j].k,querySets[0][j].queryVertices,kNNs,kNNDistances,phl,agg);
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
            ilbr.getAkNNByCOLTPHL(alt,objectIndexes[querySets[i][j].objSetFilenamePrefix],sltree,querySets[i][j].k,querySets[i][j].queryVertices,kNNs,kNNDistances,phl,agg);
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
                            querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].queryVertices.size(),querySets[i][0].objSetSize,kNNs,kNNDistances);
        std::vector<std::string> moreFields;
        moreFields.push_back(std::to_string(querySets[i][0].querySubnetSize));
        results.setAdditionalFields(moreFields);
        results.setAdditionalFields(specialFields);
        results.addSupplementaryFields("agg_func",agg.getPrintString());
        results.addSupplementaryFields("listsize",std::to_string(objectIndexes[querySets[i][0].objSetFilenamePrefix].getMaxOLSize()));
        for (std::size_t f = 0; f < supplementaryFields.size(); ++f) {
            results.addSupplementaryFields(supplementaryFields[f].first,supplementaryFields[f].second);
        }
#if defined(COLLECT_STATISTICS)
        knnStats.populateTupleFields(results,0);
#endif
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "COLT-PHL AkNN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKAggregateExperimentsCommand::runCOLTParameterTestQueries(std::string bgrFileName, std::string parameters, std::string filePathPrefix, 
                                                                 std::size_t numSets, std::string strObjectDensities, std::string objTypes, std::string objVariable, 
                                                                 std::string strFanouts, std::string strNumLandmarks, std::string strLandmarkTypes, std::string strLambdas, 
                                                                 int minLambda, AggregateFunction<EdgeWeight>& defaultAgg)
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

    // Load supporting road network index files
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);

    // Extract AkNN specific information from input parameters
    int defaultNumQueryVertices = std::stoi(parameterMap["par_aknn_qv_def_setsize"]);
    if (defaultNumQueryVertices <= 0) {
        std::cerr << "Invalid default query vertex set size " << defaultNumQueryVertices << "!\n";
        std::exit(1);    
    }
    
    int defaultQueryVertexSubnetSize = std::stoi(parameterMap["par_aknn_qv_def_setsubnetsize"]);
    int maxNumQueryVertices = std::stoi(parameterMap["par_aknn_qv_max_setsize"]);
    int numQueryVertexSets = std::stoi(parameterMap["par_aknn_num_qv_sets"]);

    if (defaultQueryVertexSubnetSize <= 0 || defaultQueryVertexSubnetSize > 100) {
        std::cerr << "Invalid default query vertex set subnet size " << defaultQueryVertexSubnetSize << "!\n";
        std::exit(1);
    }
    if (maxNumQueryVertices <= 1) {
        std::cerr << "Invalid max query vertex set size " << maxNumQueryVertices << "!\n";
        std::exit(1);
    }
    if (numQueryVertexSets <= 1) {
        std::cerr << "Invalid number of query vertex sets " << numQueryVertexSets << "!\n";
        std::exit(1);
    }

    std::string vectorType;
    unsigned int numElementsPerVec;
    
    // Load default subnetwork size query vertex sets
    std::string defSubnetQueryVertexFilePath = filePathPrefix + "/query_sets/" + utility::constructQuerySetFileName(graph.getNetworkName(),"subnet",numQueryVertexSets,maxNumQueryVertices,defaultQueryVertexSubnetSize);
    std::vector<std::vector<VertexID>> queryVertexSets = utility::getVectorSetFromFile(defSubnetQueryVertexFilePath,vectorType,numElementsPerVec);
    if (queryVertexSets.size() == 0) {
        std::cerr << "No query points were provided!\n";
        exit(1);    
    }
    
    // Create query vertex set vector for default set size
    std::vector<std::vector<VertexID>> defaultSizeQueryVertexSetsVec;
    for (std::size_t i = 0; i < queryVertexSets.size(); ++i) {
        std::vector<VertexID> tempVec;
        for (std::size_t j = 0; j < static_cast<std::size_t>(defaultNumQueryVertices); ++j) {
            tempVec.push_back(queryVertexSets[i][j]);
        }
        defaultSizeQueryVertexSetsVec.push_back(tempVec);
    }

    std::string phlIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PHL_CMD,{},{});
    PrunedHighwayLabeling phl;
    phl.LoadLabel(phlIdxFilePath.c_str());

    const std::vector<std::string> altKeys = {"landmarks","type"};
    const std::vector<std::string> altValues = {parameterMap["alt_numlandmarks"],parameterMap["alt_landmarktype"]};
    std::string altIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,altKeys,altValues);
    ALT alt = serialization::getIndexFromBinaryFile<ALT>(altIdxFilePath);

    // Query lists as runQueries2 builds them: the default queries (default k and density) and one list per density for the density experiments
    std::vector<std::vector<AknnQuery>> defaultQuerySets(1), varyingDensityQuerySets;
    for (std::size_t j = 0; j < numSets; ++j) {
        std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVariable,j);
        for (std::size_t k = 0; k < defaultSizeQueryVertexSetsVec.size(); ++k) {
            defaultQuerySets[0].push_back(AknnQuery(defaultSizeQueryVertexSetsVec[k],defaultK,defaultDensity,defaultQueryVertexSubnetSize,defaultObjType,defaultObjVariable,j,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*defaultDensity)));
        }
    }
    for (std::size_t i = 0; i < densities.size(); ++i) {
        std::vector<AknnQuery> densityQuerySet;
        for (std::size_t j = 0; j < numSets; ++j) {
            std::string objIdxFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,densities[i],defaultObjVariable,j);
            for (std::size_t k = 0; k < defaultSizeQueryVertexSetsVec.size(); ++k) {
                densityQuerySet.push_back(AknnQuery(defaultSizeQueryVertexSetsVec[k],defaultK,densities[i],defaultQueryVertexSubnetSize,defaultObjType,defaultObjVariable,j,objIdxFilenamePrefix,static_cast<int>(graph.getNumVertices()*densities[i])));
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
            this->runCOLTPHL(graph,phl,alt,sltree,defaultQuerySets,fileNameToObjIdx,statsOutputFile,"fanout vs lambda",defaultAgg,specialFields,testedParameters);
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
            this->runCOLTPHL(graph,phl,alt,sltree,defaultQuerySets,fileNameToObjIdx,statsOutputFile,"num landmarks vs lambda",defaultAgg,specialFields,testedParameters);
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
            this->runCOLTPHL(graph,phl,alt,sltree,defaultQuerySets,fileNameToObjIdx,statsOutputFile,"num landmarks vs landmark type",defaultAgg,specialFields,testedParameters);
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
            std::vector<std::vector<AknnQuery>> densityQuerySets(1,varyingDensityQuerySets[j]);
            this->runCOLTPHL(graph,phl,alt,sltree,densityQuerySets,fileNameToObjIdx,statsOutputFile,"density vs num landmarks",defaultAgg,specialFields,testedParameters);
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
            std::vector<std::vector<AknnQuery>> densityQuerySets(1,varyingDensityQuerySets[j]);
            this->runCOLTPHL(graph,phl,alt,sltree,densityQuerySets,fileNameToObjIdx,statsOutputFile,"density vs landmark type",defaultAgg,specialFields,testedParameters);
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
                std::vector<std::vector<AknnQuery>> densityQuerySets(1,varyingDensityQuerySets[j]);
                this->runCOLTPHL(graph,phl,alt,sltree,densityQuerySets,fileNameToObjIdx,statsOutputFile,"density vs lambda",defaultAgg,specialFields,testedParameters);
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
            std::vector<std::vector<AknnQuery>> densityQuerySets(1,varyingDensityQuerySets[j]);
            this->runCOLTPHL(graph,phl,alt,sltree,densityQuerySets,fileNameToObjIdx,statsOutputFile,"density vs fanout",defaultAgg,specialFields,testedParameters);
        }
    }
}
