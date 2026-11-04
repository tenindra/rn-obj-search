// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "LMKFarthestExperimentsCommand.h"

#include "../processing/INE.h"
#include "../processing/IUBR.h"
#include "../processing/SetGenerator.h"
#include "../tuple/IndexTuple.h"
#include "../tuple/ObjectIndexTuple.h"
#include "../tuple/KnnQueryTuple.h"
#include "../tuple/QuerySet.h"
#include "../common.h"
#include "../utility/StopWatch.h"
#include "../utility/Statistics.h"
#include "../utility/utility.h"
#include "../utility/serialization.h"

#include <cstdio>

void LMKFarthestExperimentsCommand::execute(int argc, char* argv[])
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
    
    /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "e:g:p:f:s:n:d:t:q:k:v:r:u:w:x:y:z:")) != -1) {
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
                minLambda = std::stoul(optarg);;
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

void LMKFarthestExperimentsCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::LMK_FARTHEST_EXPERIMENTS_CMD + " -e <experimental step>\n\n"
              << "Steps:\n"
              << utility::getFormattedUsageString(this->RUN_KNN_QUERIES,"Run all landmark FkN experiments") + "\n"
              << utility::getFormattedUsageString(this->RUN_KNN_SINGLE_QUERIES,"Run single landmark FkN experiments") + "\n"
              << utility::getFormattedUsageString(this->RUN_PARAMETER_TEST_KNN_QUERIES,"Run FkN queries for parameter experiments") + "\n";
}

void LMKFarthestExperimentsCommand::showPhaseUsage(std::string method, std::string programName)
{
    if (method == this->RUN_KNN_QUERIES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_FARTHEST_EXPERIMENTS_CMD
                  << " -e " + this->RUN_KNN_QUERIES + " -g <binary graph file>\n"
                  << "-q <query vertex file> -k <k values> -p <parameters> -n <num sets> -d <list of object set densities>\n"
                  << "-t <list of object types>  -v <list of some object variable> -f <index output file path prefix>\n"
                  << "[-r <rw object set file names>]\n";
    } else if (method == this->RUN_KNN_SINGLE_QUERIES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_FARTHEST_EXPERIMENTS_CMD
                  << " -e " + this->RUN_KNN_SINGLE_QUERIES + " -g <binary graph file>\n"
                  << "-q <query vertex file> -k <k values> -p <parameters> -n <num sets> -d <list of object set densities>\n"
                  << "-t <list of object types>  -v <list of some object variable> -f <index output file path prefix>\n"
                  << "-s <stats output file>\n";
    } else if (method == this->RUN_PARAMETER_TEST_KNN_QUERIES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::LMK_FARTHEST_EXPERIMENTS_CMD
                  << " -e " + this->RUN_PARAMETER_TEST_KNN_QUERIES + " -g <binary graph file> -p <parameters> -f <index output file path prefix>\n"
                  << "-q <query vertex file> -u <minimum lambda> -w <fanouts> -x <num landmarks> -y <landmark types> -z <lambdas>\n"
                  << "-n <num sets> -d <object densities> -t <object types> -v <object variables>";
    } else {
        std::cerr << "Invalid experiment phase!" << std::endl;
        this->showCommandUsage(programName);
    }
}

std::unordered_map<std::string, std::string> LMKFarthestExperimentsCommand::getParameters(std::string parameters)
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

void LMKFarthestExperimentsCommand::writeBlankLines(std::string fileName)
{
    if (fileName == "") {
        return; // A switched-off query list has no stats file
    }
    std::ofstream fs;
    fs.open(fileName, std::ios_base::app | std::ios_base::out);
    fs << "\n\n"; //  GNU Plot requires double line separation between data points for different plot lines
    fs.close();
}

void LMKFarthestExperimentsCommand::runQueries2(std::string bgrFileName, std::string queryVertexFile, std::string kValues, std::string parameters, std::size_t numSets, 
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
    unsigned int defaultK = std::stoul(parameterMap["fkn_def_k"]);
    double defaultDensity = std::stod(parameterMap["fkn_def_density"]);
    std::string defaultObjType = parameterMap["fkn_def_type"];
    int defaultObjVar = std::stoi(parameterMap["fkn_def_var"]);
    if (defaultObjType == "") {
        std::cerr << "Invalid default object set type of " << defaultObjType << "!\n";
        std::exit(1);    
    }
    
    /*
     * Create Query Sets for Each Experiment
     */
    
    // Each query list runs unless its switch in the parameters is off (utility::isQueryListOn). A switched-off list is not built, and
    // its stats file name is left empty, so that neither its runs nor the blank separator lines write to its stats file
    bool runVaryingK = utility::isQueryListOn(parameterMap,"fkn_varyingk");
    bool runVaryingDensity = utility::isQueryListOn(parameterMap,"fkn_varyingdensity");
    
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
    
    if (parameterMap.find("fkn_varyingk_stats_file") == parameterMap.end() || parameterMap["fkn_varyingk_stats_file"] == "") {
        std::cerr << "Varying k stats file not provided in parameters";
        std::exit(1);
    }
    std::string varyingKStatsFile = runVaryingK ? filePathPrefix + "/stats/" + parameterMap["fkn_varyingk_stats_file"] : "";
    
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
    
    if (parameterMap.find("fkn_varyingdensity_stats_file") == parameterMap.end() || parameterMap["fkn_varyingdensity_stats_file"] == "") {
        std::cerr << "Varying density stats file not provided in parameters";
        std::exit(1);
    }
    std::string varyingDensityStatsFile = runVaryingDensity ? filePathPrefix + "/stats/" + parameterMap["fkn_varyingdensity_stats_file"] : "";
    
    // Varying Real-World Object Set Experiment
    std::vector<std::vector<KnnQuery>> varyingRWObjSetQuerySets;

    std::vector<std::string> rwObjSetNamesVec = utility::splitByDelim(rwObjSetNames,',');
    bool runRealWorld = parameterMap["fkn_rw"] == "1";
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
        if (parameterMap.find("fkn_rw_qv_file") == parameterMap.end() || parameterMap["fkn_rw_qv_file"] == "") {
            std::cerr << "Real-world query vertex file not provided in parameters!\n";
            std::exit(1);
        }
        std::vector<VertexID> rwQueryVertices = utility::getPointSetFromFile(filePathPrefix + "/query_sets/" + parameterMap["fkn_rw_qv_file"]);
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

    if (parameterMap.find("fkn_varyingrwset_stats_file") == parameterMap.end() || parameterMap["fkn_varyingrwset_stats_file"] == "") {
        std::cerr << "Varying real-world stats file not provided in parameters";
        std::exit(1);
    }
    std::string varyingRWObjSetStatsFile = filePathPrefix + "/stats/" + parameterMap["fkn_varyingrwset_stats_file"];
    
    bool queryINE = parameterMap["ine"] == "1";
    bool queryLLBPHL = parameterMap["llb_phl"] == "1";
    bool queryCOLTPHL = parameterMap["colt_phl"] == "1";
    
    std::string setType;
    double setDensity;
    int setSize, setVariable;

    if (queryINE) {
        INE ine;

        // Load all object sets into memory
        std::unordered_map<std::string,std::vector<bool>> fileNameToObjIdx;
        std::unordered_map<std::string,int> fileNameToObjSetSize;
        for (std::size_t i = 0; i < objDensitiesVec.size(); ++i) {
            for (std::size_t j = 0; j < numSets; ++j) {
                std::string objSetFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,objDensitiesVec[i],defaultObjVar,j);
                std::string objSetFilePath = filePathPrefix + "/obj_indexes/" + objSetFilenamePrefix + ".txt";
                std::vector<VertexID> objSetVertices = utility::getPointSetFromFile(objSetFilePath,setType,setDensity,setSize,setVariable);
                std::vector<bool> vertexObjectStatus = ine.buildVertexObjectStatusVector(graph,objSetVertices);
                // Do we need to check it hasn't already been inserted?
                fileNameToObjIdx[objSetFilenamePrefix] = vertexObjectStatus;
                fileNameToObjSetSize[objSetFilenamePrefix] = setSize;
            }
        }
        
        if (runRealWorld) {
            for (std::size_t i = 0; i < rwObjSetNamesVec.size(); ++i) {
                std::string objSetFilenamePrefix = graph.getNetworkName() + "_" + rwObjSetNamesVec[i];
                std::string objSetFilePath = filePathPrefix + "/obj_indexes/" + objSetFilenamePrefix + ".txt";
                std::vector<VertexID> objSetVertices = utility::getPointSetFromFile(objSetFilePath,setType,setDensity,setSize,setVariable);
                std::vector<bool> vertexObjectStatus = ine.buildVertexObjectStatusVector(graph,objSetVertices);
                fileNameToObjIdx[objSetFilenamePrefix] = vertexObjectStatus;
                fileNameToObjSetSize[objSetFilenamePrefix] = setSize;
            }
        }

        this->runINEQueries2(graph,varyingKQuerySets,fileNameToObjIdx,fileNameToObjSetSize,varyingKStatsFile,"varyingk");
        this->runINEQueries2(graph,varyingDensityQuerySets,fileNameToObjIdx,fileNameToObjSetSize,varyingDensityStatsFile,"varyingdensity");
        if (runRealWorld) {
            this->runINEQueries2(graph,varyingRWObjSetQuerySets,fileNameToObjIdx,fileNameToObjSetSize,varyingRWObjSetStatsFile,"varyingrwset");
        }

        this->writeBlankLines(varyingKStatsFile);
        this->writeBlankLines(varyingDensityStatsFile);
    }
    
    if (queryLLBPHL || queryCOLTPHL) {
        std::string phlIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PHL_CMD,{},{});
        PrunedHighwayLabeling phl;
        phl.LoadLabel(phlIdxFilePath.c_str());
        
        const std::vector<std::string> altKeys = {"landmarks","type"};
        const std::vector<std::string> altValues = {parameterMap["alt_numlandmarks"],parameterMap["alt_landmarktype"]};
        std::string altIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,altKeys,altValues);
        ALT alt = serialization::getIndexFromBinaryFile<ALT>(altIdxFilePath);
        
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

void LMKFarthestExperimentsCommand::runSingleQueries2(std::string bgrFileName, std::string queryVertexFile, std::string kValues, std::string parameters, 
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
    unsigned int defaultK = std::stoul(parameterMap["fkn_def_k"]);
    double defaultDensity = std::stod(parameterMap["fkn_def_density"]);
    std::string defaultObjType = parameterMap["fkn_def_type"];
    int defaultObjVar = std::stoi(parameterMap["fkn_def_var"]);
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
    bool queryLLBPHL = parameterMap["llb_phl"] == "1";
    bool queryCOLTPHL = parameterMap["colt_phl"] == "1";
    
    std::string setType;
    double setDensity;
    int setSize, setVariable;

    if (queryINE) {
        INE ine;

        // Load all object sets into memory
        std::unordered_map<std::string,std::vector<bool>> fileNameToObjIdx;
        std::unordered_map<std::string,int> fileNameToObjSetSize;
        for (std::size_t j = 0; j < numSets; ++j) {
            std::string objSetFilenamePrefix = utility::constructObjIdxFilenamePrefix(graph.getNetworkName(),defaultObjType,defaultDensity,defaultObjVar,j);
            std::string objSetFile = filePathPrefix + "/obj_indexes/" + objSetFilenamePrefix + ".txt";
            std::vector<VertexID> objSetVertices = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
            std::vector<bool> vertexObjectStatus = ine.buildVertexObjectStatusVector(graph,objSetVertices);
            // Do we need to check it hasn't already been inserted?
            fileNameToObjIdx[objSetFilenamePrefix] = vertexObjectStatus;
            fileNameToObjSetSize[objSetFilenamePrefix] = setSize;
        }

        this->runINEQueries2(graph,singleQuerySet,fileNameToObjIdx,fileNameToObjSetSize,statsOutputFile,"single_query_set");
    }
    
    if (queryLLBPHL || queryCOLTPHL) {
        std::string phlIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PHL_CMD,{},{});
        PrunedHighwayLabeling phl;
        phl.LoadLabel(phlIdxFilePath.c_str());

        const std::vector<std::string> altKeys = {"landmarks","type"};
        const std::vector<std::string> altValues = {parameterMap["alt_numlandmarks"],parameterMap["alt_landmarktype"]};
        std::string altIdxFilePath = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,altKeys,altValues);
        ALT alt = serialization::getIndexFromBinaryFile<ALT>(altIdxFilePath);
        
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

void LMKFarthestExperimentsCommand::runINEQueries2(Graph& graph, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,std::vector<bool>>& objectIndexes, 
                                                   std::unordered_map<std::string,int>& objectSetSizes, std::string statsOutputFile, std::string desc, std::vector<std::string> specialFields) 
{
    INE ine;

    StopWatch sw;
    double totalQueryTime;
    int totalQueries;
    std::vector<VertexID> FkNs;
    std::vector<EdgeWeight> FkNDistances;
    
    // Warm-up: the first WARM_UP_QUERY_PERCENTAGE percent of the first list (at least one query) runs untimed, so that the cold start of a fresh
    // process or of freshly loaded indexes stays out of the timings (skipped in rn_obj_search_stats, whose statistics do not depend on it)
#if !defined(COLLECT_STATISTICS)
    if (querySets.size() > 0) {
        std::size_t numWarmUpQueries = utility::getNumWarmUpQueries(querySets[0].size());
        for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
            ine.getFkNs(graph,objectIndexes[querySets[0][j].objSetFilenamePrefix],objectSetSizes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertex,FkNs,FkNDistances);
            FkNs.clear();
            FkNDistances.clear();
        }
    }
#endif

    for (std::size_t i = 0; i < querySets.size(); ++i) {
        totalQueries = querySets[i].size();
        sw.reset();
        sw.start();
        for (std::size_t j = 0; j < querySets[i].size(); ++j) {
            ine.getFkNs(graph,objectIndexes[querySets[i][j].objSetFilenamePrefix],objectSetSizes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertex,FkNs,FkNDistances);
            FkNs.clear();
            FkNDistances.clear();
        }
        sw.stop();
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;
        
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::INE_KNN_QUERY,
                              querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].var,querySets[i][0].objSetSize,FkNs,FkNDistances);
        results.setAdditionalFields(specialFields);
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "INE FkN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKFarthestExperimentsCommand::runLLBPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,std::vector<VertexID>>& objectIndexes, 
                                              std::string statsOutputFile, std::string desc, std::vector<std::string> specialFields) 
{
    IUBR iubr;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries;
    std::vector<VertexID> FkNs;
    std::vector<EdgeWeight> FkNDistances;

    Statistics knnStats;
    
    // Warm-up: the first WARM_UP_QUERY_PERCENTAGE percent of the first list (at least one query) runs untimed, so that the cold start of a fresh
    // process or of freshly loaded indexes stays out of the timings (skipped in rn_obj_search_stats, whose statistics do not depend on it)
#if !defined(COLLECT_STATISTICS)
    if (querySets.size() > 0) {
        std::size_t numWarmUpQueries = utility::getNumWarmUpQueries(querySets[0].size());
        for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
            iubr.getFkNsByPHLNaive(alt,objectIndexes[querySets[0][j].objSetFilenamePrefix],querySets[0][j].k,querySets[0][j].queryVertex,FkNs,FkNDistances,phl);
            FkNs.clear();
            FkNDistances.clear();
        }
    }
#endif

    for (std::size_t i = 0; i < querySets.size(); ++i) {
        totalQueries = querySets[i].size();
#if defined(COLLECT_STATISTICS)
        // Statistics only (rn_obj_search runs the timed loop below): the queries of the list run in parallel, each thread with its
        // own IUBR and statistics, which are then summed (integer counts, so the totals do not depend on the order)
        knnStats.clear();
        sw.reset();
        sw.start();
        #pragma omp parallel
        {
            IUBR threadIUBR;
            Statistics threadStats;
            std::vector<VertexID> threadFkNs;
            std::vector<EdgeWeight> threadFkNDistances;
            #pragma omp for schedule(dynamic)
            for (std::size_t j = 0; j < querySets[i].size(); ++j) {
                std::vector<VertexID>& objects = objectIndexes.at(querySets[i][j].objSetFilenamePrefix); // at() because operator[] is not thread-safe
                threadIUBR.getFkNsByPHLNaive(alt,objects,querySets[i][j].k,querySets[i][j].queryVertex,threadFkNs,threadFkNDistances,phl);
                threadFkNs.clear();
                threadFkNDistances.clear();
                threadStats.mergeStatistics(threadIUBR.stats);
            }
            #pragma omp critical
            knnStats.mergeStatistics(threadStats);
        }
        sw.stop();
#else
        sw.reset();
        sw.start();
        for (std::size_t j = 0; j < querySets[i].size(); ++j) {
            iubr.getFkNsByPHLNaive(alt,objectIndexes[querySets[i][j].objSetFilenamePrefix],querySets[i][j].k,querySets[i][j].queryVertex,FkNs,FkNDistances,phl);
            FkNs.clear();
            FkNDistances.clear();
        }
        sw.stop();
#endif
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;  
#if defined(COLLECT_STATISTICS)
        knnStats.normalizeStatistics(totalQueries);
#endif
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::LLB_ILBR_PHL_KNN_QUERY,
                            querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].var,querySets[i][0].objSetSize,FkNs,FkNDistances);
        results.setAdditionalFields(specialFields);
#if defined(COLLECT_STATISTICS)
        knnStats.populateTupleFields(results,0);
#endif
        this->outputCommandStats(statsOutputFile,results.getTupleString());
    }
    if (querySets.size() > 0) { // Empty for a switched-off query list
        std::cout << "LLB-PHL FkN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKFarthestExperimentsCommand::runCOLTPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, SULTree& sltree, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,COLTree>& objectIndexes, 
                                       std::string statsOutputFile, std::string desc, std::vector<std::string> specialFields, std::vector<std::pair<std::string,std::string>> supplementaryFields)
{
    IUBR iubr;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries;
    std::vector<VertexID> FkNs;
    std::vector<EdgeWeight> FkNDistances;

    Statistics knnStats;
    
    // Warm-up: the first WARM_UP_QUERY_PERCENTAGE percent of the first list (at least one query) runs untimed, so that the cold start of a fresh
    // process or of freshly loaded indexes stays out of the timings (skipped in rn_obj_search_stats, whose statistics do not depend on it)
#if !defined(COLLECT_STATISTICS)
    if (querySets.size() > 0) {
        std::size_t numWarmUpQueries = utility::getNumWarmUpQueries(querySets[0].size());
        for (std::size_t j = 0; j < numWarmUpQueries && j < querySets[0].size(); ++j) {
            iubr.getFkNsByCOLTPHL(alt,objectIndexes[querySets[0][j].objSetFilenamePrefix],sltree,querySets[0][j].k,querySets[0][j].queryVertex,FkNs,FkNDistances,phl);
            FkNs.clear();
            FkNDistances.clear();
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
            iubr.getFkNsByCOLTPHL(alt,objectIndexes[querySets[i][j].objSetFilenamePrefix],sltree,querySets[i][j].k,querySets[i][j].queryVertex,FkNs,FkNDistances,phl);
            FkNs.clear();
            FkNDistances.clear();
#if defined(COLLECT_STATISTICS)
            knnStats.mergeStatistics(iubr.stats);
#endif
        }
        sw.stop();
        totalQueryTime = sw.getTimeUs();
        double queryTimeUs = totalQueryTime/totalQueries;  
#if defined(COLLECT_STATISTICS)
        knnStats.normalizeStatistics(totalQueries);
#endif
        KnnQueryTuple results(graph.getNetworkName(),graph.getNumVertices(),graph.getNumEdges(),totalQueries,constants::COLTREE_ILBR_PHL_KNN_QUERY,
                            querySets[i][0].k,queryTimeUs,querySets[i][0].type,querySets[i][0].density,querySets[i][0].var,querySets[i][0].objSetSize,FkNs,FkNDistances);
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
        std::cout << "COLT-PHL FkN queries successfully executed for " << graph.getNetworkName() << " with " << desc << std::endl;
    }
}

void LMKFarthestExperimentsCommand::runCOLTParameterTestQueries(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string queryVertexFile, 
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
    int defaultK = std::stoi(parameterMap["par_def_fkn_def_k"]);
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
