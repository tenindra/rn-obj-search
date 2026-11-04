// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _LMKEXPERIMENTSCOMMAND_H
#define _LMKEXPERIMENTSCOMMAND_H

#include "Command.h"
#include "../processing/Graph.h"
#include "../tuple/QuerySet.h"
#include "../tuple/IndexTuple.h"
#include "../processing/pruned_highway_labeling.h"
#include "../processing/ALT.h"
#include "../processing/Gtree.h"
#include "../processing/PackedRtree.h"
#include "../processing/NVDQuadtree.h"
#include "../processing/SULTree.h"

#include <unordered_map>

class LMKExperimentsCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);
        void showPhaseUsage(std::string phase, std::string programName);

        // Experimental Phases for Command
        std::string const BUILD_RN_INDEXES = "indexes";
        std::string const BUILD_PARTITIONS = "partitions";
        std::string const BUILD_OBJ_INDEXES = "objects";
        std::string const BUILD_RW_OBJ_INDEXES = "rw_obj_idxs";
        std::string const RUN_KNN_QUERIES = "query_exp";
        std::string const RUN_KNN_SINGLE_QUERIES = "squery_exp";
        std::string const RUN_SDL_EXPERIMENT = "sdl_exp";
        std::string const BUILD_PARAMETER_TEST_RN_INDEXES = "par_indexes";
        std::string const BUILD_PARAMETER_TEST_OBJ_INDEXES = "par_obj_idxs";
        std::string const RUN_PARAMETER_TEST_KNN_QUERIES = "par_query_exp";
        
    private:
        void buildIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile);
        void buildExternalIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile);
        void buildPartitions(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile);
        void runSDLExperiment(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile);
        void computeSpillFactors(SULTree& sltree, double& avg, double& min, double& max);
        void addSpillFactorFields(SULTree& sltree, IndexTuple& stats);
        void buildObjectIndexes(std::string bgrFileName, std::string parameters, std::size_t numSets, std::string objDensities, 
                                std::string objTypes, std::string objVariable, std::string filePathPrefix, std::string statsOutputFile, std::string sourceBgrFilePath = "");
        std::unordered_map<std::string,std::string> getParameters(std::string parameters);
        void writeBlankLines(std::string fileName);

        void buildGtree(Graph& graph, int fanout, std::size_t maxLeafSize, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields);
        void buildPHL(std::string networkName, int numVertices, int numEdges, std::string dataOutputFile, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields);
        void buildSULTree(Graph& graph, const HierarchicalPartition& partition, int numNodeLandmarks, std::string landmarkTypeStr, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields);
        void buildPartition(Graph& graph, int fanout, std::size_t maxLeafSize, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields, HierarchicalPartition& partition);
        void buildALT(Graph& graph, int numLandmarks, std::string selectionMethod, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields);

        void generateObjectSets(Graph& graph, std::size_t numSets, std::vector<double> objDensities, 
                                std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile,
                                unsigned int setSeed, std::string sourceBgrFilePath = "");
        void checkObjectSetHeader(std::string objSetFile, std::string objType, double objDensity, std::string setType, double setDensity);
        void buildOccurenceLists(std::string gtreeIdxFile, std::size_t numSets, std::vector<double> objDensities, 
                                 std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                                 std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string> specialFields = {});
        void buildRtrees(Graph& graph, int branchFactor, std::size_t numSets, std::vector<double> objDensities, 
                         std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string> specialFields = {});
        void buildObjectLists(Graph& graph, std::string altIdxFile, std::size_t numSets, std::vector<double> objDensities, std::vector<std::string> objTypes, 
                              std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                              std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string> specialFields = {});
        void buildNVD(Graph& graph, std::size_t numSets, std::vector<double> objDensities, std::vector<std::string> objTypes, 
                      std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string> specialFields = {});
        void buildCOLT(Graph& graph, std::string sltreeIdxFile, int objListSize, std::size_t numSets, std::vector<double> objDensities, 
                       std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                       std::unordered_map<std::string,std::string>& parameterMap, std::vector<std::string> specialFields = {});
        void buildRealWorldObjectIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix,
                                         std::string statsOutputFile, std::string rwObjSetNames);

        void runQueries2(std::string bgrFileName, std::string queryVertexFile, std::string kValues, std::string parameters, std::size_t numSets, 
                         std::string objDensities, std::string objTypes, std::string objVariable, std::string rwObjSetNames, std::string filePathPrefix);
        void runSingleQueries2(std::string bgrFileName, std::string queryVertexFile, std::string kValues, 
                               std::string parameters, std::size_t numSets, std::string objDensities, std::string objTypes, std::string objVariable, 
                               std::string filePathPrefix, std::string statsOutputFile);
        void runINEQueries2(Graph& graph, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,std::vector<bool>>& objectIndexes, 
                            std::string statsOutputFile, std::string desc, std::vector<std::string> specialFields = {});
        void runGtreeQueries2(Graph& graph, Gtree& gtree, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,OccurenceList>& objectIndexes, 
                            std::string statsOutputFile, std::string desc, 
                            std::vector<std::string> specialFields = {});
        void runIERPHL(Graph& graph, PrunedHighwayLabeling& phl, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,PackedRtree>& objectIndexes, 
                       std::string statsOutputFile, std::string desc, 
                       std::vector<std::string> specialFields = {});
        void runOLPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,ObjectList>& objectIndexes, 
                      std::string statsOutputFile, std::string desc, 
                      std::vector<std::string> specialFields = {});
        void runNVDPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,NVDQuadtree>& objectIndexes, 
                       std::string statsOutputFile, std::string desc, 
                       std::vector<std::string> specialFields = {});
        void runLLBPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,std::vector<VertexID>>& objectIndexes, 
                       std::string statsOutputFile, std::string desc, 
                       std::vector<std::string> specialFields = {});
        void runCOLTPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, SULTree& sltree, std::vector<std::vector<KnnQuery>>& querySets, std::unordered_map<std::string,COLTree>& objectIndexes, 
                        std::string statsOutputFile, std::string desc, 
                        std::vector<std::string> specialFields = {}, std::vector<std::pair<std::string,std::string>> supplementaryFields = {});

        // SULTree and COLTree parameter testing functionality
        void buildSULTreeParameterIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string strFanouts,
                                         std::string strNumLandmarks, std::string strLandmarkTypes, int minLambda);
        void buildCOLTParameterIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::size_t numSets, 
                                       std::string strObjectDensities, std::string objTypes, std::string objVariable, std::string strFanouts,
                                       std::string strNumLandmarks, std::string strLandmarkTypes, std::string strLambdas, int minLambda);
        void runCOLTParameterTestQueries(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string queryVertexFile, 
                                         std::size_t numSets, std::string strObjectDensities, std::string objTypes, std::string objVariable, 
                                         std::string strFanouts,std::string strNumLandmarks, std::string strLandmarkTypes, std::string strLambdas, 
                                         int minLambda);
};

#endif // _LMKEXPERIMENTSCOMMAND_H
