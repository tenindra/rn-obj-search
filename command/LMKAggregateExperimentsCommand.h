// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _LMKAGGEXPERIMENTSCOMMAND_H
#define _LMKAGGEXPERIMENTSCOMMAND_H

#include "Command.h"
#include "../processing/Graph.h"
#include "../tuple/QuerySet.h"
#include "../processing/pruned_highway_labeling.h"
#include "../processing/ALT.h"
#include "../processing/Gtree.h"
#include "../processing/PackedRtree.h"
#include "../processing/NVDQuadtree.h"
#include "../processing/SULTree.h"
#include "../utility/AggregateFunction.h"

#include <unordered_map>

class LMKAggregateExperimentsCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);
        void showPhaseUsage(std::string phase, std::string programName);

        // Experimental Phases for Command
        std::string const RUN_KNN_QUERIES = "query_exp";
        std::string const RUN_KNN_SINGLE_QUERIES = "squery_exp";
        std::string const RUN_PARAMETER_TEST_KNN_QUERIES = "par_query_exp";
        
    private:
        std::unordered_map<std::string,std::string> getParameters(std::string parameters);
        void writeBlankLines(std::string fileName);
        void runQueries2(std::string bgrFileName, std::string kValues, std::string parameters, std::size_t numSets, 
                         std::string objDensities, std::string objTypes, std::string objVariable, std::string rwObjSetNames, std::string queryVertexSetSizes, std::string queryVertexSetSubnetSizes, std::string filePathPrefix);
        void runSingleQueries2(std::string bgrFileName, std::string kValues, 
                               std::string parameters, std::size_t numSets, std::string objDensities, std::string objTypes, std::string objVariable, 
                               std::string filePathPrefix, std::vector<std::string>& statsOutputFiles);
        void runIERPHL(Graph& graph, PrunedHighwayLabeling& phl, std::vector<std::vector<AknnQuery>>& querySets, std::unordered_map<std::string,PackedRtree>& objectIndexes, 
                       std::string statsOutputFile, std::string desc, AggregateFunction<EdgeWeight>& agg, 
                       AggregateFunction<EuclideanDistance>& euclidAgg, std::vector<std::string> specialFields = {});
        void runNVDPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, std::vector<std::vector<AknnQuery>>& querySets, std::unordered_map<std::string,NVDQuadtree>& objectIndexes, 
                       std::string statsOutputFile, std::string desc, AggregateFunction<EdgeWeight>& agg, 
                       std::vector<std::string> specialFields = {});
        void runCOLTPHL(Graph& graph, PrunedHighwayLabeling& phl, ALT& alt, SULTree& sltree, std::vector<std::vector<AknnQuery>>& querySets, std::unordered_map<std::string,COLTree>& objectIndexes, 
                        std::string statsOutputFile, std::string desc, AggregateFunction<EdgeWeight>& agg,
                        std::vector<std::string> specialFields = {}, std::vector<std::pair<std::string,std::string>> supplementaryFields = {});

        // SULTree and COLTree parameter testing functionality
        void runCOLTParameterTestQueries(std::string bgrFileName, std::string parameters, std::string filePathPrefix, 
                                         std::size_t numSets, std::string strObjectDensities, std::string objTypes, std::string objVariable, 
                                         std::string strFanouts,std::string strNumLandmarks, std::string strLandmarkTypes, std::string strLambdas, 
                                         int minLambda, AggregateFunction<EdgeWeight>& defaultAgg);
};

#endif // _LMKAGGEXPERIMENTSCOMMAND_H
