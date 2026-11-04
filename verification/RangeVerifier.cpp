// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "RangeVerifier.h"

#include "../utility/utility.h"

#include <algorithm>
#include <sstream>

// Mirrors the synthetic query lists of LMKRangeExperimentsCommand::runQueries2. The ranges are percentages of the lower-bound
// graph diameter stored in the binary graph, and the resolved EdgeWeight range stored in each query is what the verification uses.
std::vector<QueryList<RangeQuery>> RangeVerifier::buildQueryLists(const VerifyArgs& args, Graph& graph)
{
    std::vector<QueryList<RangeQuery>> lists;
    std::vector<VertexID> queryVertices = loadQueryVertices(args.queryVertexFile);
    if (args.numSets == 0) {
        std::cerr << "No object sets were provided!" << std::endl;
        exit(1);
    }
    double defaultRangePercentage = std::stod(args.requiredParam("range_def_range"));
    double defaultDensity = std::stod(args.requiredParam("range_def_density"));
    std::string defaultObjType = args.requiredParam("range_def_type");
    int defaultObjVar = std::stoi(args.requiredParam("range_def_var"));
    std::string networkName = graph.getNetworkName();
    int numVertices = graph.getNumVertices();

    EdgeWeight lbGraphDiameter = graph.getDiameterLowerBound();
    std::cout << "INFO: Lower-bound graph diameter " << lbGraphDiameter << " used to compute ranges" << std::endl;
    EdgeWeight defaultRange = static_cast<int>(lbGraphDiameter*defaultRangePercentage/100.0);
    if (defaultRange == 0) {
        std::cerr << "Invalid default range percentage, results in default range of zero!" << std::endl;
        exit(1);
    }

    // Varying range
    if (utility::isQueryListOn(args.params,"range_varyingrange")) {
        QueryList<RangeQuery> varyingRange;
        varyingRange.name = "varyingrange";
        varyingRange.aggFunc = -1;
        std::vector<std::string> strRangeValues = utility::splitByDelim(args.rangeValues,',');
        if (strRangeValues.size() == 0) {
            std::cerr << "No range percentages provided (-a)!" << std::endl;
            exit(1);
        }
        for (std::size_t i = 0; i < strRangeValues.size(); ++i) {
            double rangePercentage = std::stod(strRangeValues[i]);
            if (rangePercentage <= 0) {
                std::cerr << "Invalid range percentage " << rangePercentage << " in list provided!" << std::endl;
                exit(1);
            }
            EdgeWeight range = static_cast<int>(lbGraphDiameter*rangePercentage/100.0);
            if (range == 0) {
                std::cerr << "Invalid varying range percentage, results in range of zero!" << std::endl;
                exit(1);
            }
            for (std::size_t j = 0; j < args.numSets; ++j) {
                std::string prefix = utility::constructObjIdxFilenamePrefix(networkName,defaultObjType,defaultDensity,defaultObjVar,j);
                for (std::size_t k = 0; k < queryVertices.size(); ++k) {
                    varyingRange.queries.push_back(RangeQuery(queryVertices[k],range,rangePercentage,defaultDensity,defaultObjType,defaultObjVar,j,prefix,static_cast<int>(numVertices*defaultDensity)));
                }
            }
        }
        lists.push_back(varyingRange);
    }

    // Varying density
    if (utility::isQueryListOn(args.params,"range_varyingdensity")) {
        QueryList<RangeQuery> varyingDensity;
        varyingDensity.name = "varyingdensity";
        varyingDensity.aggFunc = -1;
        std::vector<double> densities = parseDensities(args.objDensities);
        for (std::size_t i = 0; i < densities.size(); ++i) {
            for (std::size_t j = 0; j < args.numSets; ++j) {
                std::string prefix = utility::constructObjIdxFilenamePrefix(networkName,defaultObjType,densities[i],defaultObjVar,j);
                for (std::size_t k = 0; k < queryVertices.size(); ++k) {
                    varyingDensity.queries.push_back(RangeQuery(queryVertices[k],defaultRange,defaultRangePercentage,densities[i],defaultObjType,defaultObjVar,j,prefix,static_cast<int>(numVertices*densities[i])));
                }
            }
        }
        lists.push_back(varyingDensity);
    }

    return lists;
}

std::vector<std::string> RangeVerifier::groundTruthNames()
{
    return {"varyingrange","varyingdensity"};
}

std::vector<VertexID> RangeVerifier::queryVertices(const RangeQuery& query)
{
    return {query.queryVertex};
}

// Objects by distance from the query vertex, nearest first
void RangeVerifier::rankObjects(int aggFunc, const std::vector<VertexID>& objects, const std::vector<const std::vector<EdgeWeight>*>& distances,
                                std::vector<VertexDistancePair>& ranked)
{
    rankByDistance(*distances[0],objects,false,ranked);
}

// Every object within the range
void RangeVerifier::selectResult(const RangeQuery& query, const std::vector<VertexDistancePair>& ranked, std::vector<VertexDistancePair>& result)
{
    EdgeWeight range = query.range;
    auto end = std::partition_point(ranked.begin(),ranked.end(),[range](const VertexDistancePair& pair) { return pair.second <= range; });
    result.assign(ranked.begin(),end);
}

std::vector<std::string> RangeVerifier::methodNames()
{
    return {"ine","ier_phl","nvd_phl","colt_phl"};
}

void RangeVerifier::loadIndexes(const std::string& method, const VerifyArgs& args, Graph& graph, const std::vector<std::string>& prefixes)
{
    this->rtrees.clear();
    this->nvds.clear();
    this->colts.clear();
    if (method == "ine") {
        this->loadObjectSets(args,graph,prefixes,true);
    } else if (method == "ier_phl") {
        this->loadPHL(args,graph);
        this->loadObjectIndexes<PackedRtree>(args,prefixes,constants::OBJ_IDX_RTREE,this->rtrees);
    } else if (method == "nvd_phl") {
        this->loadPHL(args,graph);
        this->loadALT(args,graph);
        this->loadObjectIndexes<NVDQuadtree>(args,prefixes,constants::OBJ_IDX_NVD,this->nvds);
    } else if (method == "colt_phl") {
        this->loadPHL(args,graph);
        this->loadALT(args,graph);
        this->loadSULTree(args,graph);
        this->loadObjectIndexes<COLTree>(args,prefixes,constants::OBJ_IDX_COLTREE,this->colts);
    } else {
        std::cerr << "Unknown range method " << method << "!" << std::endl;
        exit(1);
    }
}

bool RangeVerifier::parallelQueries(const std::string& method)
{
    return method != "colt_phl";
}

// Same calls as the run* helpers of LMKRangeExperimentsCommand. The query objects are created per call (they hold nothing but
// statistics), so that runMethod can run in parallel
void RangeVerifier::runMethod(const std::string& method, const RangeQuery& query, int aggFunc, Graph& graph,
                              std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists)
{
    INE ine;
    IER ier;
    ILBR ilbr;
    const std::string& prefix = query.objSetFilenamePrefix;
    if (method == "ine") {
        ine.rangeSearch(graph,lookupIndex(this->objectStatus,prefix,"Object set"),lookupIndex(this->objectSetSizes,prefix,"Object set size"),
                        query.range,query.queryVertex,ids,dists);
    } else if (method == "ier_phl") {
        PackedRtree& rtree = lookupIndex(this->rtrees,prefix,"R-tree");
        if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
            ier.rangeSearchByPHL(rtree,query.range,query.queryVertex,ids,dists,graph,this->phl);
        } else {
            ier.rangeSearchByPHLTravelTimes(rtree,query.range,query.queryVertex,ids,dists,graph,this->phl);
        }
    } else if (method == "nvd_phl") {
        ilbr.rangeSearchyNVDQuadtreePHL(this->alt,lookupIndex(this->nvds,prefix,"NVD"),query.range,query.queryVertex,ids,dists,this->phl,graph);
    } else if (method == "colt_phl") {
        ilbr.rangeSearchByCOLTPHL(this->alt,lookupIndex(this->colts,prefix,"COL-Tree"),this->sltree,query.range,query.queryVertex,ids,dists,this->phl);
    } else {
        std::cerr << "Unknown range method " << method << "!" << std::endl;
        exit(1);
    }
}

bool RangeVerifier::compare(const std::vector<VertexDistancePair>& truth, std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists, std::string& message)
{
    return compareResults(truth,ids,dists,false,true,message);
}

std::string RangeVerifier::describe(const RangeQuery& query)
{
    std::stringstream ss;
    ss << "qv=" << query.queryVertex << " range=" << query.range << " (" << query.rangePerc << "%) set=" << query.objSetFilenamePrefix;
    return ss.str();
}
