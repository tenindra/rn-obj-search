// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "AknnVerifier.h"

#include "../utility/utility.h"

#include <algorithm>
#include <sstream>

// The first numVertices query vertices of every query vertex set
static std::vector<std::vector<VertexID>> truncateQueryVertexSets(const std::vector<std::vector<VertexID>>& queryVertexSets, std::size_t numVertices, const std::string& file)
{
    std::vector<std::vector<VertexID>> truncatedSets;
    for (std::size_t i = 0; i < queryVertexSets.size(); ++i) {
        if (queryVertexSets[i].size() < numVertices) {
            std::cerr << "Query vertex set " << i << " in " << file << " has " << queryVertexSets[i].size() << " query vertices but " << numVertices << " are required!" << std::endl;
            exit(1);
        }
        truncatedSets.push_back(std::vector<VertexID>(queryVertexSets[i].begin(),queryVertexSets[i].begin()+numVertices));
    }
    return truncatedSets;
}

// Mirrors the synthetic query lists of LMKAggregateExperimentsCommand::runQueries2, including their switches (a switched-off
// list is not built and its query vertex sets are not read); every list is produced for the max aggregate and for the sum aggregate
std::vector<QueryList<AknnQuery>> AknnVerifier::buildQueryLists(const VerifyArgs& args, Graph& graph)
{
    if (args.numSets == 0) {
        std::cerr << "No object sets were provided!" << std::endl;
        exit(1);
    }
    unsigned int defaultK = std::stoul(args.requiredParam("aknn_def_k"));
    double defaultDensity = std::stod(args.requiredParam("aknn_def_density"));
    std::string defaultObjType = args.requiredParam("aknn_def_type");
    int defaultObjVar = std::stoi(args.requiredParam("aknn_def_var"));
    int defaultNumQueryVertices = std::stoi(args.requiredParam("aknn_qv_def_setsize"));
    int defaultSubnetSize = std::stoi(args.requiredParam("aknn_qv_def_setsubnetsize"));
    int maxNumQueryVertices = std::stoi(args.requiredParam("aknn_qv_max_setsize"));
    int numQueryVertexSets = std::stoi(args.requiredParam("aknn_num_qv_sets"));
    if (defaultNumQueryVertices <= 0 || defaultSubnetSize <= 0 || defaultSubnetSize > 100 || maxNumQueryVertices <= 1 || numQueryVertexSets <= 1) {
        std::cerr << "Invalid AkNN query vertex set parameters!" << std::endl;
        exit(1);
    }
    std::string networkName = graph.getNetworkName();
    int numVertices = graph.getNumVertices();
    std::string vectorType;
    unsigned int numElementsPerVec;

    // Query vertex sets of the default subnetwork size, truncated to the default set size
    std::string defSubnetFile = args.filePathPrefix + "/query_sets/" + utility::constructQuerySetFileName(networkName,"subnet",numQueryVertexSets,maxNumQueryVertices,defaultSubnetSize);
    std::vector<std::vector<VertexID>> queryVertexSets = utility::getVectorSetFromFile(defSubnetFile,vectorType,numElementsPerVec);
    if (queryVertexSets.size() == 0) {
        std::cerr << "No query vertex sets were provided in " << defSubnetFile << "!" << std::endl;
        exit(1);
    }
    std::vector<std::vector<VertexID>> defaultSizeSets = truncateQueryVertexSets(queryVertexSets,defaultNumQueryVertices,defSubnetFile);

    std::vector<QueryList<AknnQuery>> baseLists;

    // Varying k
    if (utility::isQueryListOn(args.params,"aknn_varyingk")) {
        QueryList<AknnQuery> varyingK;
        varyingK.name = "varyingk";
        std::vector<unsigned int> kValues = parsePositiveIntegers(args.kValues,"k value");
        for (std::size_t i = 0; i < kValues.size(); ++i) {
            for (std::size_t j = 0; j < args.numSets; ++j) {
                std::string prefix = utility::constructObjIdxFilenamePrefix(networkName,defaultObjType,defaultDensity,defaultObjVar,j);
                for (std::size_t k = 0; k < defaultSizeSets.size(); ++k) {
                    varyingK.queries.push_back(AknnQuery(defaultSizeSets[k],kValues[i],defaultDensity,defaultSubnetSize,defaultObjType,defaultObjVar,j,prefix,static_cast<int>(numVertices*defaultDensity)));
                }
            }
        }
        baseLists.push_back(varyingK);
    }

    // Varying density
    if (utility::isQueryListOn(args.params,"aknn_varyingdensity")) {
        QueryList<AknnQuery> varyingDensity;
        varyingDensity.name = "varyingdensity";
        std::vector<double> densities = parseDensities(args.objDensities);
        for (std::size_t i = 0; i < densities.size(); ++i) {
            for (std::size_t j = 0; j < args.numSets; ++j) {
                std::string prefix = utility::constructObjIdxFilenamePrefix(networkName,defaultObjType,densities[i],defaultObjVar,j);
                for (std::size_t k = 0; k < defaultSizeSets.size(); ++k) {
                    varyingDensity.queries.push_back(AknnQuery(defaultSizeSets[k],defaultK,densities[i],defaultSubnetSize,defaultObjType,defaultObjVar,j,prefix,static_cast<int>(numVertices*densities[i])));
                }
            }
        }
        baseLists.push_back(varyingDensity);
    }

    // Varying query vertex set size (the first s query vertices of each set)
    if (utility::isQueryListOn(args.params,"aknn_varyingqvsetsize")) {
        QueryList<AknnQuery> varyingSetSize;
        varyingSetSize.name = "varyingqvsetsize";
        std::vector<unsigned int> setSizes = parsePositiveIntegers(args.qvSetSizes,"query vertex set size");
        for (std::size_t i = 0; i < setSizes.size(); ++i) {
            std::vector<std::vector<VertexID>> currentSizeSets = truncateQueryVertexSets(queryVertexSets,setSizes[i],defSubnetFile);
            for (std::size_t j = 0; j < args.numSets; ++j) {
                std::string prefix = utility::constructObjIdxFilenamePrefix(networkName,defaultObjType,defaultDensity,defaultObjVar,j);
                for (std::size_t k = 0; k < currentSizeSets.size(); ++k) {
                    varyingSetSize.queries.push_back(AknnQuery(currentSizeSets[k],defaultK,defaultDensity,defaultSubnetSize,defaultObjType,defaultObjVar,j,prefix,static_cast<int>(numVertices*defaultDensity)));
                }
            }
        }
        baseLists.push_back(varyingSetSize);
    }

    // Varying query vertex set subnetwork size (one query vertex set file per subnetwork size)
    if (utility::isQueryListOn(args.params,"aknn_varyingqvsetsubnetsize")) {
        QueryList<AknnQuery> varyingSubnetSize;
        varyingSubnetSize.name = "varyingqvsetsubnetsize";
        std::vector<unsigned int> subnetSizes = parsePositiveIntegers(args.qvSubnetSizes,"query vertex set subnet size");
        for (std::size_t i = 0; i < subnetSizes.size(); ++i) {
            if (subnetSizes[i] > 100) {
                std::cerr << "Invalid query vertex set subnet size " << subnetSizes[i] << "!" << std::endl;
                exit(1);
            }
            std::string subnetFile = args.filePathPrefix + "/query_sets/" + utility::constructQuerySetFileName(networkName,"subnet",numQueryVertexSets,maxNumQueryVertices,subnetSizes[i]);
            std::vector<std::vector<VertexID>> subnetQueryVertexSets = utility::getVectorSetFromFile(subnetFile,vectorType,numElementsPerVec);
            std::vector<std::vector<VertexID>> subnetDefaultSizeSets = truncateQueryVertexSets(subnetQueryVertexSets,defaultNumQueryVertices,subnetFile);
            for (std::size_t j = 0; j < args.numSets; ++j) {
                std::string prefix = utility::constructObjIdxFilenamePrefix(networkName,defaultObjType,defaultDensity,defaultObjVar,j);
                for (std::size_t k = 0; k < subnetDefaultSizeSets.size(); ++k) {
                    varyingSubnetSize.queries.push_back(AknnQuery(subnetDefaultSizeSets[k],defaultK,defaultDensity,subnetSizes[i],defaultObjType,defaultObjVar,j,prefix,static_cast<int>(numVertices*defaultDensity)));
                }
            }
        }
        baseLists.push_back(varyingSubnetSize);
    }

    // Every list for both aggregate functions
    std::vector<QueryList<AknnQuery>> lists;
    for (std::size_t i = 0; i < baseLists.size(); ++i) {
        baseLists[i].aggFunc = constants::AGG_MAX;
        lists.push_back(baseLists[i]);
        baseLists[i].aggFunc = constants::AGG_SUM;
        lists.push_back(baseLists[i]);
    }
    return lists;
}

std::vector<std::string> AknnVerifier::groundTruthNames()
{
    std::vector<std::string> baseNames = {"varyingk","varyingdensity","varyingqvsetsize","varyingqvsetsubnetsize"};
    std::vector<std::string> names;
    for (std::size_t i = 0; i < baseNames.size(); ++i) {
        names.push_back(baseNames[i] + "_max");
        names.push_back(baseNames[i] + "_sum");
    }
    return names;
}

std::vector<VertexID> AknnVerifier::queryVertices(const AknnQuery& query)
{
    return query.queryVertices;
}

// Objects by the aggregate (max or sum) of their distances from the query vertices, smallest first
void AknnVerifier::rankObjects(int aggFunc, const std::vector<VertexID>& objects, const std::vector<const std::vector<EdgeWeight>*>& distances,
                               std::vector<VertexDistancePair>& ranked)
{
    ranked.clear();
    ranked.reserve(objects.size());
    for (std::size_t i = 0; i < objects.size(); ++i) {
        unsigned long long aggregate = 0;
        bool reachable = true;
        for (std::size_t j = 0; j < distances.size(); ++j) {
            EdgeWeight distance = (*distances[j])[objects[i]];
            if (distance == UNREACHED) {
                reachable = false;
                break;
            } else if (aggFunc == constants::AGG_SUM) {
                aggregate += distance;
            } else if (distance > aggregate) {
                aggregate = distance;
            }
        }
        if (reachable) {
            ranked.push_back(VertexDistancePair(objects[i],static_cast<EdgeWeight>(aggregate)));
        }
    }
    std::sort(ranked.begin(),ranked.end(),nearestFirst);
}

// The k objects with the smallest aggregate distance
void AknnVerifier::selectResult(const AknnQuery& query, const std::vector<VertexDistancePair>& ranked, std::vector<VertexDistancePair>& result)
{
    firstK(ranked,query.k,result);
}

std::vector<std::string> AknnVerifier::methodNames()
{
    return {"ier_phl","nvd_phl","colt_phl"};
}

void AknnVerifier::loadIndexes(const std::string& method, const VerifyArgs& args, Graph& graph, const std::vector<std::string>& prefixes)
{
    this->rtrees.clear();
    this->nvds.clear();
    this->colts.clear();
    if (method == "ier_phl") {
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
        std::cerr << "Unknown AkNN method " << method << "!" << std::endl;
        exit(1);
    }
}

bool AknnVerifier::parallelQueries(const std::string& method)
{
    return method != "colt_phl";
}

// Same calls as the run* helpers of LMKAggregateExperimentsCommand. The query objects are created per call (they hold nothing but
// statistics), so that runMethod can run in parallel
void AknnVerifier::runMethod(const std::string& method, const AknnQuery& query, int aggFunc, Graph& graph,
                             std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists)
{
    IER ier;
    ILBR ilbr;
    const std::string& prefix = query.objSetFilenamePrefix;
    std::vector<VertexID> queryVertices = query.queryVertices;
    AggregateFunction<EdgeWeight>* agg;
    AggregateFunction<EuclideanDistance>* euclidAgg;
    if (aggFunc == constants::AGG_SUM) {
        agg = &this->aggSum;
        euclidAgg = &this->euclidAggSum;
    } else {
        agg = &this->aggMax;
        euclidAgg = &this->euclidAggMax;
    }
    if (method == "ier_phl") {
        PackedRtree& rtree = lookupIndex(this->rtrees,prefix,"R-tree");
        if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
            ier.getAkNNsByPHL(rtree,query.k,queryVertices,ids,dists,graph,this->phl,*agg,*euclidAgg);
        } else {
            ier.getAkNNsByPHLTravelTimes(rtree,query.k,queryVertices,ids,dists,graph,this->phl,*agg,*euclidAgg);
        }
    } else if (method == "nvd_phl") {
        ilbr.getAkNNByNVDQuadtreePHLByConcExps(this->alt,lookupIndex(this->nvds,prefix,"NVD"),query.k,queryVertices,ids,dists,this->phl,graph,*agg);
    } else if (method == "colt_phl") {
        ilbr.getAkNNByCOLTPHL(this->alt,lookupIndex(this->colts,prefix,"COL-Tree"),this->sltree,query.k,queryVertices,ids,dists,this->phl,*agg);
    } else {
        std::cerr << "Unknown AkNN method " << method << "!" << std::endl;
        exit(1);
    }
}

bool AknnVerifier::compare(const std::vector<VertexDistancePair>& truth, std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists, std::string& message)
{
    return compareResults(truth,ids,dists,false,false,message);
}

std::string AknnVerifier::describe(const AknnQuery& query)
{
    std::stringstream ss;
    ss << "qvs=[";
    for (std::size_t i = 0; i < query.queryVertices.size(); ++i) {
        ss << (i == 0 ? "" : ",") << query.queryVertices[i];
    }
    ss << "] k=" << query.k << " subnet=" << query.querySubnetSize << "% set=" << query.objSetFilenamePrefix;
    return ss.str();
}
