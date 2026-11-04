// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "FknVerifier.h"

#include <sstream>

std::vector<QueryList<KnnQuery>> FknVerifier::buildQueryLists(const VerifyArgs& args, Graph& graph)
{
    return buildPointQueryLists(args,graph,"fkn");
}

std::vector<std::string> FknVerifier::groundTruthNames()
{
    return {"varyingk","varyingdensity"};
}

std::vector<VertexID> FknVerifier::queryVertices(const KnnQuery& query)
{
    return {query.queryVertex};
}

// Objects by distance from the query vertex, farthest first
void FknVerifier::rankObjects(int aggFunc, const std::vector<VertexID>& objects, const std::vector<const std::vector<EdgeWeight>*>& distances,
                              std::vector<VertexDistancePair>& ranked)
{
    rankByDistance(*distances[0],objects,true,ranked);
}

// The k farthest objects
void FknVerifier::selectResult(const KnnQuery& query, const std::vector<VertexDistancePair>& ranked, std::vector<VertexDistancePair>& result)
{
    firstK(ranked,query.k,result);
}

std::vector<std::string> FknVerifier::methodNames()
{
    return {"ine","llb_phl","colt_phl"};
}

void FknVerifier::loadIndexes(const std::string& method, const VerifyArgs& args, Graph& graph, const std::vector<std::string>& prefixes)
{
    this->colts.clear();
    if (method == "ine") {
        this->loadObjectSets(args,graph,prefixes,true);
    } else if (method == "llb_phl") {
        this->loadPHL(args,graph);
        this->loadALT(args,graph);
        this->loadObjectSets(args,graph,prefixes,false);
    } else if (method == "colt_phl") {
        this->loadPHL(args,graph);
        this->loadALT(args,graph);
        this->loadSULTree(args,graph);
        this->loadObjectIndexes<COLTree>(args,prefixes,constants::OBJ_IDX_COLTREE,this->colts);
    } else {
        std::cerr << "Unknown FkN method " << method << "!" << std::endl;
        exit(1);
    }
}

bool FknVerifier::parallelQueries(const std::string& method)
{
    return method != "colt_phl";
}

// Same calls as the run* helpers of LMKFarthestExperimentsCommand. The query objects are created per call (they hold nothing but
// statistics), so that runMethod can run in parallel
void FknVerifier::runMethod(const std::string& method, const KnnQuery& query, int aggFunc, Graph& graph,
                            std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists)
{
    INE ine;
    IUBR iubr;
    const std::string& prefix = query.objSetFilenamePrefix;
    if (method == "ine") {
        ine.getFkNs(graph,lookupIndex(this->objectStatus,prefix,"Object set"),lookupIndex(this->objectSetSizes,prefix,"Object set size"),
                    query.k,query.queryVertex,ids,dists);
    } else if (method == "llb_phl") {
        iubr.getFkNsByPHLNaive(this->alt,lookupIndex(this->objectSets,prefix,"Object set"),query.k,query.queryVertex,ids,dists,this->phl);
    } else if (method == "colt_phl") {
        iubr.getFkNsByCOLTPHL(this->alt,lookupIndex(this->colts,prefix,"COL-Tree"),this->sltree,query.k,query.queryVertex,ids,dists,this->phl);
    } else {
        std::cerr << "Unknown FkN method " << method << "!" << std::endl;
        exit(1);
    }
}

bool FknVerifier::compare(const std::vector<VertexDistancePair>& truth, std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists, std::string& message)
{
    return compareResults(truth,ids,dists,true,false,message);
}

std::string FknVerifier::describe(const KnnQuery& query)
{
    std::stringstream ss;
    ss << "qv=" << query.queryVertex << " k=" << query.k << " set=" << query.objSetFilenamePrefix;
    return ss.str();
}
