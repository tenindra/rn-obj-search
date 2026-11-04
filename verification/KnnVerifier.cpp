// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "KnnVerifier.h"

#include <sstream>

std::vector<QueryList<KnnQuery>> KnnVerifier::buildQueryLists(const VerifyArgs& args, Graph& graph)
{
    return buildPointQueryLists(args,graph,"knn");
}

std::vector<std::string> KnnVerifier::groundTruthNames()
{
    return {"varyingk","varyingdensity"};
}

std::vector<VertexID> KnnVerifier::queryVertices(const KnnQuery& query)
{
    return {query.queryVertex};
}

// Objects by distance from the query vertex, nearest first
void KnnVerifier::rankObjects(int aggFunc, const std::vector<VertexID>& objects, const std::vector<const std::vector<EdgeWeight>*>& distances,
                              std::vector<VertexDistancePair>& ranked)
{
    rankByDistance(*distances[0],objects,false,ranked);
}

// The k nearest objects
void KnnVerifier::selectResult(const KnnQuery& query, const std::vector<VertexDistancePair>& ranked, std::vector<VertexDistancePair>& result)
{
    firstK(ranked,query.k,result);
}

std::vector<std::string> KnnVerifier::methodNames()
{
    return {"ine","gtree","ier_phl","ol_phl","nvd_phl","llb_phl","colt_phl"};
}

void KnnVerifier::loadIndexes(const std::string& method, const VerifyArgs& args, Graph& graph, const std::vector<std::string>& prefixes)
{
    this->occLists.clear();
    this->rtrees.clear();
    this->objLists.clear();
    this->nvds.clear();
    this->colts.clear();
    if (method == "ine") {
        this->loadObjectSets(args,graph,prefixes,true);
    } else if (method == "gtree") {
        this->loadGtree(args,graph);
        this->loadObjectIndexes<OccurenceList>(args,prefixes,constants::OBJ_IDX_GTREE,this->occLists);
    } else if (method == "ier_phl") {
        this->loadPHL(args,graph);
        this->loadObjectIndexes<PackedRtree>(args,prefixes,constants::OBJ_IDX_RTREE,this->rtrees);
    } else if (method == "ol_phl") {
        this->loadPHL(args,graph);
        this->loadALT(args,graph);
        this->loadObjectIndexes<ObjectList>(args,prefixes,constants::OBJ_IDX_LNDMRK,this->objLists);
    } else if (method == "nvd_phl") {
        this->loadPHL(args,graph);
        this->loadALT(args,graph);
        this->loadObjectIndexes<NVDQuadtree>(args,prefixes,constants::OBJ_IDX_NVD,this->nvds);
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
        std::cerr << "Unknown kNN method " << method << "!" << std::endl;
        exit(1);
    }
}

bool KnnVerifier::parallelQueries(const std::string& method)
{
    return method != "gtree" && method != "colt_phl";
}

// Same calls as the run* helpers of LMKExperimentsCommand. The query objects are created per call (they hold nothing but
// statistics), so that runMethod can run in parallel
void KnnVerifier::runMethod(const std::string& method, const KnnQuery& query, int aggFunc, Graph& graph,
                            std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists)
{
    INE ine;
    IER ier;
    ILBR ilbr;
    const std::string& prefix = query.objSetFilenamePrefix;
    if (method == "ine") {
        ine.getKNNs(graph,lookupIndex(this->objectStatus,prefix,"Object set"),query.k,query.queryVertex,ids,dists);
    } else if (method == "gtree") {
        this->gtree.getKNNs(lookupIndex(this->occLists,prefix,"Occurence list"),query.k,query.queryVertex,ids,dists,graph);
    } else if (method == "ier_phl") {
        PackedRtree& rtree = lookupIndex(this->rtrees,prefix,"R-tree");
        if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
            ier.getKNNsByPHL(rtree,query.k,query.queryVertex,ids,dists,graph,this->phl);
        } else {
            ier.getKNNsByPHLTravelTimes(rtree,query.k,query.queryVertex,ids,dists,graph,this->phl);
        }
    } else if (method == "ol_phl") {
        ilbr.getKNNsByPHLImprovedByQueryLandmark(this->alt,lookupIndex(this->objLists,prefix,"Object list"),query.k,query.queryVertex,ids,dists,this->phl);
    } else if (method == "nvd_phl") {
        ilbr.getKNNByNVDQuadtreePHL(this->alt,lookupIndex(this->nvds,prefix,"NVD"),query.k,query.queryVertex,ids,dists,this->phl,graph);
    } else if (method == "llb_phl") {
        ilbr.getKNNsByPHLNaive(this->alt,lookupIndex(this->objectSets,prefix,"Object set"),query.k,query.queryVertex,ids,dists,this->phl);
    } else if (method == "colt_phl") {
        ilbr.getKNNByCOLTPHL(this->alt,lookupIndex(this->colts,prefix,"COL-Tree"),this->sltree,query.k,query.queryVertex,ids,dists,this->phl);
    } else {
        std::cerr << "Unknown kNN method " << method << "!" << std::endl;
        exit(1);
    }
}

bool KnnVerifier::compare(const std::vector<VertexDistancePair>& truth, std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists, std::string& message)
{
    return compareResults(truth,ids,dists,false,false,message);
}

std::string KnnVerifier::describe(const KnnQuery& query)
{
    std::stringstream ss;
    ss << "qv=" << query.queryVertex << " k=" << query.k << " set=" << query.objSetFilenamePrefix;
    return ss.str();
}
