// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _AKNNVERIFIER_H
#define _AKNNVERIFIER_H

#include "Verifier.h"
#include "../processing/IER.h"
#include "../processing/ILBR.h"
#include "../processing/PackedRtree.h"
#include "../processing/NVDQuadtree.h"
#include "../utility/AggregateMax.h"
#include "../utility/AggregateSum.h"

class AknnVerifier: public Verifier<AknnQuery> {

    public:
        std::vector<QueryList<AknnQuery>> buildQueryLists(const VerifyArgs& args, Graph& graph);
        std::vector<std::string> groundTruthNames();
        std::vector<VertexID> queryVertices(const AknnQuery& query);
        void rankObjects(int aggFunc, const std::vector<VertexID>& objects, const std::vector<const std::vector<EdgeWeight>*>& distances,
                         std::vector<VertexDistancePair>& ranked);
        void selectResult(const AknnQuery& query, const std::vector<VertexDistancePair>& ranked, std::vector<VertexDistancePair>& result);
        std::vector<std::string> methodNames();
        void loadIndexes(const std::string& method, const VerifyArgs& args, Graph& graph, const std::vector<std::string>& prefixes);
        bool parallelQueries(const std::string& method);
        void runMethod(const std::string& method, const AknnQuery& query, int aggFunc, Graph& graph,
                       std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists);
        bool compare(const std::vector<VertexDistancePair>& truth, std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists, std::string& message);
        std::string describe(const AknnQuery& query);

    private:
        AggregateMax<EdgeWeight> aggMax;
        AggregateSum<EdgeWeight> aggSum;
        AggregateMax<EuclideanDistance> euclidAggMax;
        AggregateSum<EuclideanDistance> euclidAggSum;
        std::unordered_map<std::string,PackedRtree> rtrees;
        std::unordered_map<std::string,NVDQuadtree> nvds;
        std::unordered_map<std::string,COLTree> colts;

};

#endif // _AKNNVERIFIER_H
