// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _KNNVERIFIER_H
#define _KNNVERIFIER_H

#include "Verifier.h"
#include "../processing/INE.h"
#include "../processing/IER.h"
#include "../processing/ILBR.h"
#include "../processing/PackedRtree.h"
#include "../processing/NVDQuadtree.h"

class KnnVerifier: public Verifier<KnnQuery> {

    public:
        std::vector<QueryList<KnnQuery>> buildQueryLists(const VerifyArgs& args, Graph& graph);
        std::vector<std::string> groundTruthNames();
        std::vector<VertexID> queryVertices(const KnnQuery& query);
        void rankObjects(int aggFunc, const std::vector<VertexID>& objects, const std::vector<const std::vector<EdgeWeight>*>& distances,
                         std::vector<VertexDistancePair>& ranked);
        void selectResult(const KnnQuery& query, const std::vector<VertexDistancePair>& ranked, std::vector<VertexDistancePair>& result);
        std::vector<std::string> methodNames();
        void loadIndexes(const std::string& method, const VerifyArgs& args, Graph& graph, const std::vector<std::string>& prefixes);
        bool parallelQueries(const std::string& method);
        void runMethod(const std::string& method, const KnnQuery& query, int aggFunc, Graph& graph,
                       std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists);
        bool compare(const std::vector<VertexDistancePair>& truth, std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists, std::string& message);
        std::string describe(const KnnQuery& query);

    private:
        std::unordered_map<std::string,OccurenceList> occLists;
        std::unordered_map<std::string,PackedRtree> rtrees;
        std::unordered_map<std::string,ObjectList> objLists;
        std::unordered_map<std::string,NVDQuadtree> nvds;
        std::unordered_map<std::string,COLTree> colts;

};

#endif // _KNNVERIFIER_H
