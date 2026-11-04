// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _VERIFIER_H
#define _VERIFIER_H

#include "GroundTruth.h"
#include "../common.h"
#include "../processing/Graph.h"
#include "../processing/ALT.h"
#include "../processing/Gtree.h"
#include "../processing/SULTree.h"
#include "../processing/pruned_highway_labeling.h"
#include "../tuple/QuerySet.h"
#include "../utility/serialization.h"

#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

// Arguments of the verify command (the values the experiment scripts already pass to the experiment commands)
struct VerifyArgs {
    std::string queryType;       // -t: knn, aknn, fkn or range
    std::string bgrFilePath;     // -g: binary graph file
    std::string filePathPrefix;  // -f: the experiment output path
    std::string parameters;      // -p: as given
    std::unordered_map<std::string,std::string> params; // -p: parsed into key-value pairs
    std::string queryVertexFile;   // -q
    std::string kValues;         // -k
    std::string rangeValues;     // -a: percentages of the graph diameter
    std::size_t numSets = 0;     // -n: number of object sets
    std::string objDensities;    // -d
    std::string qvSetSizes;      // -z: AkNN query vertex set sizes
    std::string qvSubnetSizes;   // -y: AkNN query vertex set subnetwork sizes

    std::string objIdxDir() const { return filePathPrefix + "/verification/obj_indexes"; }
    std::string groundTruthDir() const { return filePathPrefix + "/verification/ground_truth"; }
    std::string param(const std::string& key) const;         // "" if not provided
    std::string requiredParam(const std::string& key) const; // exits if not provided
};

// Parameters, file names and query lists
std::unordered_map<std::string,std::string> parseParameters(const std::string& parameters);
std::string phlIndexPath(const VerifyArgs& args, Graph& graph);
std::string altIndexPath(const VerifyArgs& args, Graph& graph);
std::string sltreeIndexPath(const VerifyArgs& args, Graph& graph);
std::string gtreeIndexPath(const VerifyArgs& args, Graph& graph);
std::string objIndexPath(const VerifyArgs& args, const std::string& prefix, const std::string& extension);
std::vector<VertexID> loadQueryVertices(const std::string& queryVertexFile);
std::vector<unsigned int> parsePositiveIntegers(const std::string& list, const std::string& description);
std::vector<double> parseDensities(const std::string& list);
// The query lists of the kNN and FkN experiments (keyPrefix "knn" or "fkn")
std::vector<QueryList<KnnQuery>> buildPointQueryLists(const VerifyArgs& args, Graph& graph, const std::string& keyPrefix);

// Distances from a query vertex to every vertex (one full Dijkstra search); vertices not reached hold UNREACHED
const EdgeWeight UNREACHED = std::numeric_limits<EdgeWeight>::max();
void computeSSSP(Graph& graph, VertexID source, std::vector<EdgeWeight>& sssp);

// Result ordering and selection
bool nearestFirst(const VertexDistancePair& a, const VertexDistancePair& b);     // (distance asc, object asc)
bool farthestFirst(const VertexDistancePair& a, const VertexDistancePair& b);    // (distance desc, object asc)
void rankByDistance(const std::vector<EdgeWeight>& sssp, const std::vector<VertexID>& objects, bool farthest, std::vector<VertexDistancePair>& ranked);
void firstK(const std::vector<VertexDistancePair>& ranked, unsigned int k, std::vector<VertexDistancePair>& result);
void toSortedPairs(std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists, bool farthest, std::vector<VertexDistancePair>& result);

// Comparison of a method's result with the ground truth of one query (see Verifier.cpp for the rule)
bool compareResults(const std::vector<VertexDistancePair>& truth, std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists,
                    bool farthest, bool idSetOnly, std::string& message);
std::string pairsToString(const std::vector<VertexDistancePair>& pairs, std::size_t maxPairs);

// The objects step: copy the experiment object sets and their object indexes into the verification folder
void buildVerificationObjectIndexes(const VerifyArgs& args, Graph& graph);

// Index loaded for an object set
template <typename Map>
typename Map::mapped_type& lookupIndex(Map& indexes, const std::string& prefix, const std::string& description)
{
    auto it = indexes.find(prefix);
    if (it == indexes.end()) {
        std::cerr << description << " for object set " << prefix << " was not loaded!" << std::endl;
        exit(1);
    }
    return it->second;
}

// State shared by the verifiers of all query types: road network indexes and the object sets of the verification snapshot
class VerifierBase {

    public:
        virtual ~VerifierBase() {}
        // Object sets (prefix -> object vertices), their sizes and INE status vectors (for the INE method)
        std::unordered_map<std::string,std::vector<VertexID>> objectSets;
        std::unordered_map<std::string,int> objectSetSizes;
        std::unordered_map<std::string,std::vector<bool>> objectStatus;
        void loadObjectSets(const VerifyArgs& args, Graph& graph, const std::vector<std::string>& prefixes, bool buildStatus);
        // Run after a verification with failures: compares the object index of every object set that each of the methods
        // queries with the set file and lists those that do not hold exactly its objects (stale)
        void checkObjectIndexes(const VerifyArgs& args, Graph& graph, const std::vector<std::string>& methods,
                                const std::vector<std::string>& prefixes);

    protected:
        // Road network indexes, loaded once on first use
        PrunedHighwayLabeling phl;
        ALT alt;
        SULTree sltree;
        Gtree gtree;
        bool phlLoaded = false, altLoaded = false, sltreeLoaded = false, gtreeLoaded = false;
        void loadPHL(const VerifyArgs& args, Graph& graph);
        void loadALT(const VerifyArgs& args, Graph& graph);
        void loadSULTree(const VerifyArgs& args, Graph& graph);
        void loadGtree(const VerifyArgs& args, Graph& graph);

        // Object indexes of the given object sets (replaces any previously loaded)
        template <typename Index>
        void loadObjectIndexes(const VerifyArgs& args, const std::vector<std::string>& prefixes, const std::string& extension,
                               std::unordered_map<std::string,Index>& indexes)
        {
            indexes.clear();
            for (std::size_t i = 0; i < prefixes.size(); ++i) {
                indexes.emplace(prefixes[i],serialization::getIndexFromBinaryFile<Index>(objIndexPath(args,prefixes[i],extension)));
            }
        }
};

// One implementation per query type (KnnQuery for kNN and FkN, AknnQuery, RangeQuery)
template <typename Query>
class Verifier: public VerifierBase {

    public:
        // The query lists of this query type, built as the experiment command builds them
        virtual std::vector<QueryList<Query>> buildQueryLists(const VerifyArgs& args, Graph& graph) = 0;
        // Names of the ground truth files of this query type (list name, plus _max/_sum for AkNN)
        virtual std::vector<std::string> groundTruthNames() = 0;
        // The query vertices of a query: ground truth is computed with one SSSP from each of them
        virtual std::vector<VertexID> queryVertices(const Query& query) = 0;
        // Every reachable object of an object set with its (aggregate) distance from the query vertices, in result order
        // (distances[i] are the SSSP distances from queryVertices(query)[i]); done once per query vertices, object set and aggregate
        virtual void rankObjects(int aggFunc, const std::vector<VertexID>& objects, const std::vector<const std::vector<EdgeWeight>*>& distances,
                                 std::vector<VertexDistancePair>& ranked) = 0;
        // Ground truth of one query: a slice of the ranked objects (the first k, or those within the range)
        virtual void selectResult(const Query& query, const std::vector<VertexDistancePair>& ranked, std::vector<VertexDistancePair>& result) = 0;
        // Methods of this query type: their names, the indexes they need and one query
        virtual std::vector<std::string> methodNames() = 0;
        virtual void loadIndexes(const std::string& method, const VerifyArgs& args, Graph& graph, const std::vector<std::string>& prefixes) = 0;
        // Whether runMethod may be called from several threads at once (on different queries) for a method: false for the
        // methods whose indexes keep per-query state, colt_phl (COLTree::queryData) and gtree (the G-tree writes its own members)
        virtual bool parallelQueries(const std::string& method) = 0;
        virtual void runMethod(const std::string& method, const Query& query, int aggFunc, Graph& graph,
                               std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists) = 0;
        // Comparison of a method's result with the ground truth
        virtual bool compare(const std::vector<VertexDistancePair>& truth, std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists, std::string& message) = 0;
        // One-line description of a query for reports
        virtual std::string describe(const Query& query) = 0;

};

#endif // _VERIFIER_H
