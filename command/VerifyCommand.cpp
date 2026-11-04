// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "VerifyCommand.h"

#include "../common.h"
#include "../processing/Graph.h"
#include "../utility/serialization.h"
#include "../utility/utility.h"
#include "../utility/StopWatch.h"
#include "../verification/GroundTruth.h"
#include "../verification/KnnVerifier.h"
#include "../verification/FknVerifier.h"
#include "../verification/RangeVerifier.h"
#include "../verification/AknnVerifier.h"

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <utility>
#include <unistd.h>
#include <getopt.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

// Ground truth file name of a query list: the list name, plus the aggregate function for AkNN
std::string groundTruthName(const std::string& listName, int aggFunc)
{
    if (aggFunc == constants::AGG_MAX) {
        return listName + "_max";
    } else if (aggFunc == constants::AGG_SUM) {
        return listName + "_sum";
    }
    return listName;
}

std::string groundTruthPath(const VerifyArgs& args, const std::string& networkName, const std::string& name)
{
    return args.groundTruthDir() + "/" + networkName + "_" + args.queryType + "_" + name + ".gt";
}

template <typename Query>
void addObjectSetPrefixes(const std::vector<Query>& queries, std::vector<std::string>& prefixes, std::set<std::string>& seen)
{
    for (std::size_t i = 0; i < queries.size(); ++i) {
        if (seen.insert(queries[i].objSetFilenamePrefix).second) {
            prefixes.push_back(queries[i].objSetFilenamePrefix);
        }
    }
}

// A query in one of the lists
struct QueryRef {
    std::size_t list;
    std::size_t index;
};

// The queries of one query vertex group on one object set with one aggregate function: the objects are ranked once and
// every query takes its slice
struct ObjectSetQueries {
    const std::vector<VertexID>* objects;
    int aggFunc;
    std::vector<QueryRef> refs;
};

// The queries sharing the same query vertices: one SSSP per query vertex serves all of them
struct QueryVertexGroup {
    const std::vector<VertexID>* vertices;
    std::vector<ObjectSetQueries> objectSets;
    std::size_t numQueries;
};

// Computes the ground truth of every query list of the query type and saves one file per list. The experiments query
// every query vertex (or query vertex set for AkNN) against every object set and setting, so the loop is: for each query vertex
// group, one SSSP per query vertex; for each object set the group is queried against, rank its objects once; for each query
// on that object set, take its slice of the ranking. Parallel over query vertex groups.
template <typename Query>
void computeGroundTruth(Verifier<Query>& verifier, const VerifyArgs& args, Graph& graph, int numThreads)
{
    std::vector<QueryList<Query>> lists = verifier.buildQueryLists(args,graph);

    // Object sets used by the lists
    std::vector<std::string> prefixes;
    std::set<std::string> seen;
    for (std::size_t l = 0; l < lists.size(); ++l) {
        addObjectSetPrefixes(lists[l].queries,prefixes,seen);
    }
    verifier.loadObjectSets(args,graph,prefixes,false);

    // Queries grouped by query vertices, then by object set and aggregate function
    std::map<std::vector<VertexID>,std::map<std::pair<std::string,int>,std::vector<QueryRef>>> groupMap;
    std::size_t totalQueries = 0;
    for (std::size_t l = 0; l < lists.size(); ++l) {
        for (std::size_t i = 0; i < lists[l].queries.size(); ++i) {
            QueryRef ref = {l,i};
            const Query& query = lists[l].queries[i];
            groupMap[verifier.queryVertices(query)][std::make_pair(query.objSetFilenamePrefix,lists[l].aggFunc)].push_back(ref);
        }
        totalQueries += lists[l].queries.size();
    }
    std::vector<QueryVertexGroup> groups;
    for (auto git = groupMap.begin(); git != groupMap.end(); ++git) {
        QueryVertexGroup group;
        group.vertices = &git->first;
        group.numQueries = 0;
        for (auto sit = git->second.begin(); sit != git->second.end(); ++sit) {
            ObjectSetQueries objectSetQueries;
            objectSetQueries.objects = &verifier.objectSets.at(sit->first.first);
            objectSetQueries.aggFunc = sit->first.second;
            objectSetQueries.refs = sit->second;
            group.numQueries += sit->second.size();
            group.objectSets.push_back(objectSetQueries);
        }
        groups.push_back(group);
    }
    std::vector<std::vector<std::vector<VertexDistancePair>>> results(lists.size());
    for (std::size_t l = 0; l < lists.size(); ++l) {
        results[l].resize(lists[l].queries.size());
    }

#ifdef _OPENMP
    if (numThreads > 0) {
        omp_set_num_threads(numThreads);
    }
    std::cout << "Computing ground truth for " << totalQueries << " queries with " << groups.size() << " query vertex groups using "
              << omp_get_max_threads() << " threads" << std::endl;
#else
    std::cout << "Computing ground truth for " << totalQueries << " queries with " << groups.size() << " query vertex groups (built without OpenMP, single thread)" << std::endl;
#endif

    std::size_t progressStep = totalQueries >= 10 ? totalQueries/10 : 1;
    std::size_t done = 0;
    StopWatch sw;
    sw.start();
    #pragma omp parallel
    {
        std::vector<std::vector<EdgeWeight>> distances; // Per thread: one SSSP per query vertex of the current group
        std::vector<const std::vector<EdgeWeight>*> distancePtrs;
        std::vector<VertexDistancePair> ranked, result;
        #pragma omp for schedule(dynamic)
        for (long g = 0; g < static_cast<long>(groups.size()); ++g) {
            const std::vector<VertexID>& vertices = *groups[g].vertices;
            if (distances.size() < vertices.size()) {
                distances.resize(vertices.size());
            }
            distancePtrs.clear();
            for (std::size_t s = 0; s < vertices.size(); ++s) {
                computeSSSP(graph,vertices[s],distances[s]);
                distancePtrs.push_back(&distances[s]);
            }
            for (std::size_t o = 0; o < groups[g].objectSets.size(); ++o) {
                const ObjectSetQueries& objectSetQueries = groups[g].objectSets[o];
                verifier.rankObjects(objectSetQueries.aggFunc,*objectSetQueries.objects,distancePtrs,ranked);
                for (std::size_t r = 0; r < objectSetQueries.refs.size(); ++r) {
                    const QueryRef& ref = objectSetQueries.refs[r];
                    result.clear();
                    verifier.selectResult(lists[ref.list].queries[ref.index],ranked,result);
                    results[ref.list][ref.index] = result;
                }
            }
            #pragma omp critical
            {
                std::size_t before = done;
                done += groups[g].numQueries;
                if (done/progressStep != before/progressStep || done == totalQueries) {
                    std::cout << "  " << done << "/" << totalQueries << " queries" << std::endl;
                }
            }
        }
    }
    sw.stop();
    std::cout << "Ground truth computed in " << sw.getTimeMs()/1000.0 << "s" << std::endl;

    for (std::size_t l = 0; l < lists.size(); ++l) {
        std::size_t numQueries = lists[l].queries.size(), totalResults = 0;
        for (std::size_t i = 0; i < numQueries; ++i) {
            totalResults += results[l][i].size();
        }
        GroundTruthSet<Query> groundTruth;
        groundTruth.networkName = graph.getNetworkName();
        groundTruth.numVertices = graph.getNumVertices();
        groundTruth.numEdges = graph.getNumEdges();
        groundTruth.queryType = args.queryType;
        groundTruth.listName = lists[l].name;
        groundTruth.aggFunc = lists[l].aggFunc;
        groundTruth.queries = lists[l].queries;
        groundTruth.results.swap(results[l]);
        std::string path = groundTruthPath(args,graph.getNetworkName(),groundTruthName(lists[l].name,lists[l].aggFunc));
        saveGroundTruth(groundTruth,path);
        std::cout << "Saved " << path << ": " << numQueries << " queries, " << (numQueries > 0 ? totalResults/static_cast<double>(numQueries) : 0)
                  << " results per query on average" << std::endl;
    }
}

// Runs the selected methods on every saved query list and compares each result with the ground truth; returns the number of
// failures. The queries of a list run in parallel for the methods that allow it (Verifier::parallelQueries)
template <typename Query>
std::size_t verifyMethods(Verifier<Query>& verifier, const VerifyArgs& args, Graph& graph, const std::string& methods, int numThreads)
{
    std::string networkName = graph.getNetworkName();

    // Ground truth
    std::vector<GroundTruthSet<Query>> groundTruths;
    std::vector<std::string> names = verifier.groundTruthNames();
    for (std::size_t i = 0; i < names.size(); ++i) {
        // A query list switched off in the parameters is not verified, as the experiment does not run it
        std::string list = names[i].substr(0,names[i].find('_')); // e.g. varyingk of varyingk_max
        if (!utility::isQueryListOn(args.params,args.queryType + "_" + list)) {
            std::cout << "Query list " << names[i] << " is switched off, skipped" << std::endl;
            continue;
        }
        std::string path = groundTruthPath(args,networkName,names[i]);
        GroundTruthSet<Query> groundTruth;
        if (!loadGroundTruth(path,groundTruth)) {
            std::cout << "Ground truth file " << path << " not found, skipped" << std::endl;
            continue;
        }
        if (groundTruth.networkName != networkName || groundTruth.numVertices != graph.getNumVertices() || groundTruth.queryType != args.queryType) {
            std::cerr << "Ground truth file " << path << " is for " << groundTruth.queryType << " queries on " << groundTruth.networkName
                      << " (" << groundTruth.numVertices << " vertices), not " << args.queryType << " queries on " << networkName
                      << " (" << graph.getNumVertices() << " vertices)!" << std::endl;
            exit(1);
        }
        std::cout << "Loaded " << path << " (" << groundTruth.queries.size() << " queries)" << std::endl;
        groundTruths.push_back(std::move(groundTruth));
    }
    if (groundTruths.size() == 0) {
        std::cerr << "No ground truth files found in " << args.groundTruthDir() << " (run the ground_truth step first)!" << std::endl;
        exit(1);
    }

    // Methods: those enabled in the parameters by default
    std::vector<std::string> allMethods = verifier.methodNames(), selectedMethods;
    if (methods == "") {
        for (std::size_t i = 0; i < allMethods.size(); ++i) {
            if (args.param(allMethods[i]) == "1") {
                selectedMethods.push_back(allMethods[i]);
            }
        }
        if (selectedMethods.size() == 0) {
            std::cerr << "No methods enabled in -p (e.g. colt_phl=1); use -m all or -m <method,...>" << std::endl;
            exit(1);
        }
    } else if (methods == "all") {
        selectedMethods = allMethods;
    } else {
        selectedMethods = utility::splitByDelim(methods,',');
        for (std::size_t i = 0; i < selectedMethods.size(); ++i) {
            if (std::find(allMethods.begin(),allMethods.end(),selectedMethods[i]) == allMethods.end()) {
                std::cerr << "Unknown method " << selectedMethods[i] << " for " << args.queryType << " queries (methods:";
                for (std::size_t j = 0; j < allMethods.size(); ++j) {
                    std::cerr << " " << allMethods[j];
                }
                std::cerr << ")" << std::endl;
                exit(1);
            }
        }
    }

    // Object sets used by the ground truth
    std::vector<std::string> prefixes;
    std::set<std::string> seen;
    for (std::size_t i = 0; i < groundTruths.size(); ++i) {
        addObjectSetPrefixes(groundTruths[i].queries,prefixes,seen);
    }

    // Threads for the methods whose queries may run in parallel (all available unless -j was given)
    int maxThreads = 1;
#ifdef _OPENMP
    if (numThreads > 0) {
        omp_set_num_threads(numThreads);
    }
    maxThreads = omp_get_max_threads();
#endif

    std::size_t totalFailures = 0;
    std::vector<std::string> failedMethods;
    for (std::size_t m = 0; m < selectedMethods.size(); ++m) {
        const std::string& method = selectedMethods[m];
        bool threaded = maxThreads > 1 && verifier.parallelQueries(method);
        std::cout << "\nVerifying " << method << (threaded ? " using " + std::to_string(maxThreads) + " threads" : "") << "..." << std::endl;
        verifier.loadIndexes(method,args,graph,prefixes);
        StopWatch sw;
        sw.start();
        std::size_t numQueries = 0, failures = 0;
        for (std::size_t g = 0; g < groundTruths.size(); ++g) {
            GroundTruthSet<Query>& groundTruth = groundTruths[g];
            std::string name = groundTruthName(groundTruth.listName,groundTruth.aggFunc);
            // The queries of the list, in parallel if the method allows it, each thread with its own result vectors. A failed
            // query's report is kept at the query's position and printed after the list, so the output is that of a serial run
            std::vector<std::string> reports(groundTruth.queries.size());
            #pragma omp parallel if(threaded)
            {
                std::vector<VertexID> ids;
                std::vector<EdgeWeight> dists;
                std::string message;
                #pragma omp for schedule(dynamic)
                for (std::size_t i = 0; i < groundTruth.queries.size(); ++i) {
                    ids.clear();
                    dists.clear();
                    verifier.runMethod(method,groundTruth.queries[i],groundTruth.aggFunc,graph,ids,dists);
                    if (!verifier.compare(groundTruth.results[i],ids,dists,message)) {
                        std::stringstream ss;
                        ss << "FAIL " << name << " query " << i << " (" << verifier.describe(groundTruth.queries[i]) << ") method " << method << ": " << message << "\n"
                           << "  expected: " << pairsToString(groundTruth.results[i],10) << "\n"
                           << "  actual:   " << pairsToString(utility::combineVerticesAndDistanceVectors(ids,dists),10);
                        reports[i] = ss.str();
                    }
                }
            }
            numQueries += groundTruth.queries.size();
            for (std::size_t i = 0; i < reports.size(); ++i) {
                if (reports[i] != "") {
                    ++failures;
                    std::cout << reports[i] << std::endl;
                }
            }
        }
        sw.stop();
        std::cout << method << ": " << numQueries << " queries, " << failures << " failures (" << sw.getTimeMs()/1000.0 << "s)" << std::endl;
        totalFailures += failures;
        if (failures > 0) {
            failedMethods.push_back(method);
        }
    }
    std::cout << std::endl;
    if (totalFailures == 0) {
        std::cout << "All " << args.queryType << " results verified" << std::endl;
    } else {
        std::cout << "VERIFICATION FAILED: " << totalFailures << " incorrect " << args.queryType << " results" << std::endl;
        // Rules out (or finds) stale object indexes as the cause
        verifier.checkObjectIndexes(args,graph,failedMethods,prefixes);
    }
    return totalFailures;
}

template <typename Query>
void dumpGroundTruth(Verifier<Query>& verifier, const std::string& filePath)
{
    GroundTruthSet<Query> groundTruth;
    if (!loadGroundTruth(filePath,groundTruth)) {
        std::cerr << "Cannot open ground truth file " << filePath << std::endl;
        exit(1);
    }
    std::cout << "Network: " << groundTruth.networkName << " (" << groundTruth.numVertices << " vertices, " << groundTruth.numEdges << " edges)\n"
              << "Query type: " << groundTruth.queryType << "\n"
              << "List: " << groundTruthName(groundTruth.listName,groundTruth.aggFunc) << "\n"
              << "Queries: " << groundTruth.queries.size() << "\n\n";
    for (std::size_t i = 0; i < groundTruth.queries.size(); ++i) {
        std::cout << i << "\t" << verifier.describe(groundTruth.queries[i]) << "\t" << groundTruth.results[i].size() << " results\t"
                  << pairsToString(groundTruth.results[i],groundTruth.results[i].size()) << "\n";
    }
}

}

template <typename Query>
void VerifyCommand::runStep(Verifier<Query>& verifier, const std::string& step, const VerifyArgs& args, const std::string& methods,
                            const std::string& gtFile, int numThreads)
{
    if (step == this->STEP_DUMP) {
        dumpGroundTruth(verifier,gtFile);
        return;
    }
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(args.bgrFilePath);
    if (step == this->STEP_OBJECTS) {
        buildVerificationObjectIndexes(args,graph);
    } else if (step == this->STEP_GROUND_TRUTH) {
        computeGroundTruth(verifier,args,graph,numThreads);
    } else if (step == this->STEP_VERIFY) {
        if (verifyMethods(verifier,args,graph,methods,numThreads) > 0) {
            exit(1);
        }
    }
}

void VerifyCommand::execute(int argc, char* argv[])
{
    std::string step = "";
    std::string methods = "";
    std::string gtFile = "";
    int numThreads = 0;
    VerifyArgs args;

    /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "e:t:g:f:p:q:k:a:n:d:z:y:m:j:i:")) != -1) {
        switch (opt) {
            case 'e':
                step = optarg;
                break;
            case 't':
                args.queryType = optarg;
                break;
            case 'g':
                args.bgrFilePath = optarg;
                break;
            case 'f':
                args.filePathPrefix = optarg;
                break;
            case 'p':
                args.parameters = optarg;
                break;
            case 'q':
                args.queryVertexFile = optarg;
                break;
            case 'k':
                args.kValues = optarg;
                break;
            case 'a':
                args.rangeValues = optarg;
                break;
            case 'n':
                args.numSets = std::stoul(optarg);
                break;
            case 'd':
                args.objDensities = optarg;
                break;
            case 'z':
                args.qvSetSizes = optarg;
                break;
            case 'y':
                args.qvSubnetSizes = optarg;
                break;
            case 'm':
                methods = optarg;
                break;
            case 'j':
                numThreads = std::stoi(optarg);
                break;
            case 'i':
                gtFile = optarg;
                break;
            default:
                std::cerr << "Unknown option(s) provided!\n\n";
                this->showCommandUsage(argv[0]);
                exit(1);
        }
    }

    if (step == "") {
        std::cerr << "No step provided!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    if (step != this->STEP_OBJECTS && step != this->STEP_GROUND_TRUTH && step != this->STEP_VERIFY && step != this->STEP_DUMP) {
        std::cerr << "Invalid step!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    if (argc == 5) {
        // Only the step was provided
        this->showStepUsage(step,argv[0]);
        exit(1);
    }
    if (args.queryType != "knn" && args.queryType != "aknn" && args.queryType != "fkn" && args.queryType != "range") {
        std::cerr << "Invalid or missing query type (-t)!\n\n";
        this->showStepUsage(step,argv[0]);
        exit(1);
    }
    if (step == this->STEP_DUMP) {
        if (gtFile == "") {
            std::cerr << "No ground truth file provided (-i)!\n\n";
            this->showStepUsage(step,argv[0]);
            exit(1);
        }
    } else {
        if (args.bgrFilePath == "" || args.filePathPrefix == "" || args.parameters == "") {
            std::cerr << "Graph file (-g), output path (-f) and parameters (-p) are required!\n\n";
            this->showStepUsage(step,argv[0]);
            exit(1);
        }
        args.params = parseParameters(args.parameters);
    }

    if (args.queryType == "knn") {
        KnnVerifier verifier;
        this->runStep(verifier,step,args,methods,gtFile,numThreads);
    } else if (args.queryType == "aknn") {
        AknnVerifier verifier;
        this->runStep(verifier,step,args,methods,gtFile,numThreads);
    } else if (args.queryType == "fkn") {
        FknVerifier verifier;
        this->runStep(verifier,step,args,methods,gtFile,numThreads);
    } else {
        RangeVerifier verifier;
        this->runStep(verifier,step,args,methods,gtFile,numThreads);
    }
}

void VerifyCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::VERIFY_CMD + " -e <step> -t <query type>\n\n"
              << "Steps:\n"
              << utility::getFormattedUsageString(this->STEP_OBJECTS,"Copy the experiment object sets and their object indexes into <output path>/verification/obj_indexes") + "\n"
              << utility::getFormattedUsageString(this->STEP_GROUND_TRUTH,"Compute ground truth for the query lists of a query type (one SSSP per query vertex)") + "\n"
              << utility::getFormattedUsageString(this->STEP_VERIFY,"Run query methods against the saved ground truth and report mismatches") + "\n"
              << utility::getFormattedUsageString(this->STEP_DUMP,"Print a ground truth file as text") + "\n\n"
              << "Query types: knn, aknn, fkn, range\n";
}

void VerifyCommand::showStepUsage(std::string step, std::string programName)
{
    if (step == this->STEP_OBJECTS) {
        std::cerr << "Usage: " << programName << " -c " + constants::VERIFY_CMD + " -e " + this->STEP_OBJECTS
                  << " -t <query type> -g <binary graph file> -f <output path> -p <parameters> -n <num object sets> -d <object set densities>\n"
                  << "Parameters: the index parameters, the object index flags (rtree=1,nvd=1,colt=1,obj_list=0,occ_list=0) and their\n"
                  << "parameters (rtree_branchfactor, colt_listsize), and <query type>_def_type and <query type>_def_var\n";
    } else if (step == this->STEP_GROUND_TRUTH) {
        std::cerr << "Usage: " << programName << " -c " + constants::VERIFY_CMD + " -e " + this->STEP_GROUND_TRUTH
                  << " -t <query type> -g <binary graph file> -f <output path> -p <parameters>\n"
                  << "-n <num object sets> -d <object set densities> [-q <query vertex file>] [-k <k values>] [-a <range percentages>]\n"
                  << "[-z <AkNN query vertex set sizes>] [-y <AkNN query vertex set subnet sizes>] [-j <threads>]\n"
                  << "Parameters: the index parameters and the <query type>_* query parameters passed to the experiment commands\n";
    } else if (step == this->STEP_VERIFY) {
        std::cerr << "Usage: " << programName << " -c " + constants::VERIFY_CMD + " -e " + this->STEP_VERIFY
                  << " -t <query type> -g <binary graph file> -f <output path> -p <parameters> [-m all|<comma-separated methods>] [-j <threads>]\n"
                  << "Methods default to those enabled in the parameters (e.g. colt_phl=1). If any query fails, the object indexes\n"
                  << "of the failing methods are then compared with their object sets to find stale ones. The queries run in parallel\n"
                  << "(-j threads, all by default) except for the COL-Tree and G-tree methods, whose indexes keep per-query state\n";
    } else if (step == this->STEP_DUMP) {
        std::cerr << "Usage: " << programName << " -c " + constants::VERIFY_CMD + " -e " + this->STEP_DUMP
                  << " -t <query type> -i <ground truth file>\n";
    } else {
        std::cerr << "Invalid step!" << std::endl;
        this->showCommandUsage(programName);
    }
}
