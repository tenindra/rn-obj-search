// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "Verifier.h"

#include "../processing/DijkstraSearch.h"
#include "../processing/INE.h"
#include "../processing/PackedRtree.h"
#include "../processing/NVDQuadtree.h"
#include "../utility/utility.h"

#include <algorithm>
#include <cerrno>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

std::string VerifyArgs::param(const std::string& key) const
{
    auto it = this->params.find(key);
    return it == this->params.end() ? "" : it->second;
}

std::string VerifyArgs::requiredParam(const std::string& key) const
{
    std::string value = this->param(key);
    if (value == "") {
        std::cerr << "Parameter " << key << " not provided in -p!" << std::endl;
        exit(1);
    }
    return value;
}

/*
 * Parameters, file names and query lists
 */

// Same format as the experiment commands: key=value pairs separated by commas
std::unordered_map<std::string,std::string> parseParameters(const std::string& parameters)
{
    std::unordered_map<std::string,std::string> parameterMap;
    std::vector<std::string> pairs = utility::splitByDelim(parameters,',');
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        std::vector<std::string> pair = utility::splitByDelim(pairs[i],'=');
        if (pair.size() == 2) {
            parameterMap[pair[0]] = pair[1];
        } else {
            std::cerr << "Invalid key-value pair in parameter string" << std::endl;
            std::cerr << "Parameter string: " << pairs[i] << std::endl;
            exit(1);
        }
    }
    return parameterMap;
}

// Road network index file names, as the experiment commands construct them
static std::string roadIndexPath(const VerifyArgs& args, Graph& graph, const std::string& indexType,
                                 std::vector<std::string> keys, std::vector<std::string> values)
{
    return args.filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),indexType,keys,values);
}

std::string phlIndexPath(const VerifyArgs& args, Graph& graph)
{
    return roadIndexPath(args,graph,constants::IDX_PHL_CMD,{},{});
}

std::string altIndexPath(const VerifyArgs& args, Graph& graph)
{
    return roadIndexPath(args,graph,constants::IDX_ALT_CMD,{"landmarks","type"},
                         {args.requiredParam("alt_numlandmarks"),args.requiredParam("alt_landmarktype")});
}

std::string sltreeIndexPath(const VerifyArgs& args, Graph& graph)
{
    return roadIndexPath(args,graph,constants::IDX_SULTREE_CMD,{"fanout","maxleafsize","numlmks","lmktype"},
                         {args.requiredParam("sltree_fanout"),args.requiredParam("sltree_maxleafsize"),
                          args.requiredParam("sltree_numlmks"),args.requiredParam("sltree_lmktype")});
}

std::string gtreeIndexPath(const VerifyArgs& args, Graph& graph)
{
    return roadIndexPath(args,graph,constants::IDX_GTREE_CMD,{"fanout","maxleafsize"},
                         {args.requiredParam("gtree_fanout"),args.requiredParam("gtree_maxleafsize")});
}

// Object set (.txt) and object index files of the verification snapshot, named as the experiment files (with the index parameters)
std::string objIndexPath(const VerifyArgs& args, const std::string& prefix, const std::string& extension)
{
    return args.objIdxDir() + "/" + utility::constructObjectIndexFileName(prefix,extension,args.params);
}

std::vector<VertexID> loadQueryVertices(const std::string& queryVertexFile)
{
    if (queryVertexFile == "") {
        std::cerr << "No query vertex file provided (-q)!" << std::endl;
        exit(1);
    }
    std::vector<VertexID> queryVertices = utility::getPointSetFromFile(queryVertexFile);
    if (queryVertices.size() == 0) {
        std::cerr << "No query points were provided!" << std::endl;
        exit(1);
    }
    return queryVertices;
}

std::vector<unsigned int> parsePositiveIntegers(const std::string& list, const std::string& description)
{
    std::vector<std::string> strValues = utility::splitByDelim(list,',');
    std::vector<unsigned int> values;
    for (std::size_t i = 0; i < strValues.size(); ++i) {
        unsigned int value = std::stoul(strValues[i]);
        if (value > 0) {
            values.push_back(value);
        } else {
            std::cerr << "Invalid " << description << " " << strValues[i] << " in list provided!" << std::endl;
            exit(1);
        }
    }
    if (values.size() == 0) {
        std::cerr << "No " << description << " list provided!" << std::endl;
        exit(1);
    }
    return values;
}

std::vector<double> parseDensities(const std::string& list)
{
    std::vector<std::string> strDensities = utility::splitByDelim(list,',');
    std::vector<double> densities;
    for (std::size_t i = 0; i < strDensities.size(); ++i) {
        double density = std::stod(strDensities[i]);
        if (density > 0) {
            densities.push_back(density);
        } else {
            std::cerr << "Invalid density in list provided!" << std::endl;
            exit(1);
        }
    }
    if (densities.size() == 0) {
        std::cerr << "No densities provided (-d)!" << std::endl;
        exit(1);
    }
    return densities;
}

// Mirrors the synthetic query lists of LMKExperimentsCommand::runQueries2 and LMKFarthestExperimentsCommand::runQueries2,
// including their switches (a switched-off list is not built)
std::vector<QueryList<KnnQuery>> buildPointQueryLists(const VerifyArgs& args, Graph& graph, const std::string& keyPrefix)
{
    std::vector<QueryList<KnnQuery>> lists;
    std::vector<VertexID> queryVertices = loadQueryVertices(args.queryVertexFile);
    if (args.numSets == 0) {
        std::cerr << "No object sets were provided!" << std::endl;
        exit(1);
    }
    unsigned int defaultK = std::stoul(args.requiredParam(keyPrefix + "_def_k"));
    double defaultDensity = std::stod(args.requiredParam(keyPrefix + "_def_density"));
    std::string defaultObjType = args.requiredParam(keyPrefix + "_def_type");
    int defaultObjVar = std::stoi(args.requiredParam(keyPrefix + "_def_var"));
    std::string networkName = graph.getNetworkName();
    int numVertices = graph.getNumVertices();

    // Varying k
    if (utility::isQueryListOn(args.params,keyPrefix + "_varyingk")) {
        QueryList<KnnQuery> varyingK;
        varyingK.name = "varyingk";
        varyingK.aggFunc = -1;
        std::vector<unsigned int> kValues = parsePositiveIntegers(args.kValues,"k value");
        for (std::size_t i = 0; i < kValues.size(); ++i) {
            for (std::size_t j = 0; j < args.numSets; ++j) {
                std::string prefix = utility::constructObjIdxFilenamePrefix(networkName,defaultObjType,defaultDensity,defaultObjVar,j);
                for (std::size_t k = 0; k < queryVertices.size(); ++k) {
                    varyingK.queries.push_back(KnnQuery(queryVertices[k],kValues[i],defaultDensity,defaultObjType,defaultObjVar,j,prefix,static_cast<int>(numVertices*defaultDensity)));
                }
            }
        }
        lists.push_back(varyingK);
    }

    // Varying density
    if (utility::isQueryListOn(args.params,keyPrefix + "_varyingdensity")) {
        QueryList<KnnQuery> varyingDensity;
        varyingDensity.name = "varyingdensity";
        varyingDensity.aggFunc = -1;
        std::vector<double> densities = parseDensities(args.objDensities);
        for (std::size_t i = 0; i < densities.size(); ++i) {
            for (std::size_t j = 0; j < args.numSets; ++j) {
                std::string prefix = utility::constructObjIdxFilenamePrefix(networkName,defaultObjType,densities[i],defaultObjVar,j);
                for (std::size_t k = 0; k < queryVertices.size(); ++k) {
                    varyingDensity.queries.push_back(KnnQuery(queryVertices[k],defaultK,densities[i],defaultObjType,defaultObjVar,j,prefix,static_cast<int>(numVertices*densities[i])));
                }
            }
        }
        lists.push_back(varyingDensity);
    }

    return lists;
}

/*
 * Distances, result ordering and selection
 */

void computeSSSP(Graph& graph, VertexID source, std::vector<EdgeWeight>& sssp)
{
    // findSSSPDistances leaves vertices it does not reach untouched
    sssp.assign(graph.getNumVertices(),UNREACHED);
    DijkstraSearch dijkstra;
    dijkstra.findSSSPDistances(graph,source,sssp);
}

bool nearestFirst(const VertexDistancePair& a, const VertexDistancePair& b)
{
    return a.second < b.second || (a.second == b.second && a.first < b.first);
}

bool farthestFirst(const VertexDistancePair& a, const VertexDistancePair& b)
{
    return a.second > b.second || (a.second == b.second && a.first < b.first);
}

// Every reachable object with its distance, nearest (or farthest) first
void rankByDistance(const std::vector<EdgeWeight>& sssp, const std::vector<VertexID>& objects, bool farthest, std::vector<VertexDistancePair>& ranked)
{
    ranked.clear();
    ranked.reserve(objects.size());
    for (std::size_t i = 0; i < objects.size(); ++i) {
        if (sssp[objects[i]] != UNREACHED) {
            ranked.push_back(VertexDistancePair(objects[i],sssp[objects[i]]));
        }
    }
    std::sort(ranked.begin(),ranked.end(),farthest ? farthestFirst : nearestFirst);
}

// The first k ranked objects (all of them if there are fewer than k)
void firstK(const std::vector<VertexDistancePair>& ranked, unsigned int k, std::vector<VertexDistancePair>& result)
{
    std::size_t n = std::min(static_cast<std::size_t>(k),ranked.size());
    result.assign(ranked.begin(),ranked.begin()+n);
}

void toSortedPairs(std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists, bool farthest, std::vector<VertexDistancePair>& result)
{
    if (ids.size() != dists.size()) {
        std::cerr << "Number of result objects (" << ids.size() << ") does not match number of distances (" << dists.size() << ")!" << std::endl;
        exit(1);
    }
    result = utility::combineVerticesAndDistanceVectors(ids,dists);
    std::sort(result.begin(),result.end(),farthest ? farthestFirst : nearestFirst);
}

/*
 * Comparison
 */

static std::string idsToString(const std::vector<VertexID>& ids, std::size_t maxIds)
{
    std::stringstream ss;
    for (std::size_t i = 0; i < ids.size() && i < maxIds; ++i) {
        ss << (i == 0 ? "" : ",") << ids[i];
    }
    if (ids.size() > maxIds) {
        ss << ",... (" << ids.size() << " in total)";
    }
    return ss.str();
}

// Rule: same number of results; the sorted distance sequences are identical; the objects are identical at every
// position except positions whose distance equals the last (k-th) distance, where a method may return any object
// tied at that distance. With idSetOnly (range queries) only the object sets are compared, since the range
// methods do not report distances.
bool compareResults(const std::vector<VertexDistancePair>& truth, std::vector<VertexID>& ids, std::vector<EdgeWeight>& dists,
                    bool farthest, bool idSetOnly, std::string& message)
{
    std::stringstream ss;
    if (ids.size() != dists.size()) {
        ss << "method returned " << ids.size() << " objects but " << dists.size() << " distances";
        message = ss.str();
        return false;
    }
    if (ids.size() != truth.size()) {
        ss << "expected " << truth.size() << " results but got " << ids.size();
        message = ss.str();
        return false;
    }
    if (truth.size() == 0) {
        return true;
    }

    if (idSetOnly) {
        std::vector<VertexID> expected, actual = ids;
        for (std::size_t i = 0; i < truth.size(); ++i) {
            expected.push_back(truth[i].first);
        }
        std::sort(expected.begin(),expected.end());
        std::sort(actual.begin(),actual.end());
        for (std::size_t i = 1; i < actual.size(); ++i) {
            if (actual[i] == actual[i-1]) {
                ss << "object " << actual[i] << " returned more than once";
                message = ss.str();
                return false;
            }
        }
        std::vector<VertexID> missing, unexpected;
        std::set_difference(expected.begin(),expected.end(),actual.begin(),actual.end(),std::back_inserter(missing));
        std::set_difference(actual.begin(),actual.end(),expected.begin(),expected.end(),std::back_inserter(unexpected));
        if (missing.size() == 0 && unexpected.size() == 0) {
            return true;
        }
        ss << "result set differs: missing " << idsToString(missing,5) << "; unexpected " << idsToString(unexpected,5);
        message = ss.str();
        return false;
    }

    std::vector<VertexDistancePair> result;
    toSortedPairs(ids,dists,farthest,result);
    for (std::size_t i = 1; i < result.size(); ++i) {
        if (result[i].first == result[i-1].first) {
            ss << "object " << result[i].first << " returned more than once";
            message = ss.str();
            return false;
        }
    }
    for (std::size_t i = 0; i < result.size(); ++i) {
        if (result[i].second != truth[i].second) {
            ss << "distance of result " << i << " is " << result[i].second << " but should be " << truth[i].second;
            message = ss.str();
            return false;
        }
    }
    EdgeWeight boundaryDistance = truth.back().second;
    for (std::size_t i = 0; i < result.size(); ++i) {
        if (truth[i].second != boundaryDistance && result[i].first != truth[i].first) {
            ss << "result " << i << " is object " << result[i].first << " but should be object " << truth[i].first
               << " (distance " << truth[i].second << ")";
            message = ss.str();
            return false;
        }
    }
    return true;
}

std::string pairsToString(const std::vector<VertexDistancePair>& pairs, std::size_t maxPairs)
{
    std::stringstream ss;
    for (std::size_t i = 0; i < pairs.size() && i < maxPairs; ++i) {
        ss << (i == 0 ? "" : " ") << "(" << pairs[i].first << "," << pairs[i].second << ")";
    }
    if (pairs.size() > maxPairs) {
        ss << " ... [+" << pairs.size()-maxPairs << " more]";
    }
    return ss.str();
}

/*
 * The objects step
 */

static void makeDirectory(const std::string& path)
{
    if (mkdir(path.c_str(),0755) != 0 && errno != EEXIST) {
        std::cerr << "Cannot create directory " << path << std::endl;
        exit(1);
    }
}

static void copyFile(const std::string& source, const std::string& destination)
{
    std::ifstream input(source,std::ios::binary);
    if (!input.is_open()) {
        std::cerr << "Cannot open file " << source << std::endl;
        exit(1);
    }
    std::ofstream output(destination,std::ios::binary | std::ios::trunc);
    if (!output.is_open()) {
        std::cerr << "Cannot write file " << destination << std::endl;
        exit(1);
    }
    output << input.rdbuf();
}

// Copies the experiment object sets (every density and set index) into the verification folder, together with the experiment's
// object indexes of the types enabled in the parameters (built by buildDatasetIndexes and buildExperimentIndexes, same file names),
// so the verification checks the index files the experiments load. Every source file is checked first, so a missing one leaves the
// verification folder unchanged
void buildVerificationObjectIndexes(const VerifyArgs& args, Graph& graph)
{
    std::string networkName = graph.getNetworkName();
    std::string objType = args.requiredParam(args.queryType + "_def_type");
    int objVar = std::stoi(args.requiredParam(args.queryType + "_def_var"));
    std::vector<double> densities = parseDensities(args.objDensities);
    if (args.numSets == 0) {
        std::cerr << "No object sets were provided (-n)!" << std::endl;
        exit(1);
    }
    makeDirectory(args.filePathPrefix + "/verification");
    makeDirectory(args.objIdxDir());
    makeDirectory(args.groundTruthDir());

    // Object index types copied with each object set
    std::vector<std::string> indexTypes;
    const std::vector<std::pair<std::string,std::string>> typeFlags = {{"occ_list",constants::OBJ_IDX_GTREE},{"rtree",constants::OBJ_IDX_RTREE},
        {"obj_list",constants::OBJ_IDX_LNDMRK},{"nvd",constants::OBJ_IDX_NVD},{"colt",constants::OBJ_IDX_COLTREE}};
    for (std::size_t i = 0; i < typeFlags.size(); ++i) {
        if (args.param(typeFlags[i].first) == "1") {
            indexTypes.push_back(typeFlags[i].second);
        }
    }

    // Source and destination of every file: each object set and its object indexes
    std::vector<std::pair<std::string,std::string>> files;
    std::size_t numSets = 0;
    for (std::size_t i = 0; i < densities.size(); ++i) {
        const std::size_t totalSets = utility::getNumObjectSets(densities[i],args.numSets);
        for (std::size_t j = 0; j < totalSets; ++j) {
            std::string prefix = utility::constructObjIdxFilenamePrefix(networkName,objType,densities[i],objVar,j);
            files.push_back(std::make_pair(args.filePathPrefix + "/obj_indexes/" + prefix + ".txt",objIndexPath(args,prefix,"txt")));
            for (std::size_t t = 0; t < indexTypes.size(); ++t) {
                files.push_back(std::make_pair(args.filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(prefix,indexTypes[t],args.params),
                                               objIndexPath(args,prefix,indexTypes[t])));
            }
            ++numSets;
        }
    }

    std::size_t numMissing = 0;
    for (std::size_t i = 0; i < files.size(); ++i) {
        if (!std::ifstream(files[i].first).good()) {
            if (numMissing < 5) {
                std::cerr << "Missing " << files[i].first << std::endl;
            }
            ++numMissing;
        }
    }
    if (numMissing > 0) {
        std::cerr << numMissing << " object set or object index file(s) are missing; build the experiment's object indexes first "
                  << "(bash buildDatasetIndexes, then bash buildExperimentIndexes -e " << args.queryType << ")" << std::endl;
        exit(1);
    }

    for (std::size_t i = 0; i < files.size(); ++i) {
        copyFile(files[i].first,files[i].second);
    }
    std::cout << "Copied " << numSets << " object sets and " << files.size()-numSets << " object index files to " << args.objIdxDir() << std::endl;
}

/*
 * VerifierBase
 */

void VerifierBase::loadObjectSets(const VerifyArgs& args, Graph& graph, const std::vector<std::string>& prefixes, bool buildStatus)
{
    INE ine;
    std::string setType;
    double setDensity;
    int setSize, setVariable;
    for (std::size_t i = 0; i < prefixes.size(); ++i) {
        if (this->objectSets.find(prefixes[i]) == this->objectSets.end()) {
            std::string file = objIndexPath(args,prefixes[i],"txt");
            std::ifstream test(file);
            if (!test.is_open()) {
                std::cerr << "Object set file " << file << " not found (run the objects step first)!" << std::endl;
                exit(1);
            }
            std::vector<VertexID> objectVertices = utility::getPointSetFromFile(file,setType,setDensity,setSize,setVariable);
            if (static_cast<int>(objectVertices.size()) != setSize) {
                std::cerr << "Object set file " << file << " has " << objectVertices.size() << " objects but its header says " << setSize << "!" << std::endl;
                exit(1);
            }
            this->objectSets[prefixes[i]] = objectVertices;
            this->objectSetSizes[prefixes[i]] = setSize;
        }
        if (buildStatus && this->objectStatus.find(prefixes[i]) == this->objectStatus.end()) {
            this->objectStatus[prefixes[i]] = ine.buildVertexObjectStatusVector(graph,this->objectSets[prefixes[i]]);
        }
    }
}

void VerifierBase::loadPHL(const VerifyArgs& args, Graph& graph)
{
    if (!this->phlLoaded) {
        std::string file = phlIndexPath(args,graph);
        std::ifstream test(file);
        if (!test.is_open()) {
            std::cerr << "Cannot open PHL index file " << file << std::endl;
            exit(1);
        }
        this->phl.LoadLabel(file.c_str());
        this->phlLoaded = true;
    }
}

void VerifierBase::loadALT(const VerifyArgs& args, Graph& graph)
{
    if (!this->altLoaded) {
        this->alt = serialization::getIndexFromBinaryFile<ALT>(altIndexPath(args,graph));
        this->altLoaded = true;
    }
}

void VerifierBase::loadSULTree(const VerifyArgs& args, Graph& graph)
{
    if (!this->sltreeLoaded) {
        this->sltree = serialization::getIndexFromBinaryFile<SULTree>(sltreeIndexPath(args,graph));
        this->sltreeLoaded = true;
    }
}

void VerifierBase::loadGtree(const VerifyArgs& args, Graph& graph)
{
    if (!this->gtreeLoaded) {
        this->gtree = serialization::getIndexFromBinaryFile<Gtree>(gtreeIndexPath(args,graph));
        this->gtreeLoaded = true;
    }
}

/*
 * Stale object indexes
 */

// The object index a method queries, as the verifiers' loadIndexes load it ("" for the methods that read the object set files)
static std::string methodObjectIndex(const std::string& method)
{
    if (method == "gtree") {
        return constants::OBJ_IDX_GTREE;
    } else if (method == "ier_phl") {
        return constants::OBJ_IDX_RTREE;
    } else if (method == "ol_phl") {
        return constants::OBJ_IDX_LNDMRK;
    } else if (method == "nvd_phl") {
        return constants::OBJ_IDX_NVD;
    } else if (method == "colt_phl") {
        return constants::OBJ_IDX_COLTREE;
    }
    return "";
}

// The objects an object index holds, in any order
static std::vector<VertexID> indexObjects(OccurenceList& occList)
{
    // Every object is in the occurence list of its G-tree leaf
    std::vector<VertexID> objects;
    for (auto it = occList.leafOccurenceList.begin(); it != occList.leafOccurenceList.end(); ++it) {
        objects.insert(objects.end(),it->second.begin(),it->second.end());
    }
    return objects;
}

static std::vector<VertexID> indexObjects(PackedRtree& rtree)
{
    // A window covering every coordinate returns every leaf entry
    Rectangle everything;
    everything.setDimensions(std::numeric_limits<Coordinate>::min(),std::numeric_limits<Coordinate>::min(),
                             std::numeric_limits<Coordinate>::max(),std::numeric_limits<Coordinate>::max());
    return rtree.windowQuery(everything);
}

static std::vector<VertexID> indexObjects(ObjectList& objList)
{
    // Every landmark's list holds every object
    std::vector<VertexID> objects;
    if (objList.lmObjDistances.size() > 0) {
        for (std::size_t i = 0; i < objList.lmObjDistances[0].size(); ++i) {
            objects.push_back(objList.lmObjDistances[0][i].first);
        }
    }
    return objects;
}

static std::vector<VertexID> indexObjects(NVDQuadtree& nvd)
{
    return nvd.vnsObjectIDs; // One dominance region per object
}

static std::vector<VertexID> indexObjects(COLTree& colt)
{
    // Every object is in the object lists of one node, and each landmark's list of that node holds all of its objects
    std::vector<VertexID> objects;
    for (std::size_t i = 0; i < colt.treeNodes.size(); ++i) {
        COLTreeNode& node = colt.treeNodes[i];
        if (node.hasObjectList() && node.lmkObjectLists.size() > 0) {
            for (std::size_t j = 0; j < node.lmkObjectLists[0].size(); ++j) {
                objects.push_back(node.lmkObjectLists[0][j].first);
            }
        }
    }
    return objects;
}

// Loads the index of every object set, one at a time, and compares its objects with the set file; prints the indexes that differ
// and returns their number
template <typename Index>
static std::size_t checkIndexType(const VerifyArgs& args, const std::string& method, const std::string& extension, const std::vector<std::string>& prefixes,
                                  const std::unordered_map<std::string,std::vector<VertexID>>& objectSets)
{
    const std::size_t maxListed = 10;
    std::vector<std::string> stale;
    for (std::size_t i = 0; i < prefixes.size(); ++i) {
        Index index = serialization::getIndexFromBinaryFile<Index>(objIndexPath(args,prefixes[i],extension));
        std::vector<VertexID> actual = indexObjects(index), expected = objectSets.at(prefixes[i]);
        std::sort(actual.begin(),actual.end());
        std::sort(expected.begin(),expected.end());
        if (actual != expected) {
            std::vector<VertexID> missing, unexpected;
            std::set_difference(expected.begin(),expected.end(),actual.begin(),actual.end(),std::back_inserter(missing));
            std::set_difference(actual.begin(),actual.end(),expected.begin(),expected.end(),std::back_inserter(unexpected));
            std::stringstream ss;
            ss << prefixes[i] << ": " << missing.size() << " of the set's " << expected.size() << " objects missing, "
               << unexpected.size() << " objects not in the set";
            stale.push_back(ss.str());
        }
    }
    if (stale.size() == 0) {
        std::cout << method << ": the " << prefixes.size() << " object indexes (" << extension << ") hold exactly the objects of their object sets" << std::endl;
    } else {
        std::cout << method << ": " << stale.size() << " of " << prefixes.size() << " object indexes (" << extension
                  << ") do not hold exactly the objects of their object set (stale):" << std::endl;
        for (std::size_t i = 0; i < stale.size() && i < maxListed; ++i) {
            std::cout << "  " << stale[i] << std::endl;
        }
        if (stale.size() > maxListed) {
            std::cout << "  ... and " << stale.size()-maxListed << " more" << std::endl;
        }
    }
    return stale.size();
}

void VerifierBase::checkObjectIndexes(const VerifyArgs& args, Graph& graph, const std::vector<std::string>& methods,
                                      const std::vector<std::string>& prefixes)
{
    std::cout << "\nChecking the object indexes of the methods with failures against their object sets..." << std::endl;
    this->loadObjectSets(args,graph,prefixes,false);
    std::size_t numChecked = 0, numStale = 0;
    for (std::size_t m = 0; m < methods.size(); ++m) {
        std::string extension = methodObjectIndex(methods[m]);
        if (extension == constants::OBJ_IDX_GTREE) {
            numStale += checkIndexType<OccurenceList>(args,methods[m],extension,prefixes,this->objectSets);
        } else if (extension == constants::OBJ_IDX_RTREE) {
            numStale += checkIndexType<PackedRtree>(args,methods[m],extension,prefixes,this->objectSets);
        } else if (extension == constants::OBJ_IDX_LNDMRK) {
            numStale += checkIndexType<ObjectList>(args,methods[m],extension,prefixes,this->objectSets);
        } else if (extension == constants::OBJ_IDX_NVD) {
            numStale += checkIndexType<NVDQuadtree>(args,methods[m],extension,prefixes,this->objectSets);
        } else if (extension == constants::OBJ_IDX_COLTREE) {
            numStale += checkIndexType<COLTree>(args,methods[m],extension,prefixes,this->objectSets);
        } else {
            std::cout << methods[m] << ": reads the object set files, no object index to check" << std::endl;
            continue;
        }
        numChecked += prefixes.size();
    }
    if (numStale > 0) {
        std::cout << numStale << " stale object index(es): rebuild them (bash buildDatasetIndexes, then bash buildExperimentIndexes -e "
                  << args.queryType << "), then re-run the objects and verify steps" << std::endl;
    } else if (numChecked > 0) {
        std::cout << "All " << numChecked << " checked object indexes hold exactly the objects of their object sets" << std::endl;
    }
}
