// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _COMMON_H
#define _COMMON_H

#include "debug.h"

#include <string>
#include <cstdint>

typedef unsigned int EdgeWeight;
typedef unsigned int VertexID;
typedef VertexID TreeNodeID; // Same type, as the COL-Tree and G-tree queries queue tree node indexes as IDs past the last vertex
typedef int Coordinate;
typedef double EuclideanDistance;
typedef double EuclideanDistanceSqrd;
typedef unsigned long long MortonNumber;
typedef uint8_t MortonLength;
typedef unsigned long long LongPathDistance;
typedef std::pair<VertexID,EdgeWeight> VertexDistancePair;
typedef std::pair<VertexID,EdgeWeight> VertexEdgeWeightPair;
typedef std::pair<VertexID,VertexID> VertexPair;
typedef std::pair<TreeNodeID,TreeNodeID> TreeNodePair;
typedef std::pair<Coordinate,Coordinate> CoordinatePair; // x is first, y is second
typedef std::pair<int,EdgeWeight> IntDistancePair;
typedef std::pair<int,VertexID> IntVertexPair;

namespace constants
{
    // Command Names (must be less than constants::MAX_COMMAND_NAME_SIZE)
    std::string const IDX_GRAPH_CMD = "graph";
    std::string const IDX_GTREE_CMD = "gtree";
    std::string const IDX_PHL_CMD = "phl";
    std::string const IDX_ALT_CMD = "alt";
    std::string const IDX_SULTREE_CMD = "sltree";
    std::string const IDX_PARTITION_CMD = "partition"; // Hierarchical partition the SUL-Tree is built from (index file extension)
    std::string const SAMPLE_SET_CMD = "sample_set";
    std::string const TRANSFORM_INPUT_CMD = "transform";
    std::string const LMK_EXPERIMENTS_CMD = "lmk_exp";
    std::string const LMK_AGG_EXPERIMENTS_CMD = "lmk_agg_exp";
    std::string const LMK_FARTHEST_EXPERIMENTS_CMD = "lmk_farthest_exp";
    std::string const LMK_RANGE_EXPERIMENTS_CMD = "lmk_range_exp";
    std::string const VERTEX_VECTOR_CMD = "vertex_vec";
    std::string const VERIFY_CMD = "verify";

    // Point/Object Set Types
    std::string const RAND_OBJ_SET = "random";

    // Object Index Types
    std::string const OBJ_IDX_GTREE = "occ_list";
    std::string const OBJ_IDX_RTREE = "rtree";
    std::string const OBJ_IDX_LNDMRK = "obj_list";
    std::string const OBJ_IDX_NVD = "nvd";
    std::string const OBJ_IDX_COLTREE = "colt";

    // Query Methods
    std::string const INE_KNN_QUERY = "ine";
    std::string const GTREE_KNN_QUERY = "gtree";
    std::string const IER_PHL_KNN_QUERY = "ier_phl";
    
    // Landmark Methods
    std::string const OL_ILBR_QLM_PHL_KNN_QUERY = "ol_qlm_ilbr_phl";
    std::string const NVD_ILBR_PHL_KNN_QUERY = "nvd_ilbr_phl";
    std::string const LLB_ILBR_PHL_KNN_QUERY = "llb_ilbr_phl";
    std::string const COLTREE_ILBR_PHL_KNN_QUERY = "colt_ilbr_phl";
        
    // Constants
    int const DOUBLE_SWEEP_ITERATIONS = 5;
    int const MAX_COMMAND_NAME_SIZE = 20;
    int const WARM_UP_QUERY_PERCENTAGE = 10; // Share of the first query list run untimed before timing (at least one query)
    const int ROOT_PARENT_INDEX = -1;
    const VertexID UNUSED_VERTEX_ID = 100000000; // Largest VertexID is size largst graph (i.e. US with 25 million vertices)
    
    // Aggregate Functions
    const int AGG_SUM = 0;
    const int AGG_MAX = 1;

    // Other String Constants
    std::string const DISTANCE_WEIGHTS = "d";
    std::string const TIME_WEIGHTS = "t";
    std::string const VERTEX_EXT = "coordinates";
    std::string const EDGE_EXT = "graph";
}

// OpenMP pragmas compiled only into rn_obj_search_par, which is built with PARALLEL_BUILDS for the parameter-test index builds
// (CMakeLists.txt). In rn_obj_search and rn_obj_search_stats they expand to nothing, so the code those binaries time and count is unchanged
#if defined(PARALLEL_BUILDS)
#define PARALLEL_BUILDS_PRAGMA(...) _Pragma(#__VA_ARGS__)
#else
#define PARALLEL_BUILDS_PRAGMA(...)
#endif

#endif // _COMMON_H