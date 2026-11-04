// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _SETGENERATOR_H
#define _SETGENERATOR_H

#include "Graph.h"

#include <string>
#include <vector>
#include <random>

class SetGenerator {

    public:
        std::vector<VertexID> generateRandomSampleSet(Graph& graph, unsigned int numObjects, std::mt19937& engine, int status = -1, int outdegree = 0);
        std::vector<VertexID> generateRandomSampleSet(int numVertices, unsigned int numObjects, std::mt19937& engine);

        // The engine of one set file, seeded with the set seed, the kind of set and the parameters in the file's name, so that the
        // file depends only on these and the graph, not on the order of generation or the other sets generated with it
        static std::mt19937 getEngine(unsigned int setSeed, unsigned int setKind, std::vector<unsigned int> setParameters);
        // The set seed given with seedOption, which every command that generates sets requires
        static unsigned int parseSeed(std::string seedValue, std::string seedOption);

        // Kinds of set for getEngine, so that sets of different kinds with the same parameters draw different vertices
        static const unsigned int SAMPLE_SET = 1; // sample_set
        static const unsigned int VERTEX_VECTOR_SET = 2; // vertex_vec -t random
        static const unsigned int SUBNET_VECTOR_SET = 3; // vertex_vec -t subnet
        static const unsigned int OBJECT_SET = 4; // Random object sets of lmk_exp -e objects

};

#endif // _SETGENERATOR_H
