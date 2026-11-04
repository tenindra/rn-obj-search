// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _TRANSFORMINPUTCOMMAND_H
#define _TRANSFORMINPUTCOMMAND_H

#include "Command.h"
#include "../common.h"

#include <vector>
#include <unordered_set>

class TransformInputCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);
        void showMethodUsage(std::string method, std::string programName);
        
        // Input Types for Transform Command
        std::string const TRANSFORM_DIMACS = "dimacs";
        std::string const TRANSFORM_DEM = "dem";
        std::string const TRANSFORM_METIS = "metis";
    
    private:        
        void transformDIMACSFiles(std::string grFilePath, std::string coFilePath, std::string outputFilePrefix);
        // Terrain surface network from a DEM tile and a finite-element mesh from METIS files
        void transformDEMFile(std::string hgtFilePath, std::string outputFilePrefix, std::string name, std::string cropSpec, double latitude, long long unitsPerMetre);
        void transformMETISFiles(std::string grFilePath, std::string coFilePath, std::string outputFilePrefix, std::string name, std::string cropSpec, long long scale);
        void finishTransformation(std::string outputFilePrefix);
        void printGeometryStatistics();
        void addVertex(VertexID vertex, Coordinate x, Coordinate y);
        void addEdge(VertexID source, VertexID target, EdgeWeight weight);
        bool isValidEdge(VertexID source, VertexID target, int rawWeight);
        void checkValidVertex(VertexID vertex);
        EdgeWeight correctEdgeWeightByType(VertexID source, VertexID target, EdgeWeight weight);
        void fixShortestPathOverflow();
        LongPathDistance computePseudoDiameter(int numRounds);
        void computeCoordinateRanges();
        void checkNonUniqueCoordinates();
        void outputStandardFormatFiles(std::string outputFilePrefix);
        void printTransformationStatistics();
        void removeNonUniqueVertices();
        void removeIsolatedVertices();
        void removeIsolatedRegions();
        
        // Graph Information
        unsigned int numVertices = 0, numEdges = 0;
        std::string networkName = "", edgeType = "";
        
        // Error Detected/Corrected Statistics
        int numSelfLoops = 0, numParallelEdges = 0, numNegativeEdgeWeights = 0, numEuclideanCorrections = 0, numCorrectionsOverFactorTen = 0;
        int numNonUniqueVerticesRemoved = 0, numDisconnectedVerticesRemoved = 0, numNonUniqueVertexEdgesRemoved = 0, numDisconnectedVertexEdgesRemoved = 0;
        int numRegionsFound = 0, numIsolatedRegionVerticesRemoved = 0, numIsolatedRegionEdgesRemoved = 0;
        unsigned long long totalEdgeWeight = 0, weightAddedByEuclideanCorrections = 0;
        int xTranslation = 0, yTranslation = 0;
        bool wasDeflated = false;
        int deflationFactor = 0;
        unsigned long long oldTotalEdgeWeight = 0;
        double pseudoDiameterTimeMs = 0;

        // Graph
        std::vector<VertexID> vertices;
        std::vector<std::vector<VertexEdgeWeightPair>> neighbours;
        std::vector<CoordinatePair> coordinates;
        std::unordered_set<VertexID> nonUniqueVertices;
        
};

#endif // _TRANSFORMINPUTCOMMAND_H
