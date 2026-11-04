// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _GRAPH_H
#define _GRAPH_H

#include "MortonCode.h"
#include "../common.h"

#include <vector>
#include <unordered_set>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

/*
 * Directed Graph Stored In-Memory
 */

class Graph {
    
    public:
        Graph() {}
        std::string getNetworkName();
        int getNumVertices();
        int getNumEdges();
        std::string getEdgeType();
        void setGraphData(std::string name, int numVertices, int numEdges, std::string edgeType);
        EdgeWeight getDiameterLowerBound();
        void setDiameterLowerBound(EdgeWeight diameter);
        void insertVertex(VertexID vertexID, Coordinate x, Coordinate y);
        void insertEdge (VertexID sourceID, VertexID targetID, EdgeWeight edgeWeight);
        void parseGraphFile(std::string grFilePath, std::string coFilePath);
        std::vector<VertexID> getVertexIDsVector();
        std::unordered_set<VertexID> getVertexIDsUset();
        Coordinate getMaxCoordinateRange();
        void getMinMaxCoordinates(Coordinate& minX, Coordinate& maxX, Coordinate& minY, Coordinate& maxY);
        void getCoordinates(VertexID u, Coordinate& x, Coordinate& y);
        EuclideanDistance getEuclideanDistance(VertexID vertex1, VertexID vertex2);
        int getEdgeListStartIndex(VertexID vertex);
        int getEdgeListSize(VertexID vertex);
        
        // Graph Validation/Statistics Functions
        bool isConnected(VertexID vertex1, VertexID vertex2);
        bool isUndirectedGraph();
        bool isConnectedGraph();
        bool isContiguous();
        bool hasSameNeighbours(Graph& other, std::string& difference);
        std::vector<int> countOutdegreeFrequencies();
        void printGraphStats();
        double getMaxGraphSpeedByEdge();
        double getMinGraphSpeedByEdge();
        
        // Convert Graph to Other Formats and Output
        void outputToTSVFile(std::string outputFilePath);
        void outputDIMACSFiles(std::string grFilePath, std::string coFilePath);
        Graph createSubGraph(std::string subgraphName, std::unordered_set<VertexID>& subgraphVertices);
        
        // Memory Usage Functions
        double computeIndexSize();
        double computeINEIndexSize();
        
        // For each vertex in the graph (i.e. from 0 to numVertices), this vector
        // has the index of it's first edge in the edges graph. It also has one  The last
        // additional element with the size of the edges vector (so if the graph has
        // n vertices then this vector has n+1 elements).
        // Note: This is to stop us from falling off the edges vector when iterating.
        std::vector<unsigned int> firstEdgeIndex; 
        std::vector<VertexEdgeWeightPair> edges;
        std::vector<CoordinatePair> coordinates; // x is first, y is second (in pair)
        std::vector<MortonCode> mortonCodes;
        
        
    private:
        friend class boost::serialization::access;
        
        std::string networkName = "Unknown";
        unsigned int numVertices;
        unsigned int numEdges;
        std::string type;
        
        Coordinate minXCoord = 0, maxXCoord = 0, minYCoord = 0, maxYCoord = 0;
        int maxRange = 0, maxRangeExp = 0;
        VertexID minXVertex = 0, maxXVertex = 0, minYVertex = 0, maxYVertex = 0;
        double maxGraphSpeedByEdge = 0, minGraphSpeedByEdge = std::numeric_limits<double>::max();
        EdgeWeight diameterLowerBound = 0; // Set when the binary graph is built; range query radii are percentages of it

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->networkName;
            ar & this->numVertices;
            ar & this->numEdges;
            ar & this->type;
            ar & this->firstEdgeIndex;
            ar & this->edges;
            ar & this->coordinates;
            ar & this->minXCoord;
            ar & this->maxXCoord;
            ar & this->minYCoord;
            ar & this->maxYCoord;
            ar & this->minXVertex;
            ar & this->maxXVertex;
            ar & this->minYVertex;
            ar & this->maxYVertex;
            ar & this->maxGraphSpeedByEdge;
            ar & this->minGraphSpeedByEdge;
            ar & this->maxRange;
            ar & this->maxRangeExp;
            ar & this->mortonCodes;
            ar & this->diameterLowerBound;
        }    
        
};

#endif // _GRAPH_H
