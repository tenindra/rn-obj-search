// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "Graph.h"

#include "../utility/geometry.h"

#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <regex>
#include <cstdlib>
#include <limits>
#include <unordered_map>
#include <cmath>

int Graph::getNumVertices() {
    return this->numVertices;
}

int Graph::getNumEdges()
{
    return this->numEdges;
}

std::string Graph::getNetworkName()
{
    return this->networkName;
}

std::string Graph::getEdgeType()
{
    return this->type;
}

EdgeWeight Graph::getDiameterLowerBound()
{
    return this->diameterLowerBound;
}

void Graph::setDiameterLowerBound(EdgeWeight diameter)
{
    this->diameterLowerBound = diameter;
}

void Graph::setGraphData(std::string name, int numVertices, int numEdges, std::string edgeType)
{
    this->networkName = name;
    this->numVertices = numVertices;
    this->numEdges = numEdges;
    this->type = edgeType;
}

void Graph::insertVertex(VertexID vertexID, Coordinate x, Coordinate y)
{
    if (vertexID != this->coordinates.size()) {
        std::cerr << "Graph " << this->networkName << ": vertex " << vertexID << " is out of order in the coordinates file (expected vertex "
                  << this->coordinates.size() << ")" << std::endl;
        exit(1);
    }
    
    this->coordinates.push_back(CoordinatePair(x,y));
    
    // Store the minimum and maximum X and Y coordinates values
    // This is needed for quadtree contruction of this graph
    if (vertexID == 0) {
        // This is the first vertex, so we need to initialize min/max coordinates
        this->minXCoord = x;
        this->maxXCoord = x;
        this->minYCoord = y;
        this->maxYCoord = y;
        this->minXVertex = vertexID;
        this->maxXVertex = vertexID;
        this->minYVertex = vertexID;
        this->maxYVertex = vertexID;
    } else {
        if (x < this->minXCoord) {
            this->minXCoord = x;
            this->minXVertex = vertexID;
        }
        if (x > this->maxXCoord) {
            this->maxXCoord = x;
            this->maxXVertex = vertexID;
        }
        if (y < this->minYCoord) {
            this->minYCoord = y;
            this->minYVertex = vertexID;
        }
        if (y > this->maxYCoord) {
            this->maxYCoord = y;
            this->maxYVertex = vertexID;
        }
    }
}

// We assume edges are added in order (i.e. all edges for sourceID are added together
void Graph::insertEdge(VertexID sourceID, VertexID targetID, EdgeWeight edgeWeight) {
    if (sourceID == targetID) {
        std::cerr << "Graph " << this->networkName << ": edge " << sourceID << " -> " << targetID << " is a self-loop" << std::endl;
        exit(1);
    }
    if (this->type == constants::DISTANCE_WEIGHTS) {
        EdgeWeight euclidDist = std::ceil(this->getEuclideanDistance(sourceID,targetID));
        if (edgeWeight < euclidDist) {
            std::cerr << "Graph " << this->networkName << ": edge " << sourceID << " -> " << targetID << " has weight " << edgeWeight
                      << ", less than the Euclidean distance " << euclidDist << " between its vertices" << std::endl;
            exit(1);
        }
    } else if (this->type == constants::TIME_WEIGHTS) {
        EuclideanDistance euclidDist = this->getEuclideanDistance(sourceID,targetID);
        double edgeSpeed = euclidDist/edgeWeight;
        if (edgeSpeed > this->maxGraphSpeedByEdge) {
            this->maxGraphSpeedByEdge = edgeSpeed;
        }
        if (edgeSpeed < this->minGraphSpeedByEdge) {
            this->minGraphSpeedByEdge = edgeSpeed;
        }
    }
    if (edgeWeight == 0) {
        std::cerr << "Graph " << this->networkName << ": edge " << sourceID << " -> " << targetID << " has weight 0" << std::endl;
        exit(1);
    }
    
    if (sourceID == this->firstEdgeIndex.size()) {
        // Then this is the first edge encountered
        this->firstEdgeIndex.push_back(this->edges.size());
    } else {
        if (sourceID != this->firstEdgeIndex.size()-1) {
            std::cerr << "Graph " << this->networkName << ": edge " << sourceID << " -> " << targetID << " is out of order in the graph file"
                      << " (the edges must be sorted by source vertex, with at least one edge for every vertex from 0)" << std::endl;
            exit(1);
        }
        // Check if edge has already been
        int adjListStart = this->firstEdgeIndex[sourceID], adjListSize = this->edges.size();
        for (int i = adjListStart; i < adjListSize; ++i) {
            if (targetID == this->edges[i].first) {
                std::cerr << "Graph " << this->networkName << ": edge " << sourceID << " -> " << targetID << " appears twice in the graph file" << std::endl;
                exit(1);
            }
        }
    }
    this->edges.push_back(VertexEdgeWeightPair(targetID,edgeWeight));
}

void Graph::parseGraphFile(std::string grFilePath, std::string coFilePath)
{
    VertexID sourceID, targetID;
    EdgeWeight weight;
    std::string name, type;
    unsigned int vertices, edges;
    Coordinate x, y;

    // Parse Coordinates File
    // We expect all coordinates to be numbered from 0 to n-1 and are supplied in order
    std::ifstream coordFile(coFilePath, std::ios::in);
    if (coordFile.is_open()) {
        coordFile >> name >> type >> vertices;
        this->networkName = name;
        this->numVertices = vertices;
            
        std::cout << "Now parsing " << vertices << " vertices for graph " << this->networkName << std::endl;

        while (coordFile >> sourceID >> x >> y)  { // x is longitude and y is latitude (both multiplied by 1 million)
           this->insertVertex(sourceID,x,y);
        }
        if (this->coordinates.size() != static_cast<unsigned int>(this->getNumVertices())) {
            std::cerr << "Coordinates file " << coFilePath << " has " << this->coordinates.size() << " vertices, but its header says "
                      << this->getNumVertices() << std::endl;
            exit(1);
        }
    } else {
        std::cerr << "Cannot open coordinates file " << coFilePath << std::endl;
        exit(1);
    }    
    coordFile.close();
    
    // Parse Graph File
    // We expect all edges to be supplied in order of the source ID (from 0 to n-1)
    std::ifstream graphFS(grFilePath, std::ios::in);
    if (graphFS.is_open()) {
        // First line is graph information
        graphFS >> name >> type >> vertices >> edges;
        this->numEdges = edges;
        this->type = type; // Edges file sets type not coordinates
            
        std::cout << "Now parsing " << edges << " edges for graph " << this->networkName << std::endl;

        while (graphFS >> sourceID >> targetID >> weight)  {
            this->insertEdge(sourceID,targetID,weight);
        }
        // The last element firstEdgeIndex is the size of the edges vector
        this->firstEdgeIndex.push_back(this->edges.size());
        if (this->edges.size() != static_cast<unsigned int>(this->getNumEdges())) {
            std::cerr << "Graph file " << grFilePath << " has " << this->edges.size() << " edges, but its header says "
                      << this->getNumEdges() << std::endl;
            exit(1);
        }
        if (this->firstEdgeIndex.size() != static_cast<unsigned int>(this->getNumVertices()+1)) {
            std::cerr << "Graph file " << grFilePath << " has edges for " << this->firstEdgeIndex.size()-1 << " vertices, but the coordinates file has "
                      << this->getNumVertices() << std::endl;
            exit(1);
        }
    } else {
        std::cerr << "Cannot open graph file " << grFilePath << std::endl;
        exit(1);
    }
    graphFS.close();
    
    // The loop below doubles maxRange until it exceeds the coordinate range, and each coordinate fills half of a vertex's Morton
    // code, so the range must be below the largest power of two that both allow (or the loop would not end)
    long long coordinateSpan = std::max(static_cast<long long>(this->maxXCoord)-this->minXCoord,static_cast<long long>(this->maxYCoord)-this->minYCoord);
    long long coordinateSpanLimit = 1LL << std::min(std::numeric_limits<decltype(this->maxRange)>::digits-1,quadtree_constants::MAX_MORTON_LENGTH/2);
    if (coordinateSpan >= coordinateSpanLimit) {
        std::cerr << "The coordinates in " << coFilePath << " span " << coordinateSpan << " in x or y, but must span less than "
                  << coordinateSpanLimit << " for the vertices' Morton codes" << std::endl;
        exit(1);
    }
    Coordinate maxCoordinateRange = this->getMaxCoordinateRange();

    // Find number bits needed to represent largest coordinate
    this->maxRange = 1; // 2^0
    this->maxRangeExp = 0;
    // The biggest number that maxRangeExp bits can reprsent is actually
    // maxRange-1. So if maxCoordinateRange is equal to maxRange, then at
    // maxRangeExp does have no bits to represent maxCoordinateRange so
    // we must continue loop in this case and hence the <= sign
    while (this->maxRange <= maxCoordinateRange) {
        this->maxRange = this->maxRange << 1;
        this->maxRangeExp++;
    }
    this->maxRange -= 1; 
    // Worst case there will be a vertex with x and y coordinates with 
    // maxRangeExp bits. This means when we interleave x and y bits,
    // in the worst case the max morton code length is 2*maxrangeExp
    
    this->mortonCodes.reserve(this->numVertices);
    for (VertexID vertex = 0; vertex < this->numVertices; ++vertex) {
        this->getCoordinates(vertex,x,y);
        // We translate all X and Y coordinates such that (0,0) is the most
        // SW point and coordinates increase in the N and E directions
        x = x - this->minXCoord; // Translate right, e.g. -100 -(-100)
        y = y - this->minYCoord; 
        // Note: In order to have (0,0) be in the SW corner we adjust the X 
        // and Y coordinates. If coordinate is negative translate right
        // and if coordinate is position when translate left.

        // Note: We store coordinates in the quadtree so that we don't require graph
        // again. We also pre-compute morton codes to make it easier to point location.
        // Both of these help improve MortonList-based index building efficiency
        MortonCode code(x,y,maxRangeExp);
        this->mortonCodes.push_back(code);
    }
    
    std::cout << "\nSuccessfully parsed graph " << this->networkName << " with " << this->numVertices
                << " vertices and " << this->numEdges << " edges!" << std::endl;
}

std::vector<VertexID> Graph::getVertexIDsVector()
{
    std::vector<VertexID> vertexIDs;
    for (VertexID i = 0; i < this->numVertices; ++i) {
        vertexIDs.push_back(i);
    }
    return vertexIDs;
}

std::unordered_set<VertexID> Graph::getVertexIDsUset()
{
    std::unordered_set<VertexID> vertexIDs;
    for (VertexID i = 0; i < this->numVertices; ++i) {
        vertexIDs.insert(i);
    }
    return vertexIDs;
}

Coordinate Graph::getMaxCoordinateRange()
{
    // Return the largest range in coordinate values
    Coordinate xRange = this->maxXCoord-this->minXCoord, yRange = this->maxYCoord-this->minYCoord;
    return (xRange > yRange) ? xRange : yRange;
}
void Graph::getMinMaxCoordinates(Coordinate &minX, Coordinate &maxX, Coordinate &minY, Coordinate &maxY) {
    minX = this->minXCoord;
    maxX = this->maxXCoord;
    minY = this->minYCoord;
    maxY = this->maxYCoord;
}

void Graph::getCoordinates(VertexID u, Coordinate &x, Coordinate &y) {
    x = this->coordinates[u].first;
    y = this->coordinates[u].second;
}

EuclideanDistance Graph::getEuclideanDistance(VertexID vertex1, VertexID vertex2)
{
    return geometry::getEuclideanDist(this->coordinates[vertex1].first,this->coordinates[vertex1].second,
                                      this->coordinates[vertex2].first,this->coordinates[vertex2].second);
}

int Graph::getEdgeListStartIndex(VertexID vertex)
{
    return this->firstEdgeIndex[vertex];
}

int Graph::getEdgeListSize(VertexID vertex)
{
    // Since firstEdgeIndex has numVertices+1
    // element as long as vertex is a valid vertex
    // (i.e. less than numVertices) then this
    // will access a valid element of firstEdgeIndex
    return this->firstEdgeIndex[vertex+1];
}

/*
 * Graph Validation/Statistics Functions
 */

// Is vertex2 one of vertex1's neighbours
bool Graph::isConnected(VertexID vertex1, VertexID vertex2)
{
    int adjListStart, nextAdjListStart;
    adjListStart = this->getEdgeListStartIndex(vertex1);
    nextAdjListStart = this->getEdgeListSize(vertex1);
    // Note: We have to make sure we don't exceed size of graph.edges vector
    
    // Iterate over each neighbour of vertex1
    for (int i = adjListStart; i < nextAdjListStart; ++i) {
        if (this->edges[i].first == vertex2) {
            return true;
        }
    }
    // Could not find vertex2 in neighbour list
    return false;
}

bool Graph::isUndirectedGraph()
{
    bool undirected = true;
    int adjListStart, nextAdjListStart;
    for (VertexID i = 0; i < this->numVertices; ++i) {
        adjListStart = this->getEdgeListStartIndex(i);
        nextAdjListStart = this->getEdgeListSize(i);
        // Note: We have to make sure we don't exceed size of graph.edges vector
        
        // Iterate over each neighbour of current vertex i
        for (int j = adjListStart; j < nextAdjListStart; ++j) {
            if (!this->isConnected(this->edges[j].first,i)) {
                // If the reverse edge doesn't exist then graph is directed
                undirected = false;
            }
        }
    }
    return undirected;
}

bool Graph::isConnectedGraph()
{
    bool connected = true;
    for (VertexID i = 0; i < this->numVertices; ++i) {
        if (this->firstEdgeIndex[i] == this->firstEdgeIndex[i+1]) {
            connected = false;
            std::cerr << "Vertex " << i << " has no edges!" << std::endl;
        }
    }
    return connected;
}

bool Graph::isContiguous()
{
    VertexID currentVertex, neighbourVertexID;
    int adjListStart, nextAdjListStart;
    std::vector<bool> isVertexFound(this->numVertices,false);
    
    unsigned int regionSize = 0;
    std::deque<VertexID> pqueue;
    pqueue.push_back(0);
    while (pqueue.size() > 0) {
        currentVertex = pqueue.front();
        pqueue.pop_front();
        if (!isVertexFound[currentVertex]) {
            isVertexFound[currentVertex] = true;
            ++regionSize;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = this->getEdgeListStartIndex(currentVertex);
            nextAdjListStart = this->getEdgeListSize(currentVertex);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                neighbourVertexID = this->edges[i].first;
                // Only add those we haven't already seen
                if (!isVertexFound[neighbourVertexID]) {
                    pqueue.push_back(neighbourVertexID);
                }
            }

        }
    }
    
    if (regionSize == this->numVertices) {
        return true;
    } else {
        return false;
    }
}

// Does other have the same vertex IDs and each vertex the same neighbours (in any order, whatever the edge weights), so that a set
// of vertices drawn from one graph is also valid on the other; if not, difference describes the first vertex that differs
bool Graph::hasSameNeighbours(Graph& other, std::string& difference)
{
    if (this->numVertices != other.numVertices) {
        difference = this->networkName + " has " + std::to_string(this->numVertices) + " vertices, "
                     + other.networkName + " has " + std::to_string(other.numVertices);
        return false;
    }

    std::vector<VertexID> neighbours, otherNeighbours;
    for (VertexID i = 0; i < this->numVertices; ++i) {
        neighbours.clear();
        for (unsigned int j = this->firstEdgeIndex[i]; j < this->firstEdgeIndex[i+1]; ++j) {
            neighbours.push_back(this->edges[j].first);
        }
        otherNeighbours.clear();
        for (unsigned int j = other.firstEdgeIndex[i]; j < other.firstEdgeIndex[i+1]; ++j) {
            otherNeighbours.push_back(other.edges[j].first);
        }
        std::sort(neighbours.begin(),neighbours.end());
        neighbours.erase(std::unique(neighbours.begin(),neighbours.end()),neighbours.end());
        std::sort(otherNeighbours.begin(),otherNeighbours.end());
        otherNeighbours.erase(std::unique(otherNeighbours.begin(),otherNeighbours.end()),otherNeighbours.end());

        if (neighbours != otherNeighbours) {
            difference = "vertex " + std::to_string(i) + " has neighbours";
            for (VertexID neighbour : neighbours) {
                difference += " " + std::to_string(neighbour);
            }
            difference += " in " + this->networkName + " and";
            for (VertexID neighbour : otherNeighbours) {
                difference += " " + std::to_string(neighbour);
            }
            difference += " in " + other.networkName;
            return false;
        }
    }
    return true;
}

std::vector<int> Graph::countOutdegreeFrequencies()
{
    std::vector<int> frequencies;
    unsigned int outdegree;
    [[maybe_unused]] int numCountedEdges = 0;
    for (VertexID i = 0; i < this->numVertices; ++i) {
        // Outdegree for vertex i-1 (not vertex i)
        outdegree = this->firstEdgeIndex[i+1]-firstEdgeIndex[i];
        if (frequencies.size() <= outdegree) {
            // Note: resize doesn't change the values
            // of existing values (unlike assign)
            frequencies.resize(outdegree+1,0);
        }
        ++frequencies[outdegree];
        numCountedEdges += outdegree;
    }
    RN_DEBUG_ASSERT(numCountedEdges == this->getNumEdges() && "Frequency did not count all edges");
    return frequencies;
}

void Graph::printGraphStats()
{
    std::cout << "\nStatistics for Network " << this->networkName << ":" << std::endl;
    std::cout << "Edge Type = " << this->type << ", Vertices = " <<this->numVertices << ", Edges = " << this->numEdges << std::endl;
    std::cout << "Max Coordinate Range = " << this->getMaxCoordinateRange() << std::endl;
    std::cout << "Max Graph Speed by Edges = " << this->maxGraphSpeedByEdge << std::endl;
    std::cout << "Min Graph Speed by Edges = " << this->minGraphSpeedByEdge << std::endl;
    std::vector<int> outdegreeFrequencies = this->countOutdegreeFrequencies();
    std::cout << "\nOutdegree Frequencies:" << std::endl;
    for (std::size_t i = 0; i < outdegreeFrequencies.size(); ++i) {
        std::cout << "Outdegree of " << i << ": " << outdegreeFrequencies[i] << std::endl;
    }
    std::cout << "\nConnected = " << (this->isConnectedGraph() ? "Yes" : "No") << std::endl;
    std::cout << "Undirected = " << (this->isUndirectedGraph() ? "Yes" : "No") << std::endl;
    std::cout << "Contiguous = " << (this->isContiguous() ? "Yes" : "No") << std::endl;
    std::cout << "Min (X,Y) = " << "(" << minXCoord << "," << minYCoord << ")" << std::endl;
    std::cout << "Max (X,Y) = " << "(" << maxXCoord << "," << maxYCoord << ")" << std::endl;
    std::cout << "Min X Vertex= " << minXVertex << ", Max X Vertex= " << maxXVertex << std::endl;
    std::cout << "Min Y Vertex= " << minYVertex << ", Max Y Vertex= " << maxYVertex << std::endl;
}

double Graph::getMaxGraphSpeedByEdge()
{
    return this->maxGraphSpeedByEdge;
}

double Graph::getMinGraphSpeedByEdge()
{
    return this->minGraphSpeedByEdge;
}

void Graph::outputToTSVFile(std::string outputFilePath)
{
    int adjListStart, nextAdjListStart, euclideanDistance;
    std::ofstream outputTSVFile(outputFilePath, std::ios::out | std::ios::trunc);
    if (outputTSVFile.is_open()) {
        for (VertexID vertex = 0; vertex < this->numVertices; ++vertex) {
            adjListStart = this->getEdgeListStartIndex(vertex);
            nextAdjListStart = this->getEdgeListSize(vertex);

            for (int j = adjListStart; j < nextAdjListStart; ++j) {
                // Note: We lower bound above so that it never exceeds possible distance edge weight
                if (this->edges[j].first > vertex) {
                    // Only output neighbours that we have not seen yet, so that TSV file is undirected
                    euclideanDistance = std::floor(this->getEuclideanDistance(vertex,this->edges[j].first));
                    outputTSVFile << vertex << "\t" << this->edges[j].first << "\t" << this->edges[j].second << "\t" << euclideanDistance << "\n";
                    // Note: We output Euclidean distance between neighbours for travel time and distance egde weights (in addition to the edge weight)
                }
            }
        }
    } else {
        std::cerr << "Cannot open output TSV file " << outputFilePath << std::endl;
        exit(1);
    }
    outputTSVFile.close();
}

void Graph::outputDIMACSFiles(std::string grFilePath, std::string coFilePath)
{
    std::ofstream outputCoordFile(coFilePath, std::ios::out | std::ios::trunc);
    if (outputCoordFile.is_open()) {
        outputCoordFile << "c Graph Output in Format from 9th DIMACS Implementation Challenge: Shortest Paths" << "\n";
        outputCoordFile << "c Reference: http://www.dis.uniroma1.it/~challenge9" << "\n";
        outputCoordFile << "c Graph vertices coords for " << this->networkName << "\n";
        outputCoordFile << "c" << "\n";
        outputCoordFile << "p aux sp co " << this->numVertices << "\n";
        outputCoordFile << "c Graph contains " << this->numVertices << " vertices" << "\n";
        outputCoordFile << "c" << "\n";
        for (VertexID i = 0; i < this->coordinates.size(); ++i) {
            // Note: DIMACS format starts counting vertices from 1 so we offset all vertices by 1
            outputCoordFile << "v " << i+1 << " " << this->coordinates[i].first << " " << this->coordinates[i].second << "\n";
        }
    } else {
        std::cerr << "Cannot open coordinate output file " << coFilePath << std::endl;
    }
    outputCoordFile.close();
    
    std::ofstream outputGrFile(grFilePath, std::ios::out | std::ios::trunc);
    if (outputGrFile.is_open()) {
        outputGrFile << "c Graph Output in Format from 9th DIMACS Implementation Challenge: Shortest Paths" << "\n";
        outputGrFile << "c Reference: http://www.dis.uniroma1.it/~challenge9" << "\n";
        outputGrFile << "c Graph vertices coords for " << this->networkName << "\n";
        outputGrFile << "c" << "\n";
        outputGrFile << "p sp " << this->numVertices << " " << this->numEdges << "\n";
        std::string edgeType = (this->type == constants::DISTANCE_WEIGHTS ? "distance" : "travel-time");
        outputGrFile << "c Graph contains " << this->numVertices << " vertices and " << this->numEdges << " representing " << edgeType << "\n";
        outputGrFile << "c" << "\n";

        int adjListStart, nextAdjListStart;
        for (VertexID i = 0; i < this->numVertices; ++i) {
            adjListStart = this->getEdgeListStartIndex(i);
            nextAdjListStart = this->getEdgeListSize(i);
            for (int j = adjListStart; j < nextAdjListStart; ++j) {
                // Note: DIMACS format starts counting vertices from 1 so we offset all vertices by 1
                outputGrFile << "a " << i+1 << " " << this->edges[j].first+1 << " " << this->edges[j].second << "\n";
            }
        }
    } else {
        std::cerr << "Cannot open graph output file " << grFilePath << std::endl;
    }
    outputGrFile.close();
}

Graph Graph::createSubGraph(std::string subgraphName, std::unordered_set<VertexID>& subgraphVertices)
{
    int adjListStart, nextAdjListStart, numEdges = 0;
    
    // Remove any vertices that don't have edges to other subgraph vertices
    std::unordered_set<VertexID> tempSubgraphVertices;
    tempSubgraphVertices.reserve(subgraphVertices.size());
    int edgesInSubgraph;
    for (auto it = subgraphVertices.begin(); it != subgraphVertices.end(); ++it) {
        adjListStart = this->getEdgeListStartIndex(*it);
        nextAdjListStart = this->getEdgeListSize(*it);
        edgesInSubgraph = 0;
        for (int j = adjListStart; j < nextAdjListStart; ++j) {
            // Only insert edges to vertex within subgraph
            if (subgraphVertices.find(this->edges[j].first) != subgraphVertices.end()) {
                ++edgesInSubgraph;
            }
        }
        if (edgesInSubgraph != 0) {
            tempSubgraphVertices.insert(*it);
        }
    }    
    subgraphVertices.swap(tempSubgraphVertices);
    tempSubgraphVertices.clear();
    
    // Creat Graph object containing only subgraph vertices provided
    std::vector<VertexID> newVertexOrder;
    std::unordered_map<VertexID,VertexID> oldIDtoNewIDMap;
    newVertexOrder.reserve(subgraphVertices.size());
    oldIDtoNewIDMap.reserve(subgraphVertices.size());
    
    for (auto it = subgraphVertices.begin(); it != subgraphVertices.end(); ++it) {
        oldIDtoNewIDMap[*it] = newVertexOrder.size();
        newVertexOrder.push_back(*it);
    }
    
    Graph subGraph;
    // Insert Vertices
    for (VertexID vertex = 0; vertex < newVertexOrder.size(); ++vertex) {
        subGraph.insertVertex(vertex,this->coordinates[newVertexOrder[vertex]].first,this->coordinates[newVertexOrder[vertex]].second);
    }
    
    // Insert Edges
    for (VertexID vertex = 0; vertex < newVertexOrder.size(); ++vertex) {
        adjListStart = this->getEdgeListStartIndex(newVertexOrder[vertex]);
        nextAdjListStart = this->getEdgeListSize(newVertexOrder[vertex]);
        for (int j = adjListStart; j < nextAdjListStart; ++j) {
            // Only insert edges to vertex within subgraph
            if (subgraphVertices.find(this->edges[j].first) != subgraphVertices.end()) {
                ++numEdges;
                subGraph.insertEdge(vertex,oldIDtoNewIDMap[this->edges[j].first],this->edges[j].second);
            }
        }
    }
    // Need to mark end fo edges index so we don't buffer overflow
    subGraph.firstEdgeIndex.push_back(subGraph.edges.size());
    
    // Set Graph Data
    subGraph.setGraphData(subgraphName,newVertexOrder.size(),numEdges,this->type);
    
    RN_DEBUG_ASSERT(subGraph.isConnectedGraph() && "New subgraph is not a strongly connected graph!");
    RN_DEBUG_ASSERT(subGraph.getNumVertices()+1 == static_cast<int>(subGraph.firstEdgeIndex.size()) && "Subgraph has incorrect number of vertices");
    RN_DEBUG_ASSERT(subGraph.getNumEdges() == static_cast<int>(subGraph.edges.size()) && "Subgraph has incorrect number of edges");
    
    return subGraph;
}

double Graph::computeIndexSize()
{
    double memoryUsageBytes = 0;
    memoryUsageBytes += sizeof(unsigned int)*this->firstEdgeIndex.size();
    memoryUsageBytes += sizeof(VertexEdgeWeightPair)*this->edges.size();
    memoryUsageBytes += sizeof(CoordinatePair)*this->coordinates.size();
    return memoryUsageBytes/(1024*1024);
}

double Graph::computeINEIndexSize()
{
    double memoryUsageBytes = 0;
    memoryUsageBytes += sizeof(unsigned int)*this->firstEdgeIndex.size();
    memoryUsageBytes += sizeof(VertexEdgeWeightPair)*this->edges.size();
    return memoryUsageBytes/(1024*1024);
}
