// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "TransformInputCommand.h"

#include "../processing/MortonCode.h"
#include "../queue/BinaryMinHeap.h"
#include "../utility/StopWatch.h"
#include "../utility/utility.h"
#include "../utility/geometry.h"
#include "../common.h"

#include <unordered_set>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <regex>
#include <unordered_map>
#include <cstdlib>
#include <cmath>

void TransformInputCommand::execute(int argc, char* argv[])
{
    std::string method = "";
    std::string grFilePath = "";
    std::string coFilePath = "";
    std::string outputFilePrefix = "";
    std::string regionName = "";
    std::string cropSpec = "";                                   // dem, metis: -x
    double latitude = std::numeric_limits<double>::quiet_NaN();  // dem: -l
    long long units = 0;                                         // dem: units per metre (-u, default 10); metis: scale (-u, 0 = automatic)
    
   /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "m:e:n:o:r:x:l:u:")) != -1) {
        switch (opt) {
            case 'x':
                cropSpec = optarg;
                break;
            case 'l':
                latitude = std::stod(optarg);
                break;
            case 'u':
                units = std::stoll(optarg);
                break;
            case 'm':
                method = optarg;
                break;
            case 'e':
                grFilePath = optarg;
                break;
            case 'n':
                coFilePath = optarg;
                break;
            case 'o':
                outputFilePrefix = optarg;
                break;
            case 'r':
                regionName = optarg;
                break;
            default:
                std::cerr << "Unknown option(s) provided!\n\n";
                showCommandUsage(argv[0]);
                exit(1);
        }
    }  

    // Validate Command Line Arguments
    if (argc == 5 && method != "") {
        // This is 5 so that user can just enter method to find out what parameters are required for method
        this->showMethodUsage(method,argv[0]);
        exit(1);
    } 
    
    if (method == TRANSFORM_DEM) {
        // A DEM tile has no separate coordinate file; the network name is required
        if (grFilePath == "" || outputFilePrefix == "" || regionName == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
    } else if (method == TRANSFORM_METIS) {
        // A METIS graph names no network; the network name is required
        if (grFilePath == "" || coFilePath == "" || outputFilePrefix == "" || regionName == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
    } else if (grFilePath == "" || coFilePath == "" || outputFilePrefix == "") {
        std::cerr << "Invalid argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }

    if (method == TRANSFORM_DEM) {
        StopWatch sw;
        sw.start();
        this->transformDEMFile(grFilePath,outputFilePrefix,regionName,cropSpec,latitude,units == 0 ? 10 : units);
        sw.stop();
        std::cout << "\nDEM tile processed in " << sw.getTimeMs() << "ms" << std::endl;

    } else if (method == TRANSFORM_METIS) {
        StopWatch sw;
        sw.start();
        this->transformMETISFiles(grFilePath,coFilePath,outputFilePrefix,regionName,cropSpec,units);
        sw.stop();
        std::cout << "\nMETIS mesh files processed in " << sw.getTimeMs() << "ms" << std::endl;

    } else if (method == TRANSFORM_DIMACS) {
        if (argc < 11) {
            // Arguments: -e <edge file> -n <vertex file> -o <prefix for output files>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        StopWatch sw;
        
        sw.start();
        this->transformDIMACSFiles(grFilePath,coFilePath,outputFilePrefix);
        sw.stop();
        
        double processingTimeMs = sw.getTimeMs();   
    
        std::cout << "\nDIMACS input files processed in " << processingTimeMs << "ms" << std::endl;
            
    } else {
        std::cerr << "Invalid query method!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
        
}

void TransformInputCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD + " -m <query method>\n\n"
              << "Methods Options:\n"
              << utility::getFormattedUsageString(TRANSFORM_DIMACS,"Process DIMACS files") + "\n"
              << utility::getFormattedUsageString(TRANSFORM_DEM,"Terrain surface network from an SRTM .hgt tile (distance weights)") + "\n"
              << utility::getFormattedUsageString(TRANSFORM_METIS,"Mesh from a METIS graph file and its coordinates (distance weights)") + "\n";
}

void TransformInputCommand::showMethodUsage(std::string method, std::string programName)
{
    if (method == TRANSFORM_DIMACS ) {
        std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD
                  << " -m " + TRANSFORM_DIMACS + " -e <graph edge file>"
                  << " -n <graph vertex file> \n-o <prefix for output files>\n";
    } else if (method == TRANSFORM_DEM) {
        std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD
                  << " -m " + TRANSFORM_DEM + " -e <.hgt tile> -o <output directory> -r <network name>\n"
                  << "[-x <first row,first column,rows,columns>] [-l <latitude of the tile's south-west corner>] [-u <units per metre>]\n"
                  << "The tile is a square grid of big-endian 16-bit elevations in metres (3601x3601 for 1 arcsecond). Every non-void\n"
                  << "cell of the crop (default: the whole tile) is a vertex at its projected position, joined to its 8 neighbours by the\n"
                  << "3D surface length of the step. The latitude defaults to the one in the file name (e.g. N39W107); units default to 10\n";
    } else if (method == TRANSFORM_METIS) {
        std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD
                  << " -m " + TRANSFORM_METIS + " -e <METIS graph file> -n <coordinate file> -o <output directory> -r <network name>\n"
                  << "[-x <min x,min y,max x,max y>] [-u <scale>]\n"
                  << "The graph file is in METIS/Chaco format without weights (% comments, header n m, then the 1-based neighbours of each\n"
                  << "vertex); the coordinate file has one x y line per vertex. Coordinates are translated to the origin and multiplied by the\n"
                  << "scale (0 or absent: the largest power of ten keeping the extent below 2^29), the edge weight is the planar length. A\n"
                  << "bounding box keeps only the vertices inside it (induced subgraph)\n";
    } else {
        this->showCommandUsage(programName);
    }
}

void TransformInputCommand::addVertex(VertexID vertex, Coordinate x, Coordinate y)
{
    if (vertex == vertices.size()) {
        this->vertices.push_back(vertex);
        this->coordinates.push_back(CoordinatePair(x,y));
    } else {
        std::cerr << "Vertices are not in order" << std::endl;
        std::exit(1);
    }
}

void TransformInputCommand::addEdge(VertexID source, VertexID target, EdgeWeight weight)
{
    if (weight > 0) {
        this->neighbours[source].push_back(VertexEdgeWeightPair(target,weight));
    } else {
        std::cerr << "Edge weight is zero" << std::endl;
        std::exit(1);
    }
}

void TransformInputCommand::checkValidVertex(VertexID vertex)
{
    if (vertex >= this->vertices.size()) {
        std::cerr << "Vertex (" << vertex << ") does not exist (too large)" << std::endl;
        exit(1);
    }
}

bool TransformInputCommand::isValidEdge(VertexID source, VertexID target, int rawWeight)
{   
    bool isValidEdge = true;
    if (source == target) {
        // Remove self-loops
        ++this->numSelfLoops;
        isValidEdge = false;
    }
    
    if (isValidEdge && rawWeight < 0) {
        // Remove negative edges
        ++this->numNegativeEdgeWeights;
        isValidEdge = false;
    }
    
    // Remove parallel edges (i.e. multiple edges to same target)
    if (isValidEdge) {
        // Since weight it's positive if we reach here we cast to EdgeWeight
        for (std::size_t i = 0; i < this->neighbours[source].size(); ++i) {
            if (this->neighbours[source][i].first == target) {
                ++this->numParallelEdges;
                isValidEdge = false;
                break;
            }
        }
    }
    return isValidEdge;
}

EdgeWeight TransformInputCommand::correctEdgeWeightByType(VertexID source, VertexID target, EdgeWeight weight)
{
    // For distance-based edge weights, edge weight cannot be 
    // smaller than euclidean distance between source and target
    if (this->edgeType == constants::DISTANCE_WEIGHTS) {
        // We ceil to make sure edge weight is always bigger than euclidean distance
        EdgeWeight minWeight = std::ceil(geometry::getEuclideanDist(this->coordinates[source].first,this->coordinates[source].second,
                                                                    this->coordinates[target].first,this->coordinates[target].second));
        if (minWeight == 0) {
            std::cerr << "There is an edge between two vertices that have same coordinates!" << std::endl;
            std::cerr << "Vertex " << source << " has the same coordinates as vertex " << target << std::endl;
            exit(1);
        }
        if (weight < minWeight) {
            this->weightAddedByEuclideanCorrections += (minWeight-weight);
            ++this->numEuclideanCorrections;
            if (weight*10 < minWeight) {
                this->numCorrectionsOverFactorTen++;
            }
            return minWeight;
        }
    }
    return weight;
}

void TransformInputCommand::fixShortestPathOverflow()
{
    StopWatch sw;
    sw.start();
    LongPathDistance estMaxSPDist = this->computePseudoDiameter(100);
    sw.stop();
    this->pseudoDiameterTimeMs = sw.getTimeMs();
    
    // We adjust weights so that maximum possible shortest (i.e. one involving)
    // all edges fit in EdgeWeight (at it's smallest it is 32-bits)
    LongPathDistance maxPossibleSPDist = std::numeric_limits<EdgeWeight>::max();
    
    if (estMaxSPDist > maxPossibleSPDist) {
        this->wasDeflated = true;
        int deflateShift = 1;
        while ((estMaxSPDist >> deflateShift) > maxPossibleSPDist) {
            ++deflateShift;
        }
        this->deflationFactor = 1 << deflateShift;

        // We need to deflate X and Y coordinates first so Euclidean distance will match new weights
        for (std::size_t i = 0; i < this->vertices.size(); ++i) {
            this->coordinates[i].first = this->coordinates[i].first >> deflateShift;
            this->coordinates[i].second = this->coordinates[i].second>> deflateShift;
        }
        
        // Adjust all edge weights by deflation factor (but check Euclidean distances to be sure)
        unsigned long long newTotalEdgeWeight = 0;
        for (std::size_t i = 0; i < this->vertices.size(); ++i) {
            for (std::size_t j = 0; j < this->neighbours[i].size(); ++j) {
                this->neighbours[i][j].second = this->neighbours[i][j].second >> deflateShift;
                if (this->neighbours[i][j].second == 0) {
                    std::cerr << "Deflation has caused a zero edge weight!" << std::endl;
                    exit(1);
                }
                newTotalEdgeWeight += this->neighbours[i][j].second;
                EdgeWeight minWeight = std::ceil(geometry::getEuclideanDist(this->coordinates[i].first,this->coordinates[i].second,
                                                                            this->coordinates[this->neighbours[i][j].first].first,
                                                                            this->coordinates[this->neighbours[i][j].first].second));
                // Note: It's possible that if a Coordinate collision was caused by deflation the minWeight is
                // zero. That's fine, it will not affect the edge weight and we will remove non-unique Vertices later
                if (this->neighbours[i][j].second < minWeight) {
                    this->neighbours[i][j].second = minWeight;
                }
            }
        }
        
        this->oldTotalEdgeWeight = this->totalEdgeWeight;
        this->totalEdgeWeight = newTotalEdgeWeight;
    }
}

LongPathDistance TransformInputCommand::computePseudoDiameter(int numRounds)
{
    std::vector<bool> isVertexSettled(this->numVertices,false);
    BinaryMinHeap<LongPathDistance,VertexID> pqueue;
    
    LongPathDistance minDist, maxSPDist = 0, prevMaxSPDist = 0;
    VertexID minDistVertexID, adjVertex, maxSPDistVertex, maxSPDistVertexDegree = std::numeric_limits<VertexID>::max();
    
    // Initialize with priority queue with source vertex
    pqueue.insert(0,0);
    int i;
    maxSPDistVertex = 0;
    for (i = 0; i < numRounds && pqueue.size() != 0; ++i) {
        maxSPDist = 0;
        while (pqueue.size() > 0) {
            // Extract and remove vertex with smallest distance from source
            // and mark it as "settled" so we do not inspect again
            minDist = pqueue.getMinKey();
            minDistVertexID = pqueue.extractMinElement();
            if (!isVertexSettled[minDistVertexID]) {
                isVertexSettled[minDistVertexID] = true;
                
                if (maxSPDist < minDist) {
                    maxSPDist = minDist;
                    maxSPDistVertex = minDistVertexID;
                    maxSPDistVertexDegree = this->neighbours[minDistVertexID].size();
                } else if (maxSPDist == minDist && maxSPDistVertexDegree > this->neighbours[minDistVertexID].size()) {
                    // If the distances are equal we choose the vertex with the smaller degree
                    maxSPDistVertex = minDistVertexID;
                    maxSPDistVertexDegree = this->neighbours[minDistVertexID].size();
                }

                for (std::size_t i = 0; i < this->neighbours[minDistVertexID].size(); ++i) {
                    adjVertex = this->neighbours[minDistVertexID][i].first;
                    // Only update those we haven't already settled
                    if (!isVertexSettled[adjVertex]) {
                        pqueue.insert(adjVertex,minDist+this->neighbours[minDistVertexID][i].second);
                    }
                }
            }
        }    
        
        isVertexSettled.assign(this->numVertices,false);
        pqueue.clear();
        
        // Make the vertex with maximum shortest path the next vertex
        // unless we did not improve the max shortest path distance
        // (which case we stop irrespective of number of rounds)        
        if (prevMaxSPDist < maxSPDist) {
            pqueue.insert(maxSPDistVertex,0);
            prevMaxSPDist = maxSPDist;
        }
    }
    return prevMaxSPDist;
}

void TransformInputCommand::computeCoordinateRanges()
{
    Coordinate minX, minY, maxX, maxY;  
    
    minX = this->coordinates[0].first;
    maxX = this->coordinates[0].first;
    minY = this->coordinates[0].second;
    maxY = this->coordinates[0].second;    
    for (std::size_t i = 1; i < vertices.size(); ++i) {
        if (this->coordinates[i].first < minX) {
            minX = this->coordinates[i].first;
        }
        if (this->coordinates[i].first > maxX) {
            maxX = this->coordinates[i].first;
        }
        if (this->coordinates[i].second < minY) {
            minY = coordinates[i].second;
        }
        if (this->coordinates[i].second > maxY) {
            maxY = this->coordinates[i].second;
        }
    }

    // Take opportunity to note x and y translations required to make
    // origin (0,0) and hence use Morton numbers to check uniqueness later
    this->xTranslation = minX;
    this->yTranslation = minY;

}

void TransformInputCommand::checkNonUniqueCoordinates()
{
    int nonUniqueCoordinates = 0;
    std::unordered_set<MortonNumber> zNumbers;    
    for (VertexID i = 0; i < this->vertices.size(); ++i) {
        MortonCode code(this->coordinates[i].first-this->xTranslation,this->coordinates[i].second-this->yTranslation,31);
        MortonNumber zNumber = code.getZNumber();
        if (zNumbers.find(zNumber) == zNumbers.end()) {
            zNumbers.insert(code.getZNumber());
        } else {
            this->nonUniqueVertices.insert(i);
            ++nonUniqueCoordinates;
        }
    }
    if (nonUniqueCoordinates > 0) {
        std::cout << "Found " << nonUniqueCoordinates << " non-unique coordinates" << std::endl;
    }
}

void TransformInputCommand::removeNonUniqueVertices()
{
    std::cerr << "This function doesn't work - removing vertices currently may disconnect entire subgraphs but we don't remove them" << std::endl;
    if (this->nonUniqueVertices.size() == 0) {
        // Nothing to do if there are no non-unique vertices
        return;
    }
    std::unordered_map<VertexID,VertexID> newToOldVertexID;
    std::unordered_map<VertexID,VertexID> oldToNewVertexID;
    std::vector<VertexID> newVertices;
    std::vector<std::vector<VertexEdgeWeightPair>> newNeighbours;
    std::vector<CoordinatePair> newCoordinates;
    VertexID newID = 0;
    for (VertexID i = 0; i < this->vertices.size(); ++i) {
        if (this->nonUniqueVertices.find(i) == this->nonUniqueVertices.end()) {
            // This vertex should not be removed so we give it a new renumbered ID
            newToOldVertexID[newID] = i;
            oldToNewVertexID[i] = newID;
            newVertices.push_back(newID);
            newCoordinates.push_back(std::move(this->coordinates[i]));
            ++newID;
            if (newID != newVertices.size()) {
                std::cerr << "INTERNAL ERROR: renumbering failed at vertex " << i << " while removing non-unique vertices" << std::endl;
                exit(1);
            }
        } else {
            this->numNonUniqueVerticesRemoved++;
            // This vertex's edges will be implicitly removed
            this->numNonUniqueVertexEdgesRemoved += this->neighbours[i].size();
        }
    }
    this->numVertices = newVertices.size();

    // We copy the old adjacency list for each vertex, excluding
    // neighbours that are non-unique (removed above)
    
    // Note: This might result in some vertices becoming disconnected, 
    // so we later call removeIsolatedVertices() to remove them
    VertexID oldVertexID, neighbourNewID;
    this->numEdges = 0;
    for (VertexID i = 0; i < newVertices.size(); ++i) {
        if (newToOldVertexID.find(i) == newToOldVertexID.end()) {
            std::cerr << "INTERNAL ERROR: new vertex " << i << " has no old ID while removing non-unique vertices" << std::endl;
            exit(1);
        }
        oldVertexID = newToOldVertexID[i];
        std::vector<VertexEdgeWeightPair> newAdjacencyList;
        for (std::size_t j = 0; j < this->neighbours[oldVertexID].size(); ++j) {
            if (this->nonUniqueVertices.find(this->neighbours[oldVertexID][j].first) == this->nonUniqueVertices.end()) {
                if (oldToNewVertexID.find(this->neighbours[oldVertexID][j].first) == oldToNewVertexID.end()) {
                    std::cerr << "INTERNAL ERROR: old vertex " << this->neighbours[oldVertexID][j].first << " has no new ID while removing non-unique vertices" << std::endl;
                    exit(1);
                }
                neighbourNewID = oldToNewVertexID[this->neighbours[oldVertexID][j].first];
                newAdjacencyList.push_back(VertexEdgeWeightPair(neighbourNewID,this->neighbours[oldVertexID][j].second));
            } else {
                this->numNonUniqueVertexEdgesRemoved++;
            }
        }
        this->numEdges += newAdjacencyList.size();
        newNeighbours.push_back(newAdjacencyList);
    }
    this->vertices = newVertices;
    this->neighbours = newNeighbours;
    this->coordinates = newCoordinates;
    if (this->vertices.size() != this->neighbours.size() || this->neighbours.size() != this->coordinates.size()) {
        std::cerr << "INTERNAL ERROR: " << this->vertices.size() << " vertices, " << this->neighbours.size()
                  << " adjacency lists and " << this->coordinates.size() << " coordinates after removing non-unique vertices" << std::endl;
        exit(1);
    }
}

void TransformInputCommand::removeIsolatedVertices()
{
    // Note: Removing disconnected vertices will not cause more disconnected
    // to be created (unlike removing non-unique) because by definition a
    // disconnected vertex is not present in any other vertices adjacency list
    std::unordered_map<VertexID,VertexID> newToOldVertexID;
    std::unordered_map<VertexID,VertexID> oldToNewVertexID;
    std::vector<VertexID> newVertices;
    std::vector<std::vector<VertexEdgeWeightPair>> newNeighbours;
    std::vector<CoordinatePair> newCoordinates;
    std::unordered_set<VertexID> disconnectedVertices;
    
    VertexID newID = 0;
    for (VertexID i = 0; i < this->vertices.size(); ++i) {
        if (this->neighbours[i].size() != 0) {
            // This is not a disconnected vertex
            newToOldVertexID[newID] = i;
            oldToNewVertexID[i] = newID;
            newVertices.push_back(newID);
            newCoordinates.push_back(std::move(this->coordinates[i]));
            ++newID;
            if (newID != newVertices.size()) {
                std::cerr << "INTERNAL ERROR: renumbering failed at vertex " << i << " while removing isolated vertices" << std::endl;
                exit(1);
            }
        } else {
            disconnectedVertices.insert(i);
            this->numDisconnectedVerticesRemoved++;
            // Note: We do not need to count how many edges are removed
            // because we have already checked this vertex has no edges
        }
    }

    if (disconnectedVertices.size() == 0) {
        // Nothing to do, no disconnected vertices found
        return;
    }
    
    this->numVertices = newVertices.size();
    VertexID oldVertexID, neighbourNewID;
    for (VertexID i = 0; i < newVertices.size(); ++i) {
        if (newToOldVertexID.find(i) == newToOldVertexID.end()) {
            std::cerr << "INTERNAL ERROR: new vertex " << i << " has no old ID while removing isolated vertices" << std::endl;
            exit(1);
        }
        oldVertexID = newToOldVertexID[i];
        std::vector<VertexEdgeWeightPair> newAdjacencyList;
        for (std::size_t j = 0; j < this->neighbours[oldVertexID].size();  ++j) {
            if (disconnectedVertices.find(this->neighbours[oldVertexID][j].first) == disconnectedVertices.end()) {
                if (oldToNewVertexID.find(this->neighbours[oldVertexID][j].first) == oldToNewVertexID.end()) {
                    std::cerr << "INTERNAL ERROR: old vertex " << this->neighbours[oldVertexID][j].first << " has no new ID while removing isolated vertices" << std::endl;
                    exit(1);
                }
                neighbourNewID = oldToNewVertexID[this->neighbours[oldVertexID][j].first];
                newAdjacencyList.push_back(VertexEdgeWeightPair(neighbourNewID,this->neighbours[oldVertexID][j].second));
            } else {
                // An isolated vertex, by definition, is in no other vertex's adjacency list
                std::cerr << "INTERNAL ERROR: isolated vertex " << this->neighbours[oldVertexID][j].first
                          << " is a neighbour of vertex " << oldVertexID << " while removing isolated vertices" << std::endl;
                exit(1);
            }
        }
        newNeighbours.push_back(newAdjacencyList);
    }
    this->vertices = newVertices;
    this->neighbours = newNeighbours;
    this->coordinates = newCoordinates;
    if (this->vertices.size() != this->neighbours.size() || this->neighbours.size() != this->coordinates.size()) {
        std::cerr << "INTERNAL ERROR: " << this->vertices.size() << " vertices, " << this->neighbours.size()
                  << " adjacency lists and " << this->coordinates.size() << " coordinates after removing isolated vertices" << std::endl;
        exit(1);
    }
}

void TransformInputCommand::removeIsolatedRegions()
{
    std::vector<std::unordered_set<VertexID>> vertexRegions;
    VertexID currentVertex, neighbourVertexID;
    int unassignedVertices = this->numVertices;
    VertexID unassignedVertex = 0;
    
    std::vector<bool> isVertexFound(this->numVertices,false);
    std::deque<VertexID> pqueue;
    std::unordered_set<VertexID> currentRegion;
    
    std::size_t biggestRegionSize = 0;
    
    while (unassignedVertices > 0) {
        // Clear previous iteration data
        this->numRegionsFound++;
        pqueue.clear();
        currentRegion.clear();
        
        pqueue.push_back(unassignedVertex);
        while (pqueue.size() > 0) {
            currentVertex = pqueue.front();
            pqueue.pop_front();
            if (!isVertexFound[currentVertex]) {
                currentRegion.insert(currentVertex);
                isVertexFound[currentVertex] = true;
                unassignedVertices--;

                for (std::size_t i = 0; i < neighbours[currentVertex].size(); ++i) {
                    neighbourVertexID = neighbours[currentVertex][i].first;
                    // Only add those we haven't already seen
                    if (!isVertexFound[neighbourVertexID]) {
                        pqueue.push_back(neighbourVertexID);
                    }
                }
            }
        }
        
        if (currentRegion.size() > biggestRegionSize) {
            biggestRegionSize = currentRegion.size();
        }
        
        vertexRegions.push_back(std::move(currentRegion));
        
        if (unassignedVertices > 0) {
            // Find a new unassigned vertex (as root not of new BFS search)
            // (i.e. it is not present any of the region found so far)
            bool vertexFound = false;
            for (VertexID vertex = 0; vertex < this->numVertices; ++vertex) {
                if (!isVertexFound[vertex]) {
                    vertexFound = true;
                    unassignedVertex = vertex;
                    break;
                }
            }
            if (!vertexFound) {
                std::cerr << "INTERNAL ERROR: unassigned vertices remain, but none was found while removing isolated regions" << std::endl;
                exit(1);
            }
        }
    
    }
    
    // Remove all isolated regions except the biggest region
    bool biggestRegionPreserved = false; // In case there are two regions with same biggest size
    for (std::size_t i = 0; i < vertexRegions.size(); ++i) {
        if (vertexRegions[i].size() != biggestRegionSize || biggestRegionPreserved) {
            if (vertexRegions[i].size() == biggestRegionSize) {
                biggestRegionPreserved = true;
            }
            for (VertexID vertex: vertexRegions[i]) {
                // Clear it's adjacency list so it will be deleted as isolated vertex
                // This is OK because we know all of the edge are only within this region
                this->numIsolatedRegionVerticesRemoved++;
                this->numIsolatedRegionEdgesRemoved += this->neighbours[vertex].size();
                this->numEdges -= this->neighbours[vertex].size(); // this->numVertices will be decremented later
                this->neighbours[vertex].clear();
            }
        }
    }
    
}

void TransformInputCommand::outputStandardFormatFiles(std::string outputFilePrefix)
{
    std::string newCoFilePath = outputFilePrefix + "/" + this->networkName + "-" + this->edgeType + "." + constants::VERTEX_EXT;
    
    std::ofstream outputCoordFile(newCoFilePath, std::ios::out | std::ios::trunc);
    if (outputCoordFile.is_open()) {
        outputCoordFile << this->networkName << " " << this->edgeType << " " << this->numVertices << "\n";
        for (std::size_t i = 0; i < this->vertices.size(); ++i) {
            outputCoordFile << this->vertices[i] << " " << this->coordinates[i].first << " " << this->coordinates[i].second << "\n";
        }
    } else {
        std::cerr << "Cannot open coordinate output file " << newCoFilePath << std::endl;
    }
    
    std::string newGrFilePath = outputFilePrefix + "/" + networkName + "-" + edgeType + "." + constants::EDGE_EXT;
    
    std::ofstream outputGrFile(newGrFilePath, std::ios::out | std::ios::trunc);
    if (outputGrFile.is_open()) {
        outputGrFile << this->networkName << " " << this->edgeType << " " << this->numVertices << " " << this->numEdges << "\n";
        for (std::size_t i = 0; i < this->vertices.size(); ++i) {
            for (std::size_t j = 0; j < this->neighbours[i].size(); ++j) {
                outputGrFile << this->vertices[i] << " " << this->neighbours[i][j].first << " " << this->neighbours[i][j].second << "\n";
            }
        }
    } else {
        std::cerr << "Cannot open graph output file " << newGrFilePath << std::endl;
    }
}

void TransformInputCommand::printTransformationStatistics()
{
    std::cout << "\nGraph " << this->networkName << " processed with " << this->numVertices << " vertices and " << this->numEdges << " edges" << std::endl;
    std::cout << "Total edge weight = " << this->totalEdgeWeight << std::endl;
    std::cout << "Average edge weight = " << static_cast<double>(this->totalEdgeWeight)/this->numEdges << std::endl;

    std::cout << "\nStage 1: Validate Input File Graph" << std::endl;
    std::cout << "Number of self-loops removed = " << this->numSelfLoops << std::endl;
    std::cout << "Number of parallel edges removed = " << this->numParallelEdges << std::endl;        
    std::cout << "Number of negative weight edges removed = " << this->numNegativeEdgeWeights << std::endl;        
    std::cout << "Number of edges smaller than Euclidean distance corrected = " << this->numEuclideanCorrections << std::endl;        
    std::cout << "Number of edges smaller than Euclidean distance divided by 10 corrected = " << this->numCorrectionsOverFactorTen << std::endl;        
    if (this->numEuclideanCorrections != 0) {
        std::cout << "Average euclidean correction = " << static_cast<double>(this->weightAddedByEuclideanCorrections)/this->numEuclideanCorrections << std::endl;
    }

    std::cout << "\nStage 2: Correct Possible Shortest Path Distance Overflows" << std::endl;
    std::cout << "Pseudo Diameter Calculated in " << this->pseudoDiameterTimeMs << "ms" << std::endl;
    if (this->wasDeflated) {
        std::cout << "Sum of edge weights exceeded max possible shortest path distance value" << std::endl;
        std::cout << "Coordinate and EdgeWeight values were adjusted by deflation factor of " << this->deflationFactor << std::endl;
        std::cout << "Sum of edge weights was reduced from " << this->oldTotalEdgeWeight << " to " << this->totalEdgeWeight << std::endl;
    } else {
        std::cout << "No deflation performed, total of all edge weights below maximum possible shortest path distance value" << std::endl;
    }

    std::cout << "\nStage 3: Remove Non-Unique Vertices" << std::endl;
    std::cout << "Number of non-unique vertices removed = " << this->numNonUniqueVerticesRemoved << std::endl;
    std::cout << "Number of non-unique vertex edges removed = " << this->numNonUniqueVertexEdgesRemoved << std::endl;

    std::cout << "\nStage 4: Disconnect Isolated Region Vertices" << std::endl;
    std::cout << "Number of isolated regions found  = " << this->numRegionsFound-1 << std::endl; 
    std::cout << "Number of isolated region vertices made disconnected (removed later as below) = " << this->numIsolatedRegionVerticesRemoved << std::endl;
    std::cout << "Number of isolated region edges removed = " << this->numIsolatedRegionEdgesRemoved << std::endl;
    
    std::cout << "\nStage 5: Remove Disconnected Vertices" << std::endl;
    std::cout << "Number of disconnected vertices removed = " << this->numDisconnectedVerticesRemoved << std::endl; 
    std::cout << "Number of disconnected vertex edges removed = " << this->numDisconnectedVertexEdgesRemoved << std::endl; 

}

void TransformInputCommand::transformDIMACSFiles(std::string grFilePath, std::string coFilePath, std::string outputFilePrefix)
{
    std::cout << "\nProcessing Coordinates File " << coFilePath << std::endl;
    
    std::string line = "", lineType = "", dummyField;
    VertexID sourceID, targetID;
    EdgeWeight weight;
    int x, y, rawWeight;
    std::string coordEdgeType;
    
    std::ifstream coordFile(coFilePath, std::ios::in);
    if (coordFile.is_open()) {
        // Note: Boost regular expression has to match whole input string
        // partial matches will not occur
        std::regex expression(".*\\-([dt])\\.([^.]+)\\.co$");
        std::smatch matches;
        if (std::regex_match(coFilePath,matches,expression)) {
            coordEdgeType = matches[1];
            this->networkName = matches[2];
        } else {
            std::cout << "Input coordinate file not named in expected DIMACS format!" << std::endl;
            exit(1);
        }
        
        if (coordEdgeType != constants::DISTANCE_WEIGHTS && coordEdgeType != constants::TIME_WEIGHTS) {
            std::cout << "Invalid edge weight type in coordinate file name" << std::endl;
            exit(1);
        }
        while (std::getline(coordFile, line))
        {
            std::stringstream ss(line);
            ss >> lineType;
            if (lineType == "p") {
                // Header line containing graph details
                ss >> dummyField >> dummyField >> dummyField >> this->numVertices;
            } else if (lineType == "v") {
                ss >> sourceID >> x >> y;
                
                // DIMACS vertices are numbered starting from 1 and they appear in order
                --sourceID;
                this->addVertex(sourceID,x,y);
            }
        }
        
        if (this->numVertices == 0) {
            std::cerr << "File header information not found or does not contain total number of vertices" << std::endl;
            exit(1);
        }
        
        if (this->numVertices != this->vertices.size() || this->numVertices != this->coordinates.size()) {
            std::cerr << "Number of vertices read from file does not match number of vertices given in header" << std::endl;
            exit(1);
        }
        
    } else {
        std::cerr << "Cannot open coordinates file " << coFilePath << std::endl;
        exit(1);
    }    
    
    unsigned int edgesFileNumVertices = 0, edgesAdded = 0;
    this->neighbours.resize(this->numVertices);
    
    std::ifstream grFile(grFilePath, std::ios::in);
    if (grFile.is_open()) {
        std::regex expression(".*\\-([dt])\\.([^.]+)\\.gr$");
        std::smatch matches;
        if (std::regex_match(grFilePath,matches,expression)) {
            this->edgeType = matches[1];
        } else {
            std::cout << "Input graph file not named in expected DIMACS format!" << std::endl;
            exit(1);
        }
        
        if (this->edgeType != constants::DISTANCE_WEIGHTS && this->edgeType != constants::TIME_WEIGHTS) {
            std::cout << "Invalid edge weight type in graph file name" << std::endl;
            exit(1);
        }

        if (this->edgeType == constants::DISTANCE_WEIGHTS) {
            std::cout << "Processing Distance-Weight Graph File " << grFilePath << std::endl;
        } else {
            std::cout << "Processing Time-Weight Graph File " << grFilePath << std::endl;
        }

        while (std::getline(grFile, line))
        {
            std::stringstream ss(line);
            ss >> lineType;
            if (lineType == "p") {
                // Header line containing graph details
                ss >> dummyField >> edgesFileNumVertices >> this->numEdges;
            } else if (lineType == "a") {
                ss >> sourceID >> targetID >> rawWeight;
                
                // DIMACS vertices are numbered starting from 1
                --sourceID;
                this->checkValidVertex(sourceID);
                --targetID;
                this->checkValidVertex(targetID);
                
                // Validate Edges
                weight = static_cast<EdgeWeight>(rawWeight);
                if (this->isValidEdge(sourceID,targetID,rawWeight)) {
                    weight = this->correctEdgeWeightByType(sourceID,targetID,weight); // Will ensure distance weights are not less than Euclid distance
                    this->addEdge(sourceID,targetID,weight);
                    this->totalEdgeWeight += weight;
                    ++edgesAdded;
                } else {
                    // We do not add invalid edges so we need to reduce number of edges
                    --numEdges;
                    // Note: This might result in some vertices becoming isolated (no neighbours), so we
                    // later call removeIsolatedVertices() to remove them
                }
            }
        }
        
        if (edgesFileNumVertices != this->numVertices) {
            std::cerr << "Number of vertices read from coordinates file does not match number of vertices given in edge file header" << std::endl;
            exit(1);
        }
        
        if (this->numEdges != edgesAdded) {
            std::cerr << "Number of edges read from file does not match number of edges given in header" << std::endl;
            exit(1);
        }
        
    } else {
        std::cerr << "Cannot open graph file " << grFilePath << std::endl;
        exit(1);
    }
    
    this->fixShortestPathOverflow();
    
    this->computeCoordinateRanges();
    
    this->checkNonUniqueCoordinates();
    // Note: If there are non-unique vertices, you can call this->removeNonUniqueVertices()
    // after this to remove them. It doesn't call this->removeNonUniqueVertices() by
    // default because you may want to allow non-unique vertices or that it could
    // be cause by another issue, e.g., rounding coordinates too much. 
    // this->removeNonUniqueVertices() must be run after this->checkNonUniqueCoordinates()
    // because it relies on non-unique vertices being saved in a vector to process.
    
    this->removeIsolatedRegions();
    
    this->removeIsolatedVertices();
    
    this->printTransformationStatistics();
    
    this->outputStandardFormatFiles(outputFilePrefix);
}

/*
 * Common tail of the dem and metis methods (the steps of transformDIMACSFiles after the input is read)
 */
void TransformInputCommand::finishTransformation(std::string outputFilePrefix)
{
    this->fixShortestPathOverflow();
    this->computeCoordinateRanges();
    this->checkNonUniqueCoordinates();
    this->removeIsolatedRegions();
    this->removeIsolatedVertices();
    this->printTransformationStatistics();
    this->printGeometryStatistics();
    this->outputStandardFormatFiles(outputFilePrefix);
}

// Degree distribution and edge weights of the graph as it stands (after the isolated regions were removed)
void TransformInputCommand::printGeometryStatistics()
{
    std::vector<EdgeWeight> weights;
    std::vector<std::size_t> degreeCounts(1,0);
    std::size_t maxDegree = 0;
    for (std::size_t i = 0; i < this->vertices.size(); ++i) {
        std::size_t degree = this->neighbours[i].size();
        if (degree >= degreeCounts.size()) {
            degreeCounts.resize(degree+1,0);
        }
        ++degreeCounts[degree];
        if (degree > maxDegree) {
            maxDegree = degree;
        }
        for (std::size_t j = 0; j < degree; ++j) {
            weights.push_back(this->neighbours[i][j].second);
        }
    }
    std::cout << "\nGeometry" << std::endl;
    std::cout << "Maximum degree = " << maxDegree << std::endl;
    std::cout << "Degree distribution (degree: vertices) =";
    for (std::size_t d = 0; d < degreeCounts.size(); ++d) {
        if (degreeCounts[d] > 0) {
            std::cout << " " << d << ": " << degreeCounts[d];
        }
    }
    std::cout << std::endl;
    if (weights.size() > 0) {
        std::sort(weights.begin(),weights.end());
        std::cout << "Edge weights (units): shortest = " << weights.front() << ", median = " << weights[weights.size()/2]
                  << ", longest = " << weights.back() << std::endl;
    }
}

/*
 * Terrain surface network from an SRTM .hgt tile
 */
namespace {

// Length of one degree of latitude and of longitude in metres at latitude phi (degrees)
double metresPerDegreeLatitude(double phi)
{
    double r = phi*M_PI/180.0;
    return 111132.954 - 559.822*std::cos(2*r) + 1.175*std::cos(4*r);
}

double metresPerDegreeLongitude(double phi)
{
    double r = phi*M_PI/180.0;
    return 111412.84*std::cos(r) - 93.5*std::cos(3*r) + 0.118*std::cos(5*r);
}

std::vector<long long> parseIntegerList(const std::string& spec, std::size_t count, const std::string& description)
{
    std::vector<std::string> fields = utility::splitByDelim(spec,',');
    if (fields.size() != count) {
        std::cerr << description << " must have " << count << " comma-separated values, got \"" << spec << "\"" << std::endl;
        exit(1);
    }
    std::vector<long long> values;
    for (std::size_t i = 0; i < fields.size(); ++i) {
        values.push_back(std::stoll(fields[i]));
    }
    return values;
}

}

void TransformInputCommand::transformDEMFile(std::string hgtFilePath, std::string outputFilePrefix, std::string name, std::string cropSpec,
                                             double latitude, long long unitsPerMetre)
{
    const int VOID_VALUE = -32768;
    this->networkName = name;
    this->edgeType = constants::DISTANCE_WEIGHTS;
    if (unitsPerMetre <= 0) {
        std::cerr << "Units per metre must be positive" << std::endl;
        exit(1);
    }

    // The tile: side x side big-endian signed 16-bit elevations in metres, row-major from the north-west corner
    std::ifstream in(hgtFilePath, std::ios::in | std::ios::binary | std::ios::ate);
    if (!in.is_open()) {
        std::cerr << "Cannot open DEM tile " << hgtFilePath << std::endl;
        exit(1);
    }
    long long bytes = in.tellg();
    in.seekg(0);
    long long cells = bytes/2;
    long long side = std::llround(std::sqrt(static_cast<double>(cells)));
    if (bytes % 2 != 0 || side*side != cells || side < 2) {
        std::cerr << hgtFilePath << " is not a square grid of 16-bit values (" << bytes << " bytes)" << std::endl;
        exit(1);
    }
    std::vector<unsigned char> raw(bytes);
    in.read(reinterpret_cast<char*>(raw.data()),bytes);
    in.close();
    auto elevation = [&](long long row, long long col) -> int {
        std::size_t i = 2*static_cast<std::size_t>(row*side+col);
        return static_cast<short>(static_cast<unsigned short>((raw[i] << 8) | raw[i+1]));
    };
    std::cout << "\nProcessing DEM tile " << hgtFilePath << " (" << side << "x" << side << " cells)" << std::endl;

    // Latitude of the tile's south-west corner: the argument, or the file name (e.g. N39W107)
    if (std::isnan(latitude)) {
        std::regex expression(".*([NS])(\\d+)[EW](\\d+)[^/]*$");
        std::smatch matches;
        if (!std::regex_match(hgtFilePath,matches,expression)) {
            std::cerr << "Cannot read the latitude from the file name; give it with -l" << std::endl;
            exit(1);
        }
        latitude = std::stod(matches[2]) * (matches[1] == "S" ? -1 : 1);
    }

    // Crop: first row, first column, rows, columns (default: the whole tile)
    long long row0 = 0, col0 = 0, rows = side, cols = side;
    if (cropSpec != "") {
        std::vector<long long> crop = parseIntegerList(cropSpec,4,"The crop");
        row0 = crop[0]; col0 = crop[1]; rows = crop[2]; cols = crop[3];
        if (row0 < 0 || col0 < 0 || rows < 2 || cols < 2 || row0+rows > side || col0+cols > side) {
            std::cerr << "The crop " << cropSpec << " does not fit the " << side << "x" << side << " tile" << std::endl;
            exit(1);
        }
    }

    // Cell sizes in metres at the crop's centre latitude (row 0 is the tile's northern edge; 1 degree spans side-1 cells)
    double centreLatitude = latitude + 1.0 - (row0 + rows/2.0)/(side-1);
    double dy = metresPerDegreeLatitude(centreLatitude)/(side-1);
    double dx = metresPerDegreeLongitude(centreLatitude)/(side-1);
    std::cout << "Crop: rows " << row0 << " to " << row0+rows-1 << ", columns " << col0 << " to " << col0+cols-1 << "; centre latitude "
              << centreLatitude << ", cell " << dx << " m east-west by " << dy << " m north-south, " << unitsPerMetre << " units per metre" << std::endl;

    // Vertices: the non-void cells, numbered in row-major order; coordinates = projected position (y grows northwards)
    const VertexID NO_VERTEX = std::numeric_limits<VertexID>::max();
    std::vector<VertexID> cellToVertex(static_cast<std::size_t>(rows*cols),NO_VERTEX);
    std::vector<int> cellElevation(static_cast<std::size_t>(rows*cols),VOID_VALUE);
    long long voids = 0;
    int minElevation = std::numeric_limits<int>::max(), maxElevation = std::numeric_limits<int>::min();
    VertexID nextID = 0;
    for (long long r = 0; r < rows; ++r) {
        for (long long c = 0; c < cols; ++c) {
            int e = elevation(row0+r,col0+c);
            if (e == VOID_VALUE) {
                ++voids;
                continue;
            }
            cellToVertex[r*cols+c] = nextID;
            cellElevation[r*cols+c] = e;
            minElevation = std::min(minElevation,e);
            maxElevation = std::max(maxElevation,e);
            Coordinate x = static_cast<Coordinate>(std::llround(c*dx*unitsPerMetre));
            Coordinate y = static_cast<Coordinate>(std::llround((rows-1-r)*dy*unitsPerMetre));
            this->addVertex(nextID,x,y);
            ++nextID;
        }
    }
    this->numVertices = nextID;
    if (this->numVertices == 0) {
        std::cerr << "The crop has no non-void cells" << std::endl;
        exit(1);
    }
    this->neighbours.resize(this->numVertices);
    std::cout << "Cells: " << rows*cols << ", void cells dropped: " << voids << ", elevation " << minElevation << " to " << maxElevation << " m" << std::endl;

    // Edges: to the 8 neighbours inside the crop that are not void. The planar part is the distance between the rounded
    // coordinates, which is what Graph::insertEdge compares the weight with, so weight >= ceil(planar distance) holds
    const int dr[8] = {-1,-1,-1,0,0,1,1,1};
    const int dc[8] = {-1,0,1,-1,1,-1,0,1};
    std::vector<float> ratios; // weight / planar distance, every 8th edge: the ruggedness of the crop
    this->numEdges = 0;
    for (long long r = 0; r < rows; ++r) {
        for (long long c = 0; c < cols; ++c) {
            VertexID s = cellToVertex[r*cols+c];
            if (s == NO_VERTEX) {
                continue;
            }
            for (int k = 0; k < 8; ++k) {
                long long nr = r+dr[k], nc = c+dc[k];
                if (nr < 0 || nc < 0 || nr >= rows || nc >= cols) {
                    continue;
                }
                VertexID t = cellToVertex[nr*cols+nc];
                if (t == NO_VERTEX) {
                    continue;
                }
                double planar = geometry::getEuclideanDist(this->coordinates[s].first,this->coordinates[s].second,
                                                           this->coordinates[t].first,this->coordinates[t].second);
                double dz = static_cast<double>(cellElevation[nr*cols+nc]-cellElevation[r*cols+c])*unitsPerMetre;
                EdgeWeight weight = static_cast<EdgeWeight>(std::ceil(std::sqrt(planar*planar + dz*dz)));
                this->addEdge(s,t,weight);
                this->totalEdgeWeight += weight;
                ++this->numEdges;
                if ((this->numEdges & 7) == 0) {
                    ratios.push_back(static_cast<float>(weight/planar));
                }
            }
        }
    }
    if (ratios.size() > 0) {
        std::sort(ratios.begin(),ratios.end());
        std::cout << "Surface length / planar distance over the edges (sampled): median " << ratios[ratios.size()/2] << ", 90% " << ratios[ratios.size()*9/10]
                  << ", 99% " << ratios[ratios.size()*99/100] << ", maximum " << ratios.back() << std::endl;
    }

    this->finishTransformation(outputFilePrefix);
}

/*
 * Mesh from METIS/Chaco graph and coordinate files
 */
void TransformInputCommand::transformMETISFiles(std::string grFilePath, std::string coFilePath, std::string outputFilePrefix, std::string name,
                                                std::string cropSpec, long long scale)
{
    this->networkName = name;
    this->edgeType = constants::DISTANCE_WEIGHTS;

    // Coordinates: one "x y" line per vertex, % comment lines skipped
    std::cout << "\nProcessing coordinate file " << coFilePath << std::endl;
    std::vector<double> xs, ys;
    {
        std::ifstream coordFile(coFilePath, std::ios::in);
        if (!coordFile.is_open()) {
            std::cerr << "Cannot open coordinate file " << coFilePath << std::endl;
            exit(1);
        }
        std::string line;
        while (std::getline(coordFile,line)) {
            if (line.size() == 0 || line[0] == '%') {
                continue;
            }
            std::stringstream ss(line);
            double x, y;
            if (!(ss >> x >> y)) {
                std::cerr << "Cannot read two coordinates from line \"" << line << "\"" << std::endl;
                exit(1);
            }
            xs.push_back(x);
            ys.push_back(y);
        }
    }
    std::size_t numInputVertices = xs.size();
    if (numInputVertices == 0) {
        std::cerr << "No coordinates read" << std::endl;
        exit(1);
    }

    // Graph: % comments, header "n m [fmt [ncon]]" (no weights), then one line per vertex with its 1-based neighbours
    std::cout << "Processing METIS graph file " << grFilePath << std::endl;
    std::vector<std::vector<VertexID>> adjacency(numInputVertices);
    long long headerEdges = 0;
    {
        std::ifstream grFile(grFilePath, std::ios::in);
        if (!grFile.is_open()) {
            std::cerr << "Cannot open graph file " << grFilePath << std::endl;
            exit(1);
        }
        std::string line;
        bool headerRead = false;
        std::size_t vertex = 0;
        while (std::getline(grFile,line)) {
            if (line.size() > 0 && line[0] == '%') {
                continue;
            }
            std::stringstream ss(line);
            if (!headerRead) {
                long long n = 0;
                std::string format = "0";
                ss >> n >> headerEdges >> format;
                if (n != static_cast<long long>(numInputVertices)) {
                    std::cerr << "The graph has " << n << " vertices but the coordinate file " << numInputVertices << std::endl;
                    exit(1);
                }
                if (format.find_first_not_of('0') != std::string::npos) {
                    std::cerr << "METIS format " << format << " (weights) is not supported" << std::endl;
                    exit(1);
                }
                headerRead = true;
                continue;
            }
            if (vertex >= numInputVertices) {
                if (line.find_first_not_of(" \t\r") == std::string::npos) {
                    continue; // Trailing blank line
                }
                std::cerr << "More vertex lines than vertices" << std::endl;
                exit(1);
            }
            long long neighbour;
            while (ss >> neighbour) {
                if (neighbour < 1 || neighbour > static_cast<long long>(numInputVertices)) {
                    std::cerr << "Neighbour " << neighbour << " of vertex " << vertex+1 << " is out of range" << std::endl;
                    exit(1);
                }
                adjacency[vertex].push_back(static_cast<VertexID>(neighbour-1));
            }
            ++vertex;
        }
        if (vertex != numInputVertices) {
            std::cerr << "Read " << vertex << " vertex lines, expected " << numInputVertices << std::endl;
            exit(1);
        }
    }

    // Crop: the vertices inside the bounding box (induced subgraph)
    std::vector<bool> kept(numInputVertices,true);
    std::size_t numKept = numInputVertices;
    if (cropSpec != "") {
        std::vector<std::string> fields = utility::splitByDelim(cropSpec,',');
        if (fields.size() != 4) {
            std::cerr << "The bounding box must be min x,min y,max x,max y" << std::endl;
            exit(1);
        }
        double minX = std::stod(fields[0]), minY = std::stod(fields[1]), maxX = std::stod(fields[2]), maxY = std::stod(fields[3]);
        numKept = 0;
        for (std::size_t v = 0; v < numInputVertices; ++v) {
            kept[v] = xs[v] >= minX && xs[v] <= maxX && ys[v] >= minY && ys[v] <= maxY;
            numKept += kept[v];
        }
        std::cout << "Bounding box " << cropSpec << " keeps " << numKept << " of " << numInputVertices << " vertices" << std::endl;
        if (numKept == 0) {
            exit(1);
        }
    }

    // Extent and scale: coordinates are translated to the origin and multiplied by the scale. Coordinate is a 32-bit int and so
    // is the graph's coordinate range, so the extent must stay below 2^30; the automatic scale keeps it below 2^29, which also
    // keeps shortest path distances well inside EdgeWeight
    double minX = std::numeric_limits<double>::max(), minY = minX, maxX = -minX, maxY = -minX;
    for (std::size_t v = 0; v < numInputVertices; ++v) {
        if (kept[v]) {
            minX = std::min(minX,xs[v]); maxX = std::max(maxX,xs[v]);
            minY = std::min(minY,ys[v]); maxY = std::max(maxY,ys[v]);
        }
    }
    double extent = std::max(maxX-minX,maxY-minY);
    const double RANGE_LIMIT = 1073741824.0; // 2^30
    if (scale <= 0) {
        scale = 1;
        while (extent*scale*10 < RANGE_LIMIT/2) {
            scale *= 10;
        }
    }
    if (extent*scale >= RANGE_LIMIT) {
        std::cerr << "Scale " << scale << " puts the extent " << extent << " beyond the 32-bit coordinate range; use a smaller scale or a bounding box" << std::endl;
        exit(1);
    }
    std::cout << "Extent " << maxX-minX << " by " << maxY-minY << " (x " << minX << " to " << maxX << ", y " << minY << " to " << maxY << "), scale " << scale
              << ", scaled extent " << std::llround(extent*scale) << std::endl;

    // Vertices in the order of the file, renumbered over the kept ones; coincident rounded coordinates are an error, because the
    // NVD quadtree's Morton codes need distinct coordinates and removing vertices would alter the mesh
    const VertexID NO_VERTEX = std::numeric_limits<VertexID>::max();
    std::vector<VertexID> newID(numInputVertices,NO_VERTEX);
    std::unordered_set<unsigned long long> positions;
    std::size_t coincident = 0;
    VertexID nextID = 0;
    for (std::size_t v = 0; v < numInputVertices; ++v) {
        if (!kept[v]) {
            continue;
        }
        Coordinate x = static_cast<Coordinate>(std::llround((xs[v]-minX)*scale));
        Coordinate y = static_cast<Coordinate>(std::llround((ys[v]-minY)*scale));
        unsigned long long key = (static_cast<unsigned long long>(static_cast<unsigned int>(x)) << 32) | static_cast<unsigned int>(y);
        if (!positions.insert(key).second) {
            ++coincident;
        }
        newID[v] = nextID;
        this->addVertex(nextID,x,y);
        ++nextID;
    }
    this->numVertices = nextID;
    if (coincident > 0) {
        std::cerr << coincident << " vertices coincide with another after rounding at scale " << scale << "; use a larger scale or a bounding box" << std::endl;
        exit(1);
    }
    this->neighbours.resize(this->numVertices);

    // Edges: weight = ceil of the planar length between the rounded coordinates (the bound of Graph::insertEdge, so it holds with
    // equality); self-loops and duplicates are dropped and counted, and a neighbour list missing the reverse edge gets it added
    std::size_t asymmetric = 0;
    this->numEdges = 0;
    for (std::size_t v = 0; v < numInputVertices; ++v) {
        if (!kept[v]) {
            continue;
        }
        VertexID s = newID[v];
        for (std::size_t j = 0; j < adjacency[v].size(); ++j) {
            VertexID u = adjacency[v][j];
            if (!kept[u]) {
                continue;
            }
            VertexID t = newID[u];
            if (s == t) {
                ++this->numSelfLoops;
                continue;
            }
            bool duplicate = false;
            for (std::size_t k = 0; k < this->neighbours[s].size(); ++k) {
                if (this->neighbours[s][k].first == t) {
                    duplicate = true;
                    break;
                }
            }
            if (duplicate) {
                ++this->numParallelEdges;
                continue;
            }
            EdgeWeight weight = static_cast<EdgeWeight>(std::ceil(geometry::getEuclideanDist(this->coordinates[s].first,this->coordinates[s].second,
                                                                                               this->coordinates[t].first,this->coordinates[t].second)));
            this->addEdge(s,t,weight);
            this->totalEdgeWeight += weight;
            ++this->numEdges;
        }
    }
    for (VertexID s = 0; s < this->numVertices; ++s) {
        for (std::size_t j = 0; j < this->neighbours[s].size(); ++j) {
            VertexID t = this->neighbours[s][j].first;
            bool found = false;
            for (std::size_t k = 0; k < this->neighbours[t].size(); ++k) {
                if (this->neighbours[t][k].first == s) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                ++asymmetric;
                this->addEdge(t,s,this->neighbours[s][j].second);
                this->totalEdgeWeight += this->neighbours[s][j].second;
                ++this->numEdges;
            }
        }
    }
    std::cout << "Vertices kept: " << this->numVertices << " of " << numInputVertices << "; directed edges: " << this->numEdges << " (header: " << 2*headerEdges
              << "); self-loops dropped: " << this->numSelfLoops << ", duplicates dropped: " << this->numParallelEdges
              << ", reverse edges added: " << asymmetric << std::endl;

    this->finishTransformation(outputFilePrefix);
}
