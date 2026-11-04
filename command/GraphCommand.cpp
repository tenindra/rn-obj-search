// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "GraphCommand.h"

#include "../processing/Graph.h"
#include "../processing/DijkstraSearch.h"
#include "../tuple/IndexTuple.h"
#include "../common.h"
#include "../utility/StopWatch.h"
#include "../utility/serialization.h"

#include <fstream>
#include <boost/archive/binary_oarchive.hpp>

void GraphCommand::execute(int argc, char* argv[]) {
    std::string coFilePath = "";
    std::string grFilePath = "";
    std::string indexOutputFilePath = "";
    std::string statsOutputFilePath = "";
    std::string networkName = "";

    /*
     * Process Command Line Arguments
     */
    if (argc < 9) {
        // Arguments: -e <graph file> -n <coordinates file> -o <output file> -s <stats file> [-r <network name>]
        std::cerr << "Too few arguments!\n\n";
        this->showCommandUsage(argv[0]);
        return;
    }
    
    int opt;
    while ((opt = getopt (argc, argv, "e:n:o:s:r:")) != -1) {
        switch (opt) {
            case 'e':
                grFilePath = optarg;
                break;
            case 'n':
                coFilePath = optarg;
                break;
            case 'o':
                indexOutputFilePath = optarg;
                break;
            case 's':
                statsOutputFilePath = optarg;
                break;
            case 'r':
                networkName = optarg;
                break;
            default:
                std::cerr << "Unknown option(s) provided!\n\n";
                showCommandUsage(argv[0]);
                exit(1);
        }
    }    
    
    if (grFilePath == "" || indexOutputFilePath == "" 
        || coFilePath == "" || statsOutputFilePath == "") {
        std::cerr << "Invalid argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
    /*
     * Parse Graph from Text File and Serialize to Binary File
     */
    Graph graph;
    StopWatch sw;
    
    sw.start();
    graph.parseGraphFile(grFilePath,coFilePath);
    if (networkName != "") {
        // Every index and set built from the graph is named after it, so a network's text files can serve two datasets (e.g. two edge types)
        graph.setGraphData(networkName,graph.getNumVertices(),graph.getNumEdges(),graph.getEdgeType());
    }
    sw.stop();
    double processingTimeMs = sw.getTimeMs();    

    // Lower bound on the graph diameter, which range query radii are percentages of. It is computed once here (and not
    // counted in the processing time), so that every experiment and verification run uses the same ranges
    sw.reset();
    sw.start();
    DijkstraSearch dijk;
    graph.setDiameterLowerBound(dijk.getLowerboundGraphDiameterByDoubleSweep(graph,constants::DOUBLE_SWEEP_ITERATIONS));
    sw.stop();
    std::cout << "Lower-bound graph diameter " << graph.getDiameterLowerBound() << " computed in " << sw.getTimeMs() << "ms" << std::endl;

    // Print Additional Graph Statistics
    graph.printGraphStats();

    /*
     * Serialize to Binary File
     */    
    serialization::outputIndexToBinaryFile<Graph>(graph,indexOutputFilePath);
    
    /*
     * Collect Stats and Output
     */
    int vertices = graph.getNumVertices();
    int edges = graph.getNumEdges();
    double memoryUsage = graph.computeINEIndexSize(); // For experiments we report INE size (e.g., excluding coordinates information)

    IndexTuple stats(graph.getNetworkName(),vertices,edges,constants::IDX_GRAPH_CMD,processingTimeMs,memoryUsage);

    std::cout << stats.getMultilineTupleString();

    this->outputCommandStats(statsOutputFilePath,stats.getTupleString());
    
    std::cout << "Binary graph index successfully created!" << std::endl;
}

void GraphCommand::showCommandUsage(std::string programName) {
    std::cerr << "Usage: " << programName << " -c " + constants::IDX_GRAPH_CMD + " -e <text graph file>"
              << " -n <coordinates file>\n-o <binary graph output file> -s <stats file> [-r <network name>]\n"
              << "The binary graph takes its network name from the coordinates file, or from -r\n";
}
