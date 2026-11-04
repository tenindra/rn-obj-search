// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <string>

#include "utility/utility.h"
#include "common.h"
#include "command/Command.h"
#include "command/GraphCommand.h"
#include "command/GenerateSampleSetCommand.h"
#include "command/TransformInputCommand.h"
#include "command/LMKExperimentsCommand.h"
#include "command/LMKAggregateExperimentsCommand.h"
#include "command/LMKFarthestExperimentsCommand.h"
#include "command/LMKRangeExperimentsCommand.h"
#include "command/GenerateVertexVectorCommand.h"
#include "command/VerifyCommand.h"

void showUsage(std::string programName);
Command* getCommand(std::string);

int main(int argc, char* argv[])
{
    if (argc < 3) {
        // Arguments: -c <command>
        if (argc == 2) {
            std::cerr << "Too few arguments!\n\n";
        }
        showUsage(argv[0]);
        exit(1);
    }
    
    /*
     * Determine Command to Execute
     */ 
    Command* command;
    int opt = getopt(argc,argv, "c:");
    switch (opt) {
        case 'c':
            command = getCommand(optarg);
            break;
        default:
            std::cerr << "Unknown option provided!\n\n";
            showUsage(argv[0]);
            exit(1);
    }
    
    if (command != NULL) {
        if (argc == 3) {
            //  Commands are expected to have more that 2 options/arguments
            command->showCommandUsage(argv[0]);
        } else {
            command->execute(argc, argv);
        }
        delete command;
    } else {
        std::cerr << "Command not found!\n\n";
        showUsage(argv[0]);
        exit(1);
    }
    
}

void showUsage(std::string programName) {
    std::cerr << "Usage: " << programName << " -c <command>\n\n"
              << "Index Commands:\n"
              << utility::getFormattedUsageString(constants::IDX_GRAPH_CMD,"Fast in-memory graph") + "\n"
              << "Other Commands:\n"
              << utility::getFormattedUsageString(constants::SAMPLE_SET_CMD,"Generate different types of object sets") + "\n"
              << utility::getFormattedUsageString(constants::VERTEX_VECTOR_CMD,"Generate different types of query vertex sets") + "\n"
              << utility::getFormattedUsageString(constants::TRANSFORM_INPUT_CMD,"Process graphs from known input data sources") + "\n"
              << utility::getFormattedUsageString(constants::LMK_EXPERIMENTS_CMD,"Run landmark knn experiments") + "\n"
              << utility::getFormattedUsageString(constants::LMK_AGG_EXPERIMENTS_CMD,"Run landmark aknn experiments") + "\n"
              << utility::getFormattedUsageString(constants::LMK_FARTHEST_EXPERIMENTS_CMD,"Run landmark kfn experiments") + "\n"
              << utility::getFormattedUsageString(constants::LMK_RANGE_EXPERIMENTS_CMD,"Run landmark range experiments") + "\n"
              << utility::getFormattedUsageString(constants::VERIFY_CMD,"Verify query results against ground truth") + "\n\n"
              << "Note: Execute commands without parameters to view usage and all options\n";
}

Command* getCommand(std::string commandName) {
    /*
     * Determine Command
     */ 
    Command* command = NULL;
    if (commandName == constants::IDX_GRAPH_CMD) {
        command = new GraphCommand();
    } else if (commandName == constants::SAMPLE_SET_CMD) {
        command = new GenerateSampleSetCommand();
    } else if (commandName == constants::TRANSFORM_INPUT_CMD) {
        command = new TransformInputCommand();
    } else if (commandName == constants::LMK_EXPERIMENTS_CMD) {
        command = new LMKExperimentsCommand();
    } else if (commandName == constants::LMK_AGG_EXPERIMENTS_CMD) {
        command = new LMKAggregateExperimentsCommand();
    } else if (commandName == constants::LMK_FARTHEST_EXPERIMENTS_CMD) {
        command = new LMKFarthestExperimentsCommand();
    } else if (commandName == constants::LMK_RANGE_EXPERIMENTS_CMD) {
        command = new LMKRangeExperimentsCommand();
    } else if (commandName == constants::VERTEX_VECTOR_CMD) {
        command = new GenerateVertexVectorCommand();
    } else if (commandName == constants::VERIFY_CMD) {
        command = new VerifyCommand();
    }
    return command;
}
