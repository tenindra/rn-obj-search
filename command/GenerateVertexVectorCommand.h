// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _GENERATEVERTEXVECTORSCOMMAND_H
#define _GENERATEVERTEXVECTORSCOMMAND_H

#include "Command.h"

class GenerateVertexVectorCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);
        void showSetTypeUsage(std::string vectorType, std::string programName);
        
    private:
         // Set Types
        std::string const RANDOM = "random";
        std::string const SUBNET = "subnet";
       
};

#endif // _GENERATEVERTEXVECTORSCOMMAND_H
