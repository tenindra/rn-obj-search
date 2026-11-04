// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _GENERATESAMPLESETCOMMAND_H
#define _GENERATESAMPLESETCOMMAND_H

#include "Command.h"
#include "../common.h"

#include <vector>

class GenerateSampleSetCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);
        void showSetTypeUsage(std::string setType, std::string programName);

};

#endif // _GENERATESAMPLESETCOMMAND_H