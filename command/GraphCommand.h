// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _GRAPHCOMMAND_H
#define _GRAPHCOMMAND_H

#include "Command.h"

class GraphCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);

};

#endif // _GRAPHCOMMAND_H