// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _VERIFYCOMMAND_H
#define _VERIFYCOMMAND_H

#include "Command.h"
#include "../verification/Verifier.h"

#include <string>

class VerifyCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);
        void showStepUsage(std::string step, std::string programName);
        std::string const STEP_OBJECTS = "objects";
        std::string const STEP_GROUND_TRUTH = "ground_truth";
        std::string const STEP_VERIFY = "verify";
        std::string const STEP_DUMP = "dump";

    private:
        template <typename Query>
        void runStep(Verifier<Query>& verifier, const std::string& step, const VerifyArgs& args, const std::string& methods,
                     const std::string& gtFile, int numThreads);

};

#endif // _VERIFYCOMMAND_H
