// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _INDEXTUPLE_H
#define _INDEXTUPLE_H

#include "Tuple.h"

#include <string>

class IndexTuple: public Tuple {
    
public:
    std::string getTupleString();
    std::string getMultilineTupleString();
    IndexTuple(std::string networkName, int vertices, int edges, std::string indexMethod, double constructionTimeMs, double memoryUsageMB);
    
protected:
    std::string indexMethod;
    double constructionTimeMs;
    double memoryUsageMB;
};

#endif // _INDEXTUPLE_H
