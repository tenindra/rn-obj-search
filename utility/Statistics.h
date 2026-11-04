// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _STATISTICS_H
#define _STATISTICS_H

#include "../tuple/Tuple.h"

#include <string>
#include <map>

class Statistics {
    
    public:
        Statistics();
        void clear();
        void initialiseStatistic(std::string desc, double value);
        void incrementStatistic(std::string desc, double value);
        void mergeStatistics(Statistics& inputStatistics);
        void normalizeStatistics(int normalizationFactor);
        void populateTupleFields(Tuple& tuple, int decimalPlaces = 3);
        double getStatistic(std::string desc);
        
        std::map<std::string,double> statistics;
        
};

#endif // _STATISTICS_H
