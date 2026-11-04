// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "Statistics.h"
#include "../debug.h"

#include <sstream>

Statistics::Statistics() {}

void Statistics::clear()
{
    this->statistics.clear();
}

void Statistics::initialiseStatistic(std::string desc, double value)
{
    statistics[desc] = value;
}

// Assumes it has already been intialised
void Statistics::incrementStatistic(std::string desc, double value)
{
    RN_DEBUG_ASSERT(statistics.find(desc) != statistics.end() && "Statistic to be increment has not been initialised");
    statistics[desc] += value;
}

void Statistics::mergeStatistics(Statistics& inputStatistics)
{
    for (auto it = inputStatistics.statistics.begin(); it != inputStatistics.statistics.end(); ++it) {
        if (this->statistics.find(it->first) == this->statistics.end()) {
            this->statistics[it->first] = inputStatistics.statistics[it->first];
        } else {
            this->statistics[it->first] += inputStatistics.statistics[it->first];
        }
    }
}

void Statistics::normalizeStatistics(int normalizationFactor)
{
    for (auto it = this->statistics.begin(); it != this->statistics.end(); ++it) {
        this->statistics[it->first] = it->second/normalizationFactor;
    }
}

void Statistics::populateTupleFields(Tuple& tuple, int decimalPlaces)
{
    std::stringstream ss;
    for (auto it = this->statistics.begin(); it != this->statistics.end(); ++it) {
        ss.str("");
        ss.precision(decimalPlaces);
        ss << std::fixed;
        ss << it->second;
        tuple.addSupplementaryFields(it->first,ss.str());
    }
}

double Statistics::getStatistic(std::string desc)
{
    RN_DEBUG_ASSERT(statistics.find(desc) != statistics.end() && "Statistic to return has not been initialised");
    return statistics[desc];
}
