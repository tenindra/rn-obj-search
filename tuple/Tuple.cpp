// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "Tuple.h"

#include <sstream>

Tuple::Tuple(std::string networkName, int vertices, int edges)
{
    this->networkName = networkName;
    this->vertices = vertices;
    this->edges = edges;
    this->supplementaryData = "none";
}

std::string Tuple::getTupleString()
{
    std::stringstream ss;
    ss << this->networkName << "\t";
    ss << this->vertices << "\t";
    ss << this->edges << "\t";
    ss << this->supplementaryData;
    return ss.str();
}

std::string Tuple::getMultilineTupleString()
{
    std::stringstream ss;
    ss << "Road Network: " << this->networkName << " (Vertices=" << this->vertices << ",Edges=" << this->edges << ")\n";;
    ss << "Supplementary Data: " << this->supplementaryData << "\n";
    return ss.str();
}

void Tuple::addSupplementaryFields(std::string desc, std::string data)
{
    if (this->supplementaryData == "none") {
        this->supplementaryData = "";
    }
    if (this->supplementaryData != "") {
        this->supplementaryData += ",";
    }
    this->supplementaryData += desc + "=" + data;
}

std::string Tuple::getAdditionalFields()
{
    std::stringstream ss;
    for(std::size_t i = 0; i < this->additionalFields.size(); ++i) {
        ss << this->additionalFields[i];
        if (i != this->additionalFields.size()-1) {
            // If not last field we add a tab
            ss << "\t";
        }
    }
    return ss.str();
}

void Tuple::setAdditionalFields(std::vector<std::string>& fields)
{
    for(std::size_t i = 0; i < fields.size(); ++i) {
        this->additionalFields.push_back(fields[i]);
    }
}

