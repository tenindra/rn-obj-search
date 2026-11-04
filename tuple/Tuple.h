// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _TUPLE_H
#define _TUPLE_H

#include <string>
#include <vector>

class Tuple {
    
public:
    virtual std::string getTupleString();
    virtual std::string getMultilineTupleString();
    virtual ~Tuple() {};
    Tuple(std::string networkName, int vertices, int edge);
    void addSupplementaryFields(std::string desc, std::string data);
    std::string getAdditionalFields();
    void setAdditionalFields(std::vector<std::string>& fields);
    
protected:
    std::string networkName;
    int vertices;
    int edges;
    std::string supplementaryData;
    std::vector<std::string> additionalFields;
};

#endif // _TUPLE_H
