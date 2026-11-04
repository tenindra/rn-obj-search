// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _SERIALIZATION_H
#define _SERIALIZATION_H

#include <string>

namespace serialization
{
    template <typename T>
    T getIndexFromBinaryFile(std::string idxFilePath);

    template <typename T>
    void outputIndexToBinaryFile(T& index,std::string idxFilePath);    
}

#endif // _SERIALIZATION_H