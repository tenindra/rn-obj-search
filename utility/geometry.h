// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include "../common.h"

#include <cmath>

namespace geometry
{
    inline EuclideanDistance getEuclideanDist(
        Coordinate x1, Coordinate y1, Coordinate x2, Coordinate y2)
    {
        return std::sqrt(std::pow(x2-x1, 2) + std::pow(y2-y1, 2));
    }

    inline unsigned int maxPower2(unsigned int x)
    {
        RN_DEBUG_ASSERT(x > 0);
        x = x | (x >> 1);
        x = x | (x >> 2);
        x = x | (x >> 4);
        x = x | (x >> 8);
        x = x | (x >> 16);
        return x - (x >> 1);
    }
}

#endif // _GEOMETRY_H
