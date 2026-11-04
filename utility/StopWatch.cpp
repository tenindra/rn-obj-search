// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "StopWatch.h"

#include <sstream>

StopWatch::StopWatch()
{
    this->durationNs = this->durationNs.zero();
    this->pausedDurationNs = this->pausedDurationNs.zero();
}

double StopWatch::getTimeNs()
{
    // Assumes resume() has been called before stop() already
    return (this->durationNs-this->pausedDurationNs).count();
}

double StopWatch::getTimeUs()
{
    return this->getTimeNs()/1000;
}

double StopWatch::getTimeMs()
{
    return this->getTimeNs()/1000000;
}

void StopWatch::start()
{
    this->t1 = std::chrono::high_resolution_clock::now();
}

void StopWatch::stop()
{
    this->t2 = std::chrono::high_resolution_clock::now();
    this->durationNs = std::chrono::duration_cast<std::chrono::nanoseconds>(this->t2 - this->t1);
}

void StopWatch::pause()
{
    this->p_t1 = std::chrono::high_resolution_clock::now();
}

void StopWatch::resume()
{
    this->p_t2 = std::chrono::high_resolution_clock::now();
    // Accumulate amount of time that been paused 
    this->pausedDurationNs += std::chrono::duration_cast<std::chrono::nanoseconds>(this->p_t2 - this->p_t1);

}

void StopWatch::reset()
{
    this->durationNs = this->durationNs.zero();
    this->pausedDurationNs = this->pausedDurationNs.zero();
}
