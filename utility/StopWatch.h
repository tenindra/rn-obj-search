// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _STOPWATCH_H
#define _STOPWATCH_H

#include <string>
#include <chrono>

class StopWatch {
    
public:
    StopWatch();
    void start();
    void stop();
    void pause();
    void resume();
    void reset();
    double getTimeNs();
    double getTimeUs();
    double getTimeMs();
    
private:
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::high_resolution_clock::time_point p_t1;
    std::chrono::high_resolution_clock::time_point p_t2;
    std::chrono::nanoseconds durationNs;
    std::chrono::nanoseconds pausedDurationNs;
};

#endif // _STOPWATCH_H
