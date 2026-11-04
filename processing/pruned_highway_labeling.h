// Copyright (c) 2015 Yuki Kawata
// SPDX-License-Identifier: MIT

#ifndef __PRUNED_HIGHWAY_LABELING_H__
#define __PRUNED_HIGHWAY_LABELING_H__

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <queue>
#include <algorithm>
#include <malloc.h>
#include <xmmintrin.h>
#include <sys/time.h>

class PrunedHighwayLabeling {
    public:
    
    PrunedHighwayLabeling() : V(0), label(NULL), load_time(0), construct_time(0) {}
    ~PrunedHighwayLabeling() { Free(); };
    
    void ConstructLabel(const char *file);
    void LoadLabel(const char *file);
    void StoreLabel(const char *file);
    double computeIndexSize();
    double getConstructionTime();
    
    inline int Query(int v, int w) {
        if (v == w) return 0;
        if (contract[v] == w || contract[w] == v || (contract[v] != -1 && contract[v] == contract[w])) return label[v].time + label[w].time;
        int time = INF;
        unsigned *vf = label[v].path;
        int *vs = label[v].cost;
        unsigned *wf = label[w].path;
        int *ws = label[w].cost;
        
        _mm_prefetch(vf, _MM_HINT_T0);
        _mm_prefetch(vs, _MM_HINT_T0);
        _mm_prefetch(wf, _MM_HINT_T0);
        _mm_prefetch(ws, _MM_HINT_T0);
        
        while (true) {
            if ((*vf & PATH_MASK) == (*wf & PATH_MASK)) {
                if (*vf == GUARD) return time + label[v].time + label[w].time;
                for (unsigned i = 0, j = 0; ; ) {
                    if (*(vs + i) == *(ws + j)) {
                        ChangeMin(time, *(vs + i + 1) + *(ws + j + 1));
                        i += 2;
                        j += 2;
                        if (i == (*vf & NUM_MASK) || j == (*wf & NUM_MASK)) break;
                    } else if (*(vs + i) < *(ws + j)) {
                        ChangeMin(time, *(vs + i + 1) + *(ws + j + 1) + *(ws + j) - *(vs + i));
                        i += 2;
                        if (i == (*vf & NUM_MASK)) break;
                    } else {
                        ChangeMin(time, *(vs + i + 1) + *(ws + j + 1) + *(vs + i) - *(ws + j));
                        j += 2;
                        if (j == (*wf & NUM_MASK)) break;
                    }
                }
                vs += *vf & NUM_MASK;
                vf++;
                ws += *wf & NUM_MASK;
                wf++;
            } else if ((*vf & PATH_MASK) < (*wf & PATH_MASK)) {
                if (*wf == GUARD) return time + label[v].time + label[w].time;
                vs += *vf & NUM_MASK;
                vf++;
            } else {
                if (*vf == GUARD) return time + label[v].time + label[w].time;
                ws += *wf & NUM_MASK;
                wf++;
            }
        }
    }
    
    void Free(void);
    
    int NumVertices(void) { return V; }
    void Statistics(void);
    
    private:
    
    static const int LEVEL = 4;
    static const int INF = 1000000000;
    static const unsigned GUARD = 0xFFFFFFFFU;
    static const unsigned PATH_MASK = 0xFFFFFC00U;
    static const unsigned NUM_MASK = 0x000003FFU;
    
    struct road {
        int from;
        int to;
        int time;
        int dist;
    };
    
    struct edge {
        int to;
        int time;
        int level;
    };
    
    struct label_t {
        int time;
        unsigned *path;
        int *cost;
    } __attribute__((aligned(64)));
    
    int V;
    label_t *label;
    std::vector <int> contract;
    double load_time, construct_time;
    
    double GetTime(void) {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
    
    int GetLevel(int time, int dist) {
        if (time == 0) return 0;
        if (dist > time * 0.9) return 0;
        if (dist > time * 0.7) return 1;
        if (dist > time * 0.5) return 2;
        return 3;
    }
    
    inline void ChangeMin(int &now, int next) { if (next < now) now = next; }
    
    // Defined in pruned_highway_labeling.cpp, its only user
    bool Prune(std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > &lv, std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > &lw, int time);
    
#if defined(PARALLEL_BUILDS)
    // ConstructLabel's labels in batches of paths with several threads, only in rn_obj_search_par
    void BuildLabelsInBatches(std::vector <std::vector <edge> > &graph, std::vector <std::vector <int> > &order);
#endif
};

#endif
