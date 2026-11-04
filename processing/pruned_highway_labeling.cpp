// Copyright (c) 2015 Yuki Kawata
// SPDX-License-Identifier: MIT

// Changed from the original construction: both of its searches use a radix heap instead of std::priority_queue, and Prune
// (moved here from the header) tests the pairs of a shared path with AVX-512 where available (PairsWithin). The labels can
// differ from the original's (entries of equal distance leave the queues in another order), but they are exact all the same.
// It also leaves out the Prune calls whose answer is already known, prefetches the labels Prune reads, and reserves room in
// new labels, none of which changes the labels. rn_obj_search_par builds the labels in batches with several threads instead
// (BuildLabelsInBatches). A graph or label file that cannot be opened stops the program instead of leaving empty labels.

#include "pruned_highway_labeling.h"
#include "../queue/RadixHeap.h"

#if defined(__AVX512F__)
#include <immintrin.h>
#endif

#if defined(PARALLEL_BUILDS)
#include <omp.h>
#endif

const int PrunedHighwayLabeling::LEVEL;
const int PrunedHighwayLabeling::INF;
const unsigned PrunedHighwayLabeling::GUARD;
const unsigned PrunedHighwayLabeling::PATH_MASK;
const unsigned PrunedHighwayLabeling::NUM_MASK;

#if defined(__AVX512F__)
// Tests every pair rather than merging them: up to 16 pairs of the longer list at a time against each pair of the other.
// The answer is the merge's. ConstructLabel adds a vertex's pairs on a path in increasing order of position and in
// nondecreasing order of position plus distance, and skips a pair whose distance is at least the previous pair's plus the
// distance between their positions. So the distances of any two pairs of a label on a path differ by at most the distance
// between their positions, and the merge, which tests each pair against the nearest pairs of the other label on either
// side, finds the smallest sum of any two pairs.
static inline bool PairsWithinAvx512(const std::pair<int, int> *vp, unsigned vn, const std::pair<int, int> *wp, unsigned wn, int time) {
    if (vn > wn) {
        std::swap(vp, wp);
        std::swap(vn, wn);
    }
    const __m512i limit = _mm512_set1_epi32(time);
    const __m512i evenIdx = _mm512_set_epi32(30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0);
    const __m512i oddIdx = _mm512_set_epi32(31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5, 3, 1);
    for (unsigned jb = 0; jb < wn; jb += 16) {
        unsigned n = std::min(16u, wn - jb);
        __mmask16 lanes = static_cast<__mmask16>((1u << n) - 1);
        __mmask16 lowInts = static_cast<__mmask16>(n >= 8 ? 0xFFFFu : (1u << (2 * n)) - 1);
        __mmask16 highInts = static_cast<__mmask16>(n > 8 ? (1u << (2 * (n - 8))) - 1 : 0);
        const int *block = reinterpret_cast<const int *>(wp + jb);
        __m512i low = _mm512_maskz_loadu_epi32(lowInts, block);
        __m512i high = _mm512_maskz_loadu_epi32(highInts, block + 16);
        __m512i wpos = _mm512_permutex2var_epi32(low, evenIdx, high);
        __m512i wdist = _mm512_permutex2var_epi32(low, oddIdx, high);
        for (unsigned i = 0; i < vn; ++i) {
            __m512i sum = _mm512_add_epi32(_mm512_add_epi32(_mm512_set1_epi32(vp[i].second), wdist),
                                           _mm512_abs_epi32(_mm512_sub_epi32(_mm512_set1_epi32(vp[i].first), wpos)));
            if (_mm512_mask_cmple_epi32_mask(lanes, sum, limit)) return true;
        }
    }
    return false;
}
#endif

// Whether a pair of v's label and a pair of w's label on the same path (position on the path and distance, in increasing order
// of position) give a distance of at most time
static inline bool PairsWithin(const std::pair<int, int> *vp, unsigned vn, const std::pair<int, int> *wp, unsigned wn, int time) {
#if defined(__AVX512F__)
    // With fewer pairs than this the merge is faster
    if (vn + wn >= 8) return PairsWithinAvx512(vp, vn, wp, wn, time);
#endif
    for (unsigned i = 0, j = 0; ; ) {
        if (vp[i].first == wp[j].first) {
            if (vp[i].second + wp[j].second <= time) return true;
            i++;
            j++;
            if (i == vn || j == wn) return false;
        } else if (vp[i].first < wp[j].first) {
            if (vp[i].second + wp[j].second + wp[j].first - vp[i].first <= time) return true;
            i++;
            if (i == vn) return false;
        } else {
            if (vp[i].second + wp[j].second + vp[i].first - wp[j].first <= time) return true;
            j++;
            if (j == wn) return false;
        }
    }
}

bool PrunedHighwayLabeling::Prune(std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > &lv, std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > &lw, int time) {
    unsigned vfnum = 0;
    unsigned vsnum = 0;
    unsigned wfnum = 0;
    unsigned wsnum = 0;
    while (true) {
        if (vfnum == lv.first.size() || wfnum == lw.first.size()) return false;
        if ((lv.first[vfnum] & PATH_MASK) == (lw.first[wfnum] & PATH_MASK)) {
            if (PairsWithin(&lv.second[vsnum], lv.first[vfnum] & NUM_MASK, &lw.second[wsnum], lw.first[wfnum] & NUM_MASK, time)) return true;
            vsnum += lv.first[vfnum] & NUM_MASK;
            vfnum++;
            wsnum += lw.first[wfnum] & NUM_MASK;
            wfnum++;
        } else if ((lv.first[vfnum] & PATH_MASK) < (lw.first[wfnum] & PATH_MASK)) {
            vsnum += lv.first[vfnum] & NUM_MASK;
            vfnum++;
        } else {
            wsnum += lw.first[wfnum] & NUM_MASK;
            wfnum++;
        }
    }
}

void PrunedHighwayLabeling::ConstructLabel(const char *file) {
    Free();
    
    // read graph
    // each line should contain two vertices, travel time and geometrical length
    // treat a graph as an undirected graph
    load_time = -GetTime();
    std::vector <std::vector <edge> > graph;
    {
        std::vector <road> edges;
        FILE *in = fopen(file, "r");
        if (in == NULL) {
            fprintf(stderr, "Can't open the graph file %s\n", file);
            exit(1);
        }
        for (int from, to, time, dist; fscanf(in, "%d %d %d %d", &from, &to, &time, &dist) != EOF; ) {
            V = std::max(V, from + 1);
            V = std::max(V, to + 1);
            edges.push_back((road){from, to, time, dist});
        }
        fclose(in);
        
        graph.resize(V);
        for (size_t i = 0; i < edges.size(); i++) {
            int from = edges[i].from;
            int to = edges[i].to;
            int time = edges[i].time;
            int dist = edges[i].dist;
            int level = GetLevel(time, dist);
            graph[from].push_back((edge){to, time, level});
            graph[to].push_back((edge){from, time, level});
        }
    }
    load_time += GetTime();
    
    // allocate memory for labels
    label = (label_t *)memalign(64, sizeof(label_t) * V);
    if (label == NULL) {
        V = 0;
        fprintf(stderr, "Out of memory\n");
        return;
    }
    for (int v = 0; v < V; v++) {
        label[v].time = 0;
        label[v].path = NULL;
        label[v].cost = NULL;
    }
    
    // contract degree one vertices and order vertices
    construct_time = -GetTime();
    std::vector <std::vector <int> > order(LEVEL);
    {
        // contract degree one vertices
        contract.resize(V, -1);
        for (int v = 0; v < V; v++) {
            if (graph[v].size() == 1) contract[v] = graph[v][0].to;
        }
        for (int v = 0; v < V; v++) {
            if (contract[v] != -1) {
                int to = graph[v][0].to;
                for (std::vector <edge>::iterator it = graph[to].begin(); ; it++) {
                    if (it->to == v) {
                        graph[to].erase(it);
                        break;
                    }
                }
                label[v].time = graph[v][0].time;
                std::vector <edge>().swap(graph[v]);
            }
        }
        
        // order vertices
        std::vector <int> vertex_level(V, LEVEL);
        std::vector <std::vector <std::pair<double, int> > > rank(LEVEL);
        for (int v = 0; v < V; v++) {
            if (contract[v] == -1) {
                for (size_t i = 0; i < graph[v].size(); i++) {
                    vertex_level[v] = std::min(vertex_level[v], graph[v][i].level);
                }
                for (int i = vertex_level[v]; i < LEVEL; i++) rank[i].push_back(std::make_pair((double)rand() / RAND_MAX, v));
            }
        }
        for (int i = 0; i < LEVEL; i++) std::sort(rank[i].begin(), rank[i].end());
        for (int i = 0; i < LEVEL; i++) {
            if (i < LEVEL - 1) {
                for (size_t j = 0; j * 1000 < rank[i].size(); j++) order[i].push_back(rank[i][j].second);
            } else {
                for (size_t j = 0; j < rank[i].size(); j++) order[i].push_back(rank[i][j].second);
            }
        }
    }
    
    // pruned highway labeling
    // select path and add labels for the path
#if defined(PARALLEL_BUILDS)
    BuildLabelsInBatches(graph, order);
#else
    {
        unsigned path_num = 0;
        std::vector <bool> used(V, false);
        std::vector <bool> visited(V, false);
        std::vector <int> cost(V, INF);
        std::vector <int> dist(V);
        std::vector <int> parent(V, -1);
        std::vector <int> child(V, 0);
        std::vector <int> check(V);
        std::vector <int> visit(V);
        std::vector <long long> last(V, -1);
        // the path search in which Prune last returned false for a vertex, and the path and start for which it last returned
        // true for an entry of a vertex
        std::vector <int> not_pruned(V, 0);
        std::vector <long long> pruned(V, -1);
        int search_num = 0;
        std::vector <std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > > tmp_label(V);
        
        for (int i = 0; i < LEVEL; i++) {
            for (size_t j = 0; j < order[i].size(); j++) {
                // select path
                int v = order[i][j];
                if (used[v]) continue;
                int check_num = 0;
                int visit_num = 0;
                cost[v] = 0;
                check[check_num++] = v;
                search_num++;
                // A new heap for every search, as a radix heap keeps the last extracted key
                RadixHeap<int> path_que;
                path_que.push(0, v);
                while (!path_que.empty()) {
                    int time = path_que.top_key();
                    int now = path_que.top_value();
                    path_que.pop();
                    if (cost[now] < time) continue;
                    visited[now] = true;
                    visit[visit_num++] = now;
                    // the neighbours' labels, which Prune reads below
                    for (size_t k = 0; k < graph[now].size(); k++) _mm_prefetch(&tmp_label[graph[now][k].to], _MM_HINT_T0);
                    for (size_t k = 0; k < graph[now].size(); k++) {
                        int w = graph[now][k].to;
                        if (used[w]) continue;
                        if (cost[w] > time + graph[now][k].time) {
                            if (cost[w] == INF) check[check_num++] = w;
                            cost[w] = time + graph[now][k].time;
                            parent[w] = now;
                            if (graph[now][k].level <= i) {
                                child[w] = 1;
                            } else {
                                child[w] = 0;
                            }
                            // labels do not change during the search, so Prune returns false again at a lower cost
                            if (not_pruned[w] != search_num) {
                                if (Prune(tmp_label[v], tmp_label[w], cost[w])) continue;
                                not_pruned[w] = search_num;
                            }
                            path_que.push(cost[w], w);
                        }
                    }
                }
                for (int k = visit_num - 1; k > 0; k--) child[parent[visit[k]]] += child[visit[k]];
                int end = v;
                while (true) {
                    int next = -1, max_child = 0;
                    for (size_t k = 0; k < graph[end].size(); k++) {
                        int w = graph[end][k].to;
                        if (visited[w] && parent[w] == end && child[w] > max_child) {
                            next = w;
                            max_child = child[w];
                        }
                    }
                    if (next == -1) break;
                    end = next;
                }
                
                //add labels
                if (end != v) {
                    int child_num = 0;
                    RadixHeap<std::pair<int, int> > label_que;
                    while (end != -1) {
                        if (i == LEVEL - 1 || child_num * 20 <= child[end] * 19) {
                            used[end] = true;
                            label_que.push(cost[end], std::make_pair(end, end));
                        }
                        child_num = child[end];
                        end = parent[end];
                    }
                    long long num = (long long)path_num << 22;
                    while (!label_que.empty()) {
                        int time = label_que.top_key();
                        int now = label_que.top_value().first;
                        int start = label_que.top_value().second;
                        label_que.pop();
                        if (last[now] >= num + cost[start]) continue;
                        // labels only grow, so an entry pruned before (at a smaller time) is pruned again
                        if (pruned[now] == num + start) continue;
                        if ((last[now] >> 32 << 32) == num) {
                            if (dist[now] + cost[start] * 2 <= time) continue;
                            if (Prune(tmp_label[start], tmp_label[now], time - cost[start])) {
                                pruned[now] = num + start;
                                continue;
                            }
                        } else {
                            if (!used[now] && Prune(tmp_label[start], tmp_label[now], time - cost[start])) {
                                pruned[now] = num + start;
                                continue;
                            }
                        }
                        last[now] = num + cost[start];
                        dist[now] = time - cost[start] * 2;
                        // room for a few entries, which saves the first reallocations
                        if (tmp_label[now].first.capacity() == 0) {
                            tmp_label[now].first.reserve(4);
                            tmp_label[now].second.reserve(8);
                        }
                        if (tmp_label[now].first.size() == 0 || (tmp_label[now].first.back() & PATH_MASK) != path_num) {
                            tmp_label[now].first.push_back(path_num + 1);
                        } else {
                            tmp_label[now].first.back()++;
                        }
                        tmp_label[now].second.push_back(std::make_pair(cost[start], time - cost[start]));
                        for (size_t k = 0; k < graph[now].size(); k++) {
                            int w = graph[now][k].to;
                            if (!used[w] && last[w] < num + cost[start]) {
                                // the label of w, which Prune reads when the entry leaves the queue
                                _mm_prefetch(&tmp_label[w], _MM_HINT_T0);
                                label_que.push(time + graph[now][k].time, std::make_pair(w, start));
                            }
                        }
                    }
                    path_num += 1 << 10;
                }
                
                for (int k = 0; k < check_num; k++) {
                    visited[check[k]] = false;
                    cost[check[k]] = INF;
                    parent[check[k]] = -1;
                    child[check[k]] = 0;
                }
            }
        }
        
        // convert labels
        for (int v = 0; v < V; v++) {
            if (contract[v] == -1) {
                label[v].path = (unsigned *)memalign(64, sizeof(unsigned) * (tmp_label[v].first.size() + 1));
                label[v].cost = (int *)memalign(64, sizeof(int) * tmp_label[v].second.size() * 2);
                if (label[v].path == NULL || label[v].cost == NULL) {
                    Free();
                    fprintf(stderr, "Out of memory\n");
                    return;
                }
                for (size_t i = 0; i < tmp_label[v].first.size(); i++) label[v].path[i] = (tmp_label[v].first[i] & PATH_MASK) + (tmp_label[v].first[i] & NUM_MASK) * 2;
                label[v].path[tmp_label[v].first.size()] = GUARD;
                for (size_t i = 0; i < tmp_label[v].second.size(); i++) {
                    label[v].cost[i * 2] = tmp_label[v].second[i].first;
                    label[v].cost[i * 2 + 1] = tmp_label[v].second[i].second;
                }
                std::vector <unsigned>().swap(tmp_label[v].first);
                std::vector <std::pair<int, int> >().swap(tmp_label[v].second);
            }
        }
        for (int v = 0; v < V; v++) {
            if (contract[v] != -1) {
                label[v].path = label[contract[v]].path;
                label[v].cost = label[contract[v]].cost;
            }
        }
    }
#endif
    construct_time += GetTime();
}

#if defined(PARALLEL_BUILDS)
// ConstructLabel's main loop in batches of paths with several threads (OMP_NUM_THREADS sets their number), for
// rn_obj_search_par. The paths are selected one after another as in ConstructLabel, but in batches: one path at first, then
// a quarter of the paths selected so far, at most 256. The paths of a batch then add their labels in parallel, each seeing
// only the labels from before the batch, and their entries are appended in path order after the batch. A path's label
// search treats as used only the vertices of its own and earlier paths, as in ConstructLabel, where the batch's later paths
// do not exist yet. The labels are exact, but they are not ConstructLabel's: a path does not see the entries of the batch's
// earlier paths, so it keeps some entries that ConstructLabel prunes (1.9-3.8% more on the networks tested). They do not
// depend on the number of threads.
void PrunedHighwayLabeling::BuildLabelsInBatches(std::vector <std::vector <edge> > &graph, std::vector <std::vector <int> > &order) {
    const int threads = omp_get_max_threads();
    unsigned path_num = 0;
    std::vector <int> path_of(V, -1); // the number of the path a vertex was selected on (the vertex is used), -1 if none
    std::vector <int> position(V, 0); // its position on that path
    std::vector <bool> visited(V, false);
    std::vector <int> cost(V, INF);
    std::vector <int> parent(V, -1);
    std::vector <int> child(V, 0);
    std::vector <int> check(V);
    std::vector <int> visit(V);
    std::vector <int> not_pruned(V, 0);
    int search_num = 0;
    std::vector <std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > > tmp_label(V);
    // the label searches' last, dist and pruned of ConstructLabel, one each per thread; their values are keyed by the path,
    // so the values a thread left for other paths are never taken for the current path's
    std::vector <std::vector <long long> > last(threads, std::vector <long long>(V, -1));
    std::vector <std::vector <int> > dist(threads, std::vector <int>(V));
    std::vector <std::vector <long long> > pruned(threads, std::vector <long long>(V, -1));
    struct entry_t {
        int vertex;
        int position;
        int distance;
    };
    std::vector <std::pair<unsigned, std::vector <int> > > batch; // each path's path_num and its vertices that add labels
    std::vector <std::vector <entry_t> > entries; // the entries each path of the batch adds, in the order it adds them
    
    int level = 0;
    size_t next = 0;
    bool done = false;
    while (!done) {
        // select the batch's paths, with the labels from before the batch
        size_t batch_size = std::max<size_t>(1, std::min<size_t>(256, (path_num >> 10) / 4));
        batch.clear();
        while (batch.size() < batch_size) {
            while (level < LEVEL && next == order[level].size()) {
                level++;
                next = 0;
            }
            if (level == LEVEL) {
                done = true;
                break;
            }
            int i = level;
            int v = order[i][next++];
            if (path_of[v] != -1) continue;
            int check_num = 0;
            int visit_num = 0;
            cost[v] = 0;
            check[check_num++] = v;
            search_num++;
            RadixHeap<int> path_que;
            path_que.push(0, v);
            while (!path_que.empty()) {
                int time = path_que.top_key();
                int now = path_que.top_value();
                path_que.pop();
                if (cost[now] < time) continue;
                visited[now] = true;
                visit[visit_num++] = now;
                for (size_t k = 0; k < graph[now].size(); k++) _mm_prefetch(&tmp_label[graph[now][k].to], _MM_HINT_T0);
                for (size_t k = 0; k < graph[now].size(); k++) {
                    int w = graph[now][k].to;
                    if (path_of[w] != -1) continue;
                    if (cost[w] > time + graph[now][k].time) {
                        if (cost[w] == INF) check[check_num++] = w;
                        cost[w] = time + graph[now][k].time;
                        parent[w] = now;
                        if (graph[now][k].level <= i) {
                            child[w] = 1;
                        } else {
                            child[w] = 0;
                        }
                        if (not_pruned[w] != search_num) {
                            if (Prune(tmp_label[v], tmp_label[w], cost[w])) continue;
                            not_pruned[w] = search_num;
                        }
                        path_que.push(cost[w], w);
                    }
                }
            }
            for (int k = visit_num - 1; k > 0; k--) child[parent[visit[k]]] += child[visit[k]];
            int end = v;
            while (true) {
                int next_vertex = -1, max_child = 0;
                for (size_t k = 0; k < graph[end].size(); k++) {
                    int w = graph[end][k].to;
                    if (visited[w] && parent[w] == end && child[w] > max_child) {
                        next_vertex = w;
                        max_child = child[w];
                    }
                }
                if (next_vertex == -1) break;
                end = next_vertex;
            }
            if (end != v) {
                batch.push_back(std::make_pair(path_num, std::vector <int>()));
                int child_num = 0;
                while (end != -1) {
                    if (i == LEVEL - 1 || child_num * 20 <= child[end] * 19) {
                        path_of[end] = path_num >> 10;
                        position[end] = cost[end];
                        batch.back().second.push_back(end);
                    }
                    child_num = child[end];
                    end = parent[end];
                }
                path_num += 1 << 10;
            }
            for (int k = 0; k < check_num; k++) {
                visited[check[k]] = false;
                cost[check[k]] = INF;
                parent[check[k]] = -1;
                child[check[k]] = 0;
            }
        }
        if (batch.empty()) break;
        
        // add each path's labels in parallel, from the labels before the batch
        if (entries.size() < batch.size()) entries.resize(batch.size());
        #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (size_t b = 0; b < batch.size(); b++) {
            std::vector <long long> &path_last = last[omp_get_thread_num()];
            std::vector <int> &path_dist = dist[omp_get_thread_num()];
            std::vector <long long> &path_pruned = pruned[omp_get_thread_num()];
            std::vector <entry_t> &path_entries = entries[b];
            path_entries.clear();
            int path_index = batch[b].first >> 10;
            long long num = (long long)batch[b].first << 22;
            RadixHeap<std::pair<int, int> > label_que;
            for (size_t k = 0; k < batch[b].second.size(); k++) label_que.push(position[batch[b].second[k]], std::make_pair(batch[b].second[k], batch[b].second[k]));
            while (!label_que.empty()) {
                int time = label_que.top_key();
                int now = label_que.top_value().first;
                int start = label_que.top_value().second;
                label_que.pop();
                int pos = position[start];
                bool same_path = (path_last[now] >> 32 << 32) == num;
                if (same_path && path_last[now] >= num + pos) continue;
                if (path_pruned[now] == num + start) continue;
                if (same_path) {
                    if (path_dist[now] + pos * 2 <= time) continue;
                    if (Prune(tmp_label[start], tmp_label[now], time - pos)) {
                        path_pruned[now] = num + start;
                        continue;
                    }
                } else {
                    if (!(path_of[now] != -1 && path_of[now] <= path_index) && Prune(tmp_label[start], tmp_label[now], time - pos)) {
                        path_pruned[now] = num + start;
                        continue;
                    }
                }
                path_last[now] = num + pos;
                path_dist[now] = time - pos * 2;
                path_entries.push_back((entry_t){now, pos, time - pos});
                for (size_t k = 0; k < graph[now].size(); k++) {
                    int w = graph[now][k].to;
                    if (path_of[w] != -1 && path_of[w] <= path_index) continue;
                    if ((path_last[w] >> 32 << 32) == num && path_last[w] >= num + pos) continue;
                    _mm_prefetch(&tmp_label[w], _MM_HINT_T0);
                    label_que.push(time + graph[now][k].time, std::make_pair(w, start));
                }
            }
        }
        
        // append the batch's entries in path order, each thread those of its own range of vertices
        #pragma omp parallel num_threads(batch.size() > 1 ? threads : 1)
        {
            int first_vertex = (long long)V * omp_get_thread_num() / omp_get_num_threads();
            int end_vertex = (long long)V * (omp_get_thread_num() + 1) / omp_get_num_threads();
            for (size_t b = 0; b < batch.size(); b++) {
                unsigned num = batch[b].first;
                for (size_t e = 0; e < entries[b].size(); e++) {
                    const entry_t &entry = entries[b][e];
                    if (entry.vertex < first_vertex || entry.vertex >= end_vertex) continue;
                    std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > &l = tmp_label[entry.vertex];
                    if (l.first.capacity() == 0) {
                        l.first.reserve(4);
                        l.second.reserve(8);
                    }
                    if (l.first.size() == 0 || (l.first.back() & PATH_MASK) != num) {
                        l.first.push_back(num + 1);
                    } else {
                        l.first.back()++;
                    }
                    l.second.push_back(std::make_pair(entry.position, entry.distance));
                }
            }
        }
    }
    
    // convert labels (as ConstructLabel)
    for (int v = 0; v < V; v++) {
        if (contract[v] == -1) {
            label[v].path = (unsigned *)memalign(64, sizeof(unsigned) * (tmp_label[v].first.size() + 1));
            label[v].cost = (int *)memalign(64, sizeof(int) * tmp_label[v].second.size() * 2);
            if (label[v].path == NULL || label[v].cost == NULL) {
                Free();
                fprintf(stderr, "Out of memory\n");
                return;
            }
            for (size_t i = 0; i < tmp_label[v].first.size(); i++) label[v].path[i] = (tmp_label[v].first[i] & PATH_MASK) + (tmp_label[v].first[i] & NUM_MASK) * 2;
            label[v].path[tmp_label[v].first.size()] = GUARD;
            for (size_t i = 0; i < tmp_label[v].second.size(); i++) {
                label[v].cost[i * 2] = tmp_label[v].second[i].first;
                label[v].cost[i * 2 + 1] = tmp_label[v].second[i].second;
            }
            std::vector <unsigned>().swap(tmp_label[v].first);
            std::vector <std::pair<int, int> >().swap(tmp_label[v].second);
        }
    }
    for (int v = 0; v < V; v++) {
        if (contract[v] != -1) {
            label[v].path = label[contract[v]].path;
            label[v].cost = label[contract[v]].cost;
        }
    }
}
#endif

void PrunedHighwayLabeling::LoadLabel(const char *file) {
    Free();
    
    FILE *in = fopen(file, "rb");
    if (in == NULL) {
        fprintf(stderr, "Can't open the label file %s\n", file);
        exit(1);
    }
    // A short read means the file is truncated or corrupt: drop what was read (before the labels are
    // allocated, Free() cannot run, as it reads them)
    auto readFailed = [&]() {
        if (label == NULL) {
            V = 0;
        } else {
            Free();
        }
        fclose(in);
        fprintf(stderr, "Can't read the label file\n");
    };
    if (fread(&V, sizeof(int), 1, in) != 1) {
        readFailed();
        return;
    }
    contract.resize(V);
    if (fread(&contract[0], sizeof(int), V, in) != static_cast<std::size_t>(V)) {
        readFailed();
        return;
    }
    label = (label_t *)memalign(64, sizeof(label_t) * V);
    if (label == NULL) {
        V = 0;
        fprintf(stderr, "Out of memory\n");
        return;
    }
    for (int v = 0; v < V; v++) {
        label[v].time = 0;
        label[v].path = NULL;
        label[v].cost = NULL;
    }
    for (int v = 0; v < V; v++) {
        if (contract[v] != -1) {
            if (fread(&label[v].time, sizeof(int), 1, in) != 1) {
                readFailed();
                return;
            }
        } else {
            int fnum, snum;
            if (fread(&fnum, sizeof(int), 1, in) != 1 || fread(&snum, sizeof(int), 1, in) != 1) {
                readFailed();
                return;
            }
            label[v].path = (unsigned *)memalign(64, sizeof(unsigned) * fnum);
            label[v].cost = (int *)memalign(64, sizeof(int) * snum);
            if (label[v].path == NULL || label[v].cost == NULL) {
                Free();
                fprintf(stderr, "Out of memory\n");
                return;
            }
            if (fread(label[v].path, sizeof(unsigned), fnum, in) != static_cast<std::size_t>(fnum)
                || fread(label[v].cost, sizeof(int), snum, in) != static_cast<std::size_t>(snum)) {
                readFailed();
                return;
            }
        }
    }
    for (int v = 0; v < V; v++) {
        if (contract[v] != -1) {
            label[v].path = label[contract[v]].path;
            label[v].cost = label[contract[v]].cost;
        }
    }
    fclose(in);
}

void PrunedHighwayLabeling::StoreLabel(const char *file) {
    FILE *out = fopen(file, "wb");
    if (out == NULL) {
        fprintf(stderr, "Can't open the label file %s\n", file);
        exit(1);
    }
    fwrite(&V, sizeof(int), 1, out);
    fwrite(&contract[0], sizeof(int), V, out);
    for (int v = 0; v < V; v++) {
        if (contract[v] != -1) {
            fwrite(&label[v].time, sizeof(int), 1, out);
        } else {
            int fnum = 1, snum = 0;
            for (int i = 0; ; i++) {
                if (label[v].path[i] == GUARD) break;
                fnum++;
                snum += label[v].path[i] & NUM_MASK;
            }
            fwrite(&fnum, sizeof(int), 1, out);
            fwrite(&snum, sizeof(int), 1, out);
            fwrite(label[v].path, sizeof(unsigned), fnum, out);
            fwrite(label[v].cost, sizeof(int), snum, out);
        }
    }
    fclose(out);
}

void PrunedHighwayLabeling::Free(void) {
    for (int v = 0; v < V; v++) {
        if (contract[v] == -1) {
            free(label[v].path);
            free(label[v].cost);
            label[v].path = NULL;
            label[v].cost = NULL;
        }
    }
    free(label);
    V = 0;
    label = NULL;
}

void PrunedHighwayLabeling::Statistics(void) {
    long long sum_label = 0, sum_memory = 0;
    
    for (int v = 0; v < V; v++) {
        if (contract[v] != -1) continue;
        for (int i = 0; ; i++) {
            if (label[v].path[i] == GUARD) break;
            sum_label += (label[v].path[i] & NUM_MASK) / 2;
            sum_memory += sizeof(unsigned);
        }
        sum_memory += sizeof(unsigned);
    }
    sum_memory += sizeof(int) * sum_label * 2;
    
    printf("Load Time : %lf sec\n", load_time);
    printf("Construct Time : %lf sec\n", construct_time);
    printf("Average Label Size : %lld\n", sum_label / V);
    printf("Memory Size : %lld byte\n", sum_memory);
}

double PrunedHighwayLabeling::computeIndexSize() {
    long long sum_label = 0, sum_memory = 0;
    for (int v = 0; v < V; v++) {
        if (contract[v] != -1) continue;
        for (int i = 0; ; i++) {
            if (label[v].path[i] == GUARD) break;
            sum_label += (label[v].path[i] & NUM_MASK) / 2;
            sum_memory += sizeof(unsigned);
        }
        sum_memory += sizeof(unsigned);
    }
    sum_memory += sizeof(int) * sum_label * 2;
    double memoryUsage = sum_memory;
    return memoryUsage/(1024*1024); // in MB
}

double PrunedHighwayLabeling::getConstructionTime() {
    return construct_time;
}
