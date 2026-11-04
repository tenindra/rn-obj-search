// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _RADIXHEAP_H
#define _RADIXHEAP_H

#include "../debug.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>

// A monotone radix heap for unsigned 32-bit keys (EdgeWeight), the priority queue of the Dijkstra searches.

template <typename Value>
class RadixHeap {

    public:
        RadixHeap() {
            this->mins.fill(std::numeric_limits<std::uint32_t>::max());
        }
        ~RadixHeap() {
            for (Bucket& b: this->buckets) {
                delete[] b.begin;
            }
        }
        RadixHeap(const RadixHeap&) = delete;
        RadixHeap& operator=(const RadixHeap&) = delete;

        void push(std::uint32_t key, const Value& value) {
            RN_DEBUG_ASSERT(this->last <= key); // The heap is monotone: a pushed key may not be below the last key pulled
            const unsigned k = this->bucketIndex(key);
            Bucket& b = this->buckets[k];
            if (b.end == b.cap) [[unlikely]] {
                this->grow(b);
            }
            b.end->key = key;
            b.end->value = value;
            ++b.end;
            this->mins[k] = std::min(this->mins[k],key);
            this->nonEmpty |= std::uint64_t(1) << k;
            ++this->numEntries;
        }

        std::uint32_t top_key() {
            this->pull();
            return this->last;
        }

        Value& top_value() {
            this->pull();
            return (this->buckets[0].end-1)->value;
        }

        void pop() {
            this->pull();
            Bucket& b = this->buckets[0];
            --b.end;
            if (b.end == b.begin) {
                this->nonEmpty &= ~std::uint64_t(1);
            }
            --this->numEntries;
        }

        bool empty() const {
            return this->numEntries == 0;
        }

        std::size_t size() const {
            return this->numEntries;
        }

    private:
        struct Entry {
            std::uint32_t key;
            Value value;
        };
        struct Bucket {
            Entry* begin = nullptr;
            Entry* end = nullptr;
            Entry* cap = nullptr;
        };

        std::array<Bucket,33> buckets{};
        std::array<std::uint32_t,33> mins;
        std::uint64_t nonEmpty = 0; // Bit k is set if bucket k has entries
        std::uint32_t last = 0;
        std::size_t numEntries = 0;

        unsigned bucketIndex(std::uint32_t key) const {
            return 32 - std::countl_zero(key ^ this->last);
        }

        // Refills bucket 0 from the first non-empty bucket (whose entries all go to lower buckets)
        void pull() {
            RN_DEBUG_ASSERT(this->numEntries > 0);
            if (this->nonEmpty & 1) [[likely]] {
                return;
            }
            const unsigned i = std::countr_zero(this->nonEmpty);
            Bucket& src = this->buckets[i];
            this->last = this->mins[i];
            for (Entry* e = src.begin; e != src.end; ++e) {
                const unsigned k = this->bucketIndex(e->key);
                Bucket& dst = this->buckets[k];
                if (dst.end == dst.cap) [[unlikely]] {
                    this->grow(dst);
                }
                *dst.end = *e;
                ++dst.end;
                this->mins[k] = std::min(this->mins[k],e->key);
                this->nonEmpty |= std::uint64_t(1) << k;
            }
            src.end = src.begin;
            this->mins[i] = std::numeric_limits<std::uint32_t>::max();
            this->nonEmpty &= ~(std::uint64_t(1) << i);
        }

        [[gnu::noinline]] void grow(Bucket& b) {
            const std::size_t used = b.end-b.begin;
            const std::size_t capacity = std::max<std::size_t>(16,2*static_cast<std::size_t>(b.cap-b.begin));
            Entry* storage = new Entry[capacity];
            std::copy(b.begin,b.end,storage);
            delete[] b.begin;
            b.begin = storage;
            b.end = storage+used;
            b.cap = storage+capacity;
        }

};

#endif // _RADIXHEAP_H
