// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _BINARYMAXHEAP_H
#define _BINARYMAXHEAP_H

#include "../common.h"

#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>

// Can only be used for primitive types for Element and Key
// (because hash function to unordered_map is not custom)

// Note: Implemented entirely in this header so that it can be used with any
// Key/Element types without explicit instantiation and so calls can be inlined
// (member functions are marked inline because GCC inlines them more readily, which was measurably faster)

template <typename Key, typename Element>
class BinaryMaxHeap {

    public:
        void init(int n);
        void insert(Element element, Key key);
        unsigned int size();
        void clear();
        void print();
        Key getMaxKey();
        Element extractMaxElement();
        // Populate kNNs in increasing order
        void populateKNNs(std::vector<Element>& elements, std::vector<Key>& keys);

    private:
        void siftUp(int heapIndex);
        void siftDown(int heapIndex);
        struct ElementKey {
            Element element;
            Key key;
            ElementKey(Element _element, Key _key): element(_element), key(_key) {}
            // Note: Below is reversed so that larger key corresponds to higher priority
            inline bool operator<(ElementKey const& rhs) const { return key < rhs.key; }
        };

        std::vector<ElementKey> bHeap;

};

template <typename Key, typename Element>
inline void BinaryMaxHeap<Key,Element>::init(int n)
{
    this->bHeap.reserve(n);
}

template <typename Key, typename Element>
inline void BinaryMaxHeap<Key,Element>::insert(Element element, Key key)
{
    // Add element at the end of the heap and siftUp
    // until it is in the correct position
    this->bHeap.push_back(ElementKey(std::move(element),key));
    int lastHeapIdx = bHeap.size()-1;
    this->siftUp(lastHeapIdx);
}

template <typename Key, typename Element>
inline void BinaryMaxHeap<Key,Element>::clear()
{
    this->bHeap.clear();
}

template <typename Key, typename Element>
inline void BinaryMaxHeap<Key,Element>::print()
{
    std::cout << "Not implemented yet" << std::endl;
}

template <typename Key, typename Element>
inline unsigned int BinaryMaxHeap<Key,Element>::size()
{
    return this->bHeap.size();
}

template <typename Key, typename Element>
inline Key BinaryMaxHeap<Key,Element>::getMaxKey()
{
    return this->bHeap[0].key;
}

template <typename Key, typename Element>
inline Element BinaryMaxHeap<Key,Element>::extractMaxElement()
{
    Element maxElement = std::move(this->bHeap[0].element);
    // Move the last element to the root
    int lastHeapIdx = this->bHeap.size()-1;
    if (lastHeapIdx != 0) {
        // siftDown the new root element until it is in the correct position
        this->bHeap[0] = std::move(this->bHeap[lastHeapIdx]);
        this->bHeap.pop_back();
        this->siftDown(0);
    } else {
        this->bHeap.pop_back();
    }
    return maxElement;
}

// This may be destructive depending on implementation to below
template <typename Key, typename Element>
inline void BinaryMaxHeap<Key,Element>::populateKNNs(std::vector<Element>& elements, std::vector<Key>& keys)
{
    // Sort heap key value in increasing order and then copy elements one by one
    std::sort(this->bHeap.begin(),this->bHeap.end());
    for (std::size_t i = 0; i < this->bHeap.size(); ++i) {
        elements.push_back(this->bHeap[i].element);
        keys.push_back(this->bHeap[i].key);
    }
}

/*
 * Private Methods
 */

// Note: Rather than swapping at each level, siftUp and siftDown hold the sifted element aside
// and move parents (resp. children) into the vacant position, then place the sifted element
// once its position is found. The comparisons are the same as swapping, so the order in which
// elements are extracted (including elements with equal keys) is unchanged.

template <typename Key, typename Element>
inline void BinaryMaxHeap<Key,Element>::siftUp(int heapIndex)
{
    ElementKey sifted = std::move(this->bHeap[heapIndex]);
    // If parent is the smaller than the element we are sifting we move the parent down
    // Note: That if the child is larger than the parent, it must be larger
    // than the other child (the other child conforms to the heap property)
    while (heapIndex != 0) {
        int parentIndex = (heapIndex-1)/2;
        if (!(this->bHeap[parentIndex].key < sifted.key)) {
            break;
        }
        this->bHeap[heapIndex] = std::move(this->bHeap[parentIndex]);
        heapIndex = parentIndex;
    }
    this->bHeap[heapIndex] = std::move(sifted);
}

template <typename Key, typename Element>
inline void BinaryMaxHeap<Key,Element>::siftDown(int heapIndex)
{
    int heapSize = this->bHeap.size();
    ElementKey sifted = std::move(this->bHeap[heapIndex]);
    while (true) {
        int leftChildIndex = heapIndex*2 + 1;
        int rightChildIndex = heapIndex*2 + 2;
        int largestChildIdx;

        if (rightChildIndex < heapSize) {
            // If the parent has two children pick the largest (the left child on ties)
            // Note: Selected with arithmetic rather than a branch, as in BinaryMinHeap::siftDown
            largestChildIdx = leftChildIndex + (this->bHeap[rightChildIndex].key > this->bHeap[leftChildIndex].key);
        } else if (leftChildIndex < heapSize) {
            // The parent has only one child (it must be the largest)
            largestChildIdx = leftChildIndex;
        } else {
            break;
        }

        // If the largest child is larger than the sifted element we move the child up
        if (!(sifted.key < this->bHeap[largestChildIdx].key)) {
            break;
        }
        this->bHeap[heapIndex] = std::move(this->bHeap[largestChildIdx]);
        heapIndex = largestChildIdx; // Continue sifting down from the child's old position
    }
    this->bHeap[heapIndex] = std::move(sifted);
}

#endif // _BINARYMAXHEAP_H
