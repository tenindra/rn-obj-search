// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _BINARYHEAP_H
#define _BINARYHEAP_H

#include <algorithm>
#include <iostream>
#include <vector>

// Can only be used for primitive types for Key or complex types
// which have defined comparison operators

// Note: Implemented entirely in this header so that it can be used with any
// Key/Element types without explicit instantiation and so calls can be inlined
// (member functions are marked inline because GCC inlines them more readily, which was measurably faster)

template <typename Key, typename Element>
class BinaryMinHeap {

    public:
        void init(int n);
        void insert(Element element, Key key);
        unsigned int size();
        void clear();
        Key getMinKey();
        Element extractMinElement();
        void print();
        // Populate FkNs in decreasing order
        void populateFkNs(std::vector<Element>& elements, std::vector<Key>& keys);

    private:
        void siftUp(int heapIndex);
        void siftDown(int heapIndex);
        struct ElementKey {
            Element element;
            Key key;
            ElementKey(Element _element, Key _key): element(_element), key(_key) {}
            inline bool operator<(ElementKey const& rhs) const { return key > rhs.key; } // Note: This sort will be in increasing order
        };

        std::vector<ElementKey> bHeap;

};

// Note: Wherever possible we try to use move semantics on Element
// since it can be object or complex types. We expect Key to be
// primitive type so we do not force move semantics on that

template <typename Key, typename Element>
inline void BinaryMinHeap<Key,Element>::init(int n)
{
    this->bHeap.reserve(n);
}

template <typename Key, typename Element>
inline void BinaryMinHeap<Key,Element>::insert(Element element, Key key)
{
    // Add element at the end of the heap and siftUp
    // until it is in the correct position
    this->bHeap.push_back(ElementKey(std::move(element),key));
    this->siftUp(bHeap.size()-1);
}

template <typename Key, typename Element>
inline void BinaryMinHeap<Key,Element>::clear()
{
    this->bHeap.clear();
}

template <typename Key, typename Element>
inline void BinaryMinHeap<Key,Element>::print()
{
    std::cout << "Not implemented yet" << std::endl;
}

template <typename Key, typename Element>
inline unsigned int BinaryMinHeap<Key,Element>::size()
{
    return this->bHeap.size();
}

template <typename Key, typename Element>
inline Key BinaryMinHeap<Key,Element>::getMinKey()
{
    return this->bHeap[0].key;
}

template <typename Key, typename Element>
inline Element BinaryMinHeap<Key,Element>::extractMinElement()
{
    Element minElement = std::move(this->bHeap[0].element);
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
    return minElement;
}

// This may be destructive depending on implementation to below
template <typename Key, typename Element>
inline void BinaryMinHeap<Key,Element>::populateFkNs(std::vector<Element>& elements, std::vector<Key>& keys)
{
    // Return heap elements in increasing order of key
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
inline void BinaryMinHeap<Key,Element>::siftUp(int heapIndex)
{
    ElementKey sifted = std::move(this->bHeap[heapIndex]);
    // If parent is the larger than the element we are sifting we move the parent down
    // Note: That if the child is smaller than the parent, it must be smaller
    // than the other child (the other child conforms to the heap property)
    while (heapIndex != 0) {
        int parentIndex = (heapIndex-1)/2;
        if (!(this->bHeap[parentIndex].key > sifted.key)) {
            break;
        }
        this->bHeap[heapIndex] = std::move(this->bHeap[parentIndex]);
        heapIndex = parentIndex;
    }
    this->bHeap[heapIndex] = std::move(sifted);
}

template <typename Key, typename Element>
inline void BinaryMinHeap<Key,Element>::siftDown(int heapIndex)
{
    int heapSize = this->bHeap.size();
    ElementKey sifted = std::move(this->bHeap[heapIndex]);
    while (true) {
        int leftChildIndex = heapIndex*2 + 1;
        int rightChildIndex = heapIndex*2 + 2;
        int smallestChildIndex;

        // Note: In practice rightChildIndex is usually smaller than heapSize
        // because we siftDown from the start of queue (it is only larger at
        // the end of the queue). So by structuring our code in this way
        // we can ensure we take advantage of any branch prediction benefits.
        if (rightChildIndex < heapSize) {
            // If the parent has two children pick the smallest (the left child on ties)
            // Note: Which child is smaller is unpredictable, so it is selected with arithmetic rather than
            // a branch, as mispredictions roughly doubled the cost of siftDown in our benchmarks
            smallestChildIndex = leftChildIndex + (this->bHeap[rightChildIndex].key < this->bHeap[leftChildIndex].key);
        } else if (leftChildIndex < heapSize) {
            // The parent has only one child (it must be smallest)
            smallestChildIndex = leftChildIndex;
        } else {
            break;
        }

        // If the smallest child is smaller than the sifted element we move the child up
        if (!(sifted.key > this->bHeap[smallestChildIndex].key)) {
            break;
        }
        this->bHeap[heapIndex] = std::move(this->bHeap[smallestChildIndex]);
        heapIndex = smallestChildIndex; // Continue sifting down from the child's old position
    }
    this->bHeap[heapIndex] = std::move(sifted);
}

#endif // _BINARYHEAP_H
