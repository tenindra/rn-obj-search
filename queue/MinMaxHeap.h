// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _MINMAXHEAP_H
#define _MINMAXHEAP_H

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <vector>

// Double-ended priority queue (Atkinson et al., 1986): an implicit tree where nodes on even
// depths (min levels) are no larger than all of their descendants and nodes on odd depths
// (max levels) are no smaller than all of their descendants. The minimum is the root and the
// maximum is one of the root's children, so both ends are available in O(1) and can be
// extracted in O(log n). It offers the interfaces of both BinaryMinHeap and BinaryMaxHeap.

// Can only be used for primitive types for Key or complex types
// which have defined comparison operators

// Note: Implemented entirely in this header so that it can be used with any
// Key/Element types without explicit instantiation and so calls can be inlined
// (member functions are marked inline because GCC inlines them more readily)

template <typename Key, typename Element>
class MinMaxHeap {

    public:
        void init(int n);
        void insert(Element element, Key key);
        unsigned int size();
        void clear();
        void print();
        Key getMinKey();
        Key getMaxKey();
        Element extractMinElement();
        Element extractMaxElement();
        // Populate kNNs in increasing order
        void populateKNNs(std::vector<Element>& elements, std::vector<Key>& keys);
        // Populate FkNs in decreasing order
        void populateFkNs(std::vector<Element>& elements, std::vector<Key>& keys);

    private:
        struct ElementKey {
            Element element;
            Key key;
            ElementKey(Element _element, Key _key): element(std::move(_element)), key(_key) {}
        };

        std::vector<ElementKey> mmHeap;

        static bool isOnMinLevel(std::size_t heapIndex);
        std::size_t getMaxIndex();
        void siftUp(std::size_t heapIndex);
        void siftDownMin(std::size_t heapIndex);
        void siftDownMax(std::size_t heapIndex);

};

template <typename Key, typename Element>
inline void MinMaxHeap<Key,Element>::init(int n)
{
    this->mmHeap.reserve(n);
}

template <typename Key, typename Element>
inline void MinMaxHeap<Key,Element>::insert(Element element, Key key)
{
    // Add element at the end of the heap and siftUp
    // until it is in the correct position
    this->mmHeap.push_back(ElementKey(std::move(element),key));
    this->siftUp(this->mmHeap.size()-1);
}

template <typename Key, typename Element>
inline unsigned int MinMaxHeap<Key,Element>::size()
{
    return this->mmHeap.size();
}

template <typename Key, typename Element>
inline void MinMaxHeap<Key,Element>::clear()
{
    this->mmHeap.clear();
}

template <typename Key, typename Element>
inline void MinMaxHeap<Key,Element>::print()
{
    std::cout << "Not implemented yet" << std::endl;
}

template <typename Key, typename Element>
inline Key MinMaxHeap<Key,Element>::getMinKey()
{
    return this->mmHeap[0].key;
}

template <typename Key, typename Element>
inline Key MinMaxHeap<Key,Element>::getMaxKey()
{
    return this->mmHeap[this->getMaxIndex()].key;
}

template <typename Key, typename Element>
inline Element MinMaxHeap<Key,Element>::extractMinElement()
{
    Element minElement = std::move(this->mmHeap[0].element);
    // Move the last element to the root
    std::size_t lastHeapIdx = this->mmHeap.size()-1;
    if (lastHeapIdx != 0) {
        // siftDown the new root element until it is in the correct position
        this->mmHeap[0] = std::move(this->mmHeap[lastHeapIdx]);
        this->mmHeap.pop_back();
        this->siftDownMin(0);
    } else {
        this->mmHeap.pop_back();
    }
    return minElement;
}

template <typename Key, typename Element>
inline Element MinMaxHeap<Key,Element>::extractMaxElement()
{
    std::size_t maxHeapIdx = this->getMaxIndex();
    Element maxElement = std::move(this->mmHeap[maxHeapIdx].element);
    // Move the last element to where the maximum was
    std::size_t lastHeapIdx = this->mmHeap.size()-1;
    if (lastHeapIdx != maxHeapIdx) {
        // siftDown the moved element until it is in the correct position
        // Note: The maximum is only at the root (a min level) when it is the last element
        this->mmHeap[maxHeapIdx] = std::move(this->mmHeap[lastHeapIdx]);
        this->mmHeap.pop_back();
        this->siftDownMax(maxHeapIdx);
    } else {
        this->mmHeap.pop_back();
    }
    return maxElement;
}

// These are destructive (the heap is sorted in place), as for BinaryMinHeap and BinaryMaxHeap
template <typename Key, typename Element>
inline void MinMaxHeap<Key,Element>::populateKNNs(std::vector<Element>& elements, std::vector<Key>& keys)
{
    std::sort(this->mmHeap.begin(),this->mmHeap.end(),[](const ElementKey& a, const ElementKey& b) { return a.key < b.key; });
    for (std::size_t i = 0; i < this->mmHeap.size(); ++i) {
        elements.push_back(this->mmHeap[i].element);
        keys.push_back(this->mmHeap[i].key);
    }
}

template <typename Key, typename Element>
inline void MinMaxHeap<Key,Element>::populateFkNs(std::vector<Element>& elements, std::vector<Key>& keys)
{
    std::sort(this->mmHeap.begin(),this->mmHeap.end(),[](const ElementKey& a, const ElementKey& b) { return a.key > b.key; });
    for (std::size_t i = 0; i < this->mmHeap.size(); ++i) {
        elements.push_back(this->mmHeap[i].element);
        keys.push_back(this->mmHeap[i].key);
    }
}

/*
 * Private Methods
 */

template <typename Key, typename Element>
inline bool MinMaxHeap<Key,Element>::isOnMinLevel(std::size_t heapIndex)
{
    // The depth of heapIndex is floor(log2(heapIndex+1)) = 63 - clz(heapIndex+1),
    // which is even (a min level) exactly when the leading zero count is odd
    return __builtin_clzll(static_cast<unsigned long long>(heapIndex) + 1) & 1;
}

template <typename Key, typename Element>
inline std::size_t MinMaxHeap<Key,Element>::getMaxIndex()
{
    // The maximum is the larger of the root's children (or the root if it has none)
    std::size_t heapSize = this->mmHeap.size();
    if (heapSize > 2) {
        return (this->mmHeap[1].key < this->mmHeap[2].key) ? 2 : 1;
    }
    return heapSize - 1;
}

// Note: As in BinaryMinHeap, siftUp and siftDown hold the sifted element aside and move other
// elements into the vacant position, then place the sifted element once its position is found.
// Each loop iteration moves two levels (to the grandparent or a grandchild), which halves the
// number of array positions visited compared to a binary heap of the same size.

template <typename Key, typename Element>
inline void MinMaxHeap<Key,Element>::siftUp(std::size_t heapIndex)
{
    if (heapIndex == 0) {
        return;
    }
    ElementKey sifted = std::move(this->mmHeap[heapIndex]);
    std::size_t parentIndex = (heapIndex-1)/2;
    // Note: Grandparent of heapIndex is (heapIndex-3)/4, which exists when heapIndex >= 3
    if (isOnMinLevel(heapIndex)) {
        if (this->mmHeap[parentIndex].key < sifted.key) {
            // Sifted element is larger than its max-level parent, so it belongs on the max levels
            this->mmHeap[heapIndex] = std::move(this->mmHeap[parentIndex]);
            heapIndex = parentIndex;
            while (heapIndex >= 3) {
                std::size_t grandparentIndex = (heapIndex-3)/4;
                if (!(this->mmHeap[grandparentIndex].key < sifted.key)) {
                    break;
                }
                this->mmHeap[heapIndex] = std::move(this->mmHeap[grandparentIndex]);
                heapIndex = grandparentIndex;
            }
        } else {
            while (heapIndex >= 3) {
                std::size_t grandparentIndex = (heapIndex-3)/4;
                if (!(sifted.key < this->mmHeap[grandparentIndex].key)) {
                    break;
                }
                this->mmHeap[heapIndex] = std::move(this->mmHeap[grandparentIndex]);
                heapIndex = grandparentIndex;
            }
        }
    } else {
        if (sifted.key < this->mmHeap[parentIndex].key) {
            // Sifted element is smaller than its min-level parent, so it belongs on the min levels
            this->mmHeap[heapIndex] = std::move(this->mmHeap[parentIndex]);
            heapIndex = parentIndex;
            while (heapIndex >= 3) {
                std::size_t grandparentIndex = (heapIndex-3)/4;
                if (!(sifted.key < this->mmHeap[grandparentIndex].key)) {
                    break;
                }
                this->mmHeap[heapIndex] = std::move(this->mmHeap[grandparentIndex]);
                heapIndex = grandparentIndex;
            }
        } else {
            while (heapIndex >= 3) {
                std::size_t grandparentIndex = (heapIndex-3)/4;
                if (!(this->mmHeap[grandparentIndex].key < sifted.key)) {
                    break;
                }
                this->mmHeap[heapIndex] = std::move(this->mmHeap[grandparentIndex]);
                heapIndex = grandparentIndex;
            }
        }
    }
    this->mmHeap[heapIndex] = std::move(sifted);
}

// heapIndex must be on a min level
template <typename Key, typename Element>
inline void MinMaxHeap<Key,Element>::siftDownMin(std::size_t heapIndex)
{
    std::size_t heapSize = this->mmHeap.size();
    ElementKey sifted = std::move(this->mmHeap[heapIndex]);
    while (true) {
        std::size_t firstGrandchildIndex = heapIndex*4 + 3;
        if (firstGrandchildIndex+3 < heapSize) {
            // Note: In practice this is the usual case because we siftDown from the top of the heap.
            // All four grandchildren exist, so every child (a max level) has children that are no
            // larger than it, and the smallest descendant must be one of the grandchildren.
            // Which grandchild is smallest is unpredictable, so it is selected with arithmetic rather
            // than branches (mispredicted branches here doubled the cost of extractMinElement)
            const ElementKey* grandchildren = &this->mmHeap[firstGrandchildIndex];
            std::size_t leftOffset = (grandchildren[1].key < grandchildren[0].key);
            std::size_t rightOffset = 2 + (grandchildren[3].key < grandchildren[2].key);
            std::size_t rightIsSmaller = (grandchildren[rightOffset].key < grandchildren[leftOffset].key);
            std::size_t smallestIndex = firstGrandchildIndex + leftOffset + ((rightOffset - leftOffset) & (0 - rightIsSmaller));

            if (!(this->mmHeap[smallestIndex].key < sifted.key)) {
                break;
            }
            this->mmHeap[heapIndex] = std::move(this->mmHeap[smallestIndex]);
            // The sifted element must also be no larger than the max-level parent of its new position,
            // otherwise the parent moves down instead and we continue by sifting down the parent
            std::size_t parentIndex = (smallestIndex-1)/2;
            if (this->mmHeap[parentIndex].key < sifted.key) {
                std::swap(this->mmHeap[parentIndex],sifted);
            }
            heapIndex = smallestIndex; // Continue sifting down from the grandchild's old position
        } else {
            // Near the bottom of the heap some children or grandchildren are missing,
            // so compare with every child and grandchild that exists
            std::size_t leftChildIndex = heapIndex*2 + 1;
            if (leftChildIndex >= heapSize) {
                break;
            }
            std::size_t smallestIndex = leftChildIndex;
            if (leftChildIndex+1 < heapSize && this->mmHeap[leftChildIndex+1].key < this->mmHeap[smallestIndex].key) {
                smallestIndex = leftChildIndex+1;
            }
            for (std::size_t i = firstGrandchildIndex; i < heapSize; ++i) {
                if (this->mmHeap[i].key < this->mmHeap[smallestIndex].key) {
                    smallestIndex = i;
                }
            }

            if (this->mmHeap[smallestIndex].key < sifted.key) {
                this->mmHeap[heapIndex] = std::move(this->mmHeap[smallestIndex]);
                if (smallestIndex >= firstGrandchildIndex) {
                    std::size_t parentIndex = (smallestIndex-1)/2;
                    if (this->mmHeap[parentIndex].key < sifted.key) {
                        std::swap(this->mmHeap[parentIndex],sifted);
                    }
                }
                heapIndex = smallestIndex;
            }
            // Note: Positions below this node are leaves, so the sifted element has reached its position
            break;
        }
    }
    this->mmHeap[heapIndex] = std::move(sifted);
}

// heapIndex must be on a max level
template <typename Key, typename Element>
inline void MinMaxHeap<Key,Element>::siftDownMax(std::size_t heapIndex)
{
    std::size_t heapSize = this->mmHeap.size();
    ElementKey sifted = std::move(this->mmHeap[heapIndex]);
    while (true) {
        std::size_t firstGrandchildIndex = heapIndex*4 + 3;
        if (firstGrandchildIndex+3 < heapSize) {
            // All four grandchildren exist, so the largest descendant must be one of them (selected without branches)
            const ElementKey* grandchildren = &this->mmHeap[firstGrandchildIndex];
            std::size_t leftOffset = (grandchildren[0].key < grandchildren[1].key);
            std::size_t rightOffset = 2 + (grandchildren[2].key < grandchildren[3].key);
            std::size_t rightIsLarger = (grandchildren[leftOffset].key < grandchildren[rightOffset].key);
            std::size_t largestIndex = firstGrandchildIndex + leftOffset + ((rightOffset - leftOffset) & (0 - rightIsLarger));

            if (!(sifted.key < this->mmHeap[largestIndex].key)) {
                break;
            }
            this->mmHeap[heapIndex] = std::move(this->mmHeap[largestIndex]);
            // The sifted element must also be no smaller than the min-level parent of its new position
            std::size_t parentIndex = (largestIndex-1)/2;
            if (sifted.key < this->mmHeap[parentIndex].key) {
                std::swap(this->mmHeap[parentIndex],sifted);
            }
            heapIndex = largestIndex;
        } else {
            std::size_t leftChildIndex = heapIndex*2 + 1;
            if (leftChildIndex >= heapSize) {
                break;
            }
            std::size_t largestIndex = leftChildIndex;
            if (leftChildIndex+1 < heapSize && this->mmHeap[largestIndex].key < this->mmHeap[leftChildIndex+1].key) {
                largestIndex = leftChildIndex+1;
            }
            for (std::size_t i = firstGrandchildIndex; i < heapSize; ++i) {
                if (this->mmHeap[largestIndex].key < this->mmHeap[i].key) {
                    largestIndex = i;
                }
            }

            if (sifted.key < this->mmHeap[largestIndex].key) {
                this->mmHeap[heapIndex] = std::move(this->mmHeap[largestIndex]);
                if (largestIndex >= firstGrandchildIndex) {
                    std::size_t parentIndex = (largestIndex-1)/2;
                    if (sifted.key < this->mmHeap[parentIndex].key) {
                        std::swap(this->mmHeap[parentIndex],sifted);
                    }
                }
                heapIndex = largestIndex;
            }
            break;
        }
    }
    this->mmHeap[heapIndex] = std::move(sifted);
}

#endif // _MINMAXHEAP_H
