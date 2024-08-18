/*
Author: Darsh Chaurasia

Date: 2024-04-18

Description: This header file contains the implementations of MinHeap and MaxHeap data structures. 
The MinHeap is a generic minimum priority queue implementation that supports insertion and deletion in O(log n) time. 
Similarly, the MaxHeap is a generic maximum priority queue implementation with the same time complexity characteristics.
Both data structures require that the types they are instantiated with have all comparison operators overloaded to ensure proper functionality.
These heap structures are suitable for efficient priority queue operations where elements need to be constantly retrieved in sorted order.
*/

#ifndef MINMAXHEAP_CHAURASIA_HPP
#define MINMAXHEAP_CHAURASIA_HPP

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <functional>

// Custom comparator for MinHeap, if needed
template<typename T>
struct MinHeapComparator {
    bool operator()(const T& a, const T& b) const {
        return a < b;
    }
};

// Specialization for std::pair, compares only the first element
template<typename T1, typename T2>
struct MinHeapComparator<std::pair<T1, T2>> {
    bool operator()(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) const {
        return a.first < b.first;
    }
};

template<typename T, typename Compare = MinHeapComparator<T>>
class MinHeap {
private:
    std::vector<T> data;
    Compare comp;

    void bubbleUp(int index) {
        while (index > 0) {
            int parentIndex = (index - 1) / 2;
            if (comp(data[index], data[parentIndex])) { // Use comparator here
                std::swap(data[parentIndex], data[index]);
                index = parentIndex;
            } else {
                return;
            }
        }
    }

    void bubbleDown(int index) {
        int size = data.size();
        while (true) {
            int leftChildIndex = 2 * index + 1;
            int rightChildIndex = 2 * index + 2;
            int smallest = index;

            if (leftChildIndex < size && comp(data[leftChildIndex], data[smallest])) {
                smallest = leftChildIndex;
            }
            if (rightChildIndex < size && comp(data[rightChildIndex], data[smallest])) {
                smallest = rightChildIndex;
            }

            if (smallest == index) {
                break;
            }

            std::swap(data[index], data[smallest]);
            index = smallest;
        }
    }

public:
    MinHeap() {}

    void enqueue(const T& item) {
        data.push_back(item);
        bubbleUp(data.size() - 1);
    }

    void dequeue() {
        if (data.empty()) {
            throw std::out_of_range("Cannot dequeue from an empty heap.");
        }
        data[0] = data.back();
        data.pop_back();
        if (!data.empty()) {
            bubbleDown(0);
        }
    }

    const T& peek() const {
        if (data.empty()) {
            throw std::out_of_range("Heap is empty.");
        }
        return data.front();
    }

    int getSize() const {
        return data.size();
    }

    bool isEmpty() const {
        return data.empty();
    }
};

#endif // MINMAXHEAP_CHAURASIA_HPP
