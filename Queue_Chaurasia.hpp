/*

Author: Darsh Chaurasia

Date: 2024-04-18

Description: This header file declares the Queue template class, which is utilized in the MatrixGraph class for performing breadth-first search. It provides a simple queue implementation using linked lists.

*/

// Queue_Chaurasia.hpp

#ifndef QUEUE_CHAURASIA_HPP
#define QUEUE_CHAURASIA_HPP

#include <iostream>
#include <stdexcept>

template <typename T>
class Node {
public:
    T data;
    Node* next;

    Node(T val) : data(val), next(nullptr) {}
};

template <typename T>
class Queue {
private:
    Node<T>* frontNode;
    Node<T>* rearNode;
    size_t currentSize;

public:
    // Constructor
    Queue() : frontNode(nullptr), rearNode(nullptr), currentSize(0) {}

    // Destructor
    ~Queue() {
        while (frontNode != nullptr) {
            Node<T>* temp = frontNode;
            frontNode = frontNode->next;
            delete temp;
        }
    }

    // Enqueue a new item into the queue
    void enqueue(const T& item) {
        Node<T>* newNode = new Node<T>(item);
        if (rearNode != nullptr) {
            rearNode->next = newNode;
        }
        rearNode = newNode;
        if (frontNode == nullptr) {
            frontNode = rearNode;
        }
        ++currentSize;
    }

    // Dequeue the front item from the queue
    T dequeue() {
        if (isEmpty()) {
            throw std::out_of_range("Queue is empty");
        }
        Node<T>* temp = frontNode;
        T data = frontNode->data;
        frontNode = frontNode->next;
        if (frontNode == nullptr) {
            rearNode = nullptr;
        }
        delete temp;
        --currentSize;
        return data;
    }

    // Get the front item of the queue without removing it
    T front() const {
        if (isEmpty()) {
            throw std::out_of_range("Queue is empty");
        }
        return frontNode->data;
    }

    // Check if the queue is empty
    bool isEmpty() const {
        return frontNode == nullptr;
    }

    // Get the size of the queue
    size_t size() const {
        return currentSize;
    }
};

#endif // QUEUE_CHAURASIA_HPP
