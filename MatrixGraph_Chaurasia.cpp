/*

Author: Darsh Chaurasia

Date: 2024-04-18

Description: This file implements the MatrixGraph class, providing functionalities such as adding edges, checking adjacency between vertices, getting edge weights, and performing breadth-first search (BFS) to find paths between vertices in a graph represented as an adjacency matrix.

*/


// MatrixGraph_Chaurasia.cpp

#include "MatrixGraph_Chaurasia.h"
#include "Queue_Chaurasia.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <limits>
#include <algorithm>
#include "minmaxheap_chaurasia.hpp"

#include <vector>
#include <utility>  // Include for std::pair
// Returns the number of vertices in the graph
int MatrixGraph::getVertexCount() const { return vertexCount; }

// Constructor for the matrix graph: initializes the graph with a given number of vertices and directionality
MatrixGraph::MatrixGraph(int vertices, bool directed) : vertexCount(vertices), isDirected(directed) {
    // Allocate memory for the adjacency matrix and initialize all weights to 0
    adjacencyMatrix = new float*[vertexCount];
    for (int i = 0; i < vertexCount; i++) {
        adjacencyMatrix[i] = new float[vertexCount];
        for (int j = 0; j < vertexCount; j++) {
            adjacencyMatrix[i][j] = 0.0f;
        }
    }
}

// Destructor for the matrix graph: deallocates the dynamically allocated adjacency matrix
MatrixGraph::~MatrixGraph() {
    for (int i = 0; i < vertexCount; i++) {
        delete[] adjacencyMatrix[i];
    }
    delete[] adjacencyMatrix;
}

// Adds an edge between two vertices with a specified weight, updating the matrix for undirected graphs as well
void MatrixGraph::addEdge(int start, int end, float weight) {
    if (start < 0 || start >= vertexCount || end < 0 || end >= vertexCount) {
        throw std::out_of_range("Vertex index out of range.");
    }

    adjacencyMatrix[start][end] = weight;
    if (!isDirected) {
        adjacencyMatrix[end][start] = weight;
    }
}

// Removes an edge between two vertices, setting their weight to 0, again considering directionality
void MatrixGraph::removeEdge(int start, int end) {
    if (start < 0 || start >= vertexCount || end < 0 || end >= vertexCount) {
        throw std::out_of_range("Vertex index out of range.");
    }

    adjacencyMatrix[start][end] = 0.0f;
    if (!isDirected) {
        adjacencyMatrix[end][start] = 0.0f;
    }
}

// Checks if there is an edge between two vertices
bool MatrixGraph::adjacent(int v1, int v2) const {
    if (v1 < 0 || v1 >= vertexCount || v2 < 0 || v2 >= vertexCount) {
        throw std::out_of_range("Vertex index out of range.");
    }

    return adjacencyMatrix[v1][v2] != 0.0f;
}

// Returns the weight of the edge between two vertices
float MatrixGraph::getEdgeWeight(int start, int end) const {
    if (start < 0 || start >= vertexCount || end < 0 || end >= vertexCount) {
        throw std::out_of_range("Vertex index out of range.");
    }

    if (adjacencyMatrix[start][end] == 0.0f) {
        throw std::invalid_argument("Edge does not exist.");
    }

    return adjacencyMatrix[start][end];
}

// Updates the weight of an existing edge between two vertices
void MatrixGraph::setEdgeWeight(int start, int end, float weight) {
    if (start < 0 || start >= vertexCount || end < 0 || end >= vertexCount) {
        throw std::out_of_range("Vertex index out of range.");
    }

    if (adjacencyMatrix[start][end] == 0.0f) {
        throw std::invalid_argument("Edge does not exist.");
    }

    adjacencyMatrix[start][end] = weight;
    if (!isDirected) {
        adjacencyMatrix[end][start] = weight;
    }
}

// Returns a string representation of the graph, listing edges and their weights
std::string MatrixGraph::toString() const {
    std::stringstream ss;
    for (int i = 0; i < vertexCount; i++) {
        ss << "[" << std::setw(2) << i + 1 << "]:";
        for (int j = 0; j < vertexCount; j++) {
            if (adjacencyMatrix[i][j] != 0.0f) {
                ss << "-->[" << std::setw(2) << i + 1 << "," << std::setw(2) << j + 1 << "::" << std::fixed << std::setprecision(2) << std::setw(6) << adjacencyMatrix[i][j] << "]";
            }
        }
        ss << std::endl;
    }
    return ss.str();
}

// Prints the adjacency matrix to the console, including all weights
void MatrixGraph::printRaw() const {
    std::cout << "Adjacency Matrix:\n\n";
    for (int i = 0; i < vertexCount; i++) {
        for (int j = 0; j < vertexCount; j++) {
            std::cout << std::setw(7)<< std::fixed << std::setprecision(2) << adjacencyMatrix[i][j] << "";
        }
        std::cout << std::endl;
    }
}

// Checks if there is a path from start to goal using BFS
bool MatrixGraph::pathExists(int start, int goal) {
    auto path = getBFSPath(start, goal);
    return !path.empty();
}

// Calculates the total weight of a path between two vertices
float calculatePathWeight(const MatrixGraph* graph, const std::vector<int>& path) {
    float weight = 0.0f;
    for (size_t i = 1; i < path.size(); ++i) {
        weight += graph->getEdgeWeight(path[i - 1], path[i]);
    }
    return weight;
}

// Finds and returns a path from start to goal using BFS, or an empty vector if no path exists
std::vector<int> MatrixGraph::getBFSPath(int start, int goal) {
    if (start < 0 || start >= vertexCount || goal < 0 || goal >= vertexCount) {
        throw std::out_of_range("Vertex index out of range.");
    }

    std::vector<int> path;
    std::vector<bool> visited(vertexCount, false);
    std::vector<int> prev(vertexCount, -1);
    Queue<int> queue;

    visited[start] = true;
    queue.enqueue(start);

    while (!queue.isEmpty()) {
        int current = queue.dequeue();
        if (current == goal) {
            break;
        }
        for (int i = 0; i < vertexCount; i++) {
            if (adjacencyMatrix[current][i] != 0.0f && !visited[i]) {
                visited[i] = true;
                prev[i] = current;
                queue.enqueue(i);
            }
        }
    }

    for (int at = goal; at != -1; at = prev[at]) {
        path.push_back(at);
    }
    std::reverse(path.begin(), path.end());

    if (!path.empty() && path[0] == start) {
        return path;
    }

    return std::vector<int>(); // Indicates no path was found
}



std::vector<int> MatrixGraph::getDijkstraPath(int start, int goal) {
    std::vector<float> dist(vertexCount, std::numeric_limits<float>::infinity());
    std::vector<int> prev(vertexCount, -1);
    MinHeap<std::pair<float, int>> pq; // Ensure MinHeap can handle std::pair

    dist[start] = 0;
    pq.enqueue(std::make_pair(0, start)); // Use std::make_pair for clarity

    while (!pq.isEmpty()) {
        auto current = pq.peek(); // Store the top element
        float currentDist = current.first;
        int u = current.second;
        pq.dequeue();
        if (u == goal) break;
        if (currentDist > dist[u]) continue;

        for (int v = 0; v < vertexCount; v++) {
            if (adjacencyMatrix[u][v] != 0) {
                float newDist = dist[u] + adjacencyMatrix[u][v];
                if (newDist < dist[v]) {
                    dist[v] = newDist;
                    prev[v] = u;
                    pq.enqueue(std::make_pair(newDist, v)); // Use std::make_pair
                }
            }
        }
    }

    std::vector<int> path;
    for (int at = goal; at != -1; at = prev[at]) {
        path.push_back(at);
    }
    std::reverse(path.begin(), path.end());

    return (!path.empty() && path[0] == start) ? path : std::vector<int>();
}


std::vector<std::vector<int>> MatrixGraph::getDijkstraAll(int start) {
    std::vector<float> dist(vertexCount, std::numeric_limits<float>::infinity());
    std::vector<int> prev(vertexCount, -1);
    MinHeap<std::pair<float, int>> pq;

    std::vector<std::vector<int>> paths(vertexCount);
    dist[start] = 0;
    pq.enqueue({0, start});

    while (!pq.isEmpty()) {
        auto current = pq.peek();
        float currentDist = current.first;
        int u = current.second;
        pq.dequeue();

        if (currentDist > dist[u]) continue;

        for (int v = 0; v < vertexCount; v++) {
            if (adjacencyMatrix[u][v] != 0) {
                float newDist = dist[u] + adjacencyMatrix[u][v];
                if (newDist < dist[v]) {
                    dist[v] = newDist;
                    prev[v] = u;
                    pq.enqueue({newDist, v});
                }
            }
        }
    }

    // Reconstruct paths for all vertices
    for (int v = 0; v < vertexCount; v++) {
        if (dist[v] == std::numeric_limits<float>::infinity()) {
            paths[v] = std::vector<int>(); // No path to vertex v
        } else {
            // Build path for vertex v
            std::vector<int> path;
            for (int at = v; at != -1; at = prev[at]) {
                path.push_back(at);
            }
            std::reverse(path.begin(), path.end());
            paths[v] = path;
        }
    }

    return paths;
}

