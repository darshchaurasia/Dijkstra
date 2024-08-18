/*

Author: Darsh Chaurasia

Date: 2024-04-18

Description: This header file declares the MatrixGraph class, which is used to represent a graph using an adjacency matrix. It includes declarations for methods to manipulate the graph such as adding and removing edges, checking adjacency, and finding paths.

*/

#ifndef MATRIXGRAPH_CHAURASIA_H
#define MATRIXGRAPH_CHAURASIA_H

#include <vector>
#include <string>

// The MatrixGraph class represents a graph using an adjacency matrix for storage.
class MatrixGraph {
private:
    float** adjacencyMatrix; // Pointer to the adjacency matrix storing edge weights
    int vertexCount; // Number of vertices in the graph
    bool isDirected; // True if the graph is directed, false if undirected

public:
    // Constructor: Initializes a new instance of the MatrixGraph class
    MatrixGraph(int vertexCount, bool isDirected);
    // Destructor: Cleans up the resources used by the MatrixGraph instance
    ~MatrixGraph();
    // Returns the number of vertices in the graph
    int getVertexCount() const;
    //calculatepathweight
    static float calculatePathWeight(const MatrixGraph* graph, const std::vector<int>& path);
    // Adds an edge between two vertices with an optional weight (default is 1)
    void addEdge(int startVertex, int endVertex, float weight = 1);
    // Removes the edge between two vertices
    void removeEdge(int startVertex, int endVertex);
    // Returns true if there is an edge between two given vertices
    bool adjacent(int vertex1, int vertex2) const;
    // Gets the weight of the edge between two vertices
    float getEdgeWeight(int startVertex, int endVertex) const;
    // Sets the weight of the edge between two vertices
    void setEdgeWeight(int startVertex, int endVertex, float weight);
    // Returns a string representation of the graph
    std::string toString() const;
    // Prints the raw adjacency matrix to the console
    void printRaw() const;
    // Checks if there is a path between two vertices
    bool pathExists(int startVertex, int endVertex);
    // Finds and returns the shortest path between two vertices using BFS
    std::vector<int> getBFSPath(int startVertex, int endVertex);
    // Finds and returns the shortest path between two vertices using Dijkstra's algorithm
    std::vector<int> getDijkstraPath(int startVertex, int endVertex);
    // Finds and returns the shortest paths from a given start vertex to all other vertices
    std::vector<std::vector<int>> getDijkstraAll(int startVertex);
};

#endif // MATRIXGRAPH_CHAURASIA_H
