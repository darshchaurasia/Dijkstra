/*

Author: Darsh Chaurasia

Date: 2024-04-18

Description: This file contains the main function to test the functionalities of the MatrixGraph class. It includes loading graphs from files, user interaction to perform graph operations, and saving the results to files.

*/

#include "MatrixGraph_Chaurasia.h"
#include <fstream>
#include <iostream>
#include <sstream>  // This includes std::ostringstream
#include <string>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include "minmaxheap_chaurasia.hpp"


void printMenu(bool promptUser = true);
MatrixGraph* loadGraphFromFile(const std::string& filepath, bool isDirected, bool isWeighted);

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " {-u|-w} <file> [-ud]" << std::endl;
        return 1;
    }

    bool isDirected = true;
    bool isWeighted = (std::string(argv[1]) == "-w");
    std::string graphFile = argv[2];
    if (argc == 4 && std::string(argv[3]) == "-ud") {
        isDirected = false;
    }

    MatrixGraph* graph = loadGraphFromFile(graphFile, isDirected, isWeighted);
    if (graph == nullptr) {
        return 1;
    }

    std::string outputFile;
    int choice;
    do {
        printMenu(); // Display the menu options
        std::cin >> choice; // Read user choice
        switch (choice) {
            case 1: {
                // Print the graph's adjacency matrix representation
                std::cout << graph->toString();
                break;
            }
            case 2: {
                // Find and print the shortest path between two vertices
                int start, end;
                std::cin >> start >> end; // Read start and end vertices
                try {
                    std::vector<int> path = graph->getBFSPath(start - 1, end - 1); // Adjust for zero-based indexing
                    if (!path.empty()) {
                        // If a path exists, print it along with the cumulative weight
                         std::cout << "BFS path from " << start << " to " << end << " is:" << std::endl;
                        std::cout << "[" <<std::setw(2) <<path[0] + 1 << ":  0.00]==>"; // Start vertex
                        float cumulativeWeight = 0.0;
                        for (size_t i = 1; i < path.size(); ++i) {
                            cumulativeWeight += graph->getEdgeWeight(path[i - 1], path[i]);
                            std::cout << "[" <<std::setw(2)<< path[i] + 1 << ":" << std::fixed << std::setprecision(2) << std::setw(6) << cumulativeWeight << "]";
                            if (i < path.size() - 1) {
                                std::cout << "==>";
                            }
                        }
                        std::cout << std::endl;
                    } else {
                        // If no path exists, notify the user
                        std::cout << "No BFS path from " << start << " to " << end << "." << std::endl;
                    }
                } catch (const std::exception& e) {
                    // Handle exceptions, such as invalid vertex indices
                    std::cerr << "Error: " << e.what() << std::endl;
                }
                break;
            }

case 3: {  // Find a Single Dijkstra Path
    int startVertex, endVertex;
    std::cin >> startVertex >> endVertex;

    // Check if the input vertices are within the graph bounds
    if (startVertex < 1 || endVertex < 1 || startVertex > graph->getVertexCount() || endVertex > graph->getVertexCount()) {
        std::cerr << "Error: Vertex index out of range. Valid vertices are from 1 to " << graph->getVertexCount() << std::endl;
        break;
    }

    try {
        std::vector<int> path = graph->getDijkstraPath(startVertex - 1, endVertex - 1); // Convert to zero-based index
        if (!path.empty()) {
            std::cout << "DIJKSTRA path from " << startVertex << " to " << endVertex << " is:\n";
         // Output the start vertex with exact formatting
std::cout << "[" << std::setw(2) << std::right << path[0] + 1 << ":" << std::setw(6) << std::fixed << std::setprecision(2) << 0.00 << "]";

float cumulativeWeight = 0.0;
for (size_t i = 1; i < path.size(); ++i) {
    cumulativeWeight += graph->getEdgeWeight(path[i - 1], path[i]);
    // Ensure the formatting for subsequent vertices matches exactly
    std::cout << "==>[" << std::setw(2) << std::right << path[i] + 1 << ":" << std::setw(6) << cumulativeWeight << "]";
}
std::cout << std::endl;

        } else {
            std::cout << "No DIJKSTRA path from " << startVertex << " to " << endVertex << "." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    break;
}

case 4: {  // Find all Dijkstra Paths from a start
    int startVertex;
    std::cin >> startVertex;

    // Validate input to ensure the vertex is within bounds
    if (startVertex < 1 || startVertex > graph->getVertexCount()) {
        std::cerr << "Error: Vertex index out of range. Valid vertices are from 1 to " << graph->getVertexCount() << std::endl;
        break;
    }

    try {
        std::vector<std::vector<int>> paths = graph->getDijkstraAll(startVertex - 1); // Convert to zero-based index
        std::ostringstream beforeStartVertex; // Buffer for paths before the starting vertex
        std::ostringstream afterStartVertex;  // Buffer for paths after the starting vertex

        for (size_t i = 0; i < paths.size(); ++i) {
            std::ostringstream currentPath;
            if (i == startVertex - 1) {
                continue; // Skip the path from the vertex to itself
            }

            if (paths[i].empty()) {
                currentPath << "Path to " << i + 1 << ": No DIJKSTRA path from " << startVertex << " to " << i + 1 << std::endl;
            } else {
                currentPath << "Path to " << i + 1 << ": [";
                currentPath << std::setw(2) << paths[i][0] + 1 << ":  0.00]"; // Start vertex with zero cost
                float cumulativeWeight = 0.0;
                for (size_t j = 1; j < paths[i].size(); ++j) {
                    cumulativeWeight += graph->getEdgeWeight(paths[i][j - 1], paths[i][j]);
                    currentPath << "==>[" << std::setw(2) << paths[i][j] + 1 << ":" << std::setw(6) << std::fixed << std::setprecision(2) << cumulativeWeight << "]";
                }
                currentPath << std::endl;
            }

            // Decide where to buffer the output based on the vertex number
            if (i < startVertex - 1) {
                beforeStartVertex << currentPath.str();
            } else {
                afterStartVertex << currentPath.str();
            }
        }

        // Output the paths in the desired order
        std::cout << beforeStartVertex.str();
        std::cout << "DIJKSTRA Paths start at Vertex " << startVertex << "\n";
        std::cout << afterStartVertex.str();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    break;
}


            case 5: {
                // Start writing to a new file
                std::cin >> outputFile; // Read the output file name
                std::ofstream outFile(outputFile);
                if (outFile.is_open()) {
                    // If file opens successfully, write the graph's data
                    // First, count the number of edges in the graph
                    int numEdges = 0;
                    for (int i = 0; i < graph->getVertexCount(); ++i) {
                        for (int j = 0; j < graph->getVertexCount(); ++j) {
                            if (graph->adjacent(i, j)) {
                                ++numEdges;
                            }
                        }
                    }

                    // Write the number of vertices and edges to the file
                    outFile << graph->getVertexCount() << " " << numEdges << std::endl;

                    // Write each edge to the file, including the weight
                    for (int i = 0; i < graph->getVertexCount(); ++i) {
                        for (int j = 0; j < graph->getVertexCount(); ++j) {
                            if (graph->adjacent(i, j)) {
                                float weight = graph->getEdgeWeight(i, j);
                                outFile << i + 1 << " " << j + 1 << " " << std::fixed << std::setprecision(6) << weight << std::endl;
                            }
                        }
                    }

                    outFile.close(); // Close the file after writing
                } else {
                    std::cerr << "Could not open file for writing." << std::endl;
                }
                break;
            }
            case 6: {
                // Append a path to the previously created file
                if (outputFile.empty()) {
                    std::cout << "No file has been created yet." << std::endl;
                } else {
                    int start, end;
                    std::cin >> start >> end; // Read start and end vertices
                    try {
                        std::vector<int> path = graph->getBFSPath(start - 1, end - 1); // Adjust for zero-based indexing
                        std::ofstream outFile(outputFile, std::ios_base::app); // Open file in append mode
                        if (outFile.is_open()) {
                            // Write the path or a message indicating no path was found
                            if (!path.empty()) {
                                outFile << "[ " << path[0] + 1 << ":  0.00]==>"; // Start vertex
                                float cumulativeWeight = 0.0;
                                for (size_t i = 1; i < path.size(); ++i) {
                                    cumulativeWeight += graph->getEdgeWeight(path[i - 1], path[i]);
                                    outFile << "[ " << path[i] + 1 << ":" << std::fixed << std::setprecision(2) << std::setw(6) << cumulativeWeight << "]";
                                    if (i < path.size() - 1) {
                                        outFile << "==>";
                                    }
                                }
                                outFile << std::endl;
                            } else {
                                outFile << "No BFS path from " << start << " to " << end <<"."<< std::endl;
                            }
                            outFile.close(); // Close the file after appending
                        } else {
                            std::cerr << "Could not open file for appending." << std::endl;
                        }
                    } catch (const std::exception& e) {
                        std::cerr << "Error: " << e.what() << std::endl;
                    }
                }
                break;
            }
            case 7: {  // Add single Dijkstra Path to file
    if (outputFile.empty()) {
        std::cout << "No file has been created yet. Please use option 3 to start a file first." << std::endl;
        break;
    }

    int startVertex, endVertex;
    std::cout << "Enter start and end vertices for the Dijkstra path: ";
    std::cin >> startVertex >> endVertex;

    // Check if the vertices are within the valid range
    if (startVertex < 1 || endVertex < 1 || startVertex > graph->getVertexCount() || endVertex > graph->getVertexCount()) {
        std::cerr << "Error: Vertex index out of range. Valid vertices are from 1 to " << graph->getVertexCount() << std::endl;
        break;
    }

    try {
        std::vector<int> path = graph->getDijkstraPath(startVertex - 1, endVertex - 1); // Adjust for zero-based indexing
        std::ofstream outFile(outputFile, std::ios::app); // Open the output file in append mode

        if (!outFile.is_open()) {
            std::cerr << "Failed to open the file for appending." << std::endl;
            break;
        }

        if (!path.empty()) {
            outFile << "DIJKSTRA Path from " << startVertex << " to " << endVertex << ":\n";
            outFile << "[" << path[0] + 1 << ": 0.00]"; // Start vertex with zero cost
            float cumulativeWeight = 0.0;
            for (size_t i = 1; i < path.size(); ++i) {
                cumulativeWeight += graph->getEdgeWeight(path[i - 1], path[i]);
                outFile << "==>[ " << path[i] + 1 << ": " << std::fixed << std::setprecision(2) << std::setw(5) << cumulativeWeight << " ]";
            }
            outFile << "\n";
        } else {
            outFile << "No DIJKSTRA path from " << startVertex << " to " << endVertex << ".\n";
        }

        outFile.close();
    } catch (const std::exception& e) {
        std::cerr << "Error while writing to file: " << e.what() << std::endl;
    }
    break;
}
case 8: {  // Add all Dijkstra Paths from a start to the file
    if (outputFile.empty()) {
        std::cout << "No file has been created yet. Please use option 5 to start a file first." << std::endl;
        break;
    }

    int startVertex;
    std::cout << "Enter the starting vertex for Dijkstra paths: ";
    std::cin >> startVertex;

    // Check if the input vertex is within the valid range
    if (startVertex < 1 || startVertex > graph->getVertexCount()) {
        std::cerr << "Error: Vertex index out of range. Valid vertices are from 1 to " << graph->getVertexCount() << std::endl;
        break;
    }

    try {
        std::vector<std::vector<int>> paths = graph->getDijkstraAll(startVertex - 1); // Adjust for zero-based indexing
        std::ofstream outFile(outputFile, std::ios::app); // Open the output file in append mode

        if (!outFile.is_open()) {
            std::cerr << "Failed to open the file for appending." << std::endl;
            break;
        }

        outFile << "DIJKSTRA Paths from vertex " << startVertex << " to all vertices:\n";
        for (size_t i = 0; i < paths.size(); ++i) {
            if (paths[i].empty()) {
                outFile << "Path to vertex " << i + 1 << ": No path available\n";
            } else {
                outFile << "Path to vertex " << i + 1 << ": [";
                outFile << paths[i][0] + 1; // Start vertex of the path
                float cumulativeWeight = 0.0;
                for (size_t j = 1; j < paths[i].size(); ++j) {
                    cumulativeWeight += graph->getEdgeWeight(paths[i][j - 1], paths[i][j]);
                    outFile << " -> " << paths[i][j] + 1 << " (" << std::fixed << std::setprecision(2) << cumulativeWeight << ")";
                }
                outFile << "]\n";
            }
        }

        outFile.close();
    } catch (const std::exception& e) {
        std::cerr << "Error while writing to file: " << e.what() << std::endl;
    }
    break;
}

            case 9999: {
                // Debug option to print the raw adjacency matrix
                graph->printRaw();
                break;
            }
            case 0: {
                // Exit the loop and end the program
                break;
            }
            default: {
                // Handle invalid menu selections
                std::cout << "Invalid choice. Please try again." << std::endl;
                break;
            }
        }
    } while (choice != 0); // Repeat until the user chooses to quit

    delete graph; // Clean up dynamically allocated graph before exiting
    return 0;
}

void printMenu(bool promptUser) {
    std::cout << "Welcome to the Graph tester!" << std::endl;
    std::cout << "1) Print the graph" << std::endl;
    std::cout << "2) Find a BFS path" << std::endl;
    std::cout << "3) Find a Single Dijkstra Path" << std::endl;
    std::cout << "4) Find all Dijkstra Paths from a start" << std::endl;
    std::cout << "5) Start a file" << std::endl;
    std::cout << "6) Add a BFS path to the file" << std::endl;
    std::cout << "7) Add single Dijkstra Path to file" << std::endl;
    std::cout << "8) Add all Dijkstra Paths from a start" << std::endl;
    std::cout << "0) Quit" << std::endl;
    if (promptUser) {
    }
}



MatrixGraph* loadGraphFromFile(const std::string& filepath, bool isDirected, bool isWeighted) {
    std::ifstream inFile(filepath);
    if (!inFile) {
        std::cerr << "Could not open the file: " << filepath << std::endl;
        return nullptr;
    }

    int vertices, edges;
    inFile >> vertices >> edges; // Read the number of vertices and edges
    MatrixGraph* graph = new MatrixGraph(vertices, isDirected);

    int u, v;
    float weight;
    for (int i = 0; i < edges; i++) {
        inFile >> u >> v;
        weight = 1.0f; // Default weight
        if (isWeighted) {
            inFile >> weight; // Only read weights if the graph is supposed to be weighted
        }
        graph->addEdge(u - 1, v - 1, weight); // Adjusting for zero-based indexing
    }

    return graph;
}

