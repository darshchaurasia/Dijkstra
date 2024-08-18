# Author: Darsh Chaurasia
# Date: 2024-04-18
# Description: This makefile is designed to compile and link the MatrixGraph project, which includes
# the implementation of a graph using an adjacency matrix and a testing suite for the graph's functionalities.
# It defines rules for compiling C++ source files, linking the object files into an executable, cleaning up
# generated files, and running the compiled application. It utilizes flags for debugging and code analysis to
# ensure robust and error-free code.

# Compiler settings
CXX = g++
CXXFLAGS = -g3 -Wall -Wextra -Wno-unknown-pragmas -Wno-unused-parameter \
           -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract \
           -fsanitize=null -fsanitize=alignment -fsanitize=object-size \
           -fsanitize=pointer-overflow -fsanitize-recover=all -fsanitize-address-use-after-scope \
           -fstack-protector-all -fsanitize=return -fsanitize=vptr -fno-omit-frame-pointer \
           -fsanitize=undefined -fno-common -std=c++11

# Application name
APPNAME = dijkstra

# Source files
SOURCES = MatrixGraph_Chaurasia.cpp Chaurasia_TestGraph.cpp

# Object files
OBJS = $(SOURCES:.cpp=.o)

# Default rule for creating the application
all: $(APPNAME)

$(APPNAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(APPNAME) $(OBJS)

# Rule for creating object files from cpp files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up the build
clean:
	rm -f $(OBJS) $(APPNAME)

# Run the application
run: $(APPNAME)
	./$(APPNAME)

# Phony targets
.PHONY: all clean run
