#pragma once

#include <string>
#include <vector>

#include "permutation.h"

// A multiloop, defined entirely by its vertex permutation sigma.
// The other data is derived:
//   - numEdges n is the largest label appearing in sigma,
//   - epsilon = (-1,1)(-2,2)...(-n,n),
//   - phi = (sigma * epsilon)^-1, whose orbits are the faces,
//   - rho = sigma^2 * epsilon, whose orbits are the straight-ahead strands.
struct Multiloop
{
    explicit Multiloop(const Permutation& vertexPermutation);

    // Note: initialization happens in declaration order, and each member
    // below depends on the ones before it.
    int numEdges;
    Permutation sigma;
    Permutation epsilon;
    Permutation phi;
    Permutation rho;

    // Orbits of phi over all half-edges +-1..+-n. Unlike phi.Cycles(), this
    // includes fixed points, which are the 1-gon faces of kinks.
    std::vector<std::vector<int>> faces;

    // Straight-ahead Euler circuits. Each edge label 1..n belongs to exactly
    // one strand; strandOfEdge[k] is its strand index (entry 0 is unused).
    int numStrands;
    std::vector<int> strandOfEdge;
};

// Checks that sigma describes a connected 4-valent planar multiloop:
// all cycles have length 4, the labels are exactly +-1..+-n, the diagram is
// connected, and the Euler characteristic is 2. Returns an empty string if
// valid, otherwise a human-readable error message.
std::string ValidateMultiloopSigma(const Permutation& sigma);
