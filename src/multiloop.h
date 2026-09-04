#pragma once

#include <string>
#include <unordered_map>
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

enum class MultiloopPart
{
    Circle,
    Crossings
};

// One connected piece of a (possibly disconnected) diagram. Circles are not
// expressible as permutation data; crossing parts list their half-edges.
struct MultiloopComponent
{
    MultiloopPart kind = MultiloopPart::Crossings;
    std::vector<int> halfEdges; // empty for Circle
    int exteriorFace = -1;      // index into Multiloop::faces; -1 for Circle
    // Color index matching a hand-drawn curve, or -1 to assign sequentially.
    int strand = -1;
};

struct SigmaNotation
{
    Permutation sigma;
    std::vector<MultiloopComponent> components;
};

struct DenseRelabel
{
    Permutation localSigma;
    std::unordered_map<int, int> globalToLocal;
    std::unordered_map<int, int> localToGlobal;
};

// Parses cycle notation, keeping empty () as Circle slots and grouping
// 4-cycles into connected crossing components in left-to-right order.
SigmaNotation ParseSigmaNotation(const std::string& s);

// Prints components left to right: () for each circle, then the current
// sigma cycles that belong to each crossing component.
std::string FormatSigmaNotation(const Permutation& sigma,
                                const std::vector<MultiloopComponent>& components);

DenseRelabel RelabelComponent(const Permutation& sigma, const std::vector<int>& halfEdges);
int FaceIndexContainingHalfEdge(const Multiloop& loop, int halfEdge);
int FirstFaceOnComponent(const Multiloop& loop, const std::vector<int>& halfEdges);
void ApplyDefaultExteriors(const Multiloop& loop, std::vector<MultiloopComponent>& components);
std::vector<int> ExteriorFaceIndices(const std::vector<MultiloopComponent>& components);
bool ComponentOwnsFace(const MultiloopComponent& component, const Multiloop& loop, int face);

// Checks that sigma describes a 4-valent planar multiloop (possibly
// disconnected): all cycles have length 4, the labels are exactly +-1..+-n,
// and each connected component has Euler characteristic 2. Empty sigma
// (circles only) is valid. Returns an empty string if valid, otherwise a
// human-readable error message.
std::string ValidateMultiloopSigma(const Permutation& sigma);
