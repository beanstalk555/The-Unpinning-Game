#pragma once

#include "permutation.h"

struct RandomMultiloop
{
    Permutation sigma;
    // A half-edge label lying on the face of the multiloop that corresponds
    // to the unbounded region of the medial graph; the face of loop.faces
    // containing this label is the natural infinite region.
    int exteriorHalfEdge = 0;
};

// Generates a random multiloop as the medial graph of a decimated Delaunay
// triangulation:
//   1. minFaces points are sampled uniformly in the unit square (so the
//      resulting multiloop has at least minFaces faces).
//   2. The Delaunay triangulation is computed.
//   3. Edges are deleted at random, keeping the graph connected: density 1
//      deletes nothing, density 0 deletes down to a spanning tree.
//   4. The medial graph of the result is returned as a vertex permutation.
// With forceOneStrand, random bigon hooks join distinct strands at crossings
// until the multiloop is a single strand. With forceSimple, random bigon
// hooks remove strand self-intersections until every strand is a simple
// closed curve. The two options are mutually exclusive.
// Throws std::runtime_error if generation fails.
RandomMultiloop GenerateRandomMultiloop(int minFaces, double density, bool forceOneStrand,
                                        bool forceSimple, unsigned seed);
