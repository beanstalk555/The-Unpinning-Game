#pragma once

#include <memory>

#include "beads.h"
#include "multiloop.h"

// Combinatorial data for a Reidemeister 3 move on a triangular face.
//
// Indices are cyclic mod 3. Side i is the half-edge h[i], running from
// crossing X_i to X_{i+1} along the face boundary (phi order). Since sigma
// cycles are counterclockwise and interior faces are traversed
// counterclockwise, the cycle at X_i is forced to be
//
//     (h[i], -h[i-1], a[i], b[i])
//
// where a[i], b[i] are the two exterior half-edges at X_i, encountered in
// that order when walking counterclockwise around the face.
struct R3FaceInfo
{
    int face = -1;
    int h[3] = {0, 0, 0};             // phi cycle of the face
    int crossing[3] = {-1, -1, -1};   // crossing index of X_i (in the old loop)
    int vertexBead[3] = {-1, -1, -1}; // necklace vertex bead of X_i
    int a[3] = {0, 0, 0};             // exterior half-edges at X_i (ccw after -h[i-1])
    int b[3] = {0, 0, 0};
    int aTerminal[3] = {-1, -1, -1};  // edge bead adjacent to X_i toward a[i]
    int bTerminal[3] = {-1, -1, -1};  // edge bead adjacent to X_i toward b[i]
    // True when the face's phi walk is counterclockwise in grid coordinates
    // (the codebase convention); false selects the mirror-image rewrite.
    bool ccw = true;
};

// Analyzes a face for an R3 move. Returns false unless the face is a
// triangle with three distinct crossings and three distinct side edges.
// When neck is non-null, also fills vertexBead / aTerminal / bTerminal
// (returning false if a terminal is missing or not an edge bead).
bool AnalyzeR3Face(const Multiloop& loop, const Necklace* neck, int face, R3FaceInfo* out);

// The vertex permutation after the move: each triangle cycle
// (h[i], -h[i-1], a[i], b[i]) becomes (h[i], -h[i-1], b[i-1], a[i]); every
// other cycle is untouched, so no half-edge relabeling is required.
// Returns false if the rewritten sigma fails ValidateMultiloopSigma.
bool ComputeR3Sigma(const Multiloop& loop, const R3FaceInfo& info, Permutation* out);

// Rewires the necklace after the sliding animation has converged, and builds
// the new multiloop. meeting[i] is the bead on side i where the pair
// (b[i], a[i+1]) fused; it becomes the new crossing X'_{i+1} with cycle
// (h[i+1], -h[i], b[i], a[i+1]). Face beads are remapped onto the new faces
// by phi-cycle overlap. Returns nullptr (leaving the necklace untouched
// except on internal failure) if anything is inconsistent.
std::unique_ptr<Multiloop> FinalizeR3(const Multiloop& oldLoop, Necklace& neck,
                                      const R3FaceInfo& info, const int meeting[3],
                                      int oldExteriorFace, int* newExteriorFace);
