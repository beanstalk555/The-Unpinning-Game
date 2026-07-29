#pragma once

#include <raylib.h>

#include <string>
#include <vector>

#include "permutation.h"

// Result of turning a set of user-drawn closed curves into a multiloop.
struct DrawnMultiloop
{
    Permutation sigma;
    // A half-edge whose multiloop face is the unbounded region of the drawing.
    int exteriorHalfEdge = 0;
    // True when the drawing is a single crossing-free loop (a simple closed
    // curve), which has no vertex permutation and is drawn directly.
    bool isUnknot = false;
    // Nonempty on failure (e.g. a disconnected diagram, or a loop with no
    // crossings drawn alongside others).
    std::string error;
};

// Computes the vertex permutation of the multiloop traced by the given closed
// polylines, treating their transverse intersections as 4-valent crossings.
// The curves are assumed to be in general position (no three segments through
// a point, crossings spread apart), which the drawing UI enforces.
DrawnMultiloop BuildMultiloopFromCurves(const std::vector<std::vector<Vector2>>& curves);

// Proper (interior, transverse) intersection of open segments ab and cd.
// Writes the intersection point to *out when they cross.
bool SegmentsProperlyCross(Vector2 a, Vector2 b, Vector2 c, Vector2 d, Vector2* out);

// Distance from point p to segment ab.
float DistPointToSegment(Vector2 p, Vector2 a, Vector2 b);
