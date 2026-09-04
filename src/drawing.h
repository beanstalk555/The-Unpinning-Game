#pragma once

#include <raylib.h>

#include <string>
#include <vector>

#include "multiloop.h"
#include "permutation.h"

// Result of turning a set of user-drawn closed curves into a multiloop.
struct DrawnMultiloop
{
    Permutation sigma;
    std::vector<MultiloopComponent> components;
    // A half-edge whose multiloop face is the unbounded region of a crossing
    // component (0 when the drawing is circles only).
    int exteriorHalfEdge = 0;
    // True when the drawing is a single crossing-free loop (a simple closed
    // curve), which has no vertex permutation and is drawn directly.
    bool isUnknot = false;
    // Nonempty on failure (e.g. a non-4-valent crossing).
    std::string error;
    // Per-edge color index matching the source curve (entry 0 unused). Empty
    // when the drawing has no crossings.
    std::vector<int> strandOfEdge;
};

// Copies curve-index strand colors from a drawing onto a constructed multiloop
// so packed rendering matches the hand-drawn palette.
void ApplyDrawnStrandColors(Multiloop& loop, const DrawnMultiloop& drawn);

// Computes the vertex permutation of the multiloop traced by the given closed
// polylines, treating their transverse intersections as 4-valent crossings.
// Disconnected drawings are allowed: each connected crossing cluster and each
// crossing-free curve becomes a component, ordered left to right.
DrawnMultiloop BuildMultiloopFromCurves(const std::vector<std::vector<Vector2>>& curves);

// Proper (interior, transverse) intersection of open segments ab and cd.
// Writes the intersection point to *out when they cross.
bool SegmentsProperlyCross(Vector2 a, Vector2 b, Vector2 c, Vector2 d, Vector2* out);

// Distance from point p to segment ab.
float DistPointToSegment(Vector2 p, Vector2 a, Vector2 b);
