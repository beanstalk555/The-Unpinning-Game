#pragma once

#include "multiloop.h"
#include "geometry.h"

#include <vector>

struct GridPoint
{
    int x = 0;
    int y = 0;
};

// An orthogonal grid drawing of a multiloop. All coordinates are integers in
// [0, width] x [0, height], with y increasing upwards.
struct GridLayout
{
    int width = 0;
    int height = 0;

    // edgePaths[k - 1] is the polyline for edge k, running from the crossing
    // containing half-edge +k to the crossing containing -k. Consecutive
    // points differ in exactly one coordinate.
    std::vector<std::vector<GridPoint>> edgePaths;

    // Position of each crossing, indexed like loop.sigma.Cycles().
    std::vector<GridPoint> crossingPositions;
};

// Computes a grid drawing of the multiloop with the given face (an index
// into loop.faces) as the infinite region. Follows Tamassia's bend
// minimization plus the turn-regular compaction used by SnapPy's
// orthogonal.py. Throws std::runtime_error on failure.
GridLayout ComputeGridLayout(const Multiloop& loop, int exteriorFace);

// Vertex list of a face polygon traced along layout.edgePaths in phi order.
std::vector<Vec2d> GridFacePolygon(const Multiloop& loop, const GridLayout& layout, int face);

// Shoelace area of a face polygon in grid units (validates against GridFacePolygon).
double GridFaceArea(const Multiloop& loop, const GridLayout& layout, int face);
