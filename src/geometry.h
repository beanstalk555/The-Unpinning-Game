#pragma once

#include <vector>

// Cartesian point for planar geometry.
struct Vec2d
{
    double x = 0.0;
    double y = 0.0;
};

// Signed shoelace area of a simple polygon. Positive when vertices are listed
// counterclockwise. Requires at least three vertices.
double SignedShoelaceArea(const std::vector<Vec2d>& vertices);

// Absolute area enclosed by a simple (non-self-intersecting) polygon.
double ShoelaceArea(const std::vector<Vec2d>& vertices);

// Even-odd raycast: true when `point` lies inside the polygon. Vertices on
// the boundary may go either way; fewer than three vertices is outside.
bool PointInPolygon(const Vec2d& point, const std::vector<Vec2d>& vertices);
