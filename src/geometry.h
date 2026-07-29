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
