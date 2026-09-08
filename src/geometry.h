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

// Center of a largest inscribed circle (polylabel / pole of inaccessibility).
struct InscribedCircle
{
    Vec2d center;
    double radius = 0.0;
    bool found = false;
};

// Quadtree search maximizing distance to the boundary. `precision` is the
// allowed error in the radius (and roughly in the center). Non-positive
// precision picks a scale from the bounding box.
InscribedCircle LargestInscribedCircle(const std::vector<Vec2d>& vertices,
                                       double precision = 0.0);

// Same for a polygon with holes: the region inside `outer` and outside every
// hole. Distance is measured to the nearest ring (outer or hole).
InscribedCircle LargestInscribedCircle(const std::vector<Vec2d>& outer,
                                       const std::vector<std::vector<Vec2d>>& holes,
                                       double precision = 0.0);

// True when `point` is inside `outer` and inside none of the holes.
bool PointInPolygonWithHoles(const Vec2d& point, const std::vector<Vec2d>& outer,
                             const std::vector<std::vector<Vec2d>>& holes);

// Ear-clip triangulation. Writes vertex-index triples into `out` (counterclockwise
// in the given coordinate system). Returns false when n < 3 or the polygon is
// degenerate; may fall back to a fan if clipping stalls on a concave/self-touching ring.
bool EarClipTriangles(const std::vector<Vec2d>& vertices, std::vector<int>* out);

// True when chords (p0,p2) and (p1,p3) interleave around `center` — the
// geometric pairing of opposite slots at a 4-valent vertex bead. False when
// a neighbor coincides with the center or the two pairs occupy disjoint arcs.
bool OppositeSlotsLinked(const Vec2d& center, const Vec2d& p0, const Vec2d& p1,
                         const Vec2d& p2, const Vec2d& p3);
