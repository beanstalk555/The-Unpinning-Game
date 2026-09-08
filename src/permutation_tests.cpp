#include "permutation.h"
#include "multiloop.h"
#include "orthogonal.h"
#include "beads.h"
#include "physics.h"
#include "physics_params.h"
#include "move_rewrite.h"
#include "r0_move.h"
#include "r1_move.h"
#include "r2_move.h"
#include "r3_move.h"
#include "random_multiloop.h"
#include "drawing.h"
#include "geometry.h"

#include <set>
#include <unordered_set>

#include <cmath>

#include <algorithm>
#include <cstdio>
#include <exception>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace
{
    int checksRun = 0;
    int checksFailed = 0;

    void Check(bool condition, const char* description)
    {
        ++checksRun;
        if (!condition) {
            ++checksFailed;
            std::printf("[permutation] FAIL: %s\n", description);
        }
    }

    // Renders the layout as ASCII art (doubled coordinates) for eyeballing.
    void PrintAsciiLayout(const GridLayout& layout)
    {
        const int width = 2 * layout.width + 1;
        const int height = 2 * layout.height + 1;
        std::vector<std::string> rows(height, std::string(width, ' '));
        for (const std::vector<GridPoint>& path : layout.edgePaths) {
            for (size_t i = 0; i + 1 < path.size(); ++i) {
                int x = 2 * path[i].x;
                int y = 2 * path[i].y;
                const int x1 = 2 * path[i + 1].x;
                const int y1 = 2 * path[i + 1].y;
                const int dx = (x1 > x) - (x1 < x);
                const int dy = (y1 > y) - (y1 < y);
                rows[y][x] = '+';
                while (x != x1 || y != y1) {
                    x += dx;
                    y += dy;
                    rows[y][x] = (x == x1 && y == y1) ? '+' : (dx != 0 ? '-' : '|');
                }
            }
        }
        for (int y = height - 1; y >= 0; --y) {
            std::printf("    %s\n", rows[y].c_str());
        }
    }
}

bool RunPermutationTests()
{
    checksRun = 0;
    checksFailed = 0;

    const auto near = [](double a, double b) { return std::fabs(a - b) < 1e-9; };

    // ---- Shoelace area (simple polygons) ----
    {
        const std::vector<Vec2d> unitSquare = {
            {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
        Check(near(ShoelaceArea(unitSquare), 1.0), "shoelace: unit square area is 1");
        Check(SignedShoelaceArea(unitSquare) > 0.0, "shoelace: unit square is counterclockwise");

        std::vector<Vec2d> cwSquare = unitSquare;
        std::reverse(cwSquare.begin(), cwSquare.end());
        Check(near(ShoelaceArea(cwSquare), 1.0), "shoelace: clockwise square has same area");
        Check(SignedShoelaceArea(cwSquare) < 0.0, "shoelace: clockwise square is negative signed area");

        const std::vector<Vec2d> triangle = {{0.0, 0.0}, {4.0, 0.0}, {0.0, 3.0}};
        Check(near(ShoelaceArea(triangle), 6.0), "shoelace: right triangle area is 6");

        // Non-convex L-shaped hexagon, area 3.
        const std::vector<Vec2d> lShape = {
            {0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}, {1.0, 2.0}, {0.0, 2.0}};
        Check(near(ShoelaceArea(lShape), 3.0), "shoelace: non-convex L-shape area is 3");

        Check(near(ShoelaceArea({}), 0.0), "shoelace: empty polygon has area 0");
        Check(near(ShoelaceArea({{0.0, 0.0}, {1.0, 0.0}}), 0.0), "shoelace: segment has area 0");

        Check(PointInPolygon({0.5, 0.5}, unitSquare), "point-in-poly: square center is inside");
        Check(!PointInPolygon({1.5, 0.5}, unitSquare), "point-in-poly: point right of square is outside");
        Check(!PointInPolygon({-0.1, 0.5}, unitSquare), "point-in-poly: point left of square is outside");
        Check(PointInPolygon({0.5, 0.5}, lShape), "point-in-poly: L-shape stem is inside");
        Check(PointInPolygon({1.5, 0.5}, lShape), "point-in-poly: L-shape foot is inside");
        Check(PointInPolygon({0.5, 1.5}, lShape), "point-in-poly: L-shape upright is inside");
        Check(!PointInPolygon({1.5, 1.5}, lShape), "point-in-poly: L-shape notch is outside");
        Check(!PointInPolygon({0.5, 0.5}, {}), "point-in-poly: empty polygon is outside");

        {
            std::vector<int> tris;
            Check(EarClipTriangles(lShape, &tris) && tris.size() >= 9,
                  "ear-clip: L-shape yields at least 3 triangles");
            double clipArea = 0.0;
            for (size_t i = 0; i + 2 < tris.size(); i += 3) {
                const Vec2d& a = lShape[static_cast<size_t>(tris[i])];
                const Vec2d& b = lShape[static_cast<size_t>(tris[i + 1])];
                const Vec2d& c = lShape[static_cast<size_t>(tris[i + 2])];
                clipArea +=
                    std::abs((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)) * 0.5;
            }
            Check(near(clipArea, 3.0), "ear-clip: L-shape triangulation area is 3");
        }

        {
            const InscribedCircle sq = LargestInscribedCircle(unitSquare, 0.01);
            Check(sq.found, "inscribed: unit square finds a circle");
            Check(PointInPolygon(sq.center, unitSquare), "inscribed: square center is inside");
            Check(std::fabs(sq.center.x - 0.5) < 0.05 && std::fabs(sq.center.y - 0.5) < 0.05,
                  "inscribed: square center is near (0.5, 0.5)");
            Check(sq.radius > 0.45 && sq.radius < 0.55, "inscribed: square radius is near 0.5");

            const std::vector<Vec2d> rightTri = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
            const InscribedCircle tri = LargestInscribedCircle(rightTri, 0.005);
            Check(tri.found, "inscribed: right triangle finds a circle");
            Check(PointInPolygon(tri.center, rightTri), "inscribed: triangle center is inside");
            const double rExpect = (2.0 - std::sqrt(2.0)) * 0.5;
            Check(std::fabs(tri.radius - rExpect) < 0.05, "inscribed: triangle radius is the inradius");
            Check(std::fabs(tri.center.x - tri.radius) < 0.05 &&
                      std::fabs(tri.center.y - tri.radius) < 0.05,
                  "inscribed: triangle center is the incenter");

            // Largest circle sits in the inner corner against the two outer
            // walls and the reflex vertex: r = sqrt(2)/(1+sqrt(2)) = 2-sqrt(2).
            const double rL = 2.0 - std::sqrt(2.0);
            const InscribedCircle ell = LargestInscribedCircle(lShape, 0.01);
            Check(ell.found, "inscribed: L-shape finds a circle");
            Check(PointInPolygon(ell.center, lShape), "inscribed: L-shape center is inside");
            Check(std::fabs(ell.radius - rL) < 0.04, "inscribed: L-shape radius is 2-sqrt(2)");
            Check(std::fabs(ell.center.x - rL) < 0.05 && std::fabs(ell.center.y - rL) < 0.05,
                  "inscribed: L-shape center is near the inner-corner incenter");
        }
    }

    // Example multiloop from the first reference figure: 9 crossings, 18 edges.
    const Permutation sigma = Permutation::FromString(
        "(6,1,-7,-2)(11,2,-12,-3)(16,3,-17,-4)(12,7,-13,-8)(17,8,-18,-9)"
        "(4,9,-5,-10)(18,13,-1,-14)(5,14,-6,-15)(10,15,-11,-16)");
    const Permutation epsilon = MakeEdgeInvolution(18);

    // Basic application and composition order.
    Check(sigma.Apply(6) == 1, "sigma maps 6 to 1");
    Check(sigma.Apply(-2) == 6, "sigma maps -2 to 6 (cycle wraps around)");
    Check(sigma.Apply(99) == 99, "sigma fixes unmapped labels");
    Check(epsilon.Apply(7) == -7, "epsilon negates labels");
    Check((sigma * epsilon).Apply(6) == -1, "left-to-right composition: (sigma*epsilon)(6) == epsilon(sigma(6))");

    // Algebraic sanity checks.
    Check(sigma * sigma.Inverse() == Permutation(), "sigma * sigma^-1 == identity");
    Check(sigma.Inverse().Inverse() == sigma, "(sigma^-1)^-1 == sigma");
    Check(epsilon * epsilon == Permutation(), "epsilon is an involution");
    Check(Permutation::FromString(sigma.ToString()) == sigma, "FromString(ToString()) round-trips");

    // Face permutation phi = (sigma * epsilon)^-1, expected cycles from the figure.
    const Permutation phi = (sigma * epsilon).Inverse();
    const Permutation expectedPhi = Permutation::FromString(
        "(-1,6,14)(-2,11,15,-6)(-3,16,-11)(-4,-10,-16)(-5,-15,10)(-7,12,2)"
        "(-8,17,3,-12)(-9,4,-17)(-13,18,8)(-14,5,9,-18)(1,13,7)");
    Check(phi == expectedPhi, "phi == face permutation from the figure");

    // Euler characteristic on the sphere: V - E + F = 9 - 18 + 11 = 2.
    Check(sigma.Cycles().size() == 9, "sigma has 9 vertex cycles");
    Check(epsilon.Cycles().size() == 18, "epsilon has 18 edge cycles");
    Check(phi.Cycles().size() == 11, "phi has 11 face cycles");

    // Multiloop derives everything from sigma alone.
    Check(ValidateMultiloopSigma(sigma).empty(), "example 1 passes validation");
    const Multiloop loop(sigma);
    Check(loop.numEdges == 18, "Multiloop derives n = 18 from sigma");
    Check(loop.epsilon == epsilon, "Multiloop derives epsilon");
    Check(loop.phi == expectedPhi, "Multiloop derives phi");
    Check(loop.faces.size() == 11, "Multiloop lists 11 faces");
    Check(loop.rho == sigma * sigma * epsilon, "Multiloop derives rho");
    Check(loop.numStrands == 1, "example 1 is a single strand");

    // ---- Reidemeister 3 rewrite ----
    {
        int triFace = -1;
        R3FaceInfo info;
        for (int f = 0; f < static_cast<int>(loop.faces.size()); ++f) {
            if (AnalyzeR3Face(loop, nullptr, f, &info)) {
                triFace = f;
                break;
            }
        }
        Check(triFace >= 0, "R3: example 1 has a triangular face with distinct crossings");
        if (triFace >= 0) {
            Permutation movedSigma;
            Check(ComputeR3Sigma(loop, info, &movedSigma), "R3: ccw rewrite is a valid multiloop");
            const Multiloop moved(movedSigma);
            Check(moved.numEdges == loop.numEdges, "R3: edge count preserved");
            Check(moved.sigma.Cycles().size() == loop.sigma.Cycles().size(),
                  "R3: crossing count preserved");
            Check(moved.faces.size() == loop.faces.size(), "R3: face count preserved");
            Check(moved.numStrands == loop.numStrands, "R3: strand count preserved");
            // Strand colors: after rematch in FinalizeR3 the Multiloop ctor's
            // IDs may differ, but ComputeR3Sigma alone doesn't rematch. Check
            // that each old strand still has a unique best-overlap partner
            // (so rematch is well-defined).
            {
                std::vector<int> bestNew(loop.numStrands, -1);
                std::vector<int> bestOv(loop.numStrands, -1);
                for (int os = 0; os < loop.numStrands; ++os) {
                    for (int ns = 0; ns < moved.numStrands; ++ns) {
                        int ov = 0;
                        for (int k = 1; k <= loop.numEdges; ++k) {
                            if (loop.strandOfEdge[k] == os && moved.strandOfEdge[k] == ns) {
                                ++ov;
                            }
                        }
                        if (ov > bestOv[os]) {
                            bestOv[os] = ov;
                            bestNew[os] = ns;
                        }
                    }
                }
                std::set<int> used;
                bool unique = true;
                for (int os = 0; os < loop.numStrands; ++os) {
                    unique = unique && bestNew[os] >= 0 && bestOv[os] > 0 &&
                             used.insert(bestNew[os]).second;
                }
                Check(unique, "R3: strand rematch by max overlap is unique");
            }

            // The central triangle survives with the same half-edge set.
            const std::set<int> want(info.h, info.h + 3);
            bool centralSurvives = false;
            for (const std::vector<int>& face : moved.faces) {
                if (face.size() == 3 && std::set<int>(face.begin(), face.end()) == want) {
                    centralSurvives = true;
                    break;
                }
            }
            Check(centralSurvives, "R3: central triangle face survives");

            // The mirrored rewrite breaks planarity (Euler characteristic),
            // so validation rejects it: the rewrite is orientation-sensitive
            // and the ccw variant is the planar one.
            R3FaceInfo mirrored = info;
            mirrored.ccw = false;
            Permutation mirroredSigma;
            Check(!ComputeR3Sigma(loop, mirrored, &mirroredSigma),
                  "R3: mirrored rewrite is rejected as non-planar");

            // A second move on the surviving triangle is also valid (R3 moves compose).
            int triFace2 = -1;
            R3FaceInfo info2;
            for (int f = 0; f < static_cast<int>(moved.faces.size()); ++f) {
                const std::vector<int>& face = moved.faces[static_cast<size_t>(f)];
                if (face.size() == 3 && std::set<int>(face.begin(), face.end()) == want &&
                    AnalyzeR3Face(moved, nullptr, f, &info2)) {
                    triFace2 = f;
                    break;
                }
            }
            Check(triFace2 >= 0, "R3: surviving triangle qualifies for another move");
            if (triFace2 >= 0) {
                Permutation twiceSigma;
                Check(ComputeR3Sigma(moved, info2, &twiceSigma),
                      "R3: second rewrite is a valid multiloop");
                const Multiloop twice(twiceSigma);
                Check(twice.faces.size() == loop.faces.size(),
                      "R3: face count preserved after two moves");
                Check(twice.numStrands == loop.numStrands,
                      "R3: strand count preserved after two moves");
            }
        }
    }

    // Second reference figure: 10 crossings, 20 edges, two strands.
    const Permutation sigma2 = Permutation::FromString(
        "(15,14,-16,-1)(6,1,-7,-2)(11,2,-12,-3)(3,18,-4,-19)(13,8,-14,-9)"
        "(4,9,-5,-10)(19,10,-20,-11)(20,5,-15,-6)(7,16,-8,-17)(12,17,-13,-18)");
    Check(ValidateMultiloopSigma(sigma2).empty(), "example 2 passes validation");
    const Multiloop loop2(sigma2);
    Check(loop2.numEdges == 20, "example 2 has 20 edges");
    Check(loop2.faces.size() == 12, "example 2 has 12 faces (Euler characteristic 2)");
    Check(loop2.numStrands == 2, "example 2 has two strands");

    // Validation rejects bad input.
    Check(ValidateMultiloopSigma(Permutation::FromString("(1,-1,2,-2)(3,-3,4,-4)")).empty(),
          "validation accepts a disconnected diagram");
    Check(!ValidateMultiloopSigma(Permutation::FromString("(1,2,-1,-2)")).empty(),
          "validation rejects a non-planar 4-cycle");
    Check(!ValidateMultiloopSigma(Permutation::FromString("(1,2,-1)")).empty(),
          "validation rejects non-4-valent sigma");

    // Grid layout smoke test with the infinite region (1,13,7) as in the figure.
    int exterior = -1;
    for (size_t f = 0; f < loop.faces.size(); ++f) {
        const std::vector<int>& face = loop.faces[f];
        if (face.size() == 3 && face[0] == 1) {
            exterior = static_cast<int>(f);
        }
    }
    Check(exterior != -1, "found the (1,13,7) face");
    try {
        const GridLayout layout = ComputeGridLayout(loop, exterior);
        Check(layout.edgePaths.size() == 18, "layout has a path for each edge");
        Check(layout.crossingPositions.size() == 9, "layout has a position for each crossing");
        bool inBounds = true;
        for (const std::vector<GridPoint>& path : layout.edgePaths) {
            for (const GridPoint& p : path) {
                inBounds = inBounds && p.x >= 0 && p.x <= layout.width && p.y >= 0 && p.y <= layout.height;
            }
        }
        Check(inBounds, "layout points are inside the bounding box");
        std::printf("[permutation] example 1 grid drawing (%d x %d), infinite face (1,13,7):\n",
                    layout.width, layout.height);
        PrintAsciiLayout(layout);

        // Necklace structure checks.
        const Necklace necklace = BuildNecklace(loop, layout, exterior);
        int vertexBeads = 0;
        bool slotsOk = true;
        bool linksOk = true;
        for (size_t i = 0; i < necklace.beads.size(); ++i) {
            const Bead& bead = necklace.beads[i];
            if (bead.isVertex) {
                ++vertexBeads;
                for (int s = 0; s < 4; ++s) {
                    const int nb = bead.neighborOfHalfEdge[s];
                    slotsOk = slotsOk && nb >= 0 && !necklace.beads[nb].isVertex &&
                              (necklace.beads[nb].neighbors[0] == static_cast<int>(i) ||
                               necklace.beads[nb].neighbors[1] == static_cast<int>(i));
                }
            }
            else {
                for (int s = 0; s < 2; ++s) {
                    const int nb = bead.neighbors[s];
                    if (nb < 0) {
                        linksOk = false;
                        continue;
                    }
                    const Bead& other = necklace.beads[nb];
                    bool reciprocal = false;
                    if (other.isVertex) {
                        for (int t = 0; t < 4; ++t) {
                            reciprocal = reciprocal || other.neighborOfHalfEdge[t] == static_cast<int>(i);
                        }
                    }
                    else {
                        reciprocal = other.neighbors[0] == static_cast<int>(i) ||
                                     other.neighbors[1] == static_cast<int>(i);
                    }
                    linksOk = linksOk && reciprocal;
                }
            }
        }
        Check(vertexBeads == 9, "necklace has 9 vertex beads");
        Check(slotsOk, "vertex bead half-edge slots point to adjacent edge beads");
        Check(linksOk, "edge bead links are reciprocal");
        Check(necklace.faceBeads.size() == 10, "necklace has 10 face beads (finite faces)");
        bool faceBeadsOk = true;
        for (size_t i = 0; i < necklace.faceBeads.size(); ++i) {
            const FaceBead& fb = necklace.faceBeads[i];
            faceBeadsOk = faceBeadsOk && fb.label == static_cast<int>(i) + 1 &&
                          fb.vertexBeads.size() == fb.phiCycle.size() &&
                          fb.face != exterior;
        }
        Check(faceBeadsOk, "face beads are labeled 1..k and track their phi cycles");
        Check(necklace.NeighborToward(necklace.crossingBead[0], loop.sigma.Cycles()[0][0]) >= 0,
              "NeighborToward finds the edge bead for a half-edge");

        {
            PhysicsParams params;
            Necklace neckT = necklace;
            Physics phys(loop, neckT, params, exterior);
            phys.RefreshTangle();
            Check(!phys.Tangled(), "example 1 grid necklace is not tangled");
            const int vb = neckT.crossingBead[0];
            Bead& vertex = neckT.beads[static_cast<size_t>(vb)];
            const int n0 = vertex.neighborOfHalfEdge[0];
            const int n1 = vertex.neighborOfHalfEdge[1];
            Check(n0 >= 0 && n1 >= 0, "example 1 first vertex has opposite-slot neighbors");
            if (n0 >= 0 && n1 >= 0) {
                const BeadPos tmp = neckT.beads[static_cast<size_t>(n0)].pos;
                neckT.beads[static_cast<size_t>(n0)].pos = neckT.beads[static_cast<size_t>(n1)].pos;
                neckT.beads[static_cast<size_t>(n1)].pos = tmp;
                phys.RefreshTangle();
                Check(phys.Tangled(), "swapping two adjacent vertex-neighbor positions is a tangle");
            }
        }

        double minArea = 1e9;
        double maxArea = 0.0;
        int countAreaOne = 0;
        int positiveFaces = 0;
        for (size_t f = 0; f < loop.faces.size(); ++f) {
            if (static_cast<int>(f) == exterior) {
                continue;
            }
            const double a = GridFaceArea(loop, layout, static_cast<int>(f));
            if (a < 0.01) {
                continue;
            }
            ++positiveFaces;
            Check(near(ShoelaceArea(GridFacePolygon(loop, layout, static_cast<int>(f))), a),
                  "grid face area matches shoelace on grid polygon");
            // Pins the embedding convention relied on by the R3 move: interior
            // faces walk counterclockwise in grid coordinates.
            Check(SignedShoelaceArea(GridFacePolygon(loop, layout, static_cast<int>(f))) > 0.0,
                  "interior face phi walk is counterclockwise in grid coordinates");
            Check(std::fabs(a * 2.0 - std::round(a * 2.0)) < 0.02,
                  "grid face area is a half-integer in grid units");
            minArea = std::min(minArea, a);
            maxArea = std::max(maxArea, a);
            if (std::fabs(a - 1.0) < 0.01) {
                ++countAreaOne;
            }
        }
        Check(positiveFaces >= 5, "example 1 has several positive-area finite faces");
        Check(maxArea > minArea + 0.5, "example 1 faces have varying grid areas");
        Check(countAreaOne >= 1, "example 1 has at least one unit-area face");
        std::printf("[permutation] example 1 finite face grid areas (min %.0f, max %.0f):\n",
                    minArea, maxArea);
        for (size_t f = 0; f < loop.faces.size(); ++f) {
            if (static_cast<int>(f) == exterior) {
                continue;
            }
            std::printf("  face %zu: %.3f\n", f, GridFaceArea(loop, layout, static_cast<int>(f)));
        }

        // ---- End-to-end R3 move: physics slide + necklace surgery ----
        {
            Necklace neck2 = necklace;
            PhysicsParams params;
            Physics phys(loop, neck2, params, exterior);
            phys.Step();
            int expectedHoneycomb = 0;
            for (const FaceBead& fb : neck2.faceBeads) {
                if (fb.phiCycle.size() <= 3 && fb.boundaryBeads.size() == 6) {
                    ++expectedHoneycomb;
                }
            }
            Check(phys.HoneycombFortifyCount() == expectedHoneycomb,
                  "honeycomb spokes match 6-bead 1-/2-/3-gons");
            int triFace = -1;
            R3FaceInfo probe;
            for (int f = 0; f < static_cast<int>(loop.faces.size()); ++f) {
                if (f != exterior && AnalyzeR3Face(loop, &neck2, f, &probe)) {
                    triFace = f;
                    break;
                }
            }
            Check(triFace >= 0, "R3 end-to-end: found a movable triangle face");
            if (triFace >= 0) {
                Check(phys.StartR3Move(triFace), "R3 end-to-end: move starts");
                int frames = 0;
                // Pre-honeycomb + fortify settle + slide need a larger frame budget.
                while (phys.R3Active() && !phys.R3Complete() && frames < 15000) {
                    phys.R3Step();
                    ++frames;
                }
                Check(phys.R3Complete(), "R3 end-to-end: slide converges");
            }
            if (phys.R3Complete()) {
                R3FaceInfo info;
                int meeting[3] = {-1, -1, -1};
                phys.GetR3Result(&info, meeting);
                for (int i = 0; i < 3; ++i) {
                    Check(meeting[i] >= 0 && meeting[i] < static_cast<int>(neck2.beads.size()),
                          "R3 end-to-end: meeting bead is valid");
                    Check(!neck2.beads[static_cast<size_t>(meeting[i])].isVertex &&
                              !neck2.beads[static_cast<size_t>(meeting[i])].removed,
                          "R3 end-to-end: meeting is the honeycomb side bead");
                }
                int newExterior = -1;
                const std::unique_ptr<Multiloop> moved =
                    FinalizeR3(loop, neck2, info, meeting, exterior, &newExterior);
                Check(moved != nullptr, "R3 end-to-end: finalize succeeds");
                if (moved) {
                    phys.EnsureR3ExteriorHoneycombMins(*moved);
                    Check(ValidateMultiloopSigma(moved->sigma).empty(),
                          "R3 end-to-end: rewritten sigma is valid");
                    Check(newExterior >= 0, "R3 end-to-end: exterior face remapped");

                    // Every vertex bead slot must link back to an edge bead
                    // that references it.
                    bool linksOk = true;
                    for (size_t i = 0; i < neck2.beads.size(); ++i) {
                        const Bead& bead = neck2.beads[i];
                        if (bead.removed) {
                            continue;
                        }
                        if (bead.isVertex) {
                            for (int s = 0; s < 4; ++s) {
                                const int nb = bead.neighborOfHalfEdge[s];
                                linksOk = linksOk && nb >= 0 && !neck2.beads[nb].isVertex &&
                                          (neck2.beads[nb].neighbors[0] == static_cast<int>(i) ||
                                           neck2.beads[nb].neighbors[1] == static_cast<int>(i));
                            }
                        }
                        else {
                            linksOk = linksOk && bead.neighbors[0] >= 0 && bead.neighbors[1] >= 0;
                        }
                    }
                    Check(linksOk, "R3 end-to-end: necklace links are consistent");

                    // Every edge chain must run vertex-to-vertex with the
                    // correct edge label on each interior bead.
                    const std::vector<std::vector<int>> movedCycles = moved->sigma.Cycles();
                    std::map<int, int> movedCrossingOf;
                    for (int c = 0; c < static_cast<int>(movedCycles.size()); ++c) {
                        for (const int h : movedCycles[c]) {
                            movedCrossingOf[h] = c;
                        }
                    }
                    bool edgeWalksOk = true;
                    for (int k = 1; k <= moved->numEdges; ++k) {
                        const int vStart = neck2.crossingBead[movedCrossingOf[k]];
                        const int vEnd = neck2.crossingBead[movedCrossingOf[-k]];
                        int prev = vStart;
                        int cur = neck2.NeighborToward(vStart, k);
                        int guard = 0;
                        while (cur >= 0 && !neck2.beads[cur].isVertex &&
                               guard < static_cast<int>(neck2.beads.size())) {
                            edgeWalksOk = edgeWalksOk && neck2.beads[cur].edge == k;
                            const Bead& b = neck2.beads[cur];
                            const int nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
                            prev = cur;
                            cur = nxt;
                            ++guard;
                        }
                        edgeWalksOk = edgeWalksOk && cur == vEnd;
                    }
                    Check(edgeWalksOk, "R3 end-to-end: all edge chains run vertex to vertex");

                    // Physics rebuilds cleanly on the moved diagram, and the
                    // post-repair phase (gap fill + sharp-bend strip) finishes.
                    bool physicsOk = true;
                    try {
                        Physics phys2(*moved, neck2, params, newExterior);
                        phys2.Step();
                        phys2.StartR3PostRepair();
                        int postFrames = 0;
                        while (phys2.R3Active() && postFrames < 4000) {
                            phys2.R3Step();
                            ++postFrames;
                        }
                        physicsOk = !phys2.R3Active();
                    }
                    catch (const std::exception&) {
                        physicsOk = false;
                    }
                    Check(physicsOk, "R3 end-to-end: physics rebuilds on the moved diagram");

                    // Strand colors stay matched to the pre-move assignment.
                    bool colorsStable = true;
                    for (int k = 1; k <= loop.numEdges; ++k) {
                        // After rematch, each edge keeps its old strand id iff
                        // it stayed in that component; otherwise it takes the
                        // id of the component it joined. Verify rematch ran:
                        // every strand id in use is in [0, numStrands).
                        const int s = moved->strandOfEdge[k];
                        if (s < 0 || s >= moved->numStrands) {
                            colorsStable = false;
                        }
                    }
                    // Max-overlap rematch: majority of each old strand's edges
                    // should still carry that strand's id.
                    for (int os = 0; os < loop.numStrands; ++os) {
                        int keep = 0, total = 0;
                        for (int k = 1; k <= loop.numEdges; ++k) {
                            if (loop.strandOfEdge[k] != os) {
                                continue;
                            }
                            ++total;
                            if (moved->strandOfEdge[k] == os) {
                                ++keep;
                            }
                        }
                        if (keep * 2 < total) {
                            colorsStable = false;
                        }
                    }
                    Check(colorsStable, "R3 end-to-end: strand colors preserved by rematch");
                }
            }
        }
    }
    catch (const std::exception& e) {
        Check(false, (std::string("layout for example 1: ") + e.what()).c_str());
    }

    try {
        const GridLayout layout2 = ComputeGridLayout(loop2, 0);
        Check(layout2.edgePaths.size() == 20, "example 2 layout has a path for each edge");
        std::printf("[permutation] example 2 grid drawing (%d x %d):\n", layout2.width, layout2.height);
        PrintAsciiLayout(layout2);
        double minArea2 = 1e9;
        double maxArea2 = 0.0;
        int positiveFaces2 = 0;
        for (size_t f = 1; f < loop2.faces.size(); ++f) {
            const double a = GridFaceArea(loop2, layout2, static_cast<int>(f));
            if (a < 0.01) {
                continue;
            }
            ++positiveFaces2;
            Check(std::fabs(a * 2.0 - std::round(a * 2.0)) < 0.02,
                  "example 2 grid face area is a half-integer in grid units");
            minArea2 = std::min(minArea2, a);
            maxArea2 = std::max(maxArea2, a);
        }
        Check(positiveFaces2 >= 5, "example 2 has several positive-area finite faces");
        Check(maxArea2 > minArea2 + 0.5, "example 2 faces have varying grid areas");
    }
    catch (const std::exception& e) {
        Check(false, (std::string("layout for example 2: ") + e.what()).c_str());
    }

    // Random multiloop generation: dense (meshlike) case.
    try {
        const RandomMultiloop random = GenerateRandomMultiloop(20, 0.8, false, false, 12345u);
        Check(ValidateMultiloopSigma(random.sigma).empty(), "random multiloop sigma is valid");
        const Multiloop randomLoop(random.sigma);
        Check(randomLoop.faces.size() >= 20, "random multiloop has at least 20 faces");
        int exterior = -1;
        for (size_t f = 0; f < randomLoop.faces.size(); ++f) {
            for (const int h : randomLoop.faces[f]) {
                if (h == random.exteriorHalfEdge) {
                    exterior = static_cast<int>(f);
                }
            }
        }
        Check(exterior != -1, "random multiloop exterior half-edge lies on a face");
        const GridLayout randomLayout = ComputeGridLayout(randomLoop, exterior);
        Check(randomLayout.edgePaths.size() == static_cast<size_t>(randomLoop.numEdges),
              "random multiloop layout succeeds");
        std::printf("[permutation] random multiloop (20 points, density 0.8): %d crossings, "
                    "%d edges, %d faces, %d strands, grid %d x %d\n",
                    static_cast<int>(randomLoop.sigma.Cycles().size()), randomLoop.numEdges,
                    static_cast<int>(randomLoop.faces.size()), randomLoop.numStrands,
                    randomLayout.width, randomLayout.height);
    }
    catch (const std::exception& e) {
        Check(false, (std::string("random multiloop (dense): ") + e.what()).c_str());
    }

    // Random multiloop generation: treelike case (exercises kinks/monogons).
    try {
        const RandomMultiloop random = GenerateRandomMultiloop(8, 0.0, false, false, 777u);
        Check(ValidateMultiloopSigma(random.sigma).empty(), "treelike random sigma is valid");
        const Multiloop treeLoop(random.sigma);
        Check(treeLoop.sigma.Cycles().size() == 7,
              "density 0 keeps exactly a spanning tree (7 crossings from 8 points)");
        Check(treeLoop.faces.size() == 9, "tree medial has V+1 = 9 faces");
        Check(treeLoop.numStrands == 1, "tree medial is a single strand");
        int exterior = -1;
        for (size_t f = 0; f < treeLoop.faces.size(); ++f) {
            for (const int h : treeLoop.faces[f]) {
                if (h == random.exteriorHalfEdge) {
                    exterior = static_cast<int>(f);
                }
            }
        }
        Check(exterior != -1, "treelike exterior face found");
        const GridLayout treeLayout = ComputeGridLayout(treeLoop, exterior);
        Check(treeLayout.edgePaths.size() == static_cast<size_t>(treeLoop.numEdges),
              "treelike layout succeeds (1-gon faces handled)");
        std::printf("[permutation] treelike random multiloop grid drawing (%d x %d):\n",
                    treeLayout.width, treeLayout.height);
        PrintAsciiLayout(treeLayout);
    }
    catch (const std::exception& e) {
        Check(false, (std::string("random multiloop (treelike): ") + e.what()).c_str());
    }

    // Flexible parsing: square brackets, space separators, outer wrapper,
    // and PD-style duplicate labels.
    {
        const Permutation reference = Permutation::FromString("(1,2,3)(4,5)");
        Check(Permutation::FromString("[1 2 3][4 5]") == reference,
              "square brackets and spaces parse like parentheses and commas");
        Check(Permutation::FromString("[ (1, 2 3)(4,5) ]") == reference,
              "outer square bracket wrapper is stripped");
        Check(Permutation::FromString("[[1,2,3],[4,5]]") == reference,
              "fully bracketed PD-style nesting parses");
        const Permutation dup = Permutation::FromString("(1,2)(2,1)");
        Check(dup.Apply(1) == 2 && dup.Apply(2) == 1, "duplicate labels: first cycle unchanged");
        Check(dup.Apply(-2) == -1 && dup.Apply(-1) == -2,
              "duplicate labels: second occurrences negated");
        Check(Permutation::FromString("[1 2 3 4]").Apply(1) == 2,
              "single bracketed cycle is not mistaken for a wrapper");
    }

    // Bigon hooks: force one strand.
    try {
        const RandomMultiloop random = GenerateRandomMultiloop(12, 0.8, true, false, 4242u);
        Check(ValidateMultiloopSigma(random.sigma).empty(), "force-one-strand sigma is valid");
        const Multiloop hooked(random.sigma);
        Check(hooked.numStrands == 1, "force one strand yields a single strand");
        int exterior = -1;
        for (size_t f = 0; f < hooked.faces.size(); ++f) {
            for (const int h : hooked.faces[f]) {
                if (h == random.exteriorHalfEdge) {
                    exterior = static_cast<int>(f);
                }
            }
        }
        Check(exterior != -1, "force-one-strand exterior face found");
        const GridLayout hookedLayout = ComputeGridLayout(hooked, exterior);
        Check(hookedLayout.edgePaths.size() == static_cast<size_t>(hooked.numEdges),
              "force-one-strand layout succeeds");
    }
    catch (const std::exception& e) {
        Check(false, (std::string("force one strand: ") + e.what()).c_str());
    }

    // Bigon hooks: force simple (no strand self-intersections).
    try {
        const RandomMultiloop random = GenerateRandomMultiloop(12, 0.8, false, true, 4242u);
        Check(ValidateMultiloopSigma(random.sigma).empty(), "force-simple sigma is valid");
        const Multiloop simple(random.sigma);
        bool noSelfIntersections = true;
        for (const std::vector<int>& cycle : simple.sigma.Cycles()) {
            if (simple.strandOfEdge[std::abs(cycle[0])] == simple.strandOfEdge[std::abs(cycle[1])]) {
                noSelfIntersections = false;
            }
        }
        Check(noSelfIntersections, "force simple removes all strand self-intersections");
        int exterior = -1;
        for (size_t f = 0; f < simple.faces.size(); ++f) {
            for (const int h : simple.faces[f]) {
                if (h == random.exteriorHalfEdge) {
                    exterior = static_cast<int>(f);
                }
            }
        }
        Check(exterior != -1, "force-simple exterior face found");
        const GridLayout simpleLayout = ComputeGridLayout(simple, exterior);
        Check(simpleLayout.edgePaths.size() == static_cast<size_t>(simple.numEdges),
              "force-simple layout succeeds");
        std::printf("[permutation] force simple (12 points): %d crossings, %d strands\n",
                    static_cast<int>(simple.sigma.Cycles().size()), simple.numStrands);
    }
    catch (const std::exception& e) {
        Check(false, (std::string("force simple: ") + e.what()).c_str());
    }

    // Free-hand drawing: two overlapping circles (a 2-component, 2-crossing
    // link shadow).
    {
        const auto circle = [](float cx, float cy, float r, int n) {
            std::vector<Vector2> pts;
            for (int i = 0; i < n; ++i) {
                const float a = 6.2831853f * static_cast<float>(i) / static_cast<float>(n);
                pts.push_back(Vector2{cx + r * std::cos(a), cy + r * std::sin(a)});
            }
            return pts;
        };
        std::vector<std::vector<Vector2>> curves = {
            circle(300.0f, 300.0f, 120.0f, 48),
            circle(440.0f, 300.0f, 120.0f, 48),
        };
        const DrawnMultiloop drawn = BuildMultiloopFromCurves(curves);
        Check(drawn.error.empty(), "two-circle drawing builds without error");
        Check(!drawn.isUnknot, "two-circle drawing is not the unknot");
        if (drawn.error.empty()) {
            const Multiloop loop(drawn.sigma);
            Check(loop.sigma.Cycles().size() == 2, "two circles give 2 crossings");
            Check(loop.numEdges == 4, "two circles give 4 edges");
            Check(loop.numStrands == 2, "two circles are 2 strands");
            bool found = false;
            for (const std::vector<int>& face : loop.faces) {
                for (const int h : face) {
                    if (h == drawn.exteriorHalfEdge) {
                        found = true;
                    }
                }
            }
            Check(found, "two-circle exterior half-edge lies on a face");
        }
    }

    // Free-hand drawing: a single self-crossing (figure eight).
    {
        std::vector<std::vector<Vector2>> curves = {{
            Vector2{200.0f, 200.0f}, Vector2{400.0f, 400.0f},
            Vector2{400.0f, 200.0f}, Vector2{200.0f, 400.0f},
        }};
        const DrawnMultiloop drawn = BuildMultiloopFromCurves(curves);
        Check(drawn.error.empty(), "figure-eight drawing builds without error");
        if (drawn.error.empty()) {
            const Multiloop loop(drawn.sigma);
            Check(loop.sigma.Cycles().size() == 1, "figure eight has 1 crossing");
            Check(loop.numEdges == 2, "figure eight has 2 edges");
            Check(loop.numStrands == 1, "figure eight is a single strand");
            Check(loop.faces.size() == 3, "figure eight has 3 faces");
        }
    }

    // Free-hand drawing: a single crossing-free loop is the unknot.
    {
        std::vector<std::vector<Vector2>> curves = {{
            Vector2{0.0f, 0.0f}, Vector2{100.0f, 0.0f},
            Vector2{100.0f, 100.0f}, Vector2{0.0f, 100.0f},
        }};
        const DrawnMultiloop drawn = BuildMultiloopFromCurves(curves);
        Check(drawn.isUnknot, "crossing-free loop is detected as the unknot");
    }

    // Disconnected sigma, packed layout, and mixed drawings.
    {
        const auto kindsOf = [](const std::string& s) {
            return ParseSigmaNotation(s).components;
        };
        const auto kindList = [](const std::vector<MultiloopComponent>& comps) {
            std::vector<MultiloopPart> kinds;
            for (const MultiloopComponent& c : comps) {
                kinds.push_back(c.kind);
            }
            return kinds;
        };
        const auto expectKinds = [&](const std::string& s, const std::vector<MultiloopPart>& want,
                                     const char* desc) {
            const std::vector<MultiloopPart> got = kindList(kindsOf(s));
            Check(got == want, desc);
        };
        expectKinds("()(1,-1,2,-2)", {MultiloopPart::Circle, MultiloopPart::Crossings},
                    "() then figure-eight is circle then crossings");
        expectKinds("(1,-1,2,-2)()", {MultiloopPart::Crossings, MultiloopPart::Circle},
                    "figure-eight then () is crossings then circle");
        expectKinds("(1,-1,2,-2)(3,-3,4,-4)",
                    {MultiloopPart::Crossings, MultiloopPart::Crossings},
                    "two figure-eights are two crossing components");
        expectKinds("(-1,3,4,-4)()(1,2,-2,-3)",
                    {MultiloopPart::Crossings, MultiloopPart::Circle},
                    "connected 2-crossing loop then circle");
        expectKinds("()()", {MultiloopPart::Circle, MultiloopPart::Circle},
                    "()() is two circles");

        const SigmaNotation mixed = ParseSigmaNotation("()(1,-1,2,-2)");
        Check(FormatSigmaNotation(mixed.sigma, mixed.components).find("()") == 0,
              "format keeps a leading empty cycle");
        const SigmaNotation trailing = ParseSigmaNotation("(1,-1,2,-2)()");
        const std::string trailingText = FormatSigmaNotation(trailing.sigma, trailing.components);
        Check(trailingText.size() >= 2 && trailingText.substr(trailingText.size() - 2) == "()",
              "format keeps a trailing empty cycle");

        Check(ValidateMultiloopSigma(mixed.sigma).empty(), "figure-eight component validates");
        Check(ValidateMultiloopSigma(ParseSigmaNotation("()()").sigma).empty(),
              "circles-only sigma validates");

        {
            const SigmaNotation two = ParseSigmaNotation("(1,-1,2,-2)(3,-3,4,-4)");
            Check(ValidateMultiloopSigma(two.sigma).empty(), "two figure-eights validate");
            const Multiloop loop(two.sigma);
            Check(loop.numEdges == 4, "two figure-eights have 4 edges");
            Check(loop.sigma.Cycles().size() == 2, "two figure-eights have 2 crossings");
            std::vector<MultiloopComponent> comps = two.components;
            ApplyDefaultExteriors(loop, comps);
            try {
                const PackedGridLayout packed = ComputePackedGridLayout(loop, comps);
                Check(packed.layout.edgePaths.size() == 4, "packed layout has 4 edge paths");
                Check(packed.layout.crossingPositions.size() == 2,
                      "packed layout has 2 crossing positions");
                Check(packed.componentOrigin.size() == 2, "packed layout has two component slots");
                Check(near(packed.componentOrigin[0].x, 0.0),
                      "first figure-eight pack origin is 0");
                Check(near(packed.componentOrigin[1].x, 2.5),
                      "second figure-eight pack origin is width plus one bead diameter");
                const Necklace necklace = BuildPackedNecklace(loop, packed, comps);
                Check(necklace.crossingBead.size() == 2, "packed necklace has 2 vertex beads");
                int circleBeads = 0;
                for (const Bead& b : necklace.beads) {
                    if (!b.isVertex && b.edge == 0) {
                        ++circleBeads;
                    }
                }
                Check(circleBeads == 0, "two figure-eights have no circle beads");

                const int firstExterior = comps[0].exteriorFace;
                GridPoint firstCrossing = packed.layout.crossingPositions[0];
                int otherFace = -1;
                for (int f = 0; f < static_cast<int>(loop.faces.size()); ++f) {
                    if (ComponentOwnsFace(comps[1], loop, f) && f != comps[1].exteriorFace) {
                        otherFace = f;
                        break;
                    }
                }
                Check(otherFace >= 0, "second figure-eight has another face to select");
                if (otherFace >= 0) {
                    comps[1].exteriorFace = otherFace;
                    const PackedGridLayout packed2 = ComputePackedGridLayout(loop, comps);
                    Check(packed2.layout.crossingPositions[0].x == firstCrossing.x &&
                              packed2.layout.crossingPositions[0].y == firstCrossing.y,
                          "changing the other component exterior does not move the first drawing");
                    Check(comps[0].exteriorFace == firstExterior,
                          "first component exterior is unchanged");
                }
            }
            catch (const std::exception& e) {
                Check(false, (std::string("packed two figure-eights: ") + e.what()).c_str());
            }
        }

        {
            const SigmaNotation circled = ParseSigmaNotation("()(1,-1,2,-2)");
            const Multiloop loop(circled.sigma);
            std::vector<MultiloopComponent> comps = circled.components;
            ApplyDefaultExteriors(loop, comps);
            try {
                const PackedGridLayout packed = ComputePackedGridLayout(loop, comps);
                Check(packed.layout.edgePaths.size() == 2, "circle+figure-eight has 2 edge paths");
                Check(packed.layout.crossingPositions.size() == 1,
                      "circle+figure-eight has 1 crossing");
                const Necklace necklace = BuildPackedNecklace(loop, packed, comps);
                int circleBeads = 0;
                int circleFaces = 0;
                for (const Bead& b : necklace.beads) {
                    if (!b.isVertex && b.edge == 0) {
                        ++circleBeads;
                    }
                }
                for (const FaceBead& fb : necklace.faceBeads) {
                    if (fb.face < 0) {
                        ++circleFaces;
                    }
                }
                Check(circleBeads == 8, "circle component contributes 8 edge beads");
                Check(circleFaces == 1, "circle component contributes one face bead");
                {
                    Necklace neck = necklace;
                    PhysicsParams params;
                    std::unordered_set<int> exteriors;
                    for (const int f : ExteriorFaceIndices(comps)) {
                        exteriors.insert(f);
                    }
                    Physics phys(loop, neck, params, exteriors);
                    phys.Step();
                    bool sawCirclePoly = false;
                    for (const FaceBead& fb : neck.faceBeads) {
                        if (fb.face >= 0) {
                            continue;
                        }
                        std::vector<Vec2d> poly;
                        Check(phys.FaceBeadBoundaryPolygon(fb, &poly),
                              "mixed packed circle pad has a usable boundary polygon");
                        Check(poly.size() >= 3,
                              "mixed packed circle pad polygon has at least 3 vertices");
                        sawCirclePoly = true;
                    }
                    Check(sawCirclePoly, "mixed packed diagram has a circle-pad polygon");
                }
            }
            catch (const std::exception& e) {
                Check(false, (std::string("packed circle+figure-eight: ") + e.what()).c_str());
            }
        }

        {
            std::vector<std::vector<Vector2>> squares = {
                {{0.0f, 0.0f}, {80.0f, 0.0f}, {80.0f, 80.0f}, {0.0f, 80.0f}},
                {{200.0f, 0.0f}, {280.0f, 0.0f}, {280.0f, 80.0f}, {200.0f, 80.0f}},
            };
            const DrawnMultiloop drawn = BuildMultiloopFromCurves(squares);
            Check(drawn.error.empty(), "two disjoint squares build without error");
            Check(!drawn.isUnknot, "two disjoint squares are not a single unknot");
            Check(drawn.components.size() == 2, "two disjoint squares are two components");
            if (drawn.components.size() == 2) {
                Check(drawn.components[0].kind == MultiloopPart::Circle &&
                          drawn.components[1].kind == MultiloopPart::Circle,
                      "two disjoint squares are two circles");
                Check(drawn.components[0].strand == 0 && drawn.components[1].strand == 1,
                      "two disjoint squares keep left-to-right drawing colors");
            }
        }

        {
            std::vector<std::vector<Vector2>> mixedDraw = {
                {{0.0f, 0.0f}, {80.0f, 0.0f}, {80.0f, 80.0f}, {0.0f, 80.0f}},
                {{200.0f, 200.0f}, {400.0f, 400.0f}, {400.0f, 200.0f}, {200.0f, 400.0f}},
            };
            const DrawnMultiloop drawn = BuildMultiloopFromCurves(mixedDraw);
            Check(drawn.error.empty(), "square + figure-eight drawing builds without error");
            Check(!drawn.isUnknot, "square + figure-eight is not the unknot");
            bool sawCircle = false;
            bool sawCrossings = false;
            for (const MultiloopComponent& c : drawn.components) {
                if (c.kind == MultiloopPart::Circle) {
                    sawCircle = true;
                }
                if (c.kind == MultiloopPart::Crossings) {
                    sawCrossings = true;
                }
            }
            Check(sawCircle && sawCrossings, "square + figure-eight is mixed circle and crossings");
            if (drawn.error.empty() && sawCrossings) {
                const Multiloop loop(drawn.sigma);
                Check(loop.sigma.Cycles().size() == 1, "mixed drawing has 1 crossing");
                Check(loop.numEdges == 2, "mixed drawing has 2 edges");
            }
        }

        {
            const SigmaNotation twoCircles = ParseSigmaNotation("()()");
            const Permutation empty;
            const Multiloop dummy(empty);
            std::vector<MultiloopComponent> comps = twoCircles.components;
            const PackedGridLayout packed = ComputePackedGridLayout(dummy, comps);
            const Necklace necklace = BuildPackedNecklace(dummy, packed, comps);
            std::set<int> strands;
            for (const Bead& b : necklace.beads) {
                if (!b.isVertex && b.edge == 0) {
                    strands.insert(b.strand);
                }
            }
            Check(strands.count(0) == 1 && strands.count(1) == 1 && strands.size() == 2,
                  "two packed circles use distinct strand colors");
            Check(near(packed.componentOrigin[0].x, 0.0), "first packed circle origin is 0");
            Check(near(packed.componentOrigin[1].x, 1.5),
                  "second packed circle origin is 1 plus one bead diameter");
            float minCross = 1e9f;
            for (const Bead& a : necklace.beads) {
                if (a.isVertex || a.edge != 0 || a.strand != 0) {
                    continue;
                }
                for (const Bead& b : necklace.beads) {
                    if (b.isVertex || b.edge != 0 || b.strand != 1) {
                        continue;
                    }
                    const float dx = a.pos.x - b.pos.x;
                    const float dy = a.pos.y - b.pos.y;
                    minCross = std::min(minCross, std::sqrt(dx * dx + dy * dy));
                }
            }
            Check(std::fabs(minCross - 0.5f) < 1e-4f,
                  "nearest beads of two packed circles are one diameter apart");
            {
                Necklace neck = necklace;
                PhysicsParams params;
                Physics phys(dummy, neck, params, std::unordered_set<int>{});
                phys.Step();
                int circlePolys = 0;
                for (const FaceBead& fb : neck.faceBeads) {
                    if (fb.face >= 0) {
                        continue;
                    }
                    std::vector<Vec2d> poly;
                    Check(phys.FaceBeadBoundaryPolygon(fb, &poly),
                          "packed circle pad has a usable boundary polygon");
                    Check(poly.size() >= 3, "packed circle pad polygon has at least 3 vertices");
                    ++circlePolys;
                }
                Check(circlePolys == 2, "two packed circles each expose a boundary polygon");

                const auto minRingDist = [](const Necklace& n) {
                    float best = 1e9f;
                    for (const Bead& a : n.beads) {
                        if (a.removed || a.isVertex || a.edge != 0 || a.strand != 0) {
                            continue;
                        }
                        for (const Bead& b : n.beads) {
                            if (b.removed || b.isVertex || b.edge != 0 || b.strand != 1) {
                                continue;
                            }
                            const float dx = a.pos.x - b.pos.x;
                            const float dy = a.pos.y - b.pos.y;
                            best = std::min(best, std::sqrt(dx * dx + dy * dy));
                        }
                    }
                    return best;
                };
                for (Bead& b : neck.beads) {
                    if (!b.removed && !b.isVertex && b.edge == 0 && b.strand == 1) {
                        b.pos.x -= 0.35f;
                    }
                }
                for (FaceBead& fb : neck.faceBeads) {
                    if (fb.face < 0 && fb.pos.x > 0.75f) {
                        fb.pos.x -= 0.35f;
                    }
                }
                Check(minRingDist(neck) < 0.4f, "shoved packed circles start overlapping");
                for (int i = 0; i < 80; ++i) {
                    phys.Step();
                }
                Check(minRingDist(neck) + 1e-3f >= 2.0f * neck.radius,
                      "exterior collisions keep packed circle rings a diameter apart");
            }
        }

        {
            Necklace neck = BuildUnknotNecklace();
            PhysicsParams params;
            Physics phys(neck, params);
            Check(static_cast<int>(neck.beads.size()) == 8, "unknot starts with 8 edge beads");
            phys.StartAutoTighten();
            int frames = 0;
            while (phys.AutoTightening() && frames < 4000) {
                phys.AutoStep();
                ++frames;
            }
            Check(!phys.AutoTightening(), "unknot auto-tighten reaches idle");
            int live = 0;
            for (const Bead& b : neck.beads) {
                if (!b.removed && !b.isVertex) {
                    ++live;
                }
            }
            Check(live == 6, "unknot auto-tighten leaves 6 edge beads");
            Check(phys.HoneycombFortifyCount() >= 1, "unknot ring is a honeycomb pad");
        }

        {
            const SigmaNotation twoCirclesTighten = ParseSigmaNotation("()()");
            const Permutation emptyTighten;
            const Multiloop dummyTighten(emptyTighten);
            std::vector<MultiloopComponent> compsTighten = twoCirclesTighten.components;
            const PackedGridLayout packedTighten = ComputePackedGridLayout(dummyTighten, compsTighten);
            Necklace neck = BuildPackedNecklace(dummyTighten, packedTighten, compsTighten);
            PhysicsParams params;
            Physics phys(dummyTighten, neck, params, std::unordered_set<int>{});
            phys.StartAutoTighten();
            int frames = 0;
            while (phys.AutoTightening() && frames < 4000) {
                phys.AutoStep();
                ++frames;
            }
            Check(!phys.AutoTightening(), "two packed circles auto-tighten reach idle");
            int circlePads = 0;
            for (const FaceBead& fb : neck.faceBeads) {
                if (fb.face >= 0) {
                    continue;
                }
                int live = 0;
                for (const int bead : fb.boundaryBeads) {
                    if (bead >= 0 && bead < static_cast<int>(neck.beads.size()) &&
                        !neck.beads[static_cast<size_t>(bead)].removed) {
                        ++live;
                    }
                }
                Check(live == 6, "packed circle auto-tighten leaves 6 live boundary beads");
                ++circlePads;
            }
            Check(circlePads == 2, "two packed circles each have a face bead");
            Check(phys.HoneycombFortifyCount() >= 2,
                  "packed circle pads are included in honeycomb fortify count");
        }

        {
            const SigmaNotation oneCircle = ParseSigmaNotation("()");
            const Permutation emptyDrag;
            const Multiloop dummyDrag(emptyDrag);
            const PackedGridLayout packedDrag =
                ComputePackedGridLayout(dummyDrag, oneCircle.components);
            Necklace neck = BuildPackedNecklace(dummyDrag, packedDrag, oneCircle.components);
            PhysicsParams params;
            Physics phys(dummyDrag, neck, params, std::unordered_set<int>{});
            phys.Step();
            int bead = -1;
            for (int i = 0; i < static_cast<int>(neck.beads.size()); ++i) {
                if (!neck.beads[static_cast<size_t>(i)].removed &&
                    !neck.beads[static_cast<size_t>(i)].isVertex &&
                    neck.beads[static_cast<size_t>(i)].edge == 0) {
                    bead = i;
                    break;
                }
            }
            Check(bead >= 0, "packed circle has an edge bead to drag");
            if (bead >= 0) {
                const int strand = neck.beads[static_cast<size_t>(bead)].strand;
                const Vector2 p = {neck.beads[static_cast<size_t>(bead)].pos.x,
                                   neck.beads[static_cast<size_t>(bead)].pos.y};
                phys.BeginDrag(bead, p);
                neck.beads[static_cast<size_t>(bead)].pos.x += 3.0f;
                phys.Drag({p.x + 3.0f, p.y});
                for (int i = 0; i < 8; ++i) {
                    phys.Step();
                }
                phys.EndDrag();
                int live = 0;
                int wrongStrand = 0;
                for (const Bead& b : neck.beads) {
                    if (b.removed || b.isVertex || b.edge != 0) {
                        continue;
                    }
                    ++live;
                    if (b.strand != strand) {
                        ++wrongStrand;
                    }
                }
                Check(live > 8, "drag-insert grows a packed circle past 8 beads");
                Check(wrongStrand == 0, "drag-inserted circle beads keep strand color");
            }
        }

        {
            std::vector<std::vector<Vector2>> mixedDraw = {
                {{0.0f, 0.0f}, {80.0f, 0.0f}, {40.0f, 70.0f}},
                {{200.0f, 200.0f}, {400.0f, 400.0f}, {400.0f, 200.0f}, {200.0f, 400.0f}},
            };
            const DrawnMultiloop drawn = BuildMultiloopFromCurves(mixedDraw);
            Check(drawn.error.empty(), "triangle + figure-eight drawing builds without error");
            if (drawn.error.empty()) {
                Check(drawn.components.size() == 2, "triangle + figure-eight has two components");
                if (drawn.components.size() == 2) {
                    Check(drawn.components[0].kind == MultiloopPart::Circle &&
                              drawn.components[0].strand == 0,
                          "left triangle is circle strand 0");
                    Check(drawn.components[1].kind == MultiloopPart::Crossings,
                          "right figure-eight is a crossing component");
                }
                try {
                    Multiloop loop(drawn.sigma);
                    ApplyDrawnStrandColors(loop, drawn);
                    Check(loop.numStrands >= 2, "drawn mixed diagram keeps both curve colors");
                    bool eightIsStrand1 = loop.numEdges >= 1;
                    for (int k = 1; k <= loop.numEdges; ++k) {
                        eightIsStrand1 = eightIsStrand1 && loop.strandOfEdge[static_cast<size_t>(k)] == 1;
                    }
                    Check(eightIsStrand1, "figure-eight edges keep drawing color 1");
                    std::vector<MultiloopComponent> comps = drawn.components;
                    ApplyDefaultExteriors(loop, comps);
                    const PackedGridLayout packed = ComputePackedGridLayout(loop, comps);
                    Check(packed.componentOrigin.size() >= 2 &&
                              packed.componentOrigin[0].x < packed.componentOrigin[1].x,
                          "triangle packs to the left of the figure-eight");
                    Necklace neck = BuildPackedNecklace(loop, packed, comps);
                    int circleStrand = -1;
                    for (const Bead& b : neck.beads) {
                        if (!b.isVertex && b.edge == 0) {
                            circleStrand = b.strand;
                            break;
                        }
                    }
                    Check(circleStrand == 0, "packed triangle circle uses drawing color 0");
                    PhysicsParams params;
                    std::unordered_set<int> exteriors;
                    for (const int f : ExteriorFaceIndices(comps)) {
                        exteriors.insert(f);
                    }
                    Physics phys(loop, neck, params, exteriors);
                    phys.Step();
                    phys.RefreshTangle();
                    Check(true, "triangle + figure-eight packed physics runs");
                }
                catch (const std::exception& e) {
                    Check(false, (std::string("triangle + figure-eight pipeline: ") + e.what()).c_str());
                }
            }
        }

        {
            const auto polar = [](float cx, float cy, float r, float a) {
                return Vector2{cx + r * std::cos(a), cy + r * std::sin(a)};
            };
            std::vector<Vector2> pentagon;
            for (int i = 0; i < 5; ++i) {
                const float a = -1.5707963f + static_cast<float>(i) * (6.2831853f / 5.0f);
                pentagon.push_back(polar(0.0f, 0.0f, 140.0f, a));
            }
            const std::vector<Vector2> star = {pentagon[0], pentagon[2], pentagon[4], pentagon[1],
                                               pentagon[3]};
            const auto triangle = [&](float cx, float cy, float r, float rot) {
                std::vector<Vector2> pts;
                for (int i = 0; i < 3; ++i) {
                    pts.push_back(polar(cx, cy, r, rot + static_cast<float>(i) * 2.0943951f));
                }
                return pts;
            };
            std::vector<std::vector<Vector2>> curves = {
                star,
                triangle(420.0f, 0.0f, 100.0f, -1.5707963f),
                triangle(420.0f, 0.0f, 100.0f, 1.5707963f),
            };
            const DrawnMultiloop drawn = BuildMultiloopFromCurves(curves);
            Check(drawn.error.empty(), "pentagram + star of david builds without error");
            if (!drawn.error.empty()) {
                Check(false, drawn.error.c_str());
            }
            else {
                try {
                    const Multiloop loop(drawn.sigma);
                    std::vector<MultiloopComponent> comps = drawn.components;
                    ApplyDefaultExteriors(loop, comps);
                    Check(comps.size() >= 2, "pentagram + star of david has two components");
                    const PackedGridLayout packed = ComputePackedGridLayout(loop, comps);
                    Necklace neck = BuildPackedNecklace(loop, packed, comps);
                    PhysicsParams params;
                    std::unordered_set<int> exteriors;
                    for (const int f : ExteriorFaceIndices(comps)) {
                        exteriors.insert(f);
                    }
                    Physics phys(loop, neck, params, exteriors);
                    phys.Step();
                    phys.RefreshTangle();
                    Check(packed.layout.crossingPositions.size() == loop.sigma.Cycles().size(),
                          "packed pentagram + star of david has a position per crossing");
                }
                catch (const std::exception& e) {
                    Check(false,
                          (std::string("pentagram + star of david pipeline: ") + e.what()).c_str());
                }
            }
        }
    }

    // ---- Shortening Reidemeister moves (R1 / R2): rewrite, physics, nesting ----
    {
        struct Scene
        {
            std::unique_ptr<Multiloop> loop;
            std::vector<MultiloopComponent> comps;
            Necklace neck;
            std::unordered_set<int> exteriors;
        };
        // Builds a packed scene; `exterior` overrides the crossing component's
        // exterior face when >= 0.
        const auto buildScene = [&](const char* text, int exterior) {
            std::unique_ptr<Scene> scene(new Scene());
            const SigmaNotation parsed = ParseSigmaNotation(text);
            scene->loop.reset(new Multiloop(parsed.sigma));
            scene->comps = parsed.components;
            ApplyDefaultExteriors(*scene->loop, scene->comps);
            if (exterior >= 0) {
                for (MultiloopComponent& c : scene->comps) {
                    if (c.kind == MultiloopPart::Crossings) {
                        c.exteriorFace = exterior;
                    }
                }
            }
            const PackedGridLayout packed = ComputePackedGridLayout(*scene->loop, scene->comps);
            scene->neck = BuildPackedNecklace(*scene->loop, packed, scene->comps);
            for (const int f : ExteriorFaceIndices(scene->comps)) {
                scene->exteriors.insert(f);
            }
            return scene;
        };
        // Every live edge bead has two neighbors that reference it back; every
        // vertex slot links to an edge bead that references the vertex.
        const auto linksOk = [](const Necklace& neck) {
            bool ok = true;
            for (size_t i = 0; i < neck.beads.size(); ++i) {
                const Bead& bead = neck.beads[i];
                if (bead.removed) {
                    continue;
                }
                if (bead.isVertex) {
                    for (int s = 0; s < 4; ++s) {
                        const int nb = bead.neighborOfHalfEdge[s];
                        ok = ok && nb >= 0 && nb < static_cast<int>(neck.beads.size()) &&
                             !neck.beads[static_cast<size_t>(nb)].isVertex &&
                             (neck.beads[static_cast<size_t>(nb)].neighbors[0] == static_cast<int>(i) ||
                              neck.beads[static_cast<size_t>(nb)].neighbors[1] == static_cast<int>(i));
                    }
                }
                else {
                    for (int s = 0; s < 2; ++s) {
                        const int nb = bead.neighbors[s];
                        if (nb < 0 || nb >= static_cast<int>(neck.beads.size()) ||
                            neck.beads[static_cast<size_t>(nb)].removed) {
                            ok = false;
                            continue;
                        }
                        const Bead& other = neck.beads[static_cast<size_t>(nb)];
                        bool back = false;
                        if (other.isVertex) {
                            for (int t = 0; t < 4; ++t) {
                                back = back || other.neighborOfHalfEdge[t] == static_cast<int>(i);
                            }
                        }
                        else {
                            back = other.neighbors[0] == static_cast<int>(i) ||
                                   other.neighbors[1] == static_cast<int>(i);
                        }
                        ok = ok && back;
                    }
                }
            }
            return ok;
        };
        // Every edge of `loop` walks vertex to vertex with its label on each
        // interior bead; every edge-0 bead lies on a closed ring.
        const auto walksOk = [](const Multiloop& loop, const Necklace& neck) {
            bool ok = true;
            const std::vector<std::vector<int>> cycles = loop.sigma.Cycles();
            std::map<int, int> crossingOf;
            for (int c = 0; c < static_cast<int>(cycles.size()); ++c) {
                for (const int h : cycles[static_cast<size_t>(c)]) {
                    crossingOf[h] = c;
                }
            }
            for (int k = 1; k <= loop.numEdges; ++k) {
                if (crossingOf.count(k) == 0 || crossingOf.count(-k) == 0) {
                    ok = false;
                    continue;
                }
                const int vStart = neck.crossingBead[static_cast<size_t>(crossingOf[k])];
                const int vEnd = neck.crossingBead[static_cast<size_t>(crossingOf[-k])];
                int prev = vStart;
                int cur = neck.NeighborToward(vStart, k);
                int guard = 0;
                while (cur >= 0 && !neck.beads[static_cast<size_t>(cur)].isVertex &&
                       guard < static_cast<int>(neck.beads.size())) {
                    ok = ok && neck.beads[static_cast<size_t>(cur)].edge == k;
                    const Bead& b = neck.beads[static_cast<size_t>(cur)];
                    const int nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
                    prev = cur;
                    cur = nxt;
                    ++guard;
                }
                ok = ok && cur == vEnd;
            }
            for (size_t i = 0; i < neck.beads.size(); ++i) {
                const Bead& b = neck.beads[i];
                if (b.isVertex || b.removed || b.edge != 0) {
                    continue;
                }
                ok = ok && CircleRingThrough(neck, static_cast<int>(i)).size() >= 3;
            }
            return ok;
        };
        const auto runMove = [](Physics& phys) {
            int frames = 0;
            while (phys.MoveActive() && !phys.MoveComplete() && frames < 20000) {
                phys.MoveStep();
                ++frames;
            }
            return phys.MoveComplete();
        };
        const auto runPost = [](Physics& phys) {
            int frames = 0;
            while ((phys.MoveActive() || phys.FaceMergePending()) && frames < 6000) {
                phys.MoveStep();
                ++frames;
            }
            return !phys.MoveActive() && !phys.FaceMergePending();
        };
        const auto finiteFaces = [](const Multiloop& loop, const std::unordered_set<int>& ext) {
            int n = 0;
            for (int f = 0; f < static_cast<int>(loop.faces.size()); ++f) {
                if (ext.count(f) == 0) {
                    ++n;
                }
            }
            return n;
        };

        // -- Rewrite helper --
        {
            const Multiloop kink(Permutation::FromString("(1,-1,2,-2)"));
            RewriteResult rw;
            Check(DeleteCrossingsAndSplice(kink, {0}, {1}, {{2, -2}}, &rw),
                  "rewrite: deleting the only crossing of a kink succeeds");
            Check(rw.numNewEdges == 0 && rw.sigma.Cycles().empty(),
                  "rewrite: kink without its crossing has no edges");
            Check(rw.circleArcs.size() == 1 && rw.circleStrand.size() == 1 &&
                      rw.circleStrand[0] == 0,
                  "rewrite: the surviving arc closes into one circle on strand 0");

            const Multiloop twoKinks(Permutation::FromString("(1,-1,-4,2)(3,-3,-2,4)"));
            R1FaceInfo r1;
            Check(AnalyzeR1Face(twoKinks, nullptr, 1, &r1) && r1.e == 1,
                  "rewrite: face (1) of the two-kink curve is a monogon on edge 1");
            RewriteResult rw2;
            Check(DeleteCrossingsAndSplice(twoKinks, {r1.crossing}, {1}, {{r1.p, r1.q}}, &rw2),
                  "rewrite: deleting one kink of two succeeds");
            Check(rw2.numNewEdges == 2 && rw2.sigma.Cycles().size() == 1 &&
                      ValidateMultiloopSigma(rw2.sigma).empty(),
                  "rewrite: one kink left, labels compact and valid");
            Check(rw2.MapEdge(1) == 0, "rewrite: deleted edge maps to 0");
            Check(rw2.MapEdge(2) == 1 && rw2.MapEdge(4) == 1,
                  "rewrite: merged arc keeps the smallest label, shifted down");
            Check(rw2.MapEdge(3) == 2, "rewrite: untouched edge shifts down over the gap");
            Check(rw2.circleArcs.empty(), "rewrite: no circle when a crossing survives");
        }

        // -- R1 inside a larger diagram: two kinks -> one kink --
        try {
            std::unique_ptr<Scene> s = buildScene("(1,-1,-4,2)(3,-3,-2,4)", 0);
            PhysicsParams params;
            Physics phys(*s->loop, s->neck, params, s->exteriors);
            phys.Step();
            const size_t faceBeadsBefore = s->neck.faceBeads.size();
            Check(phys.StartR1Move(1), "R1: move starts on a monogon");
            Check(s->neck.faceBeads.size() + 1 == faceBeadsBefore,
                  "R1: the monogon's face bead is erased when the move starts");
            Check(runMove(phys), "R1: loop strips to two beads");
            R1FaceInfo info;
            phys.GetR1Result(&info);
            MoveComponentResult res;
            std::unique_ptr<Multiloop> moved = FinalizeR1(*s->loop, s->neck, info, s->comps, &res);
            Check(moved != nullptr, "R1: finalize succeeds");
            if (moved) {
                Check(ValidateMultiloopSigma(moved->sigma).empty() && moved->numEdges == 2 &&
                          moved->sigma.Cycles().size() == 1,
                      "R1: one crossing and two compact edges remain");
                Check(moved->numStrands == 1, "R1: strand count preserved");
                Check(res.components.size() == 1 &&
                          res.components[0].kind == MultiloopPart::Crossings &&
                          res.components[0].exteriorFace >= 0 && !res.disconnected,
                      "R1: one crossing component with a remapped exterior");
                Check(linksOk(s->neck), "R1: necklace links are consistent");
                Check(walksOk(*moved, s->neck), "R1: every edge walks vertex to vertex");
                std::unordered_set<int> ext;
                for (const int f : ExteriorFaceIndices(res.components)) {
                    ext.insert(f);
                }
                Check(static_cast<int>(s->neck.faceBeads.size()) == finiteFaces(*moved, ext),
                      "R1: one face bead per finite face after the move");
                phys.EnsureExteriorHoneycombMins(*moved, res.exteriorArcs);
                bool physicsOk = true;
                try {
                    Physics phys2(*moved, s->neck, params, ext);
                    phys2.SetNesting(res.components);
                    phys2.StartPostRepair();
                    physicsOk = runPost(phys2);
                }
                catch (const std::exception&) {
                    physicsOk = false;
                }
                Check(physicsOk, "R1: post-repair on the rewritten diagram goes idle");
            }
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R1 two kinks: ") + e.what()).c_str());
        }

        // -- R1 down to no crossings: kink -> circle --
        try {
            std::unique_ptr<Scene> s = buildScene("(1,-1,2,-2)", 0);
            PhysicsParams params;
            Physics phys(*s->loop, s->neck, params, s->exteriors);
            phys.Step();
            Check(phys.StartR1Move(1), "R1 circle: move starts on the kink's monogon");
            Check(runMove(phys), "R1 circle: loop strips to two beads");
            R1FaceInfo info;
            phys.GetR1Result(&info);
            MoveComponentResult res;
            std::unique_ptr<Multiloop> moved = FinalizeR1(*s->loop, s->neck, info, s->comps, &res);
            Check(moved != nullptr, "R1 circle: finalize succeeds");
            if (moved) {
                Check(moved->sigma.Cycles().empty() && moved->numEdges == 0,
                      "R1 circle: no crossings remain");
                Check(res.components.size() == 1 && res.components[0].kind == MultiloopPart::Circle &&
                          res.components[0].strand == 0,
                      "R1 circle: the component is a circle on strand 0");
                Check(walksOk(*moved, s->neck), "R1 circle: beads form a closed ring");
                int pads = 0;
                for (const FaceBead& fb : s->neck.faceBeads) {
                    if (fb.face < 0 && fb.boundaryBeads.size() >= 3) {
                        ++pads;
                    }
                }
                Check(pads == 1 && s->neck.faceBeads.size() == 1,
                      "R1 circle: exactly one circle pad remains");
                bool physicsOk = true;
                try {
                    Physics phys2(*moved, s->neck, params, std::unordered_set<int>{});
                    phys2.SetNesting(res.components);
                    phys2.StartPostRepair();
                    physicsOk = runPost(phys2);
                }
                catch (const std::exception&) {
                    physicsOk = false;
                }
                Check(physicsOk, "R1 circle: post-repair on the circle goes idle");
            }
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R1 kink to circle: ") + e.what()).c_str());
        }

        // -- Connected R2 on a trefoil bigon: 3 crossings -> 1 --
        try {
            std::unique_ptr<Scene> s = buildScene("(1,-3,-6,4)(2,5,-1,-4)(3,-5,-2,6)", 1);
            R2FaceInfo probe;
            Check(AnalyzeR2Face(*s->loop, nullptr, 0, &probe) && !probe.disconnects,
                  "R2 trefoil: face 0 is a connected bigon");
            PhysicsParams params;
            Physics phys(*s->loop, s->neck, params, s->exteriors);
            phys.Step();
            {
                bool hitCenter = false;
                for (size_t i = 0; i < s->neck.faceBeads.size(); ++i) {
                    const FaceBead& fb = s->neck.faceBeads[i];
                    if (!phys.FaceBeadInsideFace(fb)) {
                        continue;
                    }
                    const int hit = phys.PickRegionAt({fb.pos.x, fb.pos.y});
                    if (hit == static_cast<int>(i)) {
                        hitCenter = true;
                        break;
                    }
                }
                Check(hitCenter, "R2 trefoil: PickRegionAt hits a face bead at its center");
                Check(phys.PickRegionAt({-100.0f, -100.0f}) < 0,
                      "R2 trefoil: PickRegionAt far away is the exterior");
            }
            Check(phys.StartR2Move(0), "R2 trefoil: move starts");
            Check(runMove(phys), "R2 trefoil: sides strip and X-Y become tangent");
            R2FaceInfo info;
            phys.GetR2Result(&info);
            {
                const int xb = info.vertexBead[0];
                const int yb = info.vertexBead[1];
                Check(xb >= 0 && yb >= 0 && xb < static_cast<int>(s->neck.beads.size()) &&
                          yb < static_cast<int>(s->neck.beads.size()),
                      "R2 trefoil: vertex beads still present after the spoke");
            }
            MoveComponentResult res;
            std::unique_ptr<Multiloop> moved = FinalizeR2(*s->loop, s->neck, info, s->comps, &res);
            Check(moved != nullptr, "R2 trefoil: finalize succeeds");
            if (moved) {
                Check(ValidateMultiloopSigma(moved->sigma).empty() && moved->numEdges == 2 &&
                          moved->sigma.Cycles().size() == 1,
                      "R2 trefoil: one crossing with compact labels remains");
                Check(moved->numStrands == 1, "R2 trefoil: strand count preserved");
                Check(!res.disconnected && res.components.size() == 1,
                      "R2 trefoil: diagram stays connected");
                Check(res.faceMap[3] >= 0 && res.faceMap[3] == res.faceMap[4],
                      "R2 trefoil: the two opposing faces map to one merged face");
                Check(res.mergeKeepLabel >= 0 && res.mergeDropLabel > res.mergeKeepLabel,
                      "R2 trefoil: a merge pair keeps the smaller label");
                Check(linksOk(s->neck), "R2 trefoil: necklace links are consistent");
                Check(walksOk(*moved, s->neck), "R2 trefoil: every edge walks vertex to vertex");
                std::unordered_set<int> ext;
                for (const int f : ExteriorFaceIndices(res.components)) {
                    ext.insert(f);
                }
                Check(static_cast<int>(s->neck.faceBeads.size()) == finiteFaces(*moved, ext) + 1,
                      "R2 trefoil: merged face briefly holds two beads");
                phys.EnsureExteriorHoneycombMins(*moved, res.exteriorArcs);
                bool highlightedDrop = false;
                for (FaceBead& fb : s->neck.faceBeads) {
                    if (fb.label == res.mergeDropLabel) {
                        fb.highlighted = true;
                        highlightedDrop = true;
                    }
                }
                Check(highlightedDrop, "R2 trefoil: the drop merge bead can be highlighted");
                bool physicsOk = true;
                try {
                    Physics phys2(*moved, s->neck, params, ext);
                    phys2.SetNesting(res.components);
                    phys2.SetFaceMergePair(res.mergeKeepLabel, res.mergeDropLabel);
                    Check(phys2.FaceMergePending(), "R2 trefoil: merge pair is pending");
                    phys2.StartPostRepair();
                    physicsOk = runPost(phys2);
                }
                catch (const std::exception&) {
                    physicsOk = false;
                }
                Check(physicsOk, "R2 trefoil: post-repair and face merge finish");
                Check(static_cast<int>(s->neck.faceBeads.size()) == finiteFaces(*moved, ext),
                      "R2 trefoil: face bead count drops to one per finite face");
                bool keeperSurvives = false;
                for (const FaceBead& fb : s->neck.faceBeads) {
                    keeperSurvives = keeperSurvives || fb.label == res.mergeKeepLabel;
                }
                Check(keeperSurvives, "R2 trefoil: the smaller label survives the merge");
                bool keeperHighlighted = false;
                for (const FaceBead& fb : s->neck.faceBeads) {
                    if (fb.label == res.mergeKeepLabel) {
                        keeperHighlighted = fb.highlighted;
                    }
                }
                Check(keeperHighlighted,
                      "R2 trefoil: highlighting the dropped face leaves the survivor highlighted");
            }
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R2 trefoil: ") + e.what()).c_str());
        }

        // -- Disconnecting R2 through the exterior: Hopf shadow -> two circles --
        try {
            std::unique_ptr<Scene> s = buildScene("(1,3,-2,4)(2,-3,-1,-4)", 2);
            R2FaceInfo probe;
            Check(AnalyzeR2Face(*s->loop, nullptr, 1, &probe) && probe.disconnects &&
                      probe.faceL == 2,
                  "R2 Hopf: the lens disconnects through the exterior");
            PhysicsParams params;
            Physics phys(*s->loop, s->neck, params, s->exteriors);
            phys.Step();
            Check(phys.StartR2Move(1), "R2 Hopf exterior: move starts");
            Check(runMove(phys), "R2 Hopf exterior: animation converges");
            R2FaceInfo info;
            phys.GetR2Result(&info);
            MoveComponentResult res;
            std::unique_ptr<Multiloop> moved = FinalizeR2(*s->loop, s->neck, info, s->comps, &res);
            Check(moved != nullptr, "R2 Hopf exterior: finalize succeeds");
            if (moved) {
                Check(moved->sigma.Cycles().empty(), "R2 Hopf exterior: no crossings remain");
                Check(res.disconnected && res.components.size() == 2 &&
                          res.components[0].kind == MultiloopPart::Circle &&
                          res.components[1].kind == MultiloopPart::Circle,
                      "R2 Hopf exterior: two circle components");
                Check(!res.components[0].Nested() && !res.components[1].Nested(),
                      "R2 Hopf exterior: both circles are top level");
                Check(res.components[0].strand != res.components[1].strand,
                      "R2 Hopf exterior: circles keep distinct strands");
                Check(walksOk(*moved, s->neck), "R2 Hopf exterior: both rings close");
                Check(s->neck.faceBeads.size() == 2, "R2 Hopf exterior: one pad per circle");
                bool physicsOk = true;
                try {
                    Physics phys2(*moved, s->neck, params, std::unordered_set<int>{});
                    phys2.SetNesting(res.components);
                    phys2.StartPostRepair();
                    physicsOk = runPost(phys2);
                }
                catch (const std::exception&) {
                    physicsOk = false;
                }
                Check(physicsOk, "R2 Hopf exterior: post-repair goes idle");
            }
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R2 Hopf exterior: ") + e.what()).c_str());
        }

        // -- Disconnecting R2 through a finite face: one circle nested in the other --
        try {
            std::unique_ptr<Scene> s = buildScene("(1,3,-2,4)(2,-3,-1,-4)", 2);
            R2FaceInfo probe;
            Check(AnalyzeR2Face(*s->loop, nullptr, 0, &probe) && probe.disconnects &&
                      probe.faceL == 3,
                  "R2 Hopf nested: the crescent disconnects through a finite face");
            const int opposing = R2OpposingHostHalfEdge(*s->loop, probe, 2);
            Check(opposing != 0, "R2 Hopf nested: opposing host half-edge exists");
            PhysicsParams params;
            Physics phys(*s->loop, s->neck, params, s->exteriors);
            phys.Step();
            Check(phys.StartR2Move(0), "R2 Hopf nested: move starts");
            {
                R2FaceInfo started;
                phys.GetR2Result(&started);
                const int vb = (opposing == started.a || opposing == started.b)
                                   ? started.vertexBead[0]
                                   : started.vertexBead[1];
                int interiors = 0;
                int prev = vb;
                int cur = s->neck.NeighborToward(vb, opposing);
                int guard = 0;
                while (cur >= 0 && !s->neck.beads[static_cast<size_t>(cur)].isVertex &&
                       guard < static_cast<int>(s->neck.beads.size())) {
                    if (!s->neck.beads[static_cast<size_t>(cur)].removed) {
                        ++interiors;
                    }
                    const Bead& b = s->neck.beads[static_cast<size_t>(cur)];
                    const int nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
                    prev = cur;
                    cur = nxt;
                    ++guard;
                }
                Check(interiors >= 3, "R2 Hopf nested: opposing host edge has at least 3 interiors");
            }
            Check(runMove(phys), "R2 Hopf nested: animation converges");
            R2FaceInfo info;
            phys.GetR2Result(&info);
            MoveComponentResult res;
            std::unique_ptr<Multiloop> moved = FinalizeR2(*s->loop, s->neck, info, s->comps, &res);
            Check(moved != nullptr, "R2 Hopf nested: finalize succeeds");
            if (moved) {
                Check(res.disconnected && res.components.size() == 2,
                      "R2 Hopf nested: two components");
                int nestedCount = 0;
                int outerStrand = -1;
                int innerStrand = -1;
                for (const MultiloopComponent& c : res.components) {
                    if (c.Nested()) {
                        ++nestedCount;
                        innerStrand = c.strand;
                        Check(c.parentCircle >= 0 && c.parentFace < 0,
                              "R2 Hopf nested: the hole sits inside a circle");
                    }
                    else {
                        outerStrand = c.strand;
                    }
                }
                Check(nestedCount == 1, "R2 Hopf nested: exactly one component is nested");
                bool parentMatches = false;
                for (const MultiloopComponent& c : res.components) {
                    if (c.Nested()) {
                        parentMatches = c.parentCircle == outerStrand;
                    }
                }
                Check(parentMatches, "R2 Hopf nested: parentCircle is the outer circle's strand");
                Check(walksOk(*moved, s->neck), "R2 Hopf nested: both rings close");
                int innerBeads = 0;
                for (const Bead& b : s->neck.beads) {
                    if (!b.removed && !b.isVertex && b.strand == innerStrand) {
                        ++innerBeads;
                    }
                }
                Check(innerBeads >= 6, "R2 Hopf nested: inner ring has at least 6 beads");
                bool physicsOk = true;
                try {
                    Physics phys2(*moved, s->neck, params, std::unordered_set<int>{});
                    phys2.SetNesting(res.components);
                    phys2.Step();
                    int outerHoles = 0;
                    for (const FaceBead& fb : s->neck.faceBeads) {
                        if (fb.face < 0 && !fb.holeBeads.empty()) {
                            ++outerHoles;
                        }
                    }
                    Check(outerHoles == 1, "R2 Hopf nested: the outer pad lists one hole ring");
                    phys2.StartPostRepair();
                    physicsOk = runPost(phys2);
                    for (int i = 0; i < 200; ++i) {
                        phys2.Step();
                    }
                    // Outer pad: inside its ring, outside the hole. Inner pad:
                    // inside the inner ring. Rings keep a diameter apart.
                    bool placementOk = true;
                    for (const FaceBead& fb : s->neck.faceBeads) {
                        std::vector<Vec2d> poly;
                        std::vector<std::vector<Vec2d>> h;
                        if (!phys2.FaceBeadBoundaryPolygon(fb, &poly, &h)) {
                            placementOk = false;
                            continue;
                        }
                        const Vec2d c{static_cast<double>(fb.pos.x), static_cast<double>(fb.pos.y)};
                        placementOk = placementOk && PointInPolygonWithHoles(c, poly, h);
                    }
                    Check(placementOk, "R2 Hopf nested: pads sit inside their regions and off the hole");
                    float minRing = 1e9f;
                    for (const Bead& a : s->neck.beads) {
                        if (a.removed || a.isVertex || a.strand != outerStrand) {
                            continue;
                        }
                        for (const Bead& b : s->neck.beads) {
                            if (b.removed || b.isVertex || b.strand != innerStrand) {
                                continue;
                            }
                            const float dx = a.pos.x - b.pos.x;
                            const float dy = a.pos.y - b.pos.y;
                            minRing = std::min(minRing, std::sqrt(dx * dx + dy * dy));
                        }
                    }
                    Check(minRing + 0.15f >= 2.0f * s->neck.radius,
                          "R2 Hopf nested: hole ring and outer ring stay a diameter apart");
                }
                catch (const std::exception&) {
                    physicsOk = false;
                }
                Check(physicsOk, "R2 Hopf nested: post-repair goes idle");
            }
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R2 Hopf nested: ") + e.what()).c_str());
        }

        // -- Disconnecting R2 leaving a crossing component nested in a circle --
        // Exterior (2,3) is bounded by the circle's arc 2 and the deleted arc
        // 3, so the circle owns the exterior and the kink ends up inside it.
        try {
            std::unique_ptr<Scene> s = buildScene("(1,3,-2,-6)(2,-3,-1,4)(5,-5,-4,6)", 3);
            R2FaceInfo probe;
            Check(AnalyzeR2Face(*s->loop, nullptr, 1, &probe) && probe.disconnects &&
                      probe.faceL == 2,
                  "R2 kinked link: the lens disconnects through the finite 4-gon");
            PhysicsParams params;
            Physics phys(*s->loop, s->neck, params, s->exteriors);
            phys.Step();
            Check(phys.StartR2Move(1), "R2 kinked link: move starts");
            Check(runMove(phys), "R2 kinked link: animation converges");
            R2FaceInfo info;
            phys.GetR2Result(&info);
            MoveComponentResult res;
            std::unique_ptr<Multiloop> moved = FinalizeR2(*s->loop, s->neck, info, s->comps, &res);
            Check(moved != nullptr, "R2 kinked link: finalize succeeds");
            if (moved) {
                Check(ValidateMultiloopSigma(moved->sigma).empty() &&
                          moved->sigma.Cycles().size() == 1 && moved->numEdges == 2,
                      "R2 kinked link: the kink survives with compact labels");
                Check(res.disconnected && res.components.size() == 2,
                      "R2 kinked link: two components");
                const MultiloopComponent* circle = nullptr;
                const MultiloopComponent* kink = nullptr;
                for (const MultiloopComponent& c : res.components) {
                    if (c.kind == MultiloopPart::Circle) {
                        circle = &c;
                    }
                    else {
                        kink = &c;
                    }
                }
                Check(circle != nullptr && kink != nullptr,
                      "R2 kinked link: one circle and one crossing component");
                if (circle && kink) {
                    Check(!circle->Nested(), "R2 kinked link: the circle is top level");
                    Check(kink->parentCircle == circle->strand && kink->parentFace < 0,
                          "R2 kinked link: the kink is nested inside the circle");
                    Check(kink->exteriorFace >= 0 &&
                              kink->exteriorFace < static_cast<int>(moved->faces.size()),
                          "R2 kinked link: the nested kink has an exterior face");
                }
                Check(linksOk(s->neck) && walksOk(*moved, s->neck),
                      "R2 kinked link: necklace is consistent");
                std::unordered_set<int> ext;
                for (const int f : ExteriorFaceIndices(res.components)) {
                    ext.insert(f);
                }
                bool physicsOk = true;
                try {
                    Physics phys2(*moved, s->neck, params, ext);
                    phys2.SetNesting(res.components);
                    phys2.Step();
                    int padsWithHoles = 0;
                    for (const FaceBead& fb : s->neck.faceBeads) {
                        if (fb.face < 0 && !fb.holeBeads.empty()) {
                            ++padsWithHoles;
                        }
                    }
                    Check(padsWithHoles == 1, "R2 kinked link: the circle pad lists the kink as a hole");
                    phys2.StartPostRepair();
                    physicsOk = runPost(phys2);
                }
                catch (const std::exception&) {
                    physicsOk = false;
                }
                Check(physicsOk, "R2 kinked link: post-repair goes idle");
            }
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R2 kinked link: ") + e.what()).c_str());
        }

        // -- Same diagram, exterior (-1,-6,-4): after the move only the kink
        // bounds it, so the circle nests in a face of the kink --
        try {
            std::unique_ptr<Scene> s = buildScene("(1,3,-2,-6)(2,-3,-1,4)(5,-5,-4,6)", 0);
            PhysicsParams params;
            Physics phys(*s->loop, s->neck, params, s->exteriors);
            phys.Step();
            Check(phys.StartR2Move(1), "R2 face hole: move starts");
            Check(runMove(phys), "R2 face hole: animation converges");
            R2FaceInfo info;
            phys.GetR2Result(&info);
            MoveComponentResult res;
            std::unique_ptr<Multiloop> moved = FinalizeR2(*s->loop, s->neck, info, s->comps, &res);
            Check(moved != nullptr, "R2 face hole: finalize succeeds");
            if (moved) {
                const MultiloopComponent* circle = nullptr;
                const MultiloopComponent* kink = nullptr;
                for (const MultiloopComponent& c : res.components) {
                    if (c.kind == MultiloopPart::Circle) {
                        circle = &c;
                    }
                    else {
                        kink = &c;
                    }
                }
                Check(circle != nullptr && kink != nullptr && !kink->Nested(),
                      "R2 face hole: the kink component is top level");
                std::unordered_set<int> ext;
                for (const int f : ExteriorFaceIndices(res.components)) {
                    ext.insert(f);
                }
                if (circle && kink) {
                    Check(circle->parentFace >= 0 && circle->parentCircle < 0 &&
                              ext.count(circle->parentFace) == 0,
                          "R2 face hole: the circle sits inside a finite face of the kink");
                    Physics phys2(*moved, s->neck, params, ext);
                    phys2.SetNesting(res.components);
                    phys2.Step();
                    Check(phys2.FaceHasHoles(circle->parentFace),
                          "R2 face hole: physics reports the parent face as holed");
                    Check(!phys2.StartR1Move(circle->parentFace) &&
                              !phys2.StartR2Move(circle->parentFace) &&
                              !phys2.StartR3Move(circle->parentFace),
                          "R2 face hole: moves on a holed face are refused");
                    int holedBeads = 0;
                    for (const FaceBead& fb : s->neck.faceBeads) {
                        if (fb.face == circle->parentFace && !fb.holeBeads.empty()) {
                            ++holedBeads;
                        }
                    }
                    Check(holedBeads == 1, "R2 face hole: the parent face bead lists the ring");
                    phys2.StartPostRepair();
                    Check(runPost(phys2), "R2 face hole: post-repair goes idle");
                }
            }
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R2 face hole: ") + e.what()).c_str());
        }

        // -- Two disconnecting R2s: three concentric circles --
        try {
            const auto faceWith = [](const Multiloop& loop, int a, int b) {
                for (int f = 0; f < static_cast<int>(loop.faces.size()); ++f) {
                    const std::vector<int>& cyc = loop.faces[static_cast<size_t>(f)];
                    if (cyc.size() == 2 &&
                        ((cyc[0] == a && cyc[1] == b) || (cyc[0] == b && cyc[1] == a))) {
                        return f;
                    }
                }
                return -1;
            };
            const char* text = "(-1,6,7,8)(1,-8,-7,-4)(-2,4,-3,5)(2,-5,3,-6)";
            const Multiloop probe(Permutation::FromString(text));
            const int ext = faceWith(probe, 2, 5);
            const int firstFace = faceWith(probe, 1, 8);
            Check(ext >= 0 && firstFace >= 0, "R2 concentric: found exterior (2,5) and bigon (1,8)");
            std::unique_ptr<Scene> s = buildScene(text, ext);
            PhysicsParams params;
            {
                Physics phys(*s->loop, s->neck, params, s->exteriors);
                phys.SetNesting(s->comps);
                phys.Step();
                Check(phys.StartR2Move(firstFace), "R2 concentric: first move starts");
                Check(runMove(phys), "R2 concentric: first animation converges");
                R2FaceInfo info;
                phys.GetR2Result(&info);
                MoveComponentResult res;
                std::unique_ptr<Multiloop> moved =
                    FinalizeR2(*s->loop, s->neck, info, s->comps, &res);
                Check(moved != nullptr, "R2 concentric: first finalize succeeds");
                if (moved) {
                    s->loop = std::move(moved);
                    s->comps = std::move(res.components);
                    s->exteriors.clear();
                    for (const int f : ExteriorFaceIndices(s->comps)) {
                        s->exteriors.insert(f);
                    }
                }
            }
            int secondFace = -1;
            for (int f = 0; f < static_cast<int>(s->loop->faces.size()); ++f) {
                R2FaceInfo p;
                if (s->exteriors.count(f) != 0) {
                    continue;
                }
                if (AnalyzeR2Face(*s->loop, nullptr, f, &p) && p.disconnects) {
                    secondFace = f;
                    break;
                }
            }
            Check(secondFace >= 0, "R2 concentric: a remaining disconnecting bigon");
            if (secondFace >= 0) {
                Physics phys(*s->loop, s->neck, params, s->exteriors);
                phys.SetNesting(s->comps);
                phys.Step();
                Check(phys.StartR2Move(secondFace), "R2 concentric: second move starts");
                Check(runMove(phys), "R2 concentric: second animation converges");
                R2FaceInfo info;
                phys.GetR2Result(&info);
                MoveComponentResult res;
                std::unique_ptr<Multiloop> moved =
                    FinalizeR2(*s->loop, s->neck, info, s->comps, &res);
                Check(moved != nullptr, "R2 concentric: second finalize succeeds");
                if (moved) {
                    s->comps = std::move(res.components);
                    int circles = 0;
                    int nested = 0;
                    const MultiloopComponent* root = nullptr;
                    const MultiloopComponent* mid = nullptr;
                    const MultiloopComponent* inn = nullptr;
                    for (const MultiloopComponent& c : s->comps) {
                        if (c.kind == MultiloopPart::Circle) {
                            ++circles;
                        }
                        if (c.Nested()) {
                            ++nested;
                        }
                    }
                    Check(circles == 3 && nested == 2,
                          "R2 concentric: three circles, two nested");
                    for (const MultiloopComponent& c : s->comps) {
                        if (c.kind != MultiloopPart::Circle) {
                            continue;
                        }
                        if (!c.Nested()) {
                            root = &c;
                        }
                    }
                    if (root != nullptr) {
                        for (const MultiloopComponent& c : s->comps) {
                            if (c.kind == MultiloopPart::Circle &&
                                c.parentCircle == root->strand) {
                                mid = &c;
                            }
                        }
                    }
                    if (mid != nullptr) {
                        for (const MultiloopComponent& c : s->comps) {
                            if (c.kind == MultiloopPart::Circle &&
                                c.parentCircle == mid->strand) {
                                inn = &c;
                            }
                        }
                    }
                    Check(root != nullptr && mid != nullptr && inn != nullptr,
                          "R2 concentric: red-in-blue-in-green parent chain");
                    Physics phys2(*moved, s->neck, params, std::unordered_set<int>{});
                    phys2.SetNesting(s->comps);
                    phys2.Step();
                    int padsWithHoles = 0;
                    int doubleHoles = 0;
                    for (const FaceBead& fb : s->neck.faceBeads) {
                        if (fb.face < 0 && !fb.holeBeads.empty()) {
                            ++padsWithHoles;
                            if (fb.holeBeads.size() >= 2) {
                                ++doubleHoles;
                            }
                        }
                    }
                    Check(padsWithHoles >= 2 && doubleHoles == 0,
                          "R2 concentric: each outer pad lists one hole, not two siblings");
                }
            }
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R2 concentric: ") + e.what()).c_str());
        }

        // -- Cancel restores the erased face bead --
        try {
            std::unique_ptr<Scene> s = buildScene("(1,-3,-6,4)(2,5,-1,-4)(3,-5,-2,6)", 1);
            PhysicsParams params;
            Physics phys(*s->loop, s->neck, params, s->exteriors);
            phys.Step();
            const size_t before = s->neck.faceBeads.size();
            Check(phys.StartR2Move(0), "R2 cancel: move starts");
            Check(s->neck.faceBeads.size() + 1 == before, "R2 cancel: face bead erased");
            for (int i = 0; i < 30 && phys.MoveActive(); ++i) {
                phys.MoveStep();
            }
            phys.CancelMove();
            Check(!phys.MoveActive() && s->neck.faceBeads.size() == before,
                  "R2 cancel: face bead restored and move idle");
            bool haveFace0 = false;
            for (const FaceBead& fb : s->neck.faceBeads) {
                haveFace0 = haveFace0 || fb.face == 0;
            }
            Check(haveFace0, "R2 cancel: the restored bead is back on its face");
            phys.Step();
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R2 cancel: ") + e.what()).c_str());
        }

        const auto liveNonVertex = [](const Necklace& neck) {
            int n = 0;
            for (const Bead& b : neck.beads) {
                if (!b.removed && !b.isVertex) {
                    ++n;
                }
            }
            return n;
        };
        const auto liveOfStrand = [](const Necklace& neck, int strand) {
            int n = 0;
            for (const Bead& b : neck.beads) {
                if (!b.removed && !b.isVertex && b.strand == strand) {
                    ++n;
                }
            }
            return n;
        };
        const auto firstCirclePad = [](const Necklace& neck) {
            for (int i = 0; i < static_cast<int>(neck.faceBeads.size()); ++i) {
                if (neck.faceBeads[static_cast<size_t>(i)].face < 0) {
                    return i;
                }
            }
            return neck.faceBeads.empty() ? -1 : 0;
        };

        // -- R0 on one of two packed circles --
        try {
            const SigmaNotation two = ParseSigmaNotation("()()");
            const Permutation empty;
            const Multiloop dummy(empty);
            std::vector<MultiloopComponent> comps = two.components;
            const PackedGridLayout packed = ComputePackedGridLayout(dummy, comps);
            Necklace neck = BuildPackedNecklace(dummy, packed, comps);
            PhysicsParams params;
            Physics phys(dummy, neck, params, std::unordered_set<int>{});
            phys.Step();
            const int fi = firstCirclePad(neck);
            Check(fi >= 0 && phys.StartR0Move(fi), "R0 two circles: move starts on a pad");
            Check(runMove(phys), "R0 two circles: ring strips to three beads");
            Check(liveNonVertex(neck) >= 9, "R0 two circles: the other ring is still live");
            R0MoveInfo info;
            phys.GetR0Result(&info);
            const int popped = info.strand;
            Check(FinalizeR0(neck, info, &comps), "R0 two circles: finalize succeeds");
            Check(comps.size() == 1 && comps[0].kind == MultiloopPart::Circle,
                  "R0 two circles: one Circle component remains");
            Check(liveOfStrand(neck, popped) == 0, "R0 two circles: popped strand is gone");
            int otherLive = 0;
            for (const Bead& b : neck.beads) {
                if (!b.removed && !b.isVertex && b.strand != popped) {
                    ++otherLive;
                }
            }
            Check(otherLive >= 3, "R0 two circles: the surviving ring has beads");
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R0 two circles: ") + e.what()).c_str());
        }

        // -- R0 on the last remaining circle (unknot necklace) --
        try {
            Necklace neck = BuildUnknotNecklace();
            PhysicsParams params;
            Physics phys(neck, params);
            phys.Step();
            std::vector<MultiloopComponent> comps(1);
            comps[0].kind = MultiloopPart::Circle;
            comps[0].strand = 0;
            Check(phys.StartR0Move(0), "R0 last circle: move starts");
            Check(neck.faceBeads.empty(), "R0 last circle: face bead erased");
            Check(runMove(phys), "R0 last circle: animation completes");
            Check(liveNonVertex(neck) == 3, "R0 last circle: three beads remain before pop");
            R0MoveInfo info;
            phys.GetR0Result(&info);
            Check(FinalizeR0(neck, info, &comps), "R0 last circle: finalize succeeds");
            Check(comps.empty(), "R0 last circle: no components remain");
            Check(liveNonVertex(neck) == 0, "R0 last circle: the ring is gone");
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R0 last circle: ") + e.what()).c_str());
        }

        // -- R0 of a packed circle next to a kink: sigma unchanged --
        try {
            std::unique_ptr<Scene> s = buildScene("()(1,-1,2,-2)", 0);
            PhysicsParams params;
            Physics phys(*s->loop, s->neck, params, s->exteriors);
            phys.Step();
            const int fi = firstCirclePad(s->neck);
            Check(fi >= 0 && phys.StartR0Move(fi), "R0 mixed: move starts on the circle");
            Check(runMove(phys), "R0 mixed: animation completes");
            R0MoveInfo info;
            phys.GetR0Result(&info);
            const int nEdges = s->loop->numEdges;
            Check(FinalizeR0(s->neck, info, &s->comps), "R0 mixed: finalize succeeds");
            Check(s->loop->numEdges == nEdges && s->loop->sigma.Cycles().size() == 1,
                  "R0 mixed: kink sigma is unchanged");
            Check(s->comps.size() == 1 && s->comps[0].kind == MultiloopPart::Crossings,
                  "R0 mixed: only the crossing component remains");
            Check(linksOk(s->neck) && walksOk(*s->loop, s->neck),
                  "R0 mixed: necklace still walks the kink");
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R0 mixed: ") + e.what()).c_str());
        }

        // -- Nested Hopf circles: refuse the outer, allow the inner --
        try {
            std::unique_ptr<Scene> s = buildScene("(1,3,-2,4)(2,-3,-1,-4)", 2);
            PhysicsParams params;
            Physics phys(*s->loop, s->neck, params, s->exteriors);
            phys.Step();
            Check(phys.StartR2Move(0), "R0 nested: R2 setup starts");
            Check(runMove(phys), "R0 nested: R2 animation converges");
            R2FaceInfo r2;
            phys.GetR2Result(&r2);
            MoveComponentResult res;
            std::unique_ptr<Multiloop> moved = FinalizeR2(*s->loop, s->neck, r2, s->comps, &res);
            Check(moved != nullptr && res.disconnected && res.components.size() == 2,
                  "R0 nested: Hopf R2 yields two nested circles");
            if (moved) {
                Physics phys2(*moved, s->neck, params, std::unordered_set<int>{});
                phys2.SetNesting(res.components);
                phys2.Step();
                int outerFi = -1;
                int innerFi = -1;
                for (int i = 0; i < static_cast<int>(s->neck.faceBeads.size()); ++i) {
                    const FaceBead& fb = s->neck.faceBeads[static_cast<size_t>(i)];
                    if (fb.face >= 0) {
                        continue;
                    }
                    if (!fb.holeBeads.empty()) {
                        outerFi = i;
                    }
                    else {
                        innerFi = i;
                    }
                }
                Check(outerFi >= 0 && innerFi >= 0, "R0 nested: found outer and inner pads");
                Check(!phys2.StartR0Move(outerFi), "R0 nested: outer circle with a hole is refused");
                Check(phys2.StartR0Move(innerFi), "R0 nested: inner circle starts");
                Check(runMove(phys2), "R0 nested: inner animation completes");
                R0MoveInfo info;
                phys2.GetR0Result(&info);
                Check(FinalizeR0(s->neck, info, &res.components), "R0 nested: inner finalize succeeds");
                Check(res.components.size() == 1 && res.components[0].kind == MultiloopPart::Circle &&
                          !res.components[0].Nested(),
                      "R0 nested: the outer circle remains un-nested");
            }
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R0 nested: ") + e.what()).c_str());
        }

        // -- Cancel restores the erased circle pad --
        try {
            Necklace neck = BuildUnknotNecklace();
            PhysicsParams params;
            Physics phys(neck, params);
            phys.Step();
            const size_t before = neck.faceBeads.size();
            Check(phys.StartR0Move(0), "R0 cancel: move starts");
            Check(neck.faceBeads.size() + 1 == before, "R0 cancel: face bead erased");
            for (int i = 0; i < 30 && phys.MoveActive(); ++i) {
                phys.MoveStep();
            }
            phys.CancelMove();
            Check(!phys.MoveActive() && neck.faceBeads.size() == before,
                  "R0 cancel: face bead restored and move idle");
            phys.Step();
        }
        catch (const std::exception& e) {
            Check(false, (std::string("R0 cancel: ") + e.what()).c_str());
        }
    }

    // ---- Tangle helpers: opposite-slot linking and open segment crossings ----
    {
        const Vec2d o{0.0, 0.0};
        const Vec2d e{1.0, 0.0};
        const Vec2d n{0.0, 1.0};
        const Vec2d w{-1.0, 0.0};
        const Vec2d s{0.0, -1.0};
        Check(OppositeSlotsLinked(o, e, n, w, s),
              "tangle: CCW NESW opposite slots are linked");
        Check(OppositeSlotsLinked(o, n, w, s, e),
              "tangle: rotated opposite slots are linked");
        // Angular order 0,1,3,2: pairs {1,3} sit as consecutive samples.
        Check(!OppositeSlotsLinked(o, e, n, s, w),
              "tangle: swapped last pair is unlinked");
        const Vec2d e2{1.0, 0.1};
        const Vec2d w2{-1.0, 0.1};
        Check(!OppositeSlotsLinked(o, e, w, e2, w2),
              "tangle: two pairs on opposite sides are unlinked");
        Check(!OppositeSlotsLinked(o, o, n, w, s),
              "tangle: neighbor at the center is unlinked");

        Vector2 a{0.0f, 0.0f};
        Vector2 b{1.0f, 1.0f};
        Vector2 c{0.0f, 1.0f};
        Vector2 d{1.0f, 0.0f};
        Check(SegmentsProperlyCross(a, b, c, d, nullptr),
              "tangle: crossing diagonals of a square properly intersect");
        Check(!SegmentsProperlyCross(a, c, a, d, nullptr),
              "tangle: segments sharing an endpoint do not properly intersect");
        Check(!SegmentsProperlyCross(a, c, d, b, nullptr),
              "tangle: disjoint square sides do not intersect");
    }

    if (checksFailed == 0) {
        std::printf("[permutation] all %d checks passed\n", checksRun);
    }
    else {
        std::printf("[permutation] %d of %d checks FAILED\n", checksFailed, checksRun);
    }
    return checksFailed == 0;
}
