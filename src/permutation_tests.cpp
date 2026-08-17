#include "permutation.h"
#include "multiloop.h"
#include "orthogonal.h"
#include "beads.h"
#include "physics.h"
#include "physics_params.h"
#include "r3_move.h"
#include "random_multiloop.h"
#include "drawing.h"
#include "geometry.h"

#include <set>

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
    Check(!ValidateMultiloopSigma(Permutation::FromString("(1,2,-1,-2)(3,-3,4,-4)")).empty(),
          "validation rejects a disconnected diagram");
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
                // Pre-honeycomb settle rounds + slide need a larger frame budget.
                while (phys.R3Active() && !phys.R3Complete() && frames < 12000) {
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

    if (checksFailed == 0) {
        std::printf("[permutation] all %d checks passed\n", checksRun);
    }
    else {
        std::printf("[permutation] %d of %d checks FAILED\n", checksFailed, checksRun);
    }
    return checksFailed == 0;
}
