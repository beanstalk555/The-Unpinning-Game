#include "beads.h"

#include "geometry.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace
{
    void Require(bool condition, const char* message)
    {
        if (!condition) {
            throw std::runtime_error(std::string("Necklace failed: ") + message);
        }
    }

    float BoundaryInteriorDot(const Necklace& neck, int bead, int prev, int next)
    {
        const BeadPos& p = neck.beads[static_cast<size_t>(bead)].pos;
        const BeadPos& pa = neck.beads[static_cast<size_t>(prev)].pos;
        const BeadPos& pc = neck.beads[static_cast<size_t>(next)].pos;
        const float bax = pa.x - p.x;
        const float bay = pa.y - p.y;
        const float bcx = pc.x - p.x;
        const float bcy = pc.y - p.y;
        const float la = std::sqrt(bax * bax + bay * bay);
        const float lc = std::sqrt(bcx * bcx + bcy * bcy);
        if (la < 1e-6f || lc < 1e-6f) {
            return -1.0f;
        }
        return (bax * bcx + bay * bcy) / (la * lc);
    }

    // Points along the polyline every half grid unit, endpoints excluded.
    std::vector<BeadPos> SamplePathEveryHalfUnit(const std::vector<GridPoint>& path)
    {
        std::vector<int> cumulative = {0}; // segment ends, in half units
        for (size_t i = 0; i + 1 < path.size(); ++i) {
            const int len = std::abs(path[i + 1].x - path[i].x) + std::abs(path[i + 1].y - path[i].y);
            cumulative.push_back(cumulative.back() + 2 * len);
        }
        const int total = cumulative.back();

        std::vector<BeadPos> samples;
        size_t seg = 0;
        for (int s = 1; s < total; ++s) {
            while (seg + 1 < cumulative.size() && s > cumulative[seg + 1]) {
                ++seg;
            }
            const int span = cumulative[seg + 1] - cumulative[seg];
            if (span <= 0) {
                continue;
            }
            const float frac = static_cast<float>(s - cumulative[seg]) / static_cast<float>(span);
            const GridPoint& a = path[seg];
            const GridPoint& b = path[seg + 1];
            samples.push_back({a.x + (b.x - a.x) * frac, a.y + (b.y - a.y) * frac});
        }
        return samples;
    }
}

int Necklace::NeighborToward(int vertexBead, int halfEdge) const
{
    if (vertexBead < 0 || vertexBead >= static_cast<int>(beads.size())) {
        return -1;
    }
    const Bead& bead = beads[static_cast<size_t>(vertexBead)];
    for (size_t i = 0; i < bead.sigmaCycle.size(); ++i) {
        if (bead.sigmaCycle[i] == halfEdge) {
            return bead.neighborOfHalfEdge[i];
        }
    }
    return -1;
}

void FaceBead::WalkFaceBoundary(const Multiloop& loop, const Necklace& neck,
                                const std::unordered_map<int, int>& crossingOf) const
{
    (void)loop;
    WalkFaceBoundaryFromPhiCycle(phiCycle, neck, crossingOf, boundaryBeads, boundaryDots);
}

void FaceBead::RefreshBoundaryDots(const Necklace& neck) const
{
    const int n = static_cast<int>(boundaryBeads.size());
    if (n < 3) {
        boundaryDots.clear();
        return;
    }
    boundaryDots.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const int bead = boundaryBeads[static_cast<size_t>(i)];
        const int prev = boundaryBeads[static_cast<size_t>((i + n - 1) % n)];
        const int next = boundaryBeads[static_cast<size_t>((i + 1) % n)];
        boundaryDots[static_cast<size_t>(i)] = BoundaryInteriorDot(neck, bead, prev, next);
    }
}

void WalkFaceBoundaryFromPhiCycle(const std::vector<int>& phiCycle, const Necklace& neck,
                                  const std::unordered_map<int, int>& crossingOf,
                                  std::vector<int>& boundaryBeads,
                                  std::vector<float>& boundaryDots,
                                  std::unordered_set<int>* sharpCornerOut,
                                  float sharpDotThreshold,
                                  bool computeDots)
{
    boundaryBeads.clear();
    boundaryDots.clear();
    const bool collectSharp = sharpCornerOut != nullptr;

    const auto append = [&](int bead) {
        if (bead < 0) {
            return;
        }
        if (!neck.beads[static_cast<size_t>(bead)].isVertex &&
            neck.beads[static_cast<size_t>(bead)].removed) {
            return;
        }
        if (!boundaryBeads.empty() && boundaryBeads.back() == bead) {
            return;
        }
        boundaryBeads.push_back(bead);
    };

    const int limit = static_cast<int>(neck.beads.size()) + 8;
    for (const int h : phiCycle) {
        const auto it = crossingOf.find(h);
        if (it == crossingOf.end() || it->second < 0 ||
            it->second >= static_cast<int>(neck.crossingBead.size())) {
            continue;
        }
        const int tail = neck.crossingBead[static_cast<size_t>(it->second)];
        if (tail < 0 || tail >= static_cast<int>(neck.beads.size())) {
            continue;
        }
        append(tail);

        int prev = tail;
        int cur = neck.NeighborToward(tail, h);
        int guard = 0;
        while (cur >= 0 && !neck.beads[static_cast<size_t>(cur)].isVertex) {
            append(cur);
            const Bead& b = neck.beads[static_cast<size_t>(cur)];
            const int nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
            prev = cur;
            cur = nxt;
            if (++guard > limit) {
                break;
            }
        }
    }

    const int n = static_cast<int>(boundaryBeads.size());
    if (n < 3 || !computeDots) {
        return;
    }
    boundaryDots.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const int bead = boundaryBeads[static_cast<size_t>(i)];
        const int prev = boundaryBeads[static_cast<size_t>((i + n - 1) % n)];
        const int next = boundaryBeads[static_cast<size_t>((i + 1) % n)];
        const float dot = BoundaryInteriorDot(neck, bead, prev, next);
        boundaryDots.push_back(dot);
        if (collectSharp && !neck.beads[static_cast<size_t>(bead)].isVertex &&
            dot >= sharpDotThreshold) {
            sharpCornerOut->insert(bead);
        }
    }
}

double FaceBead::PolyArea(const Necklace& neck) const
{
    if (boundaryBeads.size() < 3) {
        return 0.0;
    }
    std::vector<Vec2d> poly;
    poly.reserve(boundaryBeads.size());
    for (const int bead : boundaryBeads) {
        const BeadPos& p = neck.beads[static_cast<size_t>(bead)].pos;
        poly.push_back(Vec2d{static_cast<double>(p.x), static_cast<double>(p.y)});
    }
    return ShoelaceArea(poly);
}

double FaceBead::ExcessArea(const Necklace& neck, float beadRadius) const
{
    if (boundaryBeads.size() < 3) {
        return 0.0;
    }
    const double area = PolyArea(neck);
    const double halfR2 = 0.5 * static_cast<double>(beadRadius) * beadRadius;
    double pie = 0.0;
    for (const float dot : boundaryDots) {
        const float clamped = std::max(-1.0f, std::min(1.0f, dot));
        pie += halfR2 * std::acos(static_cast<double>(clamped));
    }
    return area - pie;
}

double FaceBead::IsoperimetricRatio(const Necklace& neck) const
{
    const size_t n = boundaryBeads.size();
    if (n < 3) {
        return 0.0;
    }
    const double area = PolyArea(neck);
    const double perimeterBeads = static_cast<double>(n);
    return area / (perimeterBeads * perimeterBeads);
}

Necklace BuildNecklace(const Multiloop& loop, const GridLayout& layout, int exteriorFace)
{
    std::unordered_set<int> exteriors;
    exteriors.insert(exteriorFace);
    return BuildNecklace(loop, layout, exteriors);
}

Necklace BuildNecklace(const Multiloop& loop, const GridLayout& layout,
                       const std::unordered_set<int>& exteriorFaces)
{
    Necklace necklace;
    necklace.radius = 0.25f;

    const std::vector<std::vector<int>> crossings = loop.sigma.Cycles();
    const int numCrossings = static_cast<int>(crossings.size());
    std::unordered_map<int, int> crossingOf;
    for (int c = 0; c < numCrossings; ++c) {
        for (const int h : crossings[c]) {
            crossingOf[h] = c;
        }
    }

    // Vertex beads.
    necklace.crossingBead.resize(numCrossings);
    for (int c = 0; c < numCrossings; ++c) {
        Bead bead;
        bead.isVertex = true;
        bead.pos = {static_cast<float>(layout.crossingPositions[c].x),
                    static_cast<float>(layout.crossingPositions[c].y)};
        bead.crossing = c;
        bead.sigmaCycle = crossings[c];
        necklace.crossingBead[c] = static_cast<int>(necklace.beads.size());
        necklace.beads.push_back(std::move(bead));
    }

    // Edge beads, linked along the +k direction of each edge.
    const auto slotOf = [&](int vertexBead, int halfEdge) -> int& {
        Bead& bead = necklace.beads[vertexBead];
        for (size_t i = 0; i < bead.sigmaCycle.size(); ++i) {
            if (bead.sigmaCycle[i] == halfEdge) {
                return bead.neighborOfHalfEdge[i];
            }
        }
        Require(false, "half-edge missing from sigma cycle");
        return bead.neighborOfHalfEdge[0];
    };
    for (int k = 1; k <= loop.numEdges; ++k) {
        Require(k - 1 < static_cast<int>(layout.edgePaths.size()), "missing edge path");
        const std::vector<BeadPos> samples = SamplePathEveryHalfUnit(layout.edgePaths[k - 1]);
        Require(!samples.empty(), "edge with no room for beads");
        const int vbStart = necklace.crossingBead[crossingOf[k]];
        const int vbEnd = necklace.crossingBead[crossingOf[-k]];
        const int first = static_cast<int>(necklace.beads.size());
        const int count = static_cast<int>(samples.size());
        for (int i = 0; i < count; ++i) {
            Bead bead;
            bead.isVertex = false;
            bead.edge = k;
            bead.pos = samples[i];
            bead.neighbors[0] = (i == 0) ? vbStart : first + i - 1;
            bead.neighbors[1] = (i == count - 1) ? vbEnd : first + i + 1;
            necklace.beads.push_back(std::move(bead));
        }
        slotOf(vbStart, k) = first;
        slotOf(vbEnd, -k) = first + count - 1;
    }
    for (const int vb : necklace.crossingBead) {
        for (int i = 0; i < 4; ++i) {
            Require(necklace.beads[vb].neighborOfHalfEdge[i] >= 0, "unfilled vertex bead slot");
        }
    }

    // Face beads for every finite face, placed at the lexicographic top-left
    // corner of the face boundary, nudged inward by one radius diagonally.
    int label = 1;
    for (int f = 0; f < static_cast<int>(loop.faces.size()); ++f) {
        if (exteriorFaces.count(f) != 0) {
            continue;
        }
        FaceBead faceBead;
        faceBead.label = label++;
        faceBead.face = f;
        faceBead.phiCycle = loop.faces[f];
        bool first = true;
        GridPoint topLeft = {0, 0};
        for (const int h : loop.faces[f]) {
            faceBead.vertexBeads.push_back(necklace.crossingBead[crossingOf[h]]);
            for (const GridPoint& p : layout.edgePaths[std::abs(h) - 1]) {
                if (first || p.y > topLeft.y || (p.y == topLeft.y && p.x < topLeft.x)) {
                    topLeft = p;
                    first = false;
                }
            }
        }
        faceBead.pos = {topLeft.x + 2.0f*necklace.radius, topLeft.y - 2.0f*necklace.radius};
        necklace.faceBeads.push_back(std::move(faceBead));
    }

    return necklace;
}

void AppendCircleToNecklace(Necklace& necklace, BeadPos origin, int label, int strand)
{
    const BeadPos positions[8] = {
        {0.0f, 0.0f}, {0.5f, 0.0f}, {1.0f, 0.0f}, {1.0f, 0.5f},
        {1.0f, 1.0f}, {0.5f, 1.0f}, {0.0f, 1.0f}, {0.0f, 0.5f},
    };
    const int first = static_cast<int>(necklace.beads.size());
    for (int i = 0; i < 8; ++i) {
        Bead bead;
        bead.isVertex = false;
        bead.edge = 0;
        bead.strand = strand;
        bead.pos = {origin.x + positions[i].x, origin.y + positions[i].y};
        bead.neighbors[0] = first + ((i + 7) % 8);
        bead.neighbors[1] = first + ((i + 1) % 8);
        necklace.beads.push_back(std::move(bead));
    }
    FaceBead faceBead;
    faceBead.label = label;
    faceBead.face = -1;
    faceBead.pos = {origin.x + 0.5f, origin.y + 0.5f};
    for (int i = 0; i < 8; ++i) {
        faceBead.boundaryBeads.push_back(first + i);
    }
    necklace.faceBeads.push_back(std::move(faceBead));
}

Necklace BuildPackedNecklace(const Multiloop& loop, const PackedGridLayout& packed,
                             const std::vector<MultiloopComponent>& components)
{
    std::unordered_set<int> exteriors;
    for (const MultiloopComponent& component : components) {
        if (component.kind == MultiloopPart::Crossings && component.exteriorFace >= 0) {
            exteriors.insert(component.exteriorFace);
        }
    }

    Necklace necklace;
    if (loop.numEdges > 0) {
        necklace = BuildNecklace(loop, packed.layout, exteriors);
    }
    else {
        necklace.radius = 0.25f;
    }

    // Tamassia paths stay in local integer coords; shift each crossing
    // component's beads by its float pack origin.
    for (size_t i = 0; i < components.size(); ++i) {
        if (components[i].kind != MultiloopPart::Crossings ||
            i >= packed.componentOrigin.size()) {
            continue;
        }
        const PackOrigin origin = packed.componentOrigin[i];
        if (origin.x == 0.0f && origin.y == 0.0f) {
            continue;
        }
        std::unordered_set<int> owned(components[i].halfEdges.begin(),
                                      components[i].halfEdges.end());
        for (Bead& bead : necklace.beads) {
            bool mine = false;
            if (bead.isVertex) {
                for (const int h : bead.sigmaCycle) {
                    if (owned.count(h) != 0) {
                        mine = true;
                        break;
                    }
                }
            }
            else if (owned.count(bead.edge) != 0 || owned.count(-bead.edge) != 0) {
                mine = true;
            }
            if (mine) {
                bead.pos.x += origin.x;
                bead.pos.y += origin.y;
            }
        }
        for (FaceBead& fb : necklace.faceBeads) {
            if (ComponentOwnsFace(components[i], loop, fb.face)) {
                fb.pos.x += origin.x;
                fb.pos.y += origin.y;
            }
        }
    }

    int nextLabel = 1;
    for (const FaceBead& fb : necklace.faceBeads) {
        nextLabel = std::max(nextLabel, fb.label + 1);
    }
    int nextStrand = loop.numStrands;
    for (size_t i = 0; i < components.size(); ++i) {
        if (components[i].kind != MultiloopPart::Circle) {
            continue;
        }
        const PackOrigin origin =
            (i < packed.componentOrigin.size()) ? packed.componentOrigin[i] : PackOrigin{};
        const int strand =
            (components[i].strand >= 0) ? components[i].strand : nextStrand++;
        AppendCircleToNecklace(necklace, {origin.x, origin.y}, nextLabel++, strand);
    }
    return necklace;
}

Necklace BuildUnknotNecklace()
{
    Necklace necklace;
    necklace.radius = 0.25f;
    const BeadPos positions[8] = {
        {0.0f, 0.0f}, {0.5f, 0.0f}, {1.0f, 0.0f}, {1.0f, 0.5f},
        {1.0f, 1.0f}, {0.5f, 1.0f}, {0.0f, 1.0f}, {0.0f, 0.5f},
    };
    for (int i = 0; i < 8; ++i) {
        Bead bead;
        bead.isVertex = false;
        bead.edge = 1;
        bead.strand = 0;
        bead.pos = positions[i];
        bead.neighbors[0] = (i + 7) % 8;
        bead.neighbors[1] = (i + 1) % 8;
        necklace.beads.push_back(std::move(bead));
    }
    FaceBead faceBead;
    faceBead.label = 1;
    faceBead.face = 0;
    faceBead.phiCycle = {1};
    faceBead.pos = {0.5f, 0.5f};
    necklace.faceBeads.push_back(std::move(faceBead));
    return necklace;
}

GridLayout UnknotGridLayout()
{
    GridLayout layout;
    layout.width = 1;
    layout.height = 1;
    return layout;
}
