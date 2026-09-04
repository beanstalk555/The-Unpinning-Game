#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "geometry.h"
#include "multiloop.h"
#include "orthogonal.h"

struct Necklace;
struct BeadPos
{
    float x = 0.0f;
    float y = 0.0f;
};

// A bead of the necklace: either a vertex bead (at a crossing) or an edge
// bead (along an edge path, spaced half a grid unit apart). Beads form a
// doubly linked structure by index into Necklace::beads, supporting O(1)
// removal later.
struct Bead
{
    bool isVertex = false;
    // Edge beads can be removed by the auto-tighten pass. Removed beads stay
    // in the vector (so indices remain stable) but are unlinked from their
    // neighbors and skipped by physics and rendering.
    bool removed = false;
    BeadPos pos;

    // Edge beads: which edge label k (1..n) the bead lies on, and its two
    // neighbors along the +k direction (vertex or edge beads).
    int edge = -1;
    int neighbors[2] = {-1, -1};
    // Color index for crossing-free circle beads (and packed unknot rings).
    // Crossing edge beads keep -1 and use Multiloop::strandOfEdge instead.
    int strand = -1;

    // Vertex beads: the crossing index, its sigma cycle (counterclockwise
    // half-edges), and the adjacent edge bead in each half-edge direction:
    // neighborOfHalfEdge[i] lies toward sigmaCycle[i].
    int crossing = -1;
    std::vector<int> sigmaCycle;
    int neighborOfHalfEdge[4] = {-1, -1, -1, -1};
};

// One bead per finite face, placed at the lexicographic top-left corner of
// the face (offset inward by one bead radius diagonally).
struct FaceBead
{
    int label = 0; // 1..k, user-facing
    int face = -1; // index into loop.faces
    BeadPos pos;
    std::vector<int> phiCycle;
    std::vector<int> vertexBeads; // vertex beads around the face, in phi order

    // Cached boundary walk (refreshed during Physics::RebuildStructures).
    mutable std::vector<int> boundaryBeads;
    // Dot product of unit vectors to prev/next boundary beads (interior angle).
    mutable std::vector<float> boundaryDots;
    // Excess area at the last structure rebuild (topology-fresh positions).
    double cachedExcessArea = 0.0;

    // Walk boundary beads in phi order; fills boundaryBeads and boundaryDots.
    void WalkFaceBoundary(const Multiloop& loop, const Necklace& neck,
                          const std::unordered_map<int, int>& crossingOf) const;

    // Recompute boundaryDots from boundaryBeads at current bead positions.
    void RefreshBoundaryDots(const Necklace& neck) const;

    // Shoelace area of boundaryBeads (call WalkFaceBoundary first).
    double PolyArea(const Necklace& neck) const;

    // PolyArea minus pie slices (half r^2 theta) for each boundary bead.
    double ExcessArea(const Necklace& neck, float beadRadius) const;

    // PolyArea / (boundary bead count)^2.
    double IsoperimetricRatio(const Necklace& neck) const;
};

struct Necklace
{
    float radius = 0.25f; // in grid units (x/4 for grid length x)
    std::vector<Bead> beads;
    std::vector<int> crossingBead; // crossing index -> vertex bead index
    std::vector<FaceBead> faceBeads;

    // Edge bead indices on a face boundary with interior angle <= pi/3 + epsilon.
    std::vector<int> sharpCornerEdges;

    int NeighborToward(int vertexBead, int halfEdge) const;
};

// Walk a face's phi cycle; fills boundary bead indices and interior-angle dots.
// When sharpCornerOut is non-null, edge beads with dot >= sharpDotThreshold are inserted.
void WalkFaceBoundaryFromPhiCycle(const std::vector<int>& phiCycle, const Necklace& neck,
                                  const std::unordered_map<int, int>& crossingOf,
                                  std::vector<int>& boundaryBeads,
                                  std::vector<float>& boundaryDots,
                                  std::unordered_set<int>* sharpCornerOut = nullptr,
                                  float sharpDotThreshold = 2.0f,
                                  bool computeDots = true);

Necklace BuildNecklace(const Multiloop& loop, const GridLayout& layout, int exteriorFace);
Necklace BuildNecklace(const Multiloop& loop, const GridLayout& layout,
                       const std::unordered_set<int>& exteriorFaces);
Necklace BuildPackedNecklace(const Multiloop& loop, const PackedGridLayout& packed,
                             const std::vector<MultiloopComponent>& components);

// Unit-square necklace for the 1-region unknot (no crossings): 8 edge beads
// around the boundary and one interior face bead.
Necklace BuildUnknotNecklace();
GridLayout UnknotGridLayout();
void AppendCircleToNecklace(Necklace& necklace, BeadPos origin, int label, int strand);
