#pragma once

#include <vector>

#include "beads.h"
#include "multiloop.h"

// Combinatorial data for a shortening Reidemeister 0 move on a crossing-free
// circle. Physics stores this while the ring strips to three beads; finalize
// then pops those beads and drops the Circle component. Sigma is unchanged.
struct R0MoveInfo
{
    int strand = -1;
    // Index of this pad among circle pads at the moment the move started,
    // used when Circle components do not yet carry a strand id.
    int circleOrdinal = -1;
    std::vector<int> ringBeads;
};

// The k-th circle pad in `neck` corresponds to the k-th Circle component,
// unless a unique component already stores this pad's bead strand.
int CircleComponentIndexForPad(const Necklace& neck,
                               const std::vector<MultiloopComponent>& components,
                               int faceBeadIndex);

// Marks the remaining ring beads removed and erases the matching Circle
// entry. Returns false if no Circle component can be identified.
bool FinalizeR0(Necklace& neck, const R0MoveInfo& info,
                std::vector<MultiloopComponent>* components);
