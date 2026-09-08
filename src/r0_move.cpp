#include "r0_move.h"

namespace
{
    int PadBeadStrand(const Necklace& neck, const FaceBead& fb)
    {
        for (const int bead : fb.boundaryBeads) {
            if (bead < 0 || bead >= static_cast<int>(neck.beads.size())) {
                continue;
            }
            const Bead& b = neck.beads[static_cast<size_t>(bead)];
            if (!b.removed && !b.isVertex) {
                return b.strand;
            }
        }
        return -1;
    }

    bool IsCirclePad(const Necklace& neck, const FaceBead& fb)
    {
        if (fb.face < 0) {
            return true;
        }
        for (const FaceBead& other : neck.faceBeads) {
            if (other.face < 0) {
                return false;
            }
        }
        return true;
    }
}

int CircleComponentIndexForPad(const Necklace& neck,
                               const std::vector<MultiloopComponent>& components,
                               int faceBeadIndex)
{
    if (faceBeadIndex < 0 ||
        faceBeadIndex >= static_cast<int>(neck.faceBeads.size())) {
        return -1;
    }
    const FaceBead& fb = neck.faceBeads[static_cast<size_t>(faceBeadIndex)];
    const int strand = PadBeadStrand(neck, fb);
    if (strand >= 0) {
        int found = -1;
        int matches = 0;
        for (size_t i = 0; i < components.size(); ++i) {
            if (components[i].kind == MultiloopPart::Circle &&
                components[i].strand == strand) {
                found = static_cast<int>(i);
                ++matches;
            }
        }
        if (matches == 1) {
            return found;
        }
    }

    int ordinal = 0;
    for (int i = 0; i < faceBeadIndex; ++i) {
        if (IsCirclePad(neck, neck.faceBeads[static_cast<size_t>(i)])) {
            ++ordinal;
        }
    }
    int seen = 0;
    for (size_t i = 0; i < components.size(); ++i) {
        if (components[i].kind != MultiloopPart::Circle) {
            continue;
        }
        if (seen == ordinal) {
            return static_cast<int>(i);
        }
        ++seen;
    }
    return -1;
}

bool FinalizeR0(Necklace& neck, const R0MoveInfo& info,
                std::vector<MultiloopComponent>* components)
{
    if (components == nullptr) {
        return false;
    }

    const auto markBead = [&](int bead) {
        if (bead < 0 || bead >= static_cast<int>(neck.beads.size())) {
            return;
        }
        Bead& b = neck.beads[static_cast<size_t>(bead)];
        if (!b.isVertex) {
            b.removed = true;
        }
    };
    for (const int bead : info.ringBeads) {
        markBead(bead);
    }
    if (info.strand >= 0) {
        for (size_t i = 0; i < neck.beads.size(); ++i) {
            const Bead& b = neck.beads[i];
            if (!b.removed && !b.isVertex && b.strand == info.strand) {
                markBead(static_cast<int>(i));
            }
        }
    }

    int idx = -1;
    if (info.strand >= 0) {
        int found = -1;
        int matches = 0;
        for (size_t i = 0; i < components->size(); ++i) {
            if ((*components)[i].kind == MultiloopPart::Circle &&
                (*components)[i].strand == info.strand) {
                found = static_cast<int>(i);
                ++matches;
            }
        }
        if (matches == 1) {
            idx = found;
        }
    }
    if (idx < 0 && info.circleOrdinal >= 0) {
        int seen = 0;
        for (size_t i = 0; i < components->size(); ++i) {
            if ((*components)[i].kind != MultiloopPart::Circle) {
                continue;
            }
            if (seen == info.circleOrdinal) {
                idx = static_cast<int>(i);
                break;
            }
            ++seen;
        }
    }
    if (idx < 0) {
        return false;
    }
    components->erase(components->begin() + static_cast<long>(idx));
    return true;
}
