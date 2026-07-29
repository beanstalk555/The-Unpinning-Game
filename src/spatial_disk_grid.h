#pragma once

#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>

// Uniform-grid broad phase for disk-shaped bounding volumes (depth-1 BVH).
// Each body occupies one grid cell; candidates are pairs in the same or
// adjacent cells. Cell size is the bounding disk diameter (extendable to a
// deeper tree later).
class SpatialDiskGrid
{
public:
    void SetCellSize(float cellSize);
    void Clear();
    void Reserve(size_t bodyHint);

    void Insert(int body, float x, float y);

    // Invokes fn(a, b) with a < b for each unordered broad-phase pair.
    void ForEachCandidatePair(const std::function<void(int, int)>& fn) const;

private:
    static int64_t CellKey(int ix, int iy);

    float cellSize_ = 1.0f;
    struct Bucket
    {
        std::vector<int> bodies;
    };
    std::unordered_map<int64_t, Bucket> cells_;
};
