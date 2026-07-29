#include "spatial_disk_grid.h"

#include <cmath>

void SpatialDiskGrid::SetCellSize(float cellSize)
{
    if (cellSize > 1e-6f) {
        cellSize_ = cellSize;
    }
}

void SpatialDiskGrid::Clear()
{
    cells_.clear();
}

void SpatialDiskGrid::Reserve(size_t bodyHint)
{
    cells_.reserve(bodyHint / 4 + 1);
}

void SpatialDiskGrid::Insert(int body, float x, float y)
{
    const int ix = static_cast<int>(std::floor(x / cellSize_));
    const int iy = static_cast<int>(std::floor(y / cellSize_));
    cells_[CellKey(ix, iy)].bodies.push_back(body);
}

int64_t SpatialDiskGrid::CellKey(int ix, int iy)
{
    return (static_cast<int64_t>(ix) << 32) ^
           static_cast<uint32_t>(iy);
}

void SpatialDiskGrid::ForEachCandidatePair(const std::function<void(int, int)>& fn) const
{
    for (const auto& entry : cells_) {
        const int64_t key = entry.first;
        const int ix = static_cast<int>(key >> 32);
        const int iy = static_cast<int>(static_cast<uint32_t>(key));
        const std::vector<int>& here = entry.second.bodies;

        for (size_t i = 0; i < here.size(); ++i) {
            for (size_t j = i + 1; j < here.size(); ++j) {
                const int a = here[i];
                const int b = here[j];
                if (a < b) {
                    fn(a, b);
                }
                else {
                    fn(b, a);
                }
            }
        }

        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0) {
                    continue;
                }
                const auto it = cells_.find(CellKey(ix + dx, iy + dy));
                if (it == cells_.end()) {
                    continue;
                }
                const std::vector<int>& there = it->second.bodies;
                for (const int a : here) {
                    for (const int b : there) {
                        if (a < b) {
                            fn(a, b);
                        }
                        else if (b < a) {
                            fn(b, a);
                        }
                    }
                }
            }
        }
    }
}
