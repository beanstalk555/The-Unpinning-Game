#include "geometry.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <vector>

double SignedShoelaceArea(const std::vector<Vec2d>& vertices)
{
    const int n = static_cast<int>(vertices.size());
    if (n < 3) {
        return 0.0;
    }
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        const Vec2d& p = vertices[static_cast<size_t>(i)];
        const Vec2d& q = vertices[static_cast<size_t>((i + 1) % n)];
        sum += p.x * q.y - q.x * p.y;
    }
    return sum * 0.5;
}

double ShoelaceArea(const std::vector<Vec2d>& vertices)
{
    return std::abs(SignedShoelaceArea(vertices));
}

bool PointInPolygon(const Vec2d& point, const std::vector<Vec2d>& vertices)
{
    const int n = static_cast<int>(vertices.size());
    if (n < 3) {
        return false;
    }
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        const Vec2d& a = vertices[static_cast<size_t>(j)];
        const Vec2d& b = vertices[static_cast<size_t>(i)];
        // Count an edge only when it straddles the ray y = point.y, with one
        // endpoint strictly above (standard even-odd vertex handling).
        if ((a.y > point.y) == (b.y > point.y)) {
            continue;
        }
        const double dy = b.y - a.y;
        if (std::fabs(dy) < 1e-18) {
            continue;
        }
        const double xHit = (b.x - a.x) * (point.y - a.y) / dy + a.x;
        if (point.x < xHit) {
            inside = !inside;
        }
    }
    return inside;
}

namespace
{
    double DistToSegmentSq(const Vec2d& p, const Vec2d& a, const Vec2d& b)
    {
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double len2 = dx * dx + dy * dy;
        double t = 0.0;
        if (len2 > 1e-18) {
            t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2;
            t = std::max(0.0, std::min(1.0, t));
        }
        const double qx = a.x + t * dx - p.x;
        const double qy = a.y + t * dy - p.y;
        return qx * qx + qy * qy;
    }

    double MinDistSqToRing(const Vec2d& p, const std::vector<Vec2d>& poly)
    {
        const int n = static_cast<int>(poly.size());
        double minSq = 1e300;
        for (int i = 0; i < n; ++i) {
            minSq = std::min(minSq, DistToSegmentSq(p, poly[static_cast<size_t>(i)],
                                                    poly[static_cast<size_t>((i + 1) % n)]));
        }
        return minSq;
    }

    // Positive inside (outer and outside every hole), negative otherwise;
    // magnitude is distance to the nearest ring.
    double SignedBoundaryDist(const Vec2d& p, const std::vector<Vec2d>& outer,
                              const std::vector<std::vector<Vec2d>>& holes)
    {
        double minSq = MinDistSqToRing(p, outer);
        for (const std::vector<Vec2d>& hole : holes) {
            if (hole.size() >= 3) {
                minSq = std::min(minSq, MinDistSqToRing(p, hole));
            }
        }
        const double d = std::sqrt(minSq);
        return PointInPolygonWithHoles(p, outer, holes) ? d : -d;
    }

    struct PolyCell
    {
        double x = 0.0;
        double y = 0.0;
        double h = 0.0; // half-width of the square cell
        double d = 0.0; // signed distance at the cell center
        double max = 0.0; // upper bound on inscribed radius in this cell
    };

    struct PolyCellBetter
    {
        bool operator()(const PolyCell& a, const PolyCell& b) const { return a.max < b.max; }
    };

    PolyCell MakeCell(double x, double y, double h, const std::vector<Vec2d>& poly,
                      const std::vector<std::vector<Vec2d>>& holes)
    {
        PolyCell c;
        c.x = x;
        c.y = y;
        c.h = h;
        c.d = SignedBoundaryDist(Vec2d{x, y}, poly, holes);
        c.max = c.d + h * 1.4142135623730951;
        return c;
    }
}

bool PointInPolygonWithHoles(const Vec2d& point, const std::vector<Vec2d>& outer,
                             const std::vector<std::vector<Vec2d>>& holes)
{
    if (!PointInPolygon(point, outer)) {
        return false;
    }
    for (const std::vector<Vec2d>& hole : holes) {
        if (hole.size() >= 3 && PointInPolygon(point, hole)) {
            return false;
        }
    }
    return true;
}

namespace
{
    double EarClipCross(const Vec2d& a, const Vec2d& b, const Vec2d& c)
    {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    }

    // True when p is inside or on the boundary (needed on a grid, where a
    // reflex vertex can sit on a candidate ear's diagonal).
    bool PointInTriangleClosed(const Vec2d& p, const Vec2d& a, const Vec2d& b, const Vec2d& c)
    {
        const double c0 = EarClipCross(a, b, p);
        const double c1 = EarClipCross(b, c, p);
        const double c2 = EarClipCross(c, a, p);
        const double eps = 1e-12;
        return (c0 >= -eps && c1 >= -eps && c2 >= -eps) ||
               (c0 <= eps && c1 <= eps && c2 <= eps);
    }
}

bool EarClipTriangles(const std::vector<Vec2d>& vertices, std::vector<int>* out, bool* usedFan)
{
    if (usedFan != nullptr) {
        *usedFan = false;
    }
    if (out == nullptr) {
        return false;
    }
    out->clear();
    const int n0 = static_cast<int>(vertices.size());
    if (n0 < 3) {
        return false;
    }
    std::vector<int> idx(static_cast<size_t>(n0));
    for (int i = 0; i < n0; ++i) {
        idx[static_cast<size_t>(i)] = i;
    }
    const double area = SignedShoelaceArea(vertices);
    if (std::abs(area) < 1e-18) {
        return false;
    }
    const bool ccw = area > 0.0;

    auto isConvex = [&](int p, int c, int nx) {
        const double cr = EarClipCross(vertices[static_cast<size_t>(idx[static_cast<size_t>(p)])],
                                       vertices[static_cast<size_t>(idx[static_cast<size_t>(c)])],
                                       vertices[static_cast<size_t>(idx[static_cast<size_t>(nx)])]);
        return ccw ? cr > 1e-12 : cr < -1e-12;
    };
    auto isFlatOrDup = [&](int p, int c, int nx) {
        const Vec2d& A = vertices[static_cast<size_t>(idx[static_cast<size_t>(p)])];
        const Vec2d& B = vertices[static_cast<size_t>(idx[static_cast<size_t>(c)])];
        const Vec2d& C = vertices[static_cast<size_t>(idx[static_cast<size_t>(nx)])];
        const double dx = B.x - A.x;
        const double dy = B.y - A.y;
        if (dx * dx + dy * dy < 1e-18) {
            return true;
        }
        const double cr = EarClipCross(A, B, C);
        return std::fabs(cr) <= 1e-12;
    };

    int guard = n0 * n0;
    while (static_cast<int>(idx.size()) > 3 && guard-- > 0) {
        bool clipped = false;
        const int n = static_cast<int>(idx.size());
        for (int i = 0; i < n; ++i) {
            const int p = (i + n - 1) % n;
            const int nx = (i + 1) % n;
            // Collinear mid-edge beads are not convex ears; dropping them
            // keeps the outline and lets true corners clip.
            if (isFlatOrDup(p, i, nx)) {
                idx.erase(idx.begin() + i);
                clipped = true;
                break;
            }
            if (!isConvex(p, i, nx)) {
                continue;
            }
            const Vec2d& A = vertices[static_cast<size_t>(idx[static_cast<size_t>(p)])];
            const Vec2d& B = vertices[static_cast<size_t>(idx[static_cast<size_t>(i)])];
            const Vec2d& C = vertices[static_cast<size_t>(idx[static_cast<size_t>(nx)])];
            bool ear = true;
            for (int j = 0; j < n; ++j) {
                if (j == p || j == i || j == nx) {
                    continue;
                }
                if (PointInTriangleClosed(
                        vertices[static_cast<size_t>(idx[static_cast<size_t>(j)])], A, B, C)) {
                    ear = false;
                    break;
                }
            }
            if (!ear) {
                continue;
            }
            out->push_back(idx[static_cast<size_t>(p)]);
            out->push_back(idx[static_cast<size_t>(i)]);
            out->push_back(idx[static_cast<size_t>(nx)]);
            idx.erase(idx.begin() + i);
            clipped = true;
            break;
        }
        if (!clipped) {
            if (usedFan != nullptr) {
                *usedFan = true;
            }
            const int n = static_cast<int>(idx.size());
            for (int i = 1; i + 1 < n; ++i) {
                out->push_back(idx[0]);
                out->push_back(idx[static_cast<size_t>(i)]);
                out->push_back(idx[static_cast<size_t>(i + 1)]);
            }
            idx.clear();
            break;
        }
    }
    if (idx.size() == 3) {
        out->push_back(idx[0]);
        out->push_back(idx[1]);
        out->push_back(idx[2]);
    }
    if (!ccw) {
        for (size_t i = 0; i + 2 < out->size(); i += 3) {
            std::swap((*out)[i], (*out)[i + 1]);
        }
    }
    return out->size() >= 3;
}

InscribedCircle LargestInscribedCircle(const std::vector<Vec2d>& vertices, double precision)
{
    return LargestInscribedCircle(vertices, {}, precision);
}

InscribedCircle LargestInscribedCircle(const std::vector<Vec2d>& vertices,
                                       const std::vector<std::vector<Vec2d>>& holes,
                                       double precision)
{
    InscribedCircle result;
    const int n = static_cast<int>(vertices.size());
    if (n < 3) {
        return result;
    }

    double minX = vertices[0].x;
    double minY = vertices[0].y;
    double maxX = minX;
    double maxY = minY;
    double cx = 0.0;
    double cy = 0.0;
    for (const Vec2d& p : vertices) {
        minX = std::min(minX, p.x);
        minY = std::min(minY, p.y);
        maxX = std::max(maxX, p.x);
        maxY = std::max(maxY, p.y);
        cx += p.x;
        cy += p.y;
    }
    cx /= static_cast<double>(n);
    cy /= static_cast<double>(n);
    const double width = maxX - minX;
    const double height = maxY - minY;
    const double cellSize0 = std::min(width, height);
    if (cellSize0 < 1e-12) {
        return result;
    }
    if (precision <= 0.0) {
        precision = std::max(1e-4, 0.01 * cellSize0);
    }

    PolyCell best = MakeCell(cx, cy, 0.0, vertices, holes);
    const PolyCell bboxMid =
        MakeCell(minX + width * 0.5, minY + height * 0.5, 0.0, vertices, holes);
    if (bboxMid.d > best.d) {
        best = bboxMid;
    }

    std::priority_queue<PolyCell, std::vector<PolyCell>, PolyCellBetter> queue;
    // With holes the seed grid must be fine enough to land cells in the
    // annular region; a hole can swallow the single coarse cell.
    double cellSize = cellSize0;
    if (!holes.empty()) {
        cellSize = cellSize0 / 4.0;
    }
    double h = cellSize / 2.0;
    for (double x = minX; x < maxX; x += cellSize) {
        for (double y = minY; y < maxY; y += cellSize) {
            queue.push(MakeCell(x + h, y + h, h, vertices, holes));
        }
    }

    int guard = 0;
    const int guardMax = 20000;
    while (!queue.empty() && guard++ < guardMax) {
        const PolyCell cell = queue.top();
        queue.pop();
        if (cell.d > best.d) {
            best = cell;
        }
        // Max-heap on `max`: if this cell cannot beat `best`, neither can any
        // remaining cell.
        if (cell.max - best.d <= precision) {
            break;
        }
        h = cell.h / 2.0;
        queue.push(MakeCell(cell.x - h, cell.y - h, h, vertices, holes));
        queue.push(MakeCell(cell.x + h, cell.y - h, h, vertices, holes));
        queue.push(MakeCell(cell.x - h, cell.y + h, h, vertices, holes));
        queue.push(MakeCell(cell.x + h, cell.y + h, h, vertices, holes));
    }

    if (best.d > 0.0) {
        result.center = Vec2d{best.x, best.y};
        result.radius = best.d;
        result.found = true;
    }
    return result;
}

bool OppositeSlotsLinked(const Vec2d& center, const Vec2d& p0, const Vec2d& p1,
                         const Vec2d& p2, const Vec2d& p3)
{
    const Vec2d pts[4] = {p0, p1, p2, p3};
    struct Sample
    {
        double angle;
        int slot;
    };
    Sample samples[4];
    for (int i = 0; i < 4; ++i) {
        const double dx = pts[i].x - center.x;
        const double dy = pts[i].y - center.y;
        if (dx * dx + dy * dy < 1e-18) {
            return false;
        }
        samples[i].angle = std::atan2(dy, dx);
        samples[i].slot = i;
    }
    std::sort(samples, samples + 4, [](const Sample& a, const Sample& b) {
        if (a.angle != b.angle) {
            return a.angle < b.angle;
        }
        return a.slot < b.slot;
    });
    for (int i = 0; i < 4; ++i) {
        const int a = samples[i].slot;
        const int b = samples[(i + 1) % 4].slot;
        // Pairs are {0,2} (even) and {1,3} (odd). Adjacent same-pair → unlinked.
        if ((a & 1) == (b & 1)) {
            return false;
        }
    }
    return true;
}
