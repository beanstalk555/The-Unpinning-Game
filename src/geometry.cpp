#include "geometry.h"

#include <cmath>

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
