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
