#include "drawing.h"

#include "multiloop.h"

#include <algorithm>
#include <cmath>

namespace
{
    float Cross(Vector2 a, Vector2 b)
    {
        return a.x * b.y - a.y * b.x;
    }

    Vector2 Sub(Vector2 a, Vector2 b)
    {
        return Vector2{a.x - b.x, a.y - b.y};
    }

    // Proper intersection of open segments ab and cd, returning the
    // parameters t (along ab) and u (along cd) when they cross transversally.
    bool CrossParams(Vector2 a, Vector2 b, Vector2 c, Vector2 d, float* t, float* u)
    {
        const Vector2 r = Sub(b, a);
        const Vector2 s = Sub(d, c);
        const float denom = Cross(r, s);
        if (std::fabs(denom) < 1e-9f) {
            return false; // parallel or degenerate
        }
        const Vector2 ac = Sub(c, a);
        const float tt = Cross(ac, s) / denom;
        const float uu = Cross(ac, r) / denom;
        const float eps = 1e-6f;
        if (tt <= eps || tt >= 1.0f - eps || uu <= eps || uu >= 1.0f - eps) {
            return false;
        }
        *t = tt;
        *u = uu;
        return true;
    }

    struct Crossing
    {
        Vector2 pos;
    };

    // One visit to a crossing while walking along a curve.
    struct Appearance
    {
        int crossing;
        int segIndex;  // segment of the curve the crossing lies on
        float param;   // position along that segment
        Vector2 dir;   // unit travel direction along that segment
    };
}

bool SegmentsProperlyCross(Vector2 a, Vector2 b, Vector2 c, Vector2 d, Vector2* out)
{
    float t = 0.0f;
    float u = 0.0f;
    if (!CrossParams(a, b, c, d, &t, &u)) {
        return false;
    }
    if (out != nullptr) {
        *out = Vector2{a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t};
    }
    return true;
}

float DistPointToSegment(Vector2 p, Vector2 a, Vector2 b)
{
    const Vector2 ab = Sub(b, a);
    const float len2 = ab.x * ab.x + ab.y * ab.y;
    if (len2 < 1e-12f) {
        return std::sqrt((p.x - a.x) * (p.x - a.x) + (p.y - a.y) * (p.y - a.y));
    }
    float t = ((p.x - a.x) * ab.x + (p.y - a.y) * ab.y) / len2;
    t = std::max(0.0f, std::min(1.0f, t));
    const float dx = p.x - (a.x + ab.x * t);
    const float dy = p.y - (a.y + ab.y * t);
    return std::sqrt(dx * dx + dy * dy);
}

DrawnMultiloop BuildMultiloopFromCurves(const std::vector<std::vector<Vector2>>& inputCurves)
{
    DrawnMultiloop result;

    // The grid layout reflects every diagram horizontally (to enforce the
    // anticlockwise convention). Reflecting the drawn coordinates here first
    // cancels that, so the rendered multiloop matches what the user drew.
    std::vector<std::vector<Vector2>> curves = inputCurves;
    for (std::vector<Vector2>& curve : curves) {
        for (Vector2& p : curve) {
            p.x = -p.x;
        }
    }

    // Index segments globally, remembering which curve and which segment.
    struct SegRef { int curve; int seg; };
    std::vector<SegRef> segRefs;
    for (size_t c = 0; c < curves.size(); ++c) {
        for (size_t i = 0; i < curves[c].size(); ++i) {
            segRefs.push_back({static_cast<int>(c), static_cast<int>(i)});
        }
    }
    const auto segStart = [&](const SegRef& s) { return curves[s.curve][s.seg]; };
    const auto segEnd = [&](const SegRef& s) {
        const std::vector<Vector2>& curve = curves[s.curve];
        return curve[(s.seg + 1) % curve.size()];
    };

    // Find all crossings; record their parameter along each segment.
    std::vector<Crossing> crossings;
    // Per global segment, the crossings on it as (param, crossingId).
    std::vector<std::vector<std::pair<float, int>>> onSeg(segRefs.size());
    for (size_t i = 0; i < segRefs.size(); ++i) {
        for (size_t j = i + 1; j < segRefs.size(); ++j) {
            // Skip adjacent segments of the same curve (they share a vertex).
            if (segRefs[i].curve == segRefs[j].curve) {
                const int m = static_cast<int>(curves[segRefs[i].curve].size());
                const int a = segRefs[i].seg;
                const int b = segRefs[j].seg;
                const int diff = std::abs(a - b);
                if (diff == 1 || diff == m - 1) {
                    continue;
                }
            }
            float t = 0.0f;
            float u = 0.0f;
            if (CrossParams(segStart(segRefs[i]), segEnd(segRefs[i]),
                            segStart(segRefs[j]), segEnd(segRefs[j]), &t, &u)) {
                const int id = static_cast<int>(crossings.size());
                const Vector2 sa = segStart(segRefs[i]);
                const Vector2 ea = segEnd(segRefs[i]);
                crossings.push_back({Vector2{sa.x + (ea.x - sa.x) * t, sa.y + (ea.y - sa.y) * t}});
                onSeg[i].push_back({t, id});
                onSeg[j].push_back({u, id});
            }
        }
    }

    if (crossings.empty()) {
        if (curves.size() == 1) {
            result.isUnknot = true;
            return result;
        }
        result.error = "The multiloop is not connected.";
        return result;
    }

    for (std::vector<std::pair<float, int>>& list : onSeg) {
        std::sort(list.begin(), list.end());
    }

    // Walk each curve, collecting the ordered crossing appearances.
    std::vector<std::vector<Appearance>> walks(curves.size());
    int globalSeg = 0;
    for (size_t c = 0; c < curves.size(); ++c) {
        const std::vector<Vector2>& curve = curves[c];
        for (size_t i = 0; i < curve.size(); ++i, ++globalSeg) {
            const Vector2 a = curve[i];
            const Vector2 b = curve[(i + 1) % curve.size()];
            const float len = std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
            const Vector2 dir = (len > 1e-9f) ? Vector2{(b.x - a.x) / len, (b.y - a.y) / len}
                                              : Vector2{1.0f, 0.0f};
            for (const std::pair<float, int>& cr : onSeg[globalSeg]) {
                Appearance app;
                app.crossing = cr.second;
                app.segIndex = static_cast<int>(i);
                app.param = cr.first;
                app.dir = dir;
                walks[c].push_back(app);
            }
        }
        if (walks[c].empty() && curves.size() > 1) {
            result.error = "The multiloop is not connected.";
            return result;
        }
    }

    // Build edges (arcs between consecutive crossings) and darts. The dart at
    // the tail of edge k is -k (pointing forward along travel); the dart at
    // the head is +k (pointing back). epsilon swaps them.
    std::vector<std::vector<std::pair<float, int>>> dartsAt(crossings.size());
    std::vector<int> edgeTail(1, 0);
    std::vector<int> edgeHead(1, 0);
    std::vector<std::vector<Vector2>> edgePolyline(1);
    int edgeCount = 0;
    for (size_t c = 0; c < curves.size(); ++c) {
        const std::vector<Appearance>& walk = walks[c];
        const std::vector<Vector2>& curve = curves[c];
        const int m = static_cast<int>(curve.size());
        const size_t len = walk.size();
        for (size_t i = 0; i < len; ++i) {
            const Appearance& a = walk[i];
            const Appearance& b = walk[(i + 1) % len];
            const int k = ++edgeCount;

            // Tail dart at a: forward travel direction.
            dartsAt[a.crossing].push_back({std::atan2(a.dir.y, a.dir.x), -k});
            // Head dart at b: reverse of travel direction (points back).
            dartsAt[b.crossing].push_back({std::atan2(-b.dir.y, -b.dir.x), k});

            edgeTail.push_back(a.crossing);
            edgeHead.push_back(b.crossing);

            // Arc polyline from crossing a to crossing b (with bend points).
            std::vector<Vector2> poly;
            poly.push_back(crossings[a.crossing].pos);
            const int steps = ((b.segIndex - a.segIndex) % m + m) % m;
            int seg = a.segIndex;
            for (int s = 0; s < steps; ++s) {
                seg = (seg + 1) % m;
                poly.push_back(curve[seg]);
            }
            poly.push_back(crossings[b.crossing].pos);
            edgePolyline.push_back(std::move(poly));
        }
    }

    // sigma: at each crossing, the darts in counterclockwise angular order.
    std::vector<std::vector<int>> cycles;
    for (std::vector<std::pair<float, int>>& darts : dartsAt) {
        std::sort(darts.begin(), darts.end());
        std::vector<int> cycle;
        for (const std::pair<float, int>& d : darts) {
            cycle.push_back(d.second);
        }
        cycles.push_back(std::move(cycle));
    }

    try {
        result.sigma = Permutation::FromCycles(cycles);
    }
    catch (const std::exception& e) {
        result.error = std::string("Could not build the permutation: ") + e.what();
        return result;
    }
    const std::string message = ValidateMultiloopSigma(result.sigma);
    if (!message.empty()) {
        result.error = message;
        return result;
    }

    // Pick the unbounded face by signed area. The face traced from a dart
    // visits the directed arcs of its boundary; the outer face is the unique
    // one whose orientation is opposite to all the (interior) others.
    const Multiloop loop(result.sigma);
    int bestFace = 0;
    {
        std::vector<double> areas(loop.faces.size(), 0.0);
        for (size_t f = 0; f < loop.faces.size(); ++f) {
            std::vector<Vector2> poly;
            for (const int h : loop.faces[f]) {
                const std::vector<Vector2>& arc = edgePolyline[std::abs(h)];
                if (h > 0) {
                    for (size_t i = arc.size() - 1; i > 0; --i) {
                        poly.push_back(arc[i]);
                    }
                }
                else {
                    for (size_t i = 0; i + 1 < arc.size(); ++i) {
                        poly.push_back(arc[i]);
                    }
                }
            }
            double area = 0.0;
            for (size_t i = 0; i < poly.size(); ++i) {
                const Vector2& p = poly[i];
                const Vector2& q = poly[(i + 1) % poly.size()];
                area += static_cast<double>(p.x) * q.y - static_cast<double>(q.x) * p.y;
            }
            areas[f] = 0.5 * area;
        }
        int numPos = 0;
        int numNeg = 0;
        for (const double a : areas) {
            if (a > 1e-3) {
                ++numPos;
            }
            else if (a < -1e-3) {
                ++numNeg;
            }
        }
        if (numNeg == 1) {
            for (size_t f = 0; f < areas.size(); ++f) {
                if (areas[f] < -1e-3) {
                    bestFace = static_cast<int>(f);
                }
            }
        }
        else if (numPos == 1) {
            for (size_t f = 0; f < areas.size(); ++f) {
                if (areas[f] > 1e-3) {
                    bestFace = static_cast<int>(f);
                }
            }
        }
        else {
            double best = -1.0;
            for (size_t f = 0; f < areas.size(); ++f) {
                if (std::fabs(areas[f]) > best) {
                    best = std::fabs(areas[f]);
                    bestFace = static_cast<int>(f);
                }
            }
        }
    }
    result.exteriorHalfEdge = loop.faces[bestFace].front();
    return result;
}
