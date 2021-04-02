#ifndef __DUKHART_SDF_HLSL__
#define __DUKHART_SDF_HLSL__

#include "Math.cginc"

interface iSDF
{
    float getDist(float3 p);
};

/// SHAPES 2D
// Circle - exact
float sdCircle(float2 p, float r)
{
    return length(p) - r;
}
// Rounded Box- exact
float sdRoundedBox(float2 p, float2 b, float4 r)
{
    r.xy = (p.x > 0.0) ? r.xy : r.zw;
    r.x = (p.y > 0.0) ? r.x : r.y;
    float2 q = abs(p) - b + r.x;
    return min(max(q.x, q.y), 0.0) + length(max(q, 0.0)) - r.x;
}
// Box - exact
float sdBox(float2 p, float2 b)
{
    float2 d = abs(p) - b;
    return length(max(d, 0.0)) + min(max(d.x, d.y), 0.0);
}
// Oriented Box- exact
float sdOrientedBox(float2 p, float2 a, float2 b, float th)
{
    float l = length(b - a);
    float2 d = (b - a) / l;
    float2 q = (p - (a + b) * 0.5);
    q = mul(float2x2(d.x, -d.y, d.y, d.x), q);
    q = abs(q) - float2(l, th) * 0.5;
    return length(max(q, 0.0)) + min(max(q.x, q.y), 0.0);
}
// Segment - exact
float sdSegment(float2 p, float2 a, float2 b)
{
    float2 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    return length(pa - ba * h);
}
// Rhombus - exact
float sdRhombus(float2 p, float2 b)
{
    float2 q = abs(p);
    float h = clamp((-2.0 * ndot(q, b) + ndot(b, b)) / dot(b, b), -1.0, 1.0);
    float d = length(q - 0.5 * b * float2(1.0 - h, 1.0 + h));
    return d * sign(q.x * b.y + q.y * b.x - b.x * b.y);
}
//Isosceles Trapezoid- exact
float sdTrapezoid(float2 p, float r1, float r2, float he)
{
    float2 k1 = float2(r2, he);
    float2 k2 = float2(r2 - r1, 2.0 * he);
    p.x = abs(p.x);
    float2 ca = float2(p.x - min(p.x, (p.y < 0.0) ? r1 : r2), abs(p.y) - he);
    float2 cb = p - k1 + k2 * clamp(dot(k1 - p, k2) / dot2(k2), 0.0, 1.0);
    float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
    return s * sqrt(min(dot2(ca), dot2(cb)));
}
//Parallelogram - exact
float sdParallelogram(float2 p, float wi, float he, float sk)
{
    float2 e = float2(sk, he);
    p = (p.y < 0.0) ? -p : p;
    float2 w = p - e;
    w.x -= clamp(w.x, -wi, wi);
    float2 d = float2(dot(w, w), -w.y);
    float s = p.x * e.y - p.y * e.x;
    p = (s < 0.0) ? -p : p;
    float2 v = p - float2(wi, 0);
    v -= e * clamp(dot(v, e) / dot(e, e), -1.0, 1.0);
    d = min(d, float2(dot(v, v), wi * he - abs(s)));
    return sqrt(d.x) * sign(-d.y);
}
//Equilateral Triangle- exact
float sdEquilateralTriangle(float2 p)
{
    const float k = sqrt(3.0);
    p.x = abs(p.x) - 1.0;
    p.y = p.y + 1.0 / k;
    if (p.x + k * p.y > 0.0)
        p = float2(p.x - k * p.y, -k * p.x - p.y) / 2.0;
    p.x -= clamp(p.x, -2.0, 0.0);
    return -length(p) * sign(p.y);
}
// Isosceles Triangle- exact
float sdTriangleIsosceles(float2 p, float2 q)
{
    p.x = abs(p.x);
    float2 a = p - q * clamp(dot(p, q) / dot(q, q), 0.0, 1.0);
    float2 b = p - q * float2(clamp(p.x / q.x, 0.0, 1.0), 1.0);
    float s = -sign(q.y);
    float2 d = min(float2(dot(a, a), s * (p.x * q.y - p.y * q.x)),
                  float2(dot(b, b), s * (p.y - q.y)));
    return -sqrt(d.x) * sign(d.y);
}
// Triangle - exact
float sdTriangle(float2 p, float2 p0, float2 p1, float2 p2)
{
    float2 e0 = p1 - p0, e1 = p2 - p1, e2 = p0 - p2;
    float2 v0 = p - p0, v1 = p - p1, v2 = p - p2;
    float2 pq0 = v0 - e0 * clamp(dot(v0, e0) / dot(e0, e0), 0.0, 1.0);
    float2 pq1 = v1 - e1 * clamp(dot(v1, e1) / dot(e1, e1), 0.0, 1.0);
    float2 pq2 = v2 - e2 * clamp(dot(v2, e2) / dot(e2, e2), 0.0, 1.0);
    float s = sign(e0.x * e2.y - e0.y * e2.x);
    float2 d = min(min(float2(dot(pq0, pq0), s * (v0.x * e0.y - v0.y * e0.x)),
                     float2(dot(pq1, pq1), s * (v1.x * e1.y - v1.y * e1.x))),
                     float2(dot(pq2, pq2), s * (v2.x * e2.y - v2.y * e2.x)));
    return -sqrt(d.x) * sign(d.y);
}
// Uneven Capsule- exact
float sdUnevenCapsule(float2 p, float r1, float r2, float h)
{
    p.x = abs(p.x);
    float b = (r1 - r2) / h;
    float a = sqrt(1.0 - b * b);
    float k = dot(p, float2(-b, a));
    if (k < 0.0)
        return length(p) - r1;
    if (k > a * h)
        return length(p - float2(0.0, h)) - r2;
    return dot(p, float2(a, b)) - r1;
}
// Regular Pentagon- exact
float sdPentagon(float2 p, float r)
{
    const float3 k = float3(0.809016994, 0.587785252, 0.726542528);
    p.x = abs(p.x);
    p -= 2.0 * min(dot(float2(-k.x, k.y), p), 0.0) * float2(-k.x, k.y);
    p -= 2.0 * min(dot(float2(k.x, k.y), p), 0.0) * float2(k.x, k.y);
    p -= float2(clamp(p.x, -r * k.z, r * k.z), r);
    return length(p) * sign(p.y);
}
// Regular Hexagon- exact
float sdHexagon(float2 p, float r)
{
    const float3 k = float3(-0.866025404, 0.5, 0.577350269);
    p = abs(p);
    p -= 2.0 * min(dot(k.xy, p), 0.0) * k.xy;
    p -= float2(clamp(p.x, -k.z * r, k.z * r), r);
    return length(p) * sign(p.y);
}
// Regular Octogon- exact
float sdOctogon(float2 p, float r)
{
    const float3 k = float3(-0.9238795325, 0.3826834323, 0.4142135623);
    p = abs(p);
    p -= 2.0 * min(dot(float2(k.x, k.y), p), 0.0) * float2(k.x, k.y);
    p -= 2.0 * min(dot(float2(-k.x, k.y), p), 0.0) * float2(-k.x, k.y);
    p -= float2(clamp(p.x, -k.z * r, k.z * r), r);
    return length(p) * sign(p.y);
}
// Hexagram - exact
float sdHexagram(float2 p, float r)
{
    const float4 k = float4(-0.5, 0.8660254038, 0.5773502692, 1.7320508076);
    p = abs(p);
    p -= 2.0 * min(dot(k.xy, p), 0.0) * k.xy;
    p -= 2.0 * min(dot(k.yx, p), 0.0) * k.yx;
    p -= float2(clamp(p.x, r * k.z, r * k.w), r);
    return length(p) * sign(p.y);
}
// Star5 - exact
float sdStar5(float2 p, float r, float rf)
{
    const float2 k1 = float2(0.809016994375, -0.587785252292);
    const float2 k2 = float2(-k1.x, k1.y);
    p.x = abs(p.x);
    p -= 2.0 * max(dot(k1, p), 0.0) * k1;
    p -= 2.0 * max(dot(k2, p), 0.0) * k2;
    p.x = abs(p.x);
    p.y -= r;
    float2 ba = rf * float2(-k1.y, k1.x) - float2(0, 1);
    float h = clamp(dot(p, ba) / dot(ba, ba), 0.0, r);
    return length(p - ba * h) * sign(p.y * ba.x - p.x * ba.y);
}
// Regular Star- exact
float sdStar(float2 p, float r, int n, float m)
{
    // next 4 lines can be precomputed for a given shape
    float an = 3.141593 / float(n);
    float en = 3.141593 / m; // m is between 2 and n
    float2 acs = float2(cos(an), sin(an));
    float2 ecs = float2(cos(en), sin(en)); // ecs=vec2(0,1) for regular polygon

    float bn = mod(atan2(p.x, p.y), 2.0 * an) - an;
    p = length(p) * float2(cos(bn), abs(sin(bn)));
    p -= r * acs;
    p += ecs * clamp(-dot(p, ecs), 0.0, r * acs.y / ecs.y);
    return length(p) * sign(p.x);
}

// Pie - exact
float sdPie(float2 p, float2 c, float r)
{
    p.x = abs(p.x);
    float l = length(p) - r;
    float m = length(p - c * clamp(dot(p, c), 0.0, r)); // c=sin/cos of aperture
    return max(l, m * sign(c.y * p.x - c.x * p.y));
}
// Arc - exact
float sdArc(float2 p, float2 sca, float2 scb, float ra, float rb)
{
    p = mul(float2x2(sca.x, sca.y, -sca.y, sca.x), p);
    p.x = abs(p.x);
    float k = (scb.y * p.x > scb.x * p.y) ? dot(p, scb) : length(p);
    return sqrt(dot(p, p) + ra * ra - 2.0 * ra * k) - rb;
}
// Horseshoe - exact
float sdHorseshoe(float2 p, float2 c, float r, float2 w)
{
    p.x = abs(p.x);
    float l = length(p);
    p = mul(float2x2(-c.x, c.y, c.y, c.x), p);
    p = float2((p.y > 0.0) ? p.x : l * sign(-c.x),
             (p.x > 0.0) ? p.y : l);
    p = float2(p.x, abs(p.y - r)) - w;
    return length(max(p, 0.0)) + min(0.0, max(p.x, p.y));
}
// Vesica - exact
float sdVesica(float2 p, float r, float d)
{
    p = abs(p);
    float b = sqrt(r * r - d * d);
    return ((p.y - b) * d > p.x * b) ? length(p - float2(0.0, b))
                             : length(p - float2(-d, 0.0)) - r;
}
// Moon - exact
float sdMoon(float2 p, float d, float ra, float rb)
{
    p.y = abs(p.y);
    float a = (ra * ra - rb * rb + d * d) / (2.0 * d);
    float b = sqrt(max(ra * ra - a * a, 0.0));
    if (d * (p.x * b - p.y * a) > d * d * max(b - p.y, 0.0))
        return length(p - float2(a, b));
    return max((length(p) - ra),
               -(length(p - float2(d, 0)) - rb));
}
// Simple Egg- exact
float sdEgg(float2 p, float ra, float rb)
{
    const float k = sqrt(3.0);
    p.x = abs(p.x);
    float r = ra - rb;
    return ((p.y < 0.0) ? length(float2(p.x, p.y)) - r :
            (k * (p.x + r) < p.y) ? length(float2(p.x, p.y - k * r)) :
                              length(float2(p.x + r, p.y)) - 2.0 * r) - rb;
}
// Heart - exact
float sdHeart(float2 p)
{
    p.x = abs(p.x);

    if (p.y + p.x > 1.0)
        return sqrt(dot2(p - float2(0.25, 0.75))) - sqrt(2.0) / 4.0;
    return sqrt(min(dot2(p - float2(0.00, 1.00)),
                    dot2(p - 0.5 * max(p.x + p.y, 0.0)))) * sign(p.x - p.y);
}
// Cross - exact exterior, boundinterior
float sdCross(float2 p, float2 b, float r)
{
    p = abs(p);
    p = (p.y > p.x) ? p.yx : p.xy;
    float2 q = p - b;
    float k = max(q.y, q.x);
    float2 w = (k > 0.0) ? q : float2(b.y - p.x, -k);
    return sign(k) * length(max(w, 0.0)) + r;
}
// Rounded X- exact
float sdRoundedX(float2 p, float w, float r)
{
    p = abs(p);
    return length(p - min(p.x + p.y, w) * 0.5) - r;
}
// Polygon - exact
float sdPolygon(float2 p, float2 size, float sides)
{
    float d = dot(p - size, p - size);
    float s = 1.0;
    float2 thisP, lastP = p + size;
    for (int i = 0, j = sides - 1; i < sides; j = i, i++)
    {
        thisP = p + mul(size, rotMat((360 / sides)*i));
        float2 e = lastP - thisP;
        float2 w = p - thisP;
        float2 b = w - e * clamp(dot(w, e) / dot(e, e), 0.0, 1.0);
        d = min(d, dot(b, b));
        float3 c = float3(p.y >= thisP.y, p.y < lastP.y, e.x * w.y > e.y * w.x);
        if (all(c) || all(not(c)))
            s *= -1.0;
        
        lastP = thisP;
    }
    return s * sqrt(d);
}
// Ellipse - exact
float sdEllipse(float2 p, float2 ab)
{
    p = abs(p);
    if (p.x > p.y)
    {
        p = p.yx;
        ab = ab.yx;
    }
    float l = ab.y * ab.y - ab.x * ab.x;
    float m = ab.x * p.x / l;
    float m2 = m * m;
    float n = ab.y * p.y / l;
    float n2 = n * n;
    float c = (m2 + n2 - 1.0) / 3.0;
    float c3 = c * c * c;
    float q = c3 + m2 * n2 * 2.0;
    float d = c3 + m2 * n2;
    float g = m + m * n2;
    float co;
    if (d < 0.0)
    {
        float h = acos(q / c3) / 3.0;
        float s = cos(h);
        float t = sin(h) * sqrt(3.0);
        float rx = sqrt(-c * (s + t + 2.0) + m2);
        float ry = sqrt(-c * (s - t + 2.0) + m2);
        co = (ry + sign(l) * rx + abs(g) / (rx * ry) - m) / 2.0;
    }
    else
    {
        float h = 2.0 * m * n * sqrt(d);
        float s = sign(q + h) * pow(abs(q + h), 1.0 / 3.0);
        float u = sign(q - h) * pow(abs(q - h), 1.0 / 3.0);
        float rx = -s - u - c * 4.0 + 2.0 * m2;
        float ry = (s - u) * sqrt(3.0);
        float rm = sqrt(rx * rx + ry * ry);
        co = (ry / sqrt(rm - rx) + 2.0 * g / rm - m) / 2.0;
    }
    float2 r = ab * float2(co, sqrt(1.0 - co * co));
    return length(r - p) * sign(p.y - r.y);
}
// Parabola - exact

float sdParabola(float2 pos, float k)
{
    pos.x = abs(pos.x);
    float ik = 1.0 / k;
    float p = ik * (pos.y - 0.5 * ik) / 3.0;
    float q = 0.25 * ik * ik * pos.x;
    float h = q * q - p * p * p;
    float r = sqrt(abs(h));
    float x = (h > 0.0) ?
        pow(q + r, 1.0 / 3.0) - pow(abs(q - r), 1.0 / 3.0) * sign(r - q) :
        2.0 * cos(atan2(r, q) / 3.0) * sqrt(p);
    return length(pos - float2(x, k * x * x)) * sign(pos.x - x);
}

// Parabola Segment- exact
float sdParabola(float2 pos, float wi, float he)
{
    pos.x = abs(pos.x);
    float ik = wi * wi / he;
    float p = ik * (he - pos.y - 0.5 * ik) / 3.0;
    float q = pos.x * ik * ik * 0.25;
    float h = q * q - p * p * p;
    float r = sqrt(abs(h));
    float x = (h > 0.0) ?
        pow(q + r, 1.0 / 3.0) - pow(abs(q - r), 1.0 / 3.0) * sign(r - q) :
        2.0 * cos(atan(r / q) / 3.0) * sqrt(p);
    x = min(x, wi);
    return length(pos - float2(x, he - x * x / ik)) *
           sign(ik * (pos.y - he) + pos.x * pos.x);
}
// Quadratic Bezier- exact
float sdBezier(float2 pos, float2 A, float2 B, float2 C)
{
    float2 a = B - A;
    float2 b = A - 2.0 * B + C;
    float2 c = a * 2.0;
    float2 d = A - pos;
    float kk = 1.0 / dot(b, b);
    float kx = kk * dot(a, b);
    float ky = kk * (2.0 * dot(a, a) + dot(d, b)) / 3.0;
    float kz = kk * dot(d, a);
    float res = 0.0;
    float p = ky - kx * kx;
    float p3 = p * p * p;
    float q = kx * (2.0 * kx * kx - 3.0 * ky) + kz;
    float h = q * q + 4.0 * p3;
    if (h >= 0.0)
    {
        h = sqrt(h);
        float2 x = (float2(h, -h) - q) / 2.0;
        float2 uv = sign(x) * pow(abs(x), float2(1/3,1/3));
        float t = clamp(uv.x + uv.y - kx, 0.0, 1.0);
        res = dot2(d + (c + b * t) * t);
    }
    else
    {
        float z = sqrt(-p);
        float v = acos(q / (p * z * 2.0)) / 3.0;
        float m = cos(v);
        float n = sin(v) * 1.732050808;
        float3 t = clamp(float3(m + m, -n - m, n - m) * z - kx, 0.0, 1.0);
        res = min(dot2(d + (c + b * t.x) * t.x),
                   dot2(d + (c + b * t.y) * t.y));
        // the third root cannot be the closest
        // res = min(res,dot2(d+(c+b*t.z)*t.z));
    }
    return sqrt(res);
}
/// SHAPES 3D
// Sphere - exact
struct SphereSDF : iSDF
{
    float radius;
    // Sphere - exact
    static float sphereSDF(float3 p, float radius)
    {
        return length(p) - radius;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return sphereSDF(p, radius);
    }
};
// Ellipsoid - bound  - *NOT exact*
struct EllipsoidSDF : iSDF
{
    float3 radii;
    // Ellipsoid - bound  - *NOT exact*
    static float ellipsoidBndSDF(float3 p, float3 radii)
    {
        float k0 = length(p / radii);
        float k1 = length(p / (radii * radii));
        return k0 * (k0 - 1.0) / k1;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return ellipsoidBndSDF(p, radii);
    }
};

// Box - exact
struct BoxSDF : iSDF
{
    float3 size;
    // Box - exact
    static float boxSDF(float3 p, float3 s)
    {
        float3 q = abs(p) - s;
        return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return boxSDF(p, size);
    }
};
// BoxRound - exact
struct BoxRoundSDF : iSDF
{
    float3 size;
    float bevel;
    // BoxRound - exact
    static float boxRoundSDF(float3 p, float3 size, float bevel)
    {
        float3 q = abs(p) - size;
        return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0) - bevel;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return boxRoundSDF(p, size, bevel);
    }
};
// BoxFrame - exact
struct BoxFrameSDF : iSDF
{
    float3 size;
    float frame;
    // BoxFrame - exact
    static float boxFrameSDF(float3 p, float3 s, float f)
    {
        p = abs(p) - s;
        float3 q = abs(p + f) - f;
        return min(min(
    length(max(float3(p.x, q.y, q.z), 0.0)) + min(max(p.x, max(q.y, q.z)), 0.0),
    length(max(float3(q.x, p.y, q.z), 0.0)) + min(max(q.x, max(p.y, q.z)), 0.0)),
    length(max(float3(q.x, q.y, p.z), 0.0)) + min(max(q.x, max(q.y, p.z)), 0.0));
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return boxFrameSDF(p,size,frame);
    }
};

// Torus - exact
class TorusSDF : iSDF
{
    float2 radii;
    // Torus - exact
    static float torusSDF(float3 p, float2 r)
    {
        float2 q = float2(length(p.xz) - r.x, p.y);
        return length(q) - r.y;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return torusSDF(p, radii);
    }
};
// TorusCapped - exact
struct TorusCappedSDF : iSDF
{
    float2 radii;
    float2 cap;
    // TorusCapped - exact
    static float torusCappedSDF(float3 p, float2 sc, float2 r)
    {
        p.x = abs(p.x);
        float k = (sc.y * p.x > sc.x * p.y) ? dot(p.xy, sc) : length(p.xy);
        return sqrt(dot(p, p) + r.x * r.x - 2.0 * r.x * k) - r.y;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return torusCappedSDF(p, cap, radii);
    }
};
// Link - exact
struct LinkSDF : iSDF
{
    float2 radii;
    float length;
    // Link - exact
    static float linkSDF(float3 p, float l, float2 r)
    {
        float3 q = float3(p.x, max(abs(p.y) - l, 0.0), p.z);
        return length(float2(length(q.xy) - r.x, q.z)) - r.y;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return linkSDF(p, length, radii);
    }
};

// Cylinder Capped - exact
struct CylinderSDF : iSDF
{
    float height, radius;
    // Cylinder Capped - exact
    static float cylinderVerticalSDF(float3 p, float h, float r)
    {
        float2 d = abs(float2(length(p.xz), p.y)) - float2(h, r);
        return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return cylinderVerticalSDF(p, height, radius);
    }
};
// Cylinder Capped - exact
struct CylinderAdvSDF : iSDF
{
    float radius;
    float3 offsetA, offsetB;
    // Cylinder Capped - exact
    static float cylinderAdvSDF(float3 p, float3 a, float3 b, float r)
    {
        float3 ba = b - a;
        float3 pa = p - a;
        float baba = dot(ba, ba);
        float paba = dot(pa, ba);
        float x = length(pa * baba - ba * paba) - r * baba;
        float y = abs(paba - baba * 0.5) - baba * 0.5;
        float x2 = x * x;
        float y2 = y * y * baba;
        float d = (max(x, y) < 0.0) ? -min(x2, y2) : (((x > 0.0) ? x2 : 0.0) + ((y > 0.0) ? y2 : 0.0));
        return sign(d) * sqrt(abs(d)) / baba;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return cylinderAdvSDF(p, offsetA, offsetB, radius);
    }
};
// Cylinder Infinite  - exact
struct CylinderInfiniteSDF : iSDF
{
    float radius;
    float2 pos;
    // Cylinder Infinite  - exact
    static float cylinderInfiniteSDF(float3 p, float r, float2 pos = float2(0,0))
    {
        return length(p.xz - pos) - r;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return cylinderInfiniteSDF(p, radius, pos);
    }
};
// Cylinder Rounded - exact
struct CylinderRoundedSDF : iSDF
{
    float radius, height, cap;
    // Cylinder Rounded - exact
    static float cylinderRoundedSDF(float3 p, float r, float c, float h)
    {
        float2 d = float2(length(p.xz) - 2.0 * r + c, abs(p.y) - h);
        return min(max(d.x, d.y), 0.0) + length(max(d, 0.0)) - c;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return cylinderRoundedSDF(p, radius, cap, height);
    }
};

// Capsule Vertical  - exact
struct CapsuleSDF : iSDF
{
    float height, radius;
    // Capsule Vertical  - exact
    static float capsuleVerticalSDF(float3 p, float h, float r)
    {
        p.y -= clamp(p.y, 0.0, h);
        return length(p) - r;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return capsuleVerticalSDF(p, height, radius);
    }
};
// Capsule - exact
struct CapsuleAdvSDF : iSDF
{
    float radius;
    float3 offsetA, offsetB;
    // Capsule - exact
    static float capsuleSDF(float3 p, float3 a, float3 b, float r)
    {
        float3 pa = p - a, ba = b - a;
        float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
        return length(pa - ba * h) - r;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return capsuleSDF(p, offsetA, offsetB, radius);
    }
};

// Cone - exact
struct ConeSDF : iSDF
{
    float height;
    float2 radii;
    
    // Cone - exact
    static float coneSDF(in float3 p, float2 c, float h)
    {
        // c is the sin/cos of the angle, h is height
        // Alternatively pass q instead of (c,h),
        // which is the point at the base in 2D
        float2 q = h * float2(c.x / c.y, -1.0);
        float2 w = float2(length(p.xz), p.y);
        float2 a = w - q * clamp(dot(w, q) / dot(q, q), 0.0, 1.0);
        float2 b = w - q * float2(clamp(w.x / q.x, 0.0, 1.0), 1.0);
        float k = sign(q.y);
        float d = min(dot(a, a), dot(b, b));
        float s = max(k * (w.x * q.y - w.y * q.x), k * (w.y - q.y));
        return sqrt(d) * sign(s);
    }
    
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return coneSDF(p, radii, height);
    }
};
// Cone Bound - *NOT exact*
struct ConeBoundSDF : iSDF
{
    float height;
    float2 radii;
    // Cone Bound - *NOT exact*
    static float coneBoundSDF(float3 p, float2 c, float h)
    {
        float q = length(p.xz);
        return max(dot(c.xy, float2(q, p.y)), -h - p.y);
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return coneBoundSDF(p, radii, height);
    }
};
// Cone Infinite - exact
struct ConeInfiniteSDF : iSDF
{
    float2 radii;
    // Cone Infinite - exact
    static float coneInfiniteSDF(float3 p, float2 c)
    {
        // c is the sin/cos of the angle
        float2 q = float2(length(p.xz), -p.y);
        float d = length(q - c * max(dot(q, c), 0.0));
        return d * ((q.x * c.y - q.y * c.x < 0.0) ? -1.0 : 1.0);
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return coneInfiniteSDF(p, radii);
    }
};
// Cone Capped - exact
struct ConeCappedSDF : iSDF
{
    float height;
    float2 radii;
    // Cone Capped - exact
    static float coneCappedSDF(float3 p, float h, float2 r)
    {
        float2 q = float2(length(p.xz), p.y);
        float2 k1 = float2(r.y, h);
        float2 k2 = float2(r.y - r.x, 2.0 * h);
        float2 ca = float2(q.x - min(q.x, (q.y < 0.0) ? r.x : r.y), abs(q.y) - h);
        float2 cb = q - k1 + k2 * clamp(dot(k1 - q, k2) / dot2(k2), 0.0, 1.0);
        float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
        return s * sqrt(min(dot2(ca), dot2(cb)));
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return coneCappedSDF(p, height, radii);
    }
};
//Capped Cone - exact 
struct ConeCappedAdvSDF : iSDF
{
    float3 posA, posB;
    float2 radii;
    //Capped Cone - exact 
    static float coneCappedAdvSDF(float3 p, float3 a, float3 b, float2 r)
    {
        float rba = r.y - r.x;
        float baba = dot(b - a, b - a);
        float papa = dot(p - a, p - a);
        float paba = dot(p - a, b - a) / baba;
        float x = sqrt(papa - paba * paba * baba);
        float cax = max(0.0, x - ((paba < 0.5) ? r.x : r.y));
        float cay = abs(paba - 0.5) - 0.5;
        float k = rba * rba + baba;
        float f = clamp((rba * (x - r.x) + paba * baba) / k, 0.0, 1.0);
        float cbx = x - r.x - f * rba;
        float cby = paba - f;
        float s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
        return s * sqrt(min(cax * cax + cay * cay * baba,
                       cbx * cbx + cby * cby * baba));
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return coneCappedAdvSDF(p, posA, posB, radii);
    }
};
// Round Cone - exact
struct ConeRoundSDF : iSDF
{
    float height;
    float2 radii;
    
    // Round Cone - exact
    static float coneRoundSDF(float3 p, float2 r, float h)
    {
        float2 q = float2(length(p.xz), p.y);
    
        float b = (r.x - r.y) / h;
        float a = sqrt(1.0 - b * b);
        float k = dot(q, float2(-b, a));
    
        if (k < 0.0)
            return length(q) - r.x;
        if (k > a * h)
            return length(q - float2(0.0, h)) - r.y;
        
        return dot(q, float2(a, b)) - r.x;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return coneRoundSDF(p, radii, height);
    }
};
// Round Cone - exact
struct ConeRoundAdvSDF : iSDF
{
    float2 radii;
    float3 posA, posB;
    // Round Cone - exact
    static float coneRoundAdvSDF(float3 p, float3 a, float3 b, float2 r)
    {
        
        // sampling independent computations (only depend on shape)
        float3 ba = b - a;
        float l2 = dot(ba, ba);
        float rr = r.x - r.y;
        float a2 = l2 - rr * rr;
        float il2 = 1.0 / l2;
    
        // sampling dependant computations
        float3 pa = p - a;
        float y = dot(pa, ba);
        float z = y - l2;
        float x2 = dot2(pa * l2 - ba * y);
        float y2 = y * y * l2;
        float z2 = z * z * l2;

        // single square root!
        float k = sign(rr) * rr * rr * x2;
        if (sign(z) * a2 * z2 > k)
            return sqrt(x2 + z2) * il2 - r.y;
        if (sign(y) * a2 * y2 < k)
            return sqrt(x2 + y2) * il2 - r.x;
        return (sqrt(x2 * a2 * il2) + y * rr) * il2 - r.x;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return coneRoundAdvSDF(p, posA, posB, radii);
    }
};
// Solid Angle - exact
struct solidAngleSDF : iSDF
{
    float size;
    float2 angle;
    // Solid Angle - exact
    static float solidAngleSDF(float3 p, float2 c, float ra)
    {
    // c is the sin/cos of the angle
        float2 q = float2(length(p.xz), p.y);
        float l = length(q) - ra;
        float m = length(q - c * clamp(dot(q, c), 0.0, ra));
        return max(l, m * sign(c.y * q.x - c.x * q.y));
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return solidAngleSDF(p, angle, size);
    }
};

// Plane - exact
struct PlaneSDF : iSDF
{
    float h;
    float3 n;
    // Plane - exact
    static float planeSDF(float3 p, float3 n, float h)
    { 
        normalize(n); // n must be normalized
        return dot(p, n) + h;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return planeSDF(p, n, h);
    }
};
// Triangle - exact
struct TriangleSDF : iSDF
{
    float3 a, b, c;
    // Triangle - exact
    static float triangleSDF(float3 p, float3 a, float3 b, float3 c)
    {
        float3 ba = b - a;
        float3 pa = p - a;
        float3 cb = c - b;
        float3 pb = p - b;
        float3 ac = a - c;
        float3 pc = p - c;
        float3 nor = cross(ba, ac);

        return sqrt(
    (sign(dot(cross(ba, nor), pa)) +
     sign(dot(cross(cb, nor), pb)) +
     sign(dot(cross(ac, nor), pc)) < 2.0)
     ?
     min(min(
     dot2(ba * clamp(dot(ba, pa) / dot2(ba), 0.0, 1.0) - pa),
     dot2(cb * clamp(dot(cb, pb) / dot2(cb), 0.0, 1.0) - pb)),
     dot2(ac * clamp(dot(ac, pc) / dot2(ac), 0.0, 1.0) - pc))
     :
     dot(nor, pa) * dot(nor, pa) / dot2(nor));
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return triangleSDF(p, a, b, c);
    }
};
// Quad - exact
struct QuadSDF : iSDF
{
    float3 a, b, c, d;
    // Quad - exact
    static float quadSDF(float3 p, float3 a, float3 b, float3 c, float3 d)
    {
        float3 ba = b - a;
        float3 pa = p - a;
        float3 cb = c - b;
        float3 pb = p - b;
        float3 dc = d - c;
        float3 pc = p - c;
        float3 ad = a - d;
        float3 pd = p - d;
        float3 nor = cross(ba, ad);

        return sqrt(
    (sign(dot(cross(ba, nor), pa)) +
     sign(dot(cross(cb, nor), pb)) +
     sign(dot(cross(dc, nor), pc)) +
     sign(dot(cross(ad, nor), pd)) < 3.0)
     ?
     min(min(min(
     dot2(ba * clamp(dot(ba, pa) / dot2(ba), 0.0, 1.0) - pa),
     dot2(cb * clamp(dot(cb, pb) / dot2(cb), 0.0, 1.0) - pb)),
     dot2(dc * clamp(dot(dc, pc) / dot2(dc), 0.0, 1.0) - pc)),
     dot2(ad * clamp(dot(ad, pd) / dot2(ad), 0.0, 1.0) - pd))
     :
     dot(nor, pa) * dot(nor, pa) / dot2(nor));
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return quadSDF(p, a, b, c, d);
    }
};

// Prism Tri
struct TriPrismSDF : iSDF
{
    float2 size;
    // Prism Tri
    static float triPrismSDF(float3 p, float2 h)
    {
        float3 q = abs(p);
        return max(q.z - h.y, max(q.x * 0.866025 + p.y * 0.5, -p.y) - h.x * 0.5);
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return triPrismSDF(p, size);
    }
};
// Prism Hex
struct HexPrismSDF : iSDF
{
    float2 size;
    // Prism Hex
    static float hexPrismSDF(float3 p, float2 h)
    {
        const float3 k = float3(-0.8660254, 0.5, 0.57735);
        p = abs(p);
        p.xy -= 2.0 * min(dot(k.xy, p.xy), 0.0) * k.xy;
        float2 d = float2(
       length(p.xy - float2(clamp(p.x, -k.z * h.x, k.z * h.x), h.x)) * sign(p.y - h.x),
       p.z - h.y);
        return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return hexPrismSDF(p, size);
    }
};

// Octahedron - exact
struct OctahedronSDF : iSDF
{
    float radius;
    // Octahedron - exact
    static float octahedronSDF(float3 p, float s)
    {
        p = abs(p);
        float m = p.x + p.y + p.z - s;
        float3 q;
        if (3.0 * p.x < m)
            q = p.xyz;
        else if (3.0 * p.y < m)
            q = p.yzx;
        else if (3.0 * p.z < m)
            q = p.zxy;
        else
            return m * 0.57735027;
    
        float k = clamp(0.5 * (q.z - q.y + s), 0.0, s);
        return length(float3(q.x, q.y - s + k, q.z - k));
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return octahedronSDF(p, radius);
    }
};
// Octahedron - bound (not exact)
struct OctahedronBoundSDF : iSDF
{
    float radius;
    // Octahedron - bound (not exact)
    static float octahedronBoundSDF(float3 p, float s)
    {
        p = abs(p);
        return (p.x + p.y + p.z - s) * 0.57735027;
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return octahedronBoundSDF(p, radius);
    }
};
// Pyramid - exact
struct PyramidSDF : iSDF
{
    float3 size;
    // Pyramid - exact
    static float pyramidSDF(float3 p, float3 h)
    {
        float m2 = h.y * h.y + 0.25;
    
        p.xz = abs(p.xz) - h.xz;
        p.xz = (p.z > p.x) ? p.zx : p.xz;
        p.xz -= 0.5;

        float3 q = float3(p.z, h.y * p.y - 0.5 * p.x, h.y * p.x + 0.5 * p.y);
   
        float s = max(-q.x, 0.0);
        float t = clamp((q.y - 0.5 * p.z) / (m2 + 0.25), 0.0, 1.0);
    
        float a = m2 * (q.x + s) * (q.x + s) + q.y * q.y;
        float b = m2 * (q.x + 0.5 * t) * (q.x + 0.5 * t) + (q.y - m2 * t) * (q.y - m2 * t);
    
        float d2 = min(q.y, -q.x * m2 - q.y * 0.5) > 0.0 ? 0.0 : min(a, b);
    
        return sqrt((d2 + q.z * q.z) / m2) * sign(max(q.z, -p.y));
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return pyramidSDF(p, size);
    }
};
// Rhombus - exact
struct RhombusSDF : iSDF
{
    float height, radius;
    float2 length;
    // Rhombus - exact
    static float rhombusSDF(float3 p, float2 l, float h, float r)
    {
        p = abs(p);
        float2 b = float2(l.x, l.y);
        float f = clamp((ndot(b, b - 2.0 * p.xz)) / dot(b, b), -1.0, 1.0);
        float2 q = float2(length(p.xz - 0.5 * b * float2(1.0 - f, 1.0 + f)) * sign(p.x * b.y + p.z * b.x - b.x * b.y) - r, p.y - h);
        return min(max(q.x, q.y), 0.0) + length(max(q, 0.0));
    }
    // gets the sdf objects distance
    float getDist(float3 p)
    {
        return rhombusSDF(p, length, height, radius);
    }
};


// SHAPE OPERATIONS
// boolean ops
float intersectSDF(float distA, float distB)
{
    return max(distA, distB);
}
float unionSDF(float distA, float distB)
{
    return min(distA, distB);
}
float differenceSDF(float distA, float distB, int flip = 0)
{
    return max(distA, -distB);
}
float smoothUnionSDF(float distA, float distB, float amount)
{
    if (amount == 0)
        return unionSDF(distA, distB);
    float h = clamp(0.5 + 0.5 * (distB - distA) / amount, 0.0, 1.0);
    return lerp(distB, distA, h) - amount * h * (1.0 - h);
}
float smoothSubtractionSDF(float distA, float distB, float amount)
{
    if (amount == 0)
        return differenceSDF(distA, distB);
    
    float h = clamp(0.5 - 0.5 * (distB + distA) / amount, 0.0, 1.0);
    return lerp(distB, -distA, h) + amount * h * (1.0 - h);
}
float smoothIntersectionSDF(float distA, float distB, float amount)
{
    if (amount == 0)
        return intersectSDF(distA, distB);
    float h = clamp(0.5 - 0.5 * (distB - distA) / amount, 0.0, 1.0);
    return lerp(distB, distA, h) + amount * h * (1.0 - h);
}
float sdfOpScale(float3 p, float s, iSDF primitive)
{
    return primitive.getDist(p / s) * s;
}
// alteration ops
// translate / rotate
/*
vec3 opTx(in vec3 p, in transform t, float primitive)
{
    return primitive(invert(t) * p);
}
*/
float opSymX(float3 p, iSDF primitive)
{
    p.x = abs(p.x);
    return primitive.getDist(p);
}
float opSymXZ(float3 p, iSDF primitive)
{
    p.xz = abs(p.xz);
    return primitive.getDist(p);
}
float opRep(float3 p, float3 c, in iSDF primitive)
{
    float3 q = mod(p+0.5*c,c)-0.5*c;
    return primitive.getDist(q);
}
float opRepLim(float3 p, float c, float3 l, in iSDF primitive)
{
    float3 q = p - c * clamp(round(p / c), -l, l);
    return primitive.getDist(q);
}
float opDisplace(iSDF primitive, float3 p, float3 d)
{
    
    return primitive.getDist(p+d);
}
float opTwist(iSDF primitive, float3 p)
{
    const float k = 10.0; // or some other amount
    float c = cos(k*p.y);
    float s = sin(k*p.y);
    float2x2 m = float2x2(c, -s, s, c);
    float3 q = float3(mul(m, p.xz), p.y);
    return primitive.getDist(q);
}
float opCheapBend(iSDF primitive, float3 p)
{
    const float k = 10.0; // or some other amount
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    float2x2 m = float2x2(c, -s, s, c);
    float3 q = float3(mul(m, p.xy), p.z);
    return primitive.getDist(q);
}
// should be used on one axis only
float elongate1dSDF(iSDF primitive, float3 p, float3 amount)
{
    float3 q = p - clamp(p, -amount, amount);
    return primitive.getDist(q);
}
// use when modify multple axis
float elongate3dSDF(iSDF primitive, float3 p, float3 amount)
{
    float3 q = abs(p) - amount;
    return primitive.getDist(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}
float opRound(iSDF primitive, float3 p, float rad)
{
    return primitive.getDist(p) - rad;
}
float opExtrusion(iSDF primitive, float3 p, float h)
{
    float d = primitive.getDist(p);
    
    float2 w = float2(d, abs(p.z) - h);
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}
float opRevolution(iSDF primitive, float3 p, float o)
{
    float3 q = float3(length(p.xz) - o, p.y, p.z);
    return primitive.getDist(q);

}
float onionSDF(float dist, float thickness)
{
    return abs(dist) - thickness;
}

// RAY OPERATIONS
float3 GetPoint(float3 rayOrigin, float distOrigin, float3 rayDir)
{
    // get the point on the ray
    float3 p = rayOrigin + distOrigin * rayDir;
    //return
    return p;
}
float RayMarch(float3 rayOrigin, float3 rayDir, iSDF shape,float3 offset = 0, float surfDistance = 0.01f, int maxSteps = 100, int maxDistance = 1000)
{            
    float distOrigin = 0;
    float distScurface;
    for (int i = 0; i < maxSteps; i++)
    {
        // get the ray marching point
        float3 p = GetPoint(rayOrigin, distOrigin, rayDir);
        // get the distance to the surface from the ray marching point
        distScurface = shape.getDist(p + offset);

        // operation
        
        // move our origin by surface distance
        distOrigin += distScurface;
        // if distScurface < _SurfDist we hit something, 
        // if dist Origin > maxDist we reached the end of the ray and didn't hit anything
        if (distScurface < surfDistance || abs(distOrigin) > maxDistance)
        {
            break;
        }
    }
    return distOrigin;
}
float4 RayMarch2(float3 rayOrigin, float3 rayDir, iSDF shape1, iSDF shape2, float surfDistance = 0.01f, int maxSteps = 100, int maxDistance = 1000)
{
    float distOrigin = 0;
    float distScurface, distScurface1, distScurface2;
    float3 p = 0;
    for (int i = 0; i < maxSteps; i++)
    {
        // get the ray marching point
        float3 p = GetPoint(rayOrigin, distOrigin, rayDir);
        // get the distance to the surface from the ray marching point
        distScurface1 = shape1.getDist(p);
        distScurface2 = shape2.getDist(p + float3(0.1,0,0));
        
        
        /// operation
        distScurface = smoothUnionSDF(distScurface1, distScurface2, 0.05);
        
        
        // move our origin by surface distance
        distOrigin += distScurface;
        // if distScurface < _SurfDist we hit something, 
        // if dist Origin > maxDist we reached the end of the ray and didn't hit anything
        if (distScurface < surfDistance || abs(distOrigin) > maxDistance)
        {
            break;
        }
    }
    return float4(p,distOrigin);
}
#endif // __DUKHART_SDF_HLSL__
