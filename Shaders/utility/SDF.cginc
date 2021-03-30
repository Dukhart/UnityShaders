// SHAPES
float squareSDF(float3 p, float3 size) {
    return length(max(abs(p) - size,0));
    //return 0.1;
}
float squareSDF2(float3 p, float size, float wall) {
    return abs(squareSDF(p,size))- wall;
    //return 0.1;
}
float sphereSDF(float3 p, float3 size) {
    return length(p) - size;
}
float torusSDF(float3 p, float sizeA, float sizeB) {
    return length(float2(length(p.xz) - sizeA, p.y)) - sizeB;
}
float capsuleSDF(float3 p, float radius, float height) {
    float3 a = float3(0,-height * 0.5,0);
    float3 b = float3(0, 1 * height, 0);
    float3 ab = b - a;
    float3 ap = p - a;

    float t = dot(ab, ap) / dot(ab, ab);
    t = clamp(t, 0, 1);

    float3 c = a + t * ab;
    return length(p - c) * radius;
}
float cylinderSDF(float3 p, float radius, float height) {
    float3 a = float3(0, -height * 0.5, 0);
    float3 b = float3(0, 1 * height, 0);
    float3 ab = b - a;
    float3 ap = p - a;
    // cylinder
    float t = dot(ab, ap) / dot(ab, ab);
    float3 c = a + t * ab;

    float x = length(p - c) * radius;
    float y = (abs(t - 0.5) - 0.5) * length(ab);
    // exterior
    float e = length(max(float2(x, y), 0));
    // interior
    float i = min(max(x, y), 0);
    return e + i;
}
float planeSDF(float3 p, float3 dir) {
    return dot(p, normalize(dir));
}
// SHAPE OPERATIONS
float intersectSDF(float distA, float distB) {
    return max(distA, distB);
}
float unionSDF(float distA, float distB) {
    return min(distA, distB);
}
float smoothUnionSDF(float distA, float distB, float amount) {
    if (amount == 0) return min(distA, distB);
    float h = clamp(0.5 + 0.5 * (distB - distA) / amount, 0, 1);
    return lerp(distB, distA, h) - amount * h * (1 - h);
}

float differenceSDF(float distA, float distB) {
    return max(distA, -distB);
}

// get the point on the ray
float3 GetPoint(float3 rayOrigin, float distOrigin, float3 rayDir)
{
                    // get the point on the ray
    float3 p = rayOrigin + distOrigin * rayDir;
                    //return
    return p;
}