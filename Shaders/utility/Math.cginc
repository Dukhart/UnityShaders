#ifndef __DUKHART_MATH_HLSL__
#define __DUKHART_MATH_HLSL__

// for suedo random numbers
float Hash21(float2 val)
{
    val = frac(val * float2(123.345, 345.567));
    val += dot(val, val + 54.321);
    return frac(val.x * val.y);
}
// for ShaderToy mod
float mod(float x, float y)
{
    return x - y * floor(x / y);
}
float mod(float2 x, float2 y)
{
    return x - y * floor(x / y);
}
float mod(float3 x, float3 y)
{
    return x - y * floor(x / y);
}
// for 2D rotation matrix
float2x2 rotMat(float angle)
{
    float s = sin(angle);
    float c = cos(angle);
    return float2x2(c, -s, s, c);
}
float3 not(float3 v){
    v.x = v.x == 0 ? 1 : 0;
    v.y = v.y == 0 ? 1 : 0;
    v.z = v.z == 0 ? 1 : 0;
    return v;
}

float dot2(float2 v)
{
    return dot(v, v);
}
float dot2(float3 v)
{
    return dot(v, v);
}
float ndot(float2 a, float2 b)
{
    return a.x * b.x - a.y * b.y;
}

#endif // __DUKHART_MATH_HLSL__
