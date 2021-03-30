// for suedo random numbers
float Hash21(float2 val)
{
    val = frac(val * float2(123.345, 345.567));
    val += dot(val, val + 54.321);
    return frac(val.x * val.y);
}
// for mod
float mod(float x, float y)
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