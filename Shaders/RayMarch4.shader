Shader "Unlit/RayMarch4"
{
    Properties
    {
        _Color("Color", Color) = (1,1,1,1)
        _MainTex("Texture", 2D) = "white" {}
        _RayMarchColor1("RayMarchColor", Color) = (1,0,0,1)
        _RayMarchColor2("RayMarchColor2", Color) = (0,0,1,1)
        _RayMarchTex("RayMarchTex", 2D) = "white" {}
        _Operation("Mode", Range(0,4)) = 0

        _ShapeA("ShapeA", Range(0,6)) = 0
        _OffsetA("Offset A", Vector) = (0.1,0,0,0)
        _RotationA("Rotation A", Vector) = (0,0,0,0)
        _ScaleA("Scale A", Vector) = (1,1,1,1)
        _ShapeSizeA("Size A", Vector) = (.2,.1,.1,.1)

        _ShapeB("ShapeB", Range(0,6)) = 1
        _OffsetB("Offset B", Vector) = (-0.1,0,0,0)
        _RotationB("Rotation B", Vector) = (0,0,0,0)
        _ScaleB("Scale B", Vector) = (1,1,1,1)
        _ShapeSizeB("Size B", Vector) = (.2,.1,.1,.1)

        _MaxDist("Ray Distance", int) = 1000
        _MaxSteps("Ray Steps", int) = 100
        _SurfDist("Surface Distance", Float) = 0.01
    }
        SubShader
        {
            Tags { "RenderType" = "Transparent" "Queue" = "Transparent"}
            LOD 100
            GrabPass {
                "_GrabTexture"
            }
            Pass
            {
                CGPROGRAM
                #pragma vertex vert
                #pragma fragment frag
                // make fog work
                #pragma multi_compile_fog
                //#pragma surface surf Lambert alpha

                #include "UnityCG.cginc"

                struct appdata
                {
                    float4 vertex : POSITION;
                    float2 uv : TEXCOORD0;
                };

                struct v2f
                {
                    float2 uv : TEXCOORD0;
                    UNITY_FOG_COORDS(1)
                    float4 vertex : SV_POSITION;
                    float3 rayOrigin : TEXCOORD1;
                    float3 hitPos : TEXCOORD2;
                    float4 grabUV : TEXCOORD3;
                };

                sampler2D _MainTex, _GrabTexture, _RayMarchTex;
                float _SurfDist;
                float4 _ShapeSizeA, _ShapeSizeB;
                float3 _OffsetA, _RotationA, _ScaleA, _OffsetB, _RotationB, _ScaleB;
                float4 _MainTex_ST, _Color, _RayMarchColor1, _RayMarchColor2;

                int _MaxDist, _MaxSteps, _ShapeA, _ShapeB, _Operation;


                v2f vert(appdata v)
                {
                    v2f o;
                    o.vertex = UnityObjectToClipPos(v.vertex);
                    o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                    UNITY_TRANSFER_FOG(o,o.vertex);
                    o.grabUV = UNITY_PROJ_COORD(ComputeGrabScreenPos(o.vertex));
                    // get the camera
                    //o.rayOrigin = _WorldSpaceCameraPos; // world space
                    //o.hitPos = mul(unity_ObjectToWorld, v.vertex); // world space
                    o.rayOrigin = mul(unity_WorldToObject, float4(_WorldSpaceCameraPos,1)); // object space
                    o.hitPos = v.vertex; // object space
                    //o.Alpha = 
                    return o;
                }

// UTILITY FUNCTIONS
#include "utility/Math.cginc"
#include "utility/SDF.cginc"
                
                float3 ApplyPointProps(float3 p, float3 offset, float3 rotation, float3 scale) {
                    // apply offset
                    p -= offset;
                    // apply rotation
                    p.xz = mul(p.xz, rotMat(rotation.y));
                    p.xy = mul(p.xy, rotMat(rotation.z));
                    p.yz = mul(p.yz, rotMat(rotation.x));
                    // apply scale
                    p *= scale;
                    //return
                    return p;
                }
                // calculates the distance from ray point to surface
                float GetDist(float3 p, int shape, float3 scale, float4 size) {
                    float dist = 0;
                    // scene shape goes here
                    if (shape == 0) dist = sphereSDF(p, size.x);
                    else if (shape == 1) dist = boxSDF(p, size.xyz);
                    else if (shape == 2) dist = torusSDF(p, size.xy);
                    else if (shape == 3) dist = capsuleSDF(p, size.x, size.y);
                    else if (shape == 4) dist = cylinderSDF(p, size.x, size.y);
                    else if (shape == 5) dist = planeSDF(p, float3(size.x,1,size.y),size.z);
                    else if (shape == 6) dist = boxFrameSDF(p, size.xyz, size.w);
                    else dist = sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2)) - 0.5; // backup
                    // compensate for scale
                    dist = dist / max(max(scale.x, scale.y), scale.z);
                    return dist;
                }

                // calculates the normal at a point of given shape
                float3 GetNormal(float3 p, int shapeA, int shapeB) {
                    float2 e = float2(0.01f, 0);
                    float3 n;

                    float3 p1 = ApplyPointProps(p, _OffsetA, _RotationA, _ScaleA);
                    float3 p2 = ApplyPointProps(p, _OffsetB, _RotationB, _ScaleB);
                    float d1 = GetDist(p1, shapeA, _ScaleA, _ShapeSizeA);
                    float d2 = GetDist(p2, shapeB, _ScaleB, _ShapeSizeB);

                    if (min(d1, d2) == d1) {
                        n = GetDist(p1, shapeA, _ScaleA, _ShapeSizeA) - float3(
                            GetDist(p1 - e.xyy, shapeA, _ScaleA, _ShapeSizeA),
                            GetDist(p1 - e.yxy, shapeA, _ScaleA, _ShapeSizeA),
                            GetDist(p1 - e.yyx, shapeA, _ScaleA, _ShapeSizeA)
                            );
                    }
                    else {
                        n = GetDist(p2, shapeB, _ScaleB, _ShapeSizeB) - float3(
                            GetDist(p2 - e.xyy, shapeB, _ScaleB, _ShapeSizeB),
                            GetDist(p2 - e.yxy, shapeB, _ScaleB, _ShapeSizeB),
                            GetDist(p2 - e.yyx, shapeB, _ScaleB, _ShapeSizeB)
                            );
                    }

                    return normalize(n);
                }



                float2 RayMarch(float3 rayOrigin, float3 rayDir, int shapeA, int shapeB) {
                    float distOrigin = 0;
                    float distSurface, distSurface1, distSurface2, colLerp;
                    float time = mod(_Time.y, 72000);
                    for (int i = 0; i < _MaxSteps; i++)
                    {
                        // get the ray marching point
                        float3 p = GetPoint(rayOrigin, distOrigin, rayDir);
                        float3 p1 = ApplyPointProps(p,_OffsetA, _RotationA, _ScaleA);
                        float3 p2 = ApplyPointProps(p, _OffsetB, _RotationB, _ScaleB);
                        p2.z += sin(p2.x * 10 + time) * 0.1;
                        // get the distance to the surface from the ray marching point
                        distSurface1 = GetDist(p1, shapeA, _ScaleA, _ShapeSizeA);// -sin(p1.x * 20 + time) * .1) * 0.2;
                        distSurface2 = GetDist(p2, shapeB, _ScaleB, _ShapeSizeB);// -cos(p2.x * 20 + time) * .1) * 0.2;;


                        if (_Operation == 0) {
                            distSurface = unionSDF(distSurface1, distSurface2);
                            colLerp = distSurface == distSurface1 ? 0 : 1;
                        }
                        else if (_Operation == 1) {
                            distSurface = intersectSDF(distSurface1, distSurface2);
                            colLerp = distSurface == distSurface1 ? 0 : 1;
                        }
                        else if (_Operation == 2) {
                            distSurface = differenceSDF(distSurface1, distSurface2);
                            colLerp = distSurface == distSurface1 ? 0 : 1;
                        }
                        else if (_Operation == 3) {
                            float t = (_ShapeSizeA.z + _ShapeSizeB.z) * 0.05;
                            float v = 0.01;
                            distSurface = smoothUnionSDF(distSurface1, distSurface2, t);

                            float2 d = float2(distSurface + v, distSurface - v);

                            if (distSurface1 > d.y && distSurface1 < d.x) colLerp = 0;
                            else if (distSurface2 > d.y && distSurface2 < d.x) colLerp = 1;
                            else {
                                float d1 = distSurface1 - distSurface;
                                float d2 = distSurface2 - distSurface;
                                if (d1 < d2)
                                    colLerp = 0 + d1;
                                else
                                    colLerp = 1 - d2;
                            }
                        }
                        else if (_Operation == 4) {
                            float t = sin(mod(_Time.y * ((_ShapeSizeA.z + _ShapeSizeB.z) * .5), 7200)) * .5 + .5;
                            distSurface = lerp(distSurface1, distSurface2, t);
                            if (distSurface == distSurface1) colLerp = 0;
                            else if (distSurface == distSurface2) colLerp = 1;
                            else colLerp = t;
                        }
                        else {
                            distSurface = min(distSurface1, distSurface2);
                            colLerp = distSurface == distSurface1 ? 0 : 1;
                        }

                        // move our origin by surface distance
                        distOrigin += distSurface;
                        // if distSurface < _SurfDist we hit something, 
                        // if dist Origin > maxDist we reached the end of the ray and didn't hit anything
                        if (distSurface < _SurfDist || abs(distOrigin) > _MaxDist) {
                            break;
                        }
                    }
                    return float2(distOrigin, colLerp);
                }

                fixed4 frag(v2f i) : SV_Target
                {
                    fixed4 tex = tex2D(_MainTex, i.uv) * _Color;
                    fixed4 col = tex2D(_GrabTexture, i.grabUV.xy / i.grabUV.w);

                    // -0.5 to center uv on each face
                    float2 uv = i.uv - 0.5;
                    float mask = dot(uv, uv);
                    // set the ray to be just behind the camera
                    float3 rayOrigin = i.rayOrigin;
                    // 
                    float3 rayDir = normalize(i.hitPos - rayOrigin);
                    float2 rm = RayMarch(rayOrigin, rayDir, _ShapeA, _ShapeB);
                    float dist = rm.x;


                    if (dist < _MaxDist) {
                        float3 p = GetPoint(rayOrigin, dist, rayDir);
                        //float3 n = GetNormal(p, _ShapeA, _ShapeB);
                        col = tex2D(_RayMarchTex, i.uv) * lerp(_RayMarchColor1, _RayMarchColor2, rm.y);
                    }

                    //col = lerp(col, tex, smoothstep(.1, .2, mask));

                    //col.rg = uv;
                    return col;
                }
                ENDCG
            }
        }
}
