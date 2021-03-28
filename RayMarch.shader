Shader "Unlit/RayMarch"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _MaxDist("Ray Distance", int) = 1000
        _MaxSteps("Ray Steps", int) = 100
        _SurfDist("Surface Distance", Float) = 0.01
        _Shape("Shape", Range(0,4)) = 0
        _ShapeSize1("Size", Float) = 0.3
        _ShapeSize2("Size", Float) = 0.2
    }
    SubShader
    {
        Tags { "RenderType"="Transparent" "Queue" = "Transparent"}
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

            sampler2D _MainTex, _GrabTexture;
            float4 _MainTex_ST;
            int _MaxDist, _MaxSteps, _Shape;
            float _SurfDist, _ShapeSize1, _ShapeSize2;

            v2f vert (appdata v)
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

            // shapes
            float squareSDF(float3 p, float size) {
                return length(max(abs(p) - size,0));
                //return 0.1;
            }
            float sphereSDF(float3 p, float size) {
                return length(p) - size;
            }
            float torusSDF(float3 p, float sizeA, float sizeB) {
                return length(float2(length(p.xz) - sizeA, p.y)) - sizeB;
            }
            float capsuleSDF(float3 p, float radius, float height) {
                float3 a = float3(0,-height*0.5,0);
                float3 b = float3(0, 1*height, 0);
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
            // operations
            float intersectSDF(float distA, float distB) {
                return max(distA, distB);
            }
            float unionSDF(float distA, float distB) {
                return min(distA, distB);
            }
            float differenceSDF(float distA, float distB) {
                return max(distA, -distB);
            }

            // calculates the distance from ray point to surface
            float GetDist(float3 p, int shape) {
                // scene shape goes here
                if (shape == 0) return sphereSDF(p, _ShapeSize1);
                else if (shape == 1) return squareSDF(p, _ShapeSize1);
                else if (shape == 2) return torusSDF(p, _ShapeSize1, _ShapeSize2);
                else if (shape == 3) return capsuleSDF(p, _ShapeSize1, _ShapeSize2);
                else if (shape == 4) return cylinderSDF(p, _ShapeSize1, _ShapeSize2);
                
                // backup
                float dist = sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2)) - 0.5; // = circle of radius 0.5 uv
                return dist;
            }

            // calculates the normal at a point of given shape
            float3 GetNormal(float3 p, int shape) {
                float2 e = float2(0.01f, 0);
                float3 n = GetDist(p,shape) - float3(
                    GetDist(p - e.xyy, shape),
                    GetDist(p - e.yxy, shape),
                    GetDist(p - e.yyx, shape)
                    );
                return normalize(n);
            }

            // get the point on the ray
            float3 GetPoint(float3 rayOrigin, float distOrigin, float3 rayDir) {
                return rayOrigin + distOrigin * rayDir;
            }

            float RayMarch(float3 rayOrigin, float3 rayDir, int shape) {
                float distOrigin = 0;
                float distScurface;
                for (int i = 0; i < _MaxSteps; i++)
                {
                    // get the ray marching point
                    float3 p = GetPoint(rayOrigin, distOrigin, rayDir);
                    // get the distance to the surface from the ray marching point
                    distScurface = GetDist(p, shape);
                    // move our origin by surface distance
                    distOrigin += distScurface;
                    // if distScurface < _SurfDist we hit something, 
                    // if dist Origin > maxDist we reached the end of the ray and didn't hit anything
                    if (distScurface < _SurfDist || distOrigin > _MaxDist) {
                        break;
                    }
                }
                return distOrigin;
            }
            
            fixed4 frag(v2f i) : SV_Target
            {
                fixed4 tex = tex2D(_MainTex, i.uv);
                fixed4 col = tex2D(_GrabTexture, i.grabUV.xy / i.grabUV.w);
                
                // -0.5 to center uv on each face
                float2 uv = i.uv - 0.5;
                float mask = dot(uv, uv);
                // set the ray to be just behind the camera
                float3 rayOrigin = i.rayOrigin;
                // 
                float3 rayDir = normalize(i.hitPos - rayOrigin);
                float dist = RayMarch(rayOrigin, rayDir, _Shape);

                if (dist < _MaxDist) {
                    float3 p = GetPoint(rayOrigin, dist, rayDir);
                    float3 n = GetNormal(p, _Shape);
                    col.rgb = n;
                }

                col = lerp(col, tex, smoothstep(.1, .2, mask));

                //col.rg = uv;
                return col;
            }
            ENDCG
        }
    }
}
