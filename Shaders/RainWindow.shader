Shader "Unlit/RainWindow"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _Size("Resolution", Float) = 10

        _DropSize("_DropSize1", Float) = 10
        _StreakSize("_StreakSize", Float) = 0.001

        _Distortion("Drop Distortion", Range(-5,5)) = -2
        _Blur("Glass Blur", Range(0,1)) = 0.1
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" "Queue"="Transparent"}
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

            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 grabUV : TEXCOORD1;
                UNITY_FOG_COORDS(1)
                float4 vertex : SV_POSITION;
            };

            sampler2D _MainTex, _GrabTexture;
            float4 _MainTex_ST;
            float _Size, _Distortion, _Blur, _StreakSize, _DropSize;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                // center the grid
                o.uv -= 0.1 / _Size / 2;
                o.grabUV = UNITY_PROJ_COORD(ComputeGrabScreenPos(o.vertex));
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }

            float RandomNum(float2 p) {
                p = frac(p * float2(123.345, 345.567));
                p += dot(p, p + 54.321);
                return frac(p.x * p.y);
            }
            float3 DropLayer(float2 uvLayer, float time) {
                float tOffset = .18;
                float2 aspect = float2(2, 1);
                // get uv from the vertex shader and multiply by grid size
                float2 uv = uvLayer * _Size * aspect;
                uv.y += time * tOffset;
                // create the grid
                float2 gv = frac(uv) - 0.5;
                // give each grid a unique id
                float2 id = floor(uv);
                float n = RandomNum(id); // 0 - 1
                time += n * 6.2831;

                float w = uvLayer.y * 10;
                w += n * 6.2831;
                // apply the formula for our x movement curve (drop side to side movement)
                float x = sin(3 * w) * pow(sin(w * n), 6) * tOffset * 2;
                // apply the formula for our y movement curve (drop vertical speed)
                float y = -sin(time + sin(time + sin(time) * 0.5)) * tOffset * 2;
                y += -(gv.x - x) * (gv.x - x);

                float size = _DropSize * (n + 0.1);

                // drops
                float2 dropPos = (gv - float2(x, y)) / aspect;
                float drop = smoothstep(size, size *0.5, length(dropPos));
                // drop trails
                float2 trailPos = (gv - float2(x, tOffset)) / aspect;
                trailPos.y = (frac(trailPos.y * 4) - 0.5) / 4;
                float trail = smoothstep(size, size *0.5, length(trailPos));
                // mask using main drop y position
                float fogTrail = smoothstep(-.005, .005, dropPos.y);
                //fade out trial
                fogTrail *= smoothstep(0.5, y, gv.y);
                //fogTrail *= smoothstep(10.45, y, gv.y);
                trail *= fogTrail;
                //float fogTrail = smoothstep(-0.05, 0.05, dropPos.y);
                fogTrail *= smoothstep(size, size *0.5, abs(dropPos.x));

                float2 offset = drop * dropPos + trail * trailPos;
                return float3(offset, fogTrail);
            }
            float3 StreakLayer(float2 uvLayer, float time) {
                float tOffset = .18;
                float2 aspect = float2(2, 1);
                // get uv from the vertex shader and multiply by grid size
                float2 uv = uvLayer * _Size * aspect;
                uv.y += time * tOffset;
                // create the grid
                float2 gv = frac(uv) - 0.5;
                // give each grid a unique id
                float2 id = floor(uv);
                float n = RandomNum(id); // 0 - 1
                time += n * 6.2831;

                float w = uvLayer.y * 10;
                w += n * 6.2831;
                // apply the formula for our x movement curve (drop side to side movement)
                float x = sin(3 * w) * pow(sin(w * n), 6) * tOffset * 2;
                // apply the formula for our y movement curve (drop vertical speed)
                float y = -sin(time + sin(time + sin(time) * 0.5)) * tOffset * 2;
                y += -(gv.x - x) * (gv.x - x);


                //shapes
                // drops
                float2 dropPos = (gv - float2(x, y)) / aspect;
                float drop = smoothstep(_StreakSize, _StreakSize * 0.5, length(dropPos)); // dot size
                
                // drop trails
                float2 trailPos = (gv - float2(x, tOffset)) / aspect;
                //trailPos.y = (frac(trailPos.y * _StreakVar4) - _StreakVar5) / _StreakVar4;
                float trail = smoothstep(_StreakSize, _StreakSize * 0.5, length(trailPos)); // streak size
                // mask using main drop y position
                float fogTrail = smoothstep(-.005, .005, dropPos.y);
                //fade out trial
                fogTrail *= smoothstep(0.5, y, gv.y);
                //fogTrail *= smoothstep(10.45, y, gv.y);
                trail *= fogTrail;
                //float fogTrail = smoothstep(-0.05, 0.05, dropPos.y);
                fogTrail *= smoothstep(_StreakSize, _StreakSize*0.5, abs(dropPos.x));

                float2 offset = drop * dropPos + trail * trailPos;
                return float3(offset, fogTrail);
                //return float3(offset, trail);

            }
            fixed4 frag(v2f i) : SV_Target
            {
                float4 col = 0;
                // prevent time from overflowing
                float time = fmod(_Time.y, 7200);
                float fade = 1 - saturate(fwidth(i.uv) * 50);

                // grid outline
                //if (gv.x > .5 - 0.1 || gv.y > .5 - 0.1) col = float4(1, 0, 0, 1);

                
                float3 drops = DropLayer(i.uv, time*2);
                drops += DropLayer(i.uv * 1.5 + 3.3, time);
                drops += DropLayer(i.uv * 1.7 - 3.6, time);
                drops += StreakLayer(i.uv - 3.1, time * 10);
                drops += StreakLayer(i.uv * 1.5 + 3.7, time*3);

                float blur = _Blur - drops.z*fade;
                float2 projUV = i.grabUV.xy / i.grabUV.w;
                projUV += drops.xy * _Distortion;// *fade;

                blur *= .1;
                const float numSamples = 32;
                float a = RandomNum(i.uv) * 6.2831;
                for (int i = 0; i < numSamples; i++) {
                    float2 offset = float2(sin(a), cos(a)) * blur;
                    float d = frac(sin((i + 1) * 546) * 8694);
                    offset *= sqrt(d);
                    col += tex2D(_GrabTexture, projUV+offset);
                    ++a;
                }
                col /= numSamples;
                
                return col;
            }
            ENDCG
        }
    }
}
