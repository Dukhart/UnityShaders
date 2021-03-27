Shader "Unlit/TruchetTiling2"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _Width("Line Width", Float) = .1
        _Fade("Line Fade", Float) = .1
        _Size("Size", Float) = 10
        _Speed("Speed", Float) = 0.01
        _TimeOffset("Time", Float) = 0
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

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
                UNITY_FOG_COORDS(1)
                float4 vertex : SV_POSITION;
            };

            sampler2D _MainTex;
            float4 _MainTex_ST;
            float _Width, _Size, _Speed, _Fade, _TimeOffset;

            float Hash21(float2 val) {
                val = frac(val * float2(123.345, 345.567));
                val += dot(val, val + 54.321);
                return frac(val.x * val.y);
            }

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                // prevent time from overflowing
                float time = fmod(_Time.y, 7200);
                time += _TimeOffset;
                // sample the texture
                fixed4 col = tex2D(_MainTex, i.uv);
                col = 0; // init black
                float2 uv = i.uv - 0.5;
                float gridSize = _Size;
                float2 gv = frac(uv * gridSize) - 0.5;
                // get a unique id for each face
                float2 id = floor(uv * gridSize);
                // get a random number using uinique id as seed
                float n = Hash21(id * _TimeOffset);
                // draws face edges
                //if (gv.x > .48 || gv.y > .48) col += fixed4(1, 0, 0, 1);
                //col.rg += gv;
                float width = _Width;
                n = sin(n * 6.2831 + time * _Speed);
                //flip direction based on random number
                if (n < .3) gv.x *= -1;
                else if (n > 6) gv.y *= -1;
                // create lines
                float dist = length(gv - sign(gv.x+gv.y)*.5) - .5;
                dist = abs(dist);
                // mask for our lines
                float mask = smoothstep(_Fade, -_Fade, dist - width);

                col += mask;
                //col.rg += id*0.3;
                //col += n;
                // apply fog
                //UNITY_APPLY_FOG(i.fogCoord, col);
                return col;
            }
            ENDCG
        }
    }
}
