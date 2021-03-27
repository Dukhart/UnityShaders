Shader "Unlit/TruchetTiling3"
{
    Properties
    {
        _MainColor ("Color", Color) = (1,1,1,1)
        _MainTex("Texture", 2D) = "white" {}
        _BGColor("BG Color", Color) = (1,1,1,1)
        _BGTex("Background Texture", 2D) = "white" {}
        _Width("Line Width", Float) = .1
        _Fade("Line Fade", Float) = .1
        _Size("Size", Float) = 10
        _Speed("Speed", Float) = 0.01
        _TimeOffset("Time", Float) = 0
    }
        SubShader
        {
            Tags { "RenderType" = "Opaque" }
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

                sampler2D _MainTex, _BGTex;
                float4 _MainTex_ST, _MainColor, _BGColor;
                float _Width, _Size, _Speed, _Fade, _TimeOffset;

                float Hash21(float2 val) {
                    val = frac(val * float2(123.345, 345.567));
                    val += dot(val, val + 54.321);
                    return frac(val.x * val.y);
                }
                float mod(float x, float y)
                {
                    return x - y * floor(x / y);
                }
                v2f vert(appdata v)
                {
                    v2f o;
                    o.vertex = UnityObjectToClipPos(v.vertex);
                    o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                    UNITY_TRANSFER_FOG(o,o.vertex);
                    return o;
                }

                fixed4 frag(v2f i) : SV_Target
                {
                    // prevent time from overflowing
                    float time = fmod(_Time.y, 7200);
                    time += _TimeOffset;
                    // sample the texture
                    fixed4 col = tex2D(_BGTex, i.uv) * _BGColor;
                    //col = 0; // init black
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
                    if (n < .5) gv.x *= -1;
                    //else if (n > 6) gv.y *= -1;
                    float2 cUV = gv - sign(gv.x + gv.y) * .5;
                    // create lines
                    float dist = length(cUV);
                    // mask for our lines
                    float mask = smoothstep(_Fade, -_Fade, abs(dist - .5) - width);

                    float angle = atan2(cUV.x, cUV.y);
                    float checker = mod(id.x + id.y, 2)*2-1;
                    float flow = sin(checker * angle * 10 + time);

                    float x = frac(checker * angle / 1.57 + time); // 1/4 pi
                    float y = (dist-(.5-width))/(2*width);
                    y = abs(y - .5) * 2;
                    float2 tUV = float2(x, y);
                    //float y = dist;
                    col *= 1 - mask;
                    col.rgb += tex2D(_MainTex, tUV) * mask * _MainColor;

                    //col *= 1 - tUV.y;
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
