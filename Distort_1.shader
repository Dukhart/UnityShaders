Shader "Unlit/Distort_1"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _Distort("Distort", Float) = 1
        _DistortSpeedX("DistortSpeedX", Float) = 0.5
        _DistortSpeedY("DistortSpeedY", Float) = 0.5
        _DistortRangeX("DistortRangeX", Range(0,0.6)) = 0.5
        _DistortRangeY("DistortRangeY", Range(0,0.6)) = 0.5
        _DistortXOffset("DistortXOffset", Range(0,1)) = 0.5
        _DistortYOffset("DistortYOffset", Range(0,1)) = 0.5
        _DistortSize("DistortSize", Range(0,1)) = 0.1
        _DistortSizeInternal("DistortSizeInternal", Range(0,1)) = 0
        
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
            float _DistortSpeedX;
            float _DistortSpeedY;
            float _DistortRangeX;
            float _DistortRangeY;
            float _DistortXOffset;
            float _DistortYOffset;
            float _DistortSize;
            float _DistortSizeInternal;
            float _Distort;

            v2f vert (appdata v)
            {
                v2f o;
                //v.vertex.y += sin(v.vertex.x + _Time.y) * 1;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                // sample the texture
                //float2 uv = i.uv - _DistortVarA;
                float2 uv = float2(i.uv.x - _DistortXOffset, i.uv.y - _DistortYOffset);
                float a = _Time.y * _DistortSpeedX;
                float b = _Time.y * _DistortSpeedY;
                //float a = 0.1;
                //float2 p = float2(sin(a), cos(b)) * _DistortRange;
                float2 p = float2(sin(a) * _DistortRangeX, cos(b) * _DistortRangeY);
                float2 distort = uv - p;
                float d = length(distort);
                float m = smoothstep(_DistortSize, _DistortSizeInternal, d);
                distort = distort * m * _Distort;
                fixed4 col = tex2D(_MainTex, i.uv + distort);

                return col;
            }
            ENDCG
        }
    }
}
