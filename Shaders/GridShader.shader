Shader "Unlit/Grid"
{
    Properties
    {
        [HideInInspector]_MainTex("Texture", 2D) = "white" {}
        _Color("Color", Color) = (0,0,0,1)
        _StrokeColor("Stroke Color", Color) = (1,1,1,1)
        _GridSize("Grid Size", Float) = 3
        _GridAspectX("Grid Aspect X", Float) = 1
        _GridAspectY("Grid Aspect Y", Float) = 1
        _StrokeSize("Stroke Size", Range(0, 1)) = 1
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
            float _GridSize;
            float _GridAspectX;
            float _GridAspectY;
            float _StrokeSize;
            float4 _StrokeColor;
            float4 _Color;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                // center the grid
                o.uv -= _StrokeSize / _GridSize /2;
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                float4 col = _Color;
                float2 aspect = float2(_GridAspectX, _GridAspectY);
                // get uv from the vertex shader and multiply by grid size
                float2 uv = i.uv * _GridSize * aspect;
                float2 gv = frac(uv) - 0.5;
                // give each grid a unique id
                //float2 id = floor(uv);

                // grid outline
                if (gv.x > .5 - _StrokeSize || gv.y > .5 - _StrokeSize) col = _StrokeColor;

                //col.rg = gv;
                return col;
            }
            ENDCG
        }
    }
}
