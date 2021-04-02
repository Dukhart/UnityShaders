Shader "Unlit/RayMarch4"
{
	Properties
	{
		_Color("Color", Color) = (1,1,1,1)
		_MainTex("Texture", 2D) = "white" {}
		_Color2("Color 2", Color) = (1,1,1,1)
		_Tex2("Tex 2", 2D) = "white" {}
		_Size("Size", Vector) = (.1,.1,.1,.1)
		_Offset("Offset", Vector) = (0,0,0,0)
		_Rotation("Rotation", Vector) = (0,0,0,0)
		_Scale("Scale", Vector) = (1,1,1,1)
		_MaxDist("Ray Distance", int) = 1000
		_MaxSteps("Ray Steps", int) = 100
		_SurfDist("Surface Distance", Float) = 0.01
		_Shape("Shape", Range(0,5)) = 0
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
				#include "utility/SDF.cginc"

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
				//Buffer<SphereSDF> mybuffer : register(t0);

				sampler2D _MainTex, _GrabTexture, _Tex2;
				float _SurfDist;
				float3 _Offset, _Rotation, _Scale;
				float4 _MainTex_ST, _Color, _Color2, _Size;
				int _MaxDist, _MaxSteps, _Shape;

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
					return o;
				}

				fixed4 frag(v2f i) : SV_Target
				{

					 //fixed4 tex = tex2D(_MainTex, i.uv) * _Color;
					 fixed4 col = tex2D(_GrabTexture, i.grabUV.xy / i.grabUV.w);

					 // -0.5 to center uv on each face
					 float2 uv = i.uv - 0.5;
					 //float mask = dot(uv, uv);

					 // set the ray to be just behind the camera
					 float3 rayOrigin = i.rayOrigin;
					 float3 rayDir = normalize(i.hitPos - rayOrigin);
					 float4 rm = 0;
					 float dist = 0;

					 SphereSDF testSDF;
					 testSDF.radius = _Scale;

					 BoxSDF testSDF2;
					 testSDF2.size = _Size;


					 dist = RayMarch(rayOrigin, rayDir, testSDF);

					 //rm = RayMarch2(rayOrigin, rayDir, testSDF, testSDF2);
					 //dist = rm.w;

					 if (dist < _MaxDist) {
						 //float3 p = GetPoint(rayOrigin, dist, rayDir);
						 //float3 n = GetNormal(p, _Shape);
						 float4 c = float4(0, 0, 0, 1);
						 float4 c1 = tex2D(_MainTex, uv) * _Color;
						 //float4 c2 = tex2D(_Tex2, uv) * _Color2;
						 //float dist2 = RayMarch(rayOrigin, rayDir, testSDF2, float3(.1,0,0));
						 col = c1;
					 }
					 //col = lerp(col, tex, smoothstep(.1, .2, mask));

					 //col.rg = uv;
					 return col;
				
				}
				ENDCG
			 }
		}
}
