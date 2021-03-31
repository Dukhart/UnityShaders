Shader "Unlit/RayMarch4"
{
	Properties
	{
		_Color("Color", Color) = (1,1,1,1)
		_MainTex("Texture", 2D) = "white" {}
		_RayMarchColor("RayMarchColor", Color) = (1,1,1,1)
		_RayMarchTex("RayMarchTex", 2D) = "white" {}
		_Offset("Offset", Vector) = (0,0,0,0)
		_Rotation("Rotation", Vector) = (0,0,0,0)
		_Scale("Scale", Vector) = (1,1,1,1)
		_MaxDist("Ray Distance", int) = 1000
		_MaxSteps("Ray Steps", int) = 100
		_SurfDist("Surface Distance", Float) = 0.01
		_Shape("Shape", Range(0,5)) = 0
		_ShapeSize1("Size", Float) = 0.3
		_ShapeSize2("Size", Float) = 0.2
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

				cbuffer cbOnce
				{
			/*
					static SphereSDF g_SphereSDF
					{
						radius = 0.5f;
					};
					
					static BoxSDF g_BoxSDF
					{
						static x = 0.5f;
					};
					*/
				};
		/*
				static const int MAXLIGHTS = 64;
				CBUFFER_START(Once);
				int layer;
				static SphereSDF sphereSDF[MAXLIGHTS];
				int activeLights;
				CBUFFER_END
*/

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

				sampler2D _MainTex, _GrabTexture, _RayMarchTex;
				float _SurfDist, _ShapeSize1, _ShapeSize2;
				float3 _Offset, _Rotation, _Scale;
				float4 _MainTex_ST, _Color, _RayMarchColor;
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

				// calculates the distance from ray point to surface
				/*
				float GetDist(float3 p, int shape) {
					float dist = 0;
					// scene shape goes here
					if (shape == 0) dist = SphereSDF::sphereSDF(p, _ShapeSize1);
					else dist = sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2)) - 0.5; // backup

					return dist;
				}
				*/
				/*
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
				*/
				/*
				 float RayMarch(float3 rayOrigin, float3 rayDir, int shape) {

					 float distOrigin = 0;
					 float distScurface;
					 for (int i = 0; i < _MaxSteps; i++)
					 {
						 // get the ray marching point
						 float3 p = GetPoint(rayOrigin, distOrigin, rayDir);
						 // get the distance to the surface from the ray marching point

						 //distScurface = g_SphereSDF.getDist(p);
						 distScurface = 0;
						 //distScurface = GetDist(p, shape);
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
				 */
				 /*
				 SphereSDF g_SphereSDF
				 {
					 radius = 0.5f;
				 };
				 */
				 fixed4 frag(v2f i) : SV_Target
				 {

					 fixed4 tex = tex2D(_MainTex, i.uv) * _Color;
					 fixed4 col = tex2D(_GrabTexture, i.grabUV.xy / i.grabUV.w);

					 // -0.5 to center uv on each face
					 float2 uv = i.uv - 0.5;
					 float mask = dot(uv, uv);

					 // set the ray to be just behind the camera
					 float3 rayOrigin = i.rayOrigin;
					 float3 rayDir = normalize(i.hitPos - rayOrigin);
					 float dist = 0;

					 SphereSDF testSDF;
					 testSDF.radius = 0.1;
					 dist = RayMarch(rayOrigin, rayDir, testSDF);

					 BoxSDF testSDF2;
					 testSDF2.size = float3(0.1, 0.1, 0.1);
					 dist = RayMarch(rayOrigin, rayDir, testSDF2);
					 

					 if (dist < _MaxDist) {
						 float3 p = GetPoint(rayOrigin, dist, rayDir);
						 //float3 n = GetNormal(p, _Shape);
						 col.rgb = tex2D(_RayMarchTex, i.uv) * _RayMarchColor;
					 }
					 col = lerp(col, tex, smoothstep(.1, .2, mask));

					 //col.rg = uv;
					 return col;
				 }
				 ENDCG
			 }
		}
}
