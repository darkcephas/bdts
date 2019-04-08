#pragma once
#include <math.h>
#include "btds_math.h"
#include <algorithm>
#if BDTS_IS_CPP_COMPILE 
namespace bdts
{
#endif
		inline bool Equal(float2 v1, float2 v2)
		{
			return  v1.x == v2.x && v1.y == v2.y;
		}


		inline bool Quadratic(float a, float b, float c, float& x0, float& x1)
		{
			x0 = 0;
			x1 = 0;
			float valSq = (b * b) - (4.0f *a*c);
			if (valSq < 0.0f)
			{
				return false;// no imagination
			}

			x0 = (-b + sqrtf(valSq)) / (2.0f * a);
			x1 = (-b - sqrtf(valSq)) / (2.0f * a);
			return true;
		}



		inline float Dot(float2 v1, float2 v2) { return v1.x * v2.x + v1.y *v2.y; }
		inline float Dot(float3 v1, float3 v2) { return v1.x * v2.x + v1.y *v2.y + v1.z*v2.z; }
		inline float Dot(float4 v1, float4 v2) { return v1.x * v2.x + v1.y *v2.y + v1.z*v2.z + v1.w * v2.w; }

		inline float Length(float2 v) { return sqrtf(Dot(v, v)); }
		inline float Length(float3 v) { return sqrtf(Dot(v, v)); }
		inline float Length(float4 v) { return sqrtf(Dot(v, v)); }

		inline float2 Normalize(float2 v) { return  v / Length(v); }
		inline float3 Normalize(float3 v) { return  v / Length(v); }
		inline float4 Normalize(float4 v) { return  v / Length(v); }

		inline Matrix22 Invert(const Matrix22& mat)
		{
			float a = mat.Get().x;
			float b = mat.Get().y;
			float c = mat.Get().z;
			float d = mat.Get().w;

			float invDet = 1.0f / (a*d - b * c);
			float4 asVec = float4(d, -b, -c, a) * invDet;
			return Matrix22(asVec);
		}

		// takes a vector and a start
		// converts to untranslated/rot space with dir pointing up in the pos Y
		class CononicalTrans2d
		{
		public:
			CononicalTrans2d(float2 dir, float2 point)
			{
				transXY = -point;
				rowY = dir;
				rowX = -dir.GetCross();
			}

			inline float2 transDir(float2 dir) const
			{
				return float2(Dot(dir, rowX), Dot(dir, rowY));
			}

			inline float2 transPos(float2 point) const
			{
				float2 temp = point + transXY;
				return float2(Dot(temp, rowX), Dot(temp, rowY));
			}

			float2 rowX;
			float2 rowY;
			float2 transXY;
		};


	



		inline static bool FiniteLineCutLine(const Line2d& segToCut, const Line2d& cuttingPlaneLine, float2& posOut, float2& tOut)
		{
			float2 li = segToCut.pi;
			float2 lf = segToCut.pf;
			float2 vi = cuttingPlaneLine.pi;
			float2 vf = cuttingPlaneLine.pf;
			double x3mx4 = ((double)vi.x) - ((double)vf.x);
			double x1mx3 = ((double)li.x) - ((double)vi.x);
			double x1mx2 = ((double)li.x) - ((double)lf.x);
			double y1my2 = ((double)li.y) - ((double)lf.y);
			double y1my3 = ((double)li.y) - ((double)vi.y);
			double y3my4 = ((double)vi.y) - ((double)vf.y);

			double demom = x1mx2 * y3my4 - y1my2 * x3mx4;
			double t;
			double u;
			if (abs(demom) < 0.000001)
			{
				return false;
			}
			else
			{
				t = (x1mx3 * y3my4 - y1my3 * x3mx4) / demom; // det
				u = -(x1mx2* y1my3 - y1my2 * x1mx3) / demom;
			}

			tOut = { (float)t, (float)u };
			posOut = (lf - li) * t + li;
			if (t == 1.0)
			{
				posOut = lf;
			}

			if (t == 0.0)
			{
				posOut = li;
			}

			if (t < 0.0 || u < 0.0)
			{
				return false;
			}

			if (t > 1.0 || u > 1.0f)
			{
				return false;
			}



			return  true;
		}



		inline static float4 toFalseColor(float x)
		{
			x = std::min(std::max(x, 0.0f), 1.0f);
			static const int numFalseColor = 7;
			float4 listCol[numFalseColor] =
			{ { 0.2f,0.2f,.2f,1 } , { 0,0.3f,1,1 } ,{ 0,0,1,1 } ,{ 0,1,1,1 } ,{ 0,1,0,1 } ,{ 1,1,0,1 } ,{ 1,0,0,1 } };

			int indexlow = std::min(std::max((int)(x * (numFalseColor - 1.0f)), 0), numFalseColor - 1);
			int indexHigh = std::min(std::max(indexlow + 1, 1), numFalseColor - 1);

			float lerp = x * (numFalseColor - 1.0f) - indexlow;

			return listCol[indexlow] * (1.0f - lerp) + listCol[indexHigh] * lerp;
		}

#if BDTS_IS_CPP_COMPILE 
}
#endif