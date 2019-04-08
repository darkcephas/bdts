#pragma once
#include "btds_math.h"
#include "btds_func.h"

#if BDTS_IS_CPP_COMPILE 
namespace bdts
{
#endif

	bool Triangle2d::isLeft() const
	{
		float2 vec = v[1] - v[0];
		float2 asCross = vec.GetCross();
		return Dot(asCross, v[2] - v[1]) <= 0.0f;
	}

	bool Triangle2d::hasVolume() const
	{
		return getArea() > 0.0f;
	}

	float Triangle2d::getArea() const
	{
		float2 vec = v[1] - v[0];
		float2 asCross = vec.GetCross();
		return fabsf(Dot(asCross, v[2] - v[1])) * 0.5f;
	}

	float Triangle2d::robustness() const
	{

		for (int i = 0; i < 3; i++)
		{
			float length = Length(v[i] - v[(i + 1) % 3]);
			if (length == 0.0f)
				return 0.0f;
		}

		float2 vec0 = Normalize(v[1] - v[0]);
		float2 vec1 = Normalize(v[2] - v[1]);
		float2 vec2 = Normalize(v[0] - v[2]);

		float dot0 = fabsf(Dot(vec0, vec1));
		float dot1 = fabsf(Dot(vec2, vec1));
		float dot2 = fabsf(Dot(vec0, vec2));

		float maxDot = bdts_max(dot0, bdts_max(dot1, dot2));
		return 1.0f - maxDot;
	}


	float2 Triangle2d::getCentroid() const
	{
		float2 avg = v[0] + v[1] + v[2];
		return avg / 3.0f;
	}

	float Triangle2d::getMaxAngle() const
	{
		float maxSoFar = -1.0f;

		for (int i = 0; i < 3; i++)
		{
			float2 forwardVec = Normalize(v[(i + 1) % 3] - v[i]);
			float2 backwardVec = Normalize(v[(i + 2) % 3] - v[i]);
			float currDotSame = 1.0f - Dot(forwardVec, backwardVec);
			maxSoFar = bdts_max(maxSoFar, currDotSame);
		}

		return maxSoFar;
	}

	Triangle2d Triangle2d::toRight() const
	{
		if (isLeft())
		{
			return Triangle2d({ v[0] , v[2], v[1] });
		}
		else
		{
			return *this;
		}
	}

	Triangle2d Triangle2d::toLeft() const
	{
		if (!isLeft())
		{
			return Triangle2d({ v[0] , v[2], v[1] });
		}
		else
		{
			return *this;
		}
	}

	float2 Line2d::getVec() const
	{
		return pf - pi;
	}



	float Line2d::dirXToCoordX(float dirX)
	{
		float dirY = sqrt(1.0f - dirX * dirX);
		float a = getCompA();
		float b = getCompB();

		float t = b / (dirY - a * dirX);
		const float2 vf = pf - pi;
		float2 	outStart = float2(dirX, dirY) * t;

		return outStart.x;
	}

	bool Line2d::getStartStopForDistance(float distance, float2& outStart, float2& outStop, bool limited)
	{
		const float r = distance;
		const float pix = pi.x;
		const float piy = pi.y;
		const float2 vf = pf - pi;
		const float vfx = vf.x;
		const float vfy = vf.y;

		// t*t( vfx*vfx + vfy *vfy)+    2.0f t*( vix*vfx + viy*vfy)      + (vix*vix + viy *viy - r * r)  = 0
		const float a = (vfx*vfx + vfy * vfy);
		const float b = (pix*vfx + piy * vfy)* 2.0f;
		const float c = (pix*pix + piy * piy) - (r * r);
		float t0;
		float t1;
		if (Quadratic(a, b, c, t0, t1))
		{
			float tStart = bdts_min(t0, t1);
			float tStop = bdts_max(t0, t1);

			if (limited)
			{
				tStart = bdts_max(0.0f, bdts_min(1.0f, tStart));
				tStop = bdts_max(0.0f, bdts_min(1.0f, tStop));
			}
			outStart = pi + vf * tStart;
			outStop = pi + vf * tStop;
			return true;
		}

		return false;
	}

	bool  Line2d::hasLength() const
	{
		return Length(pf - pi) > 0.001f;
	}


	// the pi to pf forms a vector which spits the plane. The cross is the up direction
	bool Line2d::isUpRelative(const float2& ptest) const
	{
		float2 lineCross = (pf - pi).GetCross();
		float2 vecTest = (ptest - pi);
		float dotRes = Dot(lineCross, vecTest);

		return dotRes >= 0.0f;// even for degenerate testing it will show as above.
	}

	float Line2d::length() const
	{
		return Length(getVec());
	}

	float Line2d::projPointLineFloat(const float2& ptest) const
	{
		float2 vecTest = (ptest - pi);
		float dotRes = Dot(Normalize(getVec()), vecTest);

		return dotRes;
	}

	// basically projection will be in the range of the original line
	bool Line2d::isPointAlongLine(const float2& ptest) const
	{
		float dotRes = projPointLineFloat(ptest);

		return dotRes >= 0.0f && dotRes <= Length(getVec());// even for degenerate testing it will show as above.
	}


	float Line2d::distToLine(const float2& ptest) const
	{
		float2 lineCross = Normalize((pf - pi).GetCross());
		float2 vecTest = (ptest - pi);
		float dotRes = Dot(lineCross, vecTest);

		return dotRes;
	}
#if BDTS_IS_CPP_COMPILE 
}
#endif