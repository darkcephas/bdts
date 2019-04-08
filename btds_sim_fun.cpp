#pragma once
#include "btds_sim_fun.h"

#if BDTS_IS_CPP_COMPILE 
namespace bdts
{
#endif
	// DO NOT CHANGE UNLES YOU CHANGE THE INTENGRAL IN  findEnergyIntegral
	float ContactEnergy(float dist, float area, float elemMultEnergy)
	{
		// This function was designed to have a zero energy and zero force at dist = minContactDistance
		// then to gradually increase to infinity as dist goes to zero.


		if (dist >= minContactDistance) return 0.0f;
		float factor = area * elemMultEnergy *contactEnergyMult;
		float k = minContactDistance;
		float x = dist;
		float energy = (1.0f / x) + (x / (k * k)) - (2.0f / k);
		// force = (1/ (k *k )) - (1/ (x*x))
		energy = energy * factor;
		return bdts_max(energy, 0.0f);
	}


	float estimateForceFromContactEnergy(float contactEnergy, float area, float elemMultEnergy)
	{
		// found by derivative of 
		// required for current friction mechanism
		float k = minContactDistance;
		float factor = area * elemMultEnergy *contactEnergyMult;
		float y = contactEnergy / factor;
		float f0 = 1.0f / (k *k);
		float f1 = -4 / (k *k *pow(k* y - sqrt(k) *sqrt(y) *sqrt(k *y + 4.0f) + 2.0f, 2.0f));

		float force = (f0 + f1) * factor;
		return bdts_max(-force, 0.0f);
	}

	float pureIntegratePointEnergyFunction(float xIn, Line2d asLine)
	{
		float k = minContactDistance;
		float a = asLine.getCompA();
		float b = asLine.getCompB();
		float x = xIn;
		float e0 = k * (-a * k * cos(x) - 2 * b * x + k * sin(x));
		float e1 = -(2 * b *b* atanh((a - tan(x / 2)) / sqrt(a *a + 1))) / sqrt(a *a + 1);

		float divFull = (b * k *k);

		float eout = (e0 + e1) / divFull;
		return eout;
	}

	float pureIntegrateLineEnergyFunction(float xInput, Line2d lineIn)
	{
		//return (logN(lineIn.getCompA()* xInput + lineIn.getCompB()) / lineIn.getCompA()) - (xInput / minContactDistance);
		float x = xInput;
		float k = minContactDistance;
		float a = lineIn.getCompA();
		float b = lineIn.getCompB();
		float eout = (((b + a * x) * (b - 4 * k + a * x)) / (k * k) + 2 * logN(b + a * x)) / (2 * a);
		return eout;
	};

	float integratePointEnergy(Line2d asLine, float xi, float xf, float elemMultEnergy)
	{
		float2 vt0;
		float2 vt1;
		if (asLine.getStartStopForDistance(minContactDistance, vt0, vt1, true))
		{
			vt0 = Normalize(vt0);
			vt1 = Normalize(vt1);

			//tracingArrows.push_back({ float2(0,0),  vt0, float4(1,0.5f,1,1) });
			//tracingArrows.push_back({ float2(0,0), vt1, float4(1,1,0.5f,1) });
			vt0 = horzNormDir(vt0);
			vt1 = horzNormDir(vt1);
			xi = bdts_min(vt1.x, bdts_max(xi, vt0.x));
			xf = bdts_min(vt1.x, bdts_max(xf, vt0.x));
		}
		else
		{
			return 0.0f;// entire line is too far away.
		}


		// conversion from direction to actual x coord
		float anglei = float2(xi, sqrtf(1.0f - xi * xi)).GetAngle();
		float anglej = float2(xf, sqrtf(1.0f - xf * xf)).GetAngle();


		if (anglei - anglej < 0.0001f) return 0.0f;

#if 0
		tracingArrows.push_back({ float2(0,0), float2(cos(anglei), sin(anglei)), float4(1,0.5f,1,1) });
		tracingArrows.push_back({ float2(0,0), float2(cos(anglej), sin(anglej)), float4(1,1,0.5f,1) });
#endif

		// ########### we reorientated the angles such that y is zero rad up. 
		// this gives y = cos(@) and x = -sin(@)

		anglei = anglei - (bdts_pi / 2);
		anglej = anglej - (bdts_pi / 2);

		//BDTS_ASSERT(fabsf(xi - xf) > 0.001f, 0.0f);

		float energy = contactEnergyMult * elemMultEnergy*(pureIntegratePointEnergyFunction(anglei, asLine) - pureIntegratePointEnergyFunction(anglej, asLine));

		//BDTS_ASSERT(energy >= -0.005f , 0.0f);
		return bdts_max(0.0f, energy) / bdts_pi;
	}


	float4 UnitCellFromStrainToStress(const float4  &strain, float ym, float pr)
	{
		float multFactor = ym / (1.0f - pr * pr);
		float stressX = (strain.x + pr * strain.y) * multFactor;
		float stressY = (strain.x * pr + strain.y) * multFactor;
		float stressShear0 = (1 - pr) * strain.z * multFactor;
		float stressShear1 = (1 - pr) * strain.w * multFactor;

		return float4(stressX, stressY, stressShear0, stressShear1);
	}

	float4 UnitCellFromStressToStrain(const float4 & stress, float ym, float pr)
	{
		float multFactor = 1.0f / ym;
		float strainX = (stress.x - pr * stress.y) * multFactor;
		float strainY = (-stress.x * pr + stress.y) * multFactor;
		float strainShear0 = (1 + pr) * stress.z * multFactor;
		float strainShear1 = (1 + pr) * stress.w * multFactor;
		return float4(strainX, strainY, strainShear0, strainShear1);
	}


	float UnitCellStrainEnergyRelation(const float4  &strain, float ym, float pr)
	{
		float4 stress = UnitCellFromStrainToStress(strain, ym, pr);
		float energyPerUnit = Dot(stress, strain);
		return energyPerUnit * 0.5f;
	}



	float findEnergyIntegral(Line2d lineIn, float2 winxixf, float elemMultEnergy)
	{
		float xi = winxixf.x;
		float xf = winxixf.y;


		if (lineIn.isNearlyVertical())// very sharp
			return 0.f;// effective zero energy

		float di = lineIn.fromXToY(xi, minContactDistance);
		float df = lineIn.fromXToY(xf, minContactDistance);

		if (lineIn.isNearlyHorizontal())
		{
			// very flat line. 
			// this number was discovered by the fact its not numerically stable for an (a)
			float avgEnergy = (ContactEnergy(di, (xf - xi), elemMultEnergy) + ContactEnergy(df, (xf - xi), elemMultEnergy)) *0.5f;
			return avgEnergy;
		}

		if (df >= minContactDistance && di > minContactDistance)
			return 0.0f;// behind

		if (di >= minContactDistance)
		{
			xi = lineIn.fromYToX(minContactDistance, xi);
			di = minContactDistance;
		}


		if (df >= minContactDistance)
		{
			xf = lineIn.fromYToX(minContactDistance, xf);
			df = minContactDistance;
		}

		//BDTS_ASSERT(df <= minContactDistance, 0.0f);
		//BDTS_ASSERT(di <= minContactDistance, 0.0f);

		float intRes = (pureIntegrateLineEnergyFunction(xf, lineIn) - pureIntegrateLineEnergyFunction(xi, lineIn))*contactEnergyMult*elemMultEnergy;

		if (intRes < -0.0001f)
		{
			//	BDTS_ASSERT( false, 0.0f);
		}

		return    bdts_max(0.0f, intRes);
	}


	float4 GetTriShapedParams(const float2& vert0, const float2& vert1, const float2& vert2)
	{
		float2 base = vert1 - vert0;
		float2 bar = vert2 - vert0;
		float2 baseNorm = Normalize(base);
		float2 crossNorm = baseNorm.GetCross();

		float barLength = Length(base);
		float barPerpCross = Dot(crossNorm, bar);
		float barPerpBase = Dot(baseNorm, bar);

		float barPerpCrossSign = barPerpCross >= 0.0f ? 1.0f : -1.0f;

		barPerpCross = barPerpCross * barPerpCrossSign; // forces the cross to be positive
		return { barLength, barPerpCross, barPerpBase,  barPerpCrossSign };
	}



	Matrix22 rotRemoval(const Matrix22& matIn)
	{
		//   [a b]
		//   [c d]   for any a,b,c,d
		float a = matIn.Get().x;
		float b = matIn.Get().y;
		float c = matIn.Get().z;
		float d = matIn.Get().w;

		float asX = (b - c) / (a + d);

		float cosAtan = 1 / sqrtf(asX*asX + 1.0f);  //cos(atan(x)) from wolfram
		float sinAtan = asX / sqrtf(asX*asX + 1.0f);  //sin(atan(x)) from wolfram

		// we desire a rotation free deformation matrix
		//   [matH matS]
		//   [matS matV]    as deformation matrix without rotation
		return Matrix22({ cosAtan, -sinAtan, sinAtan, cosAtan });
	}


	Line2d transLine(const Line2d& line, const CononicalTrans2d& trans)
	{
		float2 pi = trans.transPos(line.pi);
		float2 pf = trans.transPos(line.pf);
		return Line2d(pi, pf);
	}

#if BDTS_IS_CPP_COMPILE 
}
#endif