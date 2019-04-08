#ifndef BDTS_SIM_FUN__INCLUDED
#define BDTS_SIM_FUN__INCLUDED

#include "bdts_base.h"
#include "btds_math.h"
#include "btds_func.h"

#if BDTS_IS_CPP_COMPILE 
namespace bdts
{
#endif
	float ContactEnergy(float dist, float area, float elemMultEnergy);

	float estimateForceFromContactEnergy(float contactEnergy, float area, float elemMultEnergy);

	float pureIntegratePointEnergyFunction(float xIn, Line2d asLine);

	float pureIntegrateLineEnergyFunction(float xInput, Line2d lineIn);

	float integratePointEnergy(Line2d asLine, float xi, float xf, float elemMultEnergy);


	float4 UnitCellFromStrainToStress(const float4  &strain, float ym, float pr);

	float4 UnitCellFromStressToStrain(const float4 & stress, float ym, float pr);

	float UnitCellStrainEnergyRelation(const float4  &strain, float ym, float pr);


	float findEnergyIntegral(Line2d lineIn, float2 winxixf, float elemMultEnergy);

	float4 GetTriShapedParams(const float2 & vert0, const float2 & vert1, const float2 & vert2);

	Matrix22 rotRemoval(const Matrix22& matIn);

	Line2d transLine(const Line2d& line, const CononicalTrans2d& trans);
#if BDTS_IS_CPP_COMPILE 
}
#endif

#endif