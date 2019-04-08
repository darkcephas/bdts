#ifndef BDTS_SIM_PHY__INCLUDED
#define BDTS_SIM_PHY__INCLUDED

#include "bdts_base.h"
#include "btds_math.h"
#include "btds_func.h"
#include "btds_sim_obj.h"

#if BDTS_IS_CPP_COMPILE 
namespace bdts
{
#endif
	const static float2 gravity = float2(0.0f, -10.0f);// ;

	float EnergyForContactSurfaceElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai, eEnergyComputePhase ePhase);
	int GetSurfaceNeighIndexForContact(SimObj* o, int elementIndex);
	float getElementSideDistance(SimObj* o, ReplacePos* rp, int elementIndex, int triIndex, int ri);
	float EnergyForContactFrictionElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai, eEnergyComputePhase ePhase, float energySurface);
	void InitElementsSideDistances(SimObj* o, int elementIndex, int ri);
	float EnergyForContactPointElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai, eEnergyComputePhase ePhase);
	float EnergyForContactElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai, eEnergyComputePhase ePhase);
	float EnergyForElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai, eEnergyComputePhase ePhase);
	float4 GetShapeSkewedMat(SimObj* o, int elementIndex, const float4& currShapedParams);
	float4 rotFreeDeformForElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai);
	float DiffusionForElement(SimObj* o, int elementIndex, int ri, int vi, int ai, float deltaT);
	bool tryAdaptive(SimObj* o, int elemIndexi, int elemIndexj, int neighVertIndexi, int neighVertIndexj);
	void ClearForce(SimObj* o, int ai);
	void ForceForElement(SimObj* o, int elementIndex, int ri, int vi, int ai);
	Triangle2d toTriangle(SimObj* o, int* indexList);

#if BDTS_IS_CPP_COMPILE 
}
#endif

#endif //BDTS_SIM_PHY__INCLUDED