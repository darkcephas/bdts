#pragma once
#include "btds_sim_obj.h"

#if BDTS_IS_CPP_COMPILE 
namespace bdts
{
#endif
	float nced_evalEnergyState(NonConEnergyDissipator * obj, eEnergyComputePhase ePhase, float resitutitveMult, float energyIn, float& lossOut)
	{
		float currEnergy = bdts_max(nced_getCurrEnergy(obj, resitutitveMult, energyIn), 0.0f);
		if (ePhase == eEnergyComputePhase::PREPERATION_COMPUTE_PHASE)
		{
			nced_startSample(obj, resitutitveMult, energyIn);
		}
		else if (ePhase == eEnergyComputePhase::LOSS_COMPUTE_PHASE)
		{
			lossOut = nced_endLoss(obj, resitutitveMult, energyIn);
			currEnergy = bdts_max(nced_getCurrEnergy(obj, resitutitveMult, energyIn), 0.0f);
			obj->prevEnergy = currEnergy;
		}
		else
		{
			obj->frameEnergyVariation = bdts_max( obj->frameEnergyVariation, abs(obj->prevEnergy - currEnergy));
		}

		return currEnergy;
	}


	float nced_getCurrEnergy(NonConEnergyDissipator * obj, float resitutitveMult, float energyIn) 
	{
		return resitutitveMult * (energyIn - obj->nominalEnergy);
	}

	void nced_startSample(NonConEnergyDissipator * obj, float resitutitveMult, float energyIn)
	{

	}

	float nced_endLoss(NonConEnergyDissipator * obj, float resitutitveMult, float energyIn)
	{
		float outEnergyLoss = 0.0f;
		float currE = nced_getCurrEnergy(obj, resitutitveMult, energyIn);


		float maxDiffAllowed = obj->frameEnergyVariation * 4.0f;

		if (fabsf(currE) > maxDiffAllowed)
		{
			if (currE < 0.0f)
			{
				if (energyIn > obj->nominalEnergy)
				{
					obj->nominalEnergy = energyIn - (maxDiffAllowed / resitutitveMult);
				}
				else
				{
					obj->nominalEnergy = energyIn + (maxDiffAllowed / resitutitveMult);
				}

			}
			else
			{
				outEnergyLoss = (currE - maxDiffAllowed);
				obj->nominalEnergy = energyIn - (maxDiffAllowed / resitutitveMult);
			}
			obj->nominalEnergy = bdts_max(0.0f, obj->nominalEnergy);
		}
		return outEnergyLoss;
	}


	EMat::EMat(float youngsMod, float maxElasticPerArea, float frictionCof, float frictionSurface, int colIndex)
	{
		em_youngsMod = youngsMod;
		em_maxElasticPerArea = maxElasticPerArea;
		em_frictionCof = frictionCof;
		em_frictionSurface = frictionSurface;


		debugColor = { rndFromIndex(colIndex, 0), rndFromIndex(colIndex, 1), rndFromIndex(colIndex, 2), 1.0f };
	}

	bool microset_contains(tSetObj* obj, int checkid)
	{
		for (int i = 0; i < obj->count; i++)
		{
			if (obj->ids[i] == checkid)
			{
				return true;
			}
		}

		return false;
	}

	void microset_init(tSetObj* obj)
	{
		obj->count = 0;
	}

	bool microset_insert(tSetObj* obj, int checkid)
	{

		if (microset_contains(obj, checkid))
		{
			return false;
		}

		obj->ids[obj->count++] = checkid;
		return true;
	}


	float2 rcheck_replacePos(ReplacePos* rp, int vertIndex, Pos p)
	{
		if (vertIndex == rp->vertIndex)
		{
			return rp->pos;
		}
		return p.v;
	}

#if BDTS_IS_CPP_COMPILE 
}
#endif