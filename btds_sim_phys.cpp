#pragma once
#include "btds_sim_phys.h"
#include "btds_sim_fun.h"
#define BDTS_ASSERT( expression, rtnValue)
#define ADAPTIVE_ENERGY_FIX_ENABLED 1

#if BDTS_IS_CPP_COMPILE 
namespace bdts
{
#endif


	bool tryAdaptive(SimObj* o, int elemIndexi, int elemIndexj, int neighVertIndexi, int neighVertIndexj)
	{
		auto& elemi = o->elements[elemIndexi];
		auto& elemj = o->elements[elemIndexj];

		ReplacePos rp;
		rp.vertIndex = -1;

		if (elemi.elementType != elemj.elementType) return false;
		// before
		Triangle2d origTri_i = toTriangle(o, elemi.vertindex);
		Triangle2d origTri_j = toTriangle(o, elemj.vertindex);

		float orig_max_angle = bdts_max(origTri_i.getMaxAngle(), origTri_j.getMaxAngle());

		int oi0 = neighVertIndexi;
		int oi1 = (neighVertIndexi + 1) % 3;
		int oi2 = (neighVertIndexi + 2) % 3;

		int oj0 = neighVertIndexj;
		int oj1 = (neighVertIndexj + 1) % 3;
		int oj2 = (neighVertIndexj + 2) % 3;

		BDTS_ASSERT(elemi.vertindex[oi0] == elemj.vertindex[oj1], false);
		BDTS_ASSERT(elemi.vertindex[oi1] == elemj.vertindex[oj0], false);

		int newTriVerti[3] = { elemi.vertindex[oi0], elemj.vertindex[oj2],  elemi.vertindex[oi2] };
		int newTriVertj[3] = { elemj.vertindex[oj0], elemi.vertindex[oi2],  elemj.vertindex[oj2] };

		Triangle2d newTri_i = toTriangle(o, newTriVerti);
		Triangle2d newTri_j = toTriangle(o, newTriVertj);
		float new_max_angle = bdts_max(newTri_i.getMaxAngle(), newTri_j.getMaxAngle());

		float histScale = 1.0f;
		if (newTri_i.isLeft() && newTri_j.isLeft())
		{
			if (orig_max_angle > (new_max_angle*histScale))
			{
#if 0			// adapting non contact. reevaluate this idea later

				// grab this before we mess it up
				auto propdatai = GetElmentUnitAreaProp(elemIndexi);
				auto propdataj = GetElmentUnitAreaProp(elemIndexj);
				auto areai = elemi.origArea;
				auto areaj = elemj.origArea;

				VertMassChangeForElement(elemi, -1.0f, areai);
				VertMassChangeForElement(elemj, -1.0f, areaj);
#endif


				float energyiOrig = EnergyForElement(o, &rp, elemIndexi, 0, 0, eEnergyComputePhase::VARIATION_COMPUTE_PHASE);
				float energyjOrig = EnergyForElement(o, &rp, elemIndexj, 0, 0, eEnergyComputePhase::VARIATION_COMPUTE_PHASE);
				float eOrig = energyiOrig + energyjOrig;

				int oldNeighi[3];
				int oldNeighj[3];
				for (int ee = 0; ee < 3; ee++)
				{
					elemi.vertindex[ee] = newTriVerti[ee];
					elemj.vertindex[ee] = newTriVertj[ee];
					oldNeighi[ee] = elemi.neighElem[ee];
					oldNeighj[ee] = elemj.neighElem[ee];
				}


				// reassign neigh
				elemi.neighElem[0] = oldNeighj[oj1];
				elemi.neighElem[1] = oldNeighi[oi0];
				elemi.neighElem[2] = oldNeighi[oi2];

				elemj.neighElem[0] = oldNeighi[oi1];
				elemj.neighElem[1] = oldNeighj[oj0];
				elemj.neighElem[2] = oldNeighj[oj2];


				// neigh neigh replacement
				// There will be 2 neigh that will need to flip neighs
				int fixOrigNeighi = oldNeighj[oj1];
				if (fixOrigNeighi != -1)
				{
					for (int ee = 0; ee < 3; ee++)
					{
						if (o->elements[fixOrigNeighi].neighElem[ee] == elemIndexj)
							o->elements[fixOrigNeighi].neighElem[ee] = elemIndexi;
					}
				}


				int fixOrigNeighj = oldNeighi[oi1];
				if (fixOrigNeighj != -1)
				{
					for (int ee = 0; ee < 3; ee++)
					{
						if (o->elements[fixOrigNeighj].neighElem[ee] == elemIndexi)
							o->elements[fixOrigNeighj].neighElem[ee] = elemIndexj;
					}
				}

#if 0	 // adapting noncontact objs
				if (elemi.elementType != Element::CONTACT_ELEMENT || elemj.elementType != Element::CONTACT_ELEMENT)
				{
					// blend
					float areaTotal = areai + areaj;
					auto newMat = (propdataj.worldDefMat * areaj + propdatai.worldDefMat * areai) *(1.0f / areaTotal);
					auto tempNew = propdataj;
					tempNew.worldDefMat = newMat;

					volatile float disti[3];
					volatile float distj[3];
					for (int i = 0; i < 3; i++)
					{
						distj[i] = Length(propdataj.expandedVert[i] - propdataj.expandedVert[(i + 1) % 3]);
						disti[i] = Length(propdatai.expandedVert[i] - propdatai.expandedVert[(i + 1) % 3]);
					}


					InitElementsAsCurrent(elemIndexi, &tempNew);
					InitElementsAsCurrent(elemIndexj, &tempNew);
				}
				else
#endif
				{
#if 0 // not actually required for contact since non of the solid element data is used
					InitElementsAsCurrent(elemIndexi);
					InitElementsAsCurrent(elemIndexj);
#endif
				}


				float energyiNew = EnergyForElement(o, &rp, elemIndexi, 0, 0, eEnergyComputePhase::VARIATION_COMPUTE_PHASE);
				float energyjNew = EnergyForElement(o, &rp, elemIndexj, 0, 0, eEnergyComputePhase::VARIATION_COMPUTE_PHASE);

				*o->energyOffset += (energyiOrig + energyjOrig) - (energyiNew + energyjNew);


				if (eOrig < 0.01f)
				{
					return true;
				}

				if (energyiNew + energyjNew < 0.01f)
				{
					return true;
				}

				elemi.multEnergy = 1.0f;
				elemj.multEnergy = 1.0f;
				energyiNew = EnergyForElement(o, &rp, elemIndexi, 0, 0, eEnergyComputePhase::VARIATION_COMPUTE_PHASE);
				energyjNew = EnergyForElement(o, &rp, elemIndexj, 0, 0, eEnergyComputePhase::VARIATION_COMPUTE_PHASE);

#if ADAPTIVE_ENERGY_FIX_ENABLED
				float energyRatio = eOrig / (energyiNew + energyjNew);

				elemi.multEnergy = 1.0f*energyRatio;
				elemj.multEnergy = 1.0f*energyRatio;
#endif
				return true;
			}
		}



		return false;
	}

	// converts 3 index to a full fleshed triangle
	Triangle2d toTriangle(SimObj* o, int* indexList)
	{
		static const int ri = 0;// assumed init
		return Triangle2d(o->verts[indexList[0]].r[ri].v,
			o->verts[indexList[1]].r[ri].v,
			o->verts[indexList[2]].r[ri].v);
	}

	float DiffusionForElement(SimObj* o, int elemIndex, int ri, int vi, int ai, float deltaT)
	{
		auto& element = o->elements[elemIndex];

		if (element.elementType == Element::CONTACT_ELEMENT)
		{
			return  0.0f;
		}

		auto& emat = o->materials[element.materialId];

		float materialDensity = 1.0f;
		float elemMass = element.origArea  * materialDensity;
		float contVertMass = elemMass * 0.3333333333f;
		float2 vert0 = o->verts[element.vertindex[0]].r[ri].v;
		float2 vert1 = o->verts[element.vertindex[1]].r[ri].v;
		float2 vert2 = o->verts[element.vertindex[2]].r[ri].v;

		float2 vertV0 = o->verts[element.vertindex[0]].v[vi];
		float2 vertV1 = o->verts[element.vertindex[1]].v[vi];
		float2 vertV2 = o->verts[element.vertindex[2]].v[vi];

		float v0tov1 = Length(vert0 - vert1);
		float v1tov2 = Length(vert1 - vert2);
		float v2tov0 = Length(vert2 - vert0);
		// momentum is conserved but energy is not (heat due to viscosity)
		float2 m0Transfer = deltaT * emat.em_viscosity * vertV0 * contVertMass;
		float2 m1Transfer = deltaT * emat.em_viscosity * vertV1 * contVertMass;
		float2 m2Transfer = deltaT * emat.em_viscosity * vertV2 * contVertMass;

		float2 m0TransferTom1 = m0Transfer / v0tov1;
		float2 m1TransferTom0 = m1Transfer / v0tov1;

		float2 m1TransferTom2 = m1Transfer / v1tov2;
		float2 m2TransferTom1 = m2Transfer / v1tov2;

		float2 m2TransferTom0 = m2Transfer / v2tov0;
		float2 m0TransferTom2 = m0Transfer / v2tov0;

		// dd/ddx error here ?????? why are we not dividing
		// by distance twice
		float2 v0MomentumChange = (m1TransferTom0 - m0TransferTom1) + (m2TransferTom0 - m0TransferTom2);
		float2 v1MomentumChange = (m2TransferTom1 - m1TransferTom2) + (m0TransferTom1 - m1TransferTom0);
		float2 v2MomentumChange = (m0TransferTom2 - m2TransferTom0) + (m1TransferTom2 - m2TransferTom1);

		vertV0 = ((o->verts[element.vertindex[0]].mMassCount * vertV0) + v0MomentumChange) / o->verts[element.vertindex[0]].mMassCount;
		vertV1 = ((o->verts[element.vertindex[1]].mMassCount * vertV1) + v1MomentumChange) / o->verts[element.vertindex[1]].mMassCount;
		vertV2 = ((o->verts[element.vertindex[2]].mMassCount * vertV2) + v2MomentumChange) / o->verts[element.vertindex[2]].mMassCount;


		o->verts[element.vertindex[0]].v[vi] = vertV0;
		o->verts[element.vertindex[1]].v[vi] = vertV1;
		o->verts[element.vertindex[2]].v[vi] = vertV2;

		return 0.0f;
	}


	float EnergyForContactPointElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai, eEnergyComputePhase ePhase)
	{
		auto& element = o->elements[elementIndex];
		int elementLayer = element.layerIndex;
		float energy = 0.0f;
		for (int eint = 0; eint < 3; eint++)
		{
			int sv0 = eint;
			int sv1 = (eint + 1) % 3;
			int sv2 = (eint + 2) % 3;

			int v0Index = element.vertindex[sv0];

			float2 v1 = rcheck_replacePos(rp, element.vertindex[sv1], o->verts[element.vertindex[sv1]].r[ri]);
			float2 v2 = rcheck_replacePos(rp, element.vertindex[sv2], o->verts[element.vertindex[sv2]].r[ri]);

			float2 v0 = rcheck_replacePos(rp, v0Index, o->verts[v0Index].r[ri]);
			BDTS_ASSERT(o->verts[v0Index].neighSurfLayer[elementLayer].neighSurfaceIn != -1, 0.0f);
			BDTS_ASSERT(o->verts[v0Index].neighSurfLayer[elementLayer].neighSurfaceOut != -1, 0.0f);
			float2 v0In = rcheck_replacePos(rp, o->verts[v0Index].neighSurfLayer[elementLayer].neighSurfaceIn, o->verts[o->verts[v0Index].neighSurfLayer[elementLayer].neighSurfaceIn].r[ri]);
			float2 v0Out = rcheck_replacePos(rp, o->verts[v0Index].neighSurfLayer[elementLayer].neighSurfaceOut, o->verts[o->verts[v0Index].neighSurfLayer[elementLayer].neighSurfaceOut].r[ri]);

			// point to surface radial energy
			float2 inVec(v0 - v0In);
			float2 outVec(v0Out - v0);
			float2 surfDirR = -Normalize(inVec.GetCross());
			float2 surfDirL = -Normalize(outVec.GetCross());

			float2 dirM = Normalize(surfDirR + surfDirL);
			float2 posM = v0;
			CononicalTrans2d ctrans(dirM, posM);

			surfDirR = ctrans.transDir(surfDirR);
			surfDirL = ctrans.transDir(surfDirL);

			float2 pV1 = ctrans.transPos(v1);
			float2 pV2 = ctrans.transPos(v2);

			// ok at this stage we have everything retransformed into the space of the points normal
			float2 evDirL = Normalize(pV1);
			float2 evDirR = Normalize(pV2);
			evDirL = horzNormDir(evDirL);
			evDirR = horzNormDir(evDirR);

			float2 cutDirL = ((float)evDirL.x > (float)surfDirL.x) ? evDirL : surfDirL;
			float2 cutDirR = ((float)evDirR.x > (float)surfDirR.x) ? surfDirR : evDirR;


			float avgArea = (Length(inVec) + Length(outVec)) * 0.5f;

			//now do EACH surface
			float runCutXMin = cutDirL.x;
			float runCutXMax = cutDirR.x;
			for (int i2 = 1; i2 < 3; i2++)
			{
				int oppCountIndex = (i2 + eint) % 3;
				int oppVertIndex = element.vertindex[oppCountIndex];

				int vertIndexSurf[3] = { o->verts[oppVertIndex].neighSurfLayer[elementLayer].neighSurfaceIn ,oppVertIndex ,o->verts[oppVertIndex].neighSurfLayer[elementLayer].neighSurfaceOut };
				for (int i3 = 0; i3 < 2; i3++)
				{
					if (vertIndexSurf[i3] == -1 || vertIndexSurf[i3 + 1] == -1) continue;

					float2 vi = rcheck_replacePos(rp, vertIndexSurf[i3], o->verts[vertIndexSurf[i3]].r[ri]);
					float2 vf = rcheck_replacePos(rp, vertIndexSurf[i3 + 1], o->verts[vertIndexSurf[i3 + 1]].r[ri]);
					float2 dirVi = Normalize(ctrans.transPos(vi));
					float2 dirVf = Normalize(ctrans.transPos(vf));



					if (dirVi.y <= 0.0f && dirVf.y <= 0.0f)
						continue; // completely behind

					dirVi = horzNormDir(dirVi);
					dirVf = horzNormDir(dirVf);
					float xi = dirVi.x;
					float xf = dirVf.x;



					if (xi >= xf)
						continue; // wrong pointing direction


					Line2d asLine(ctrans.transPos(vi), ctrans.transPos(vf));

					xi = bdts_min(runCutXMax, bdts_max(xi, runCutXMin));
					xf = bdts_min(runCutXMax, bdts_max(xf, runCutXMin));

					if (fabsf(xi - xf) < 0.00001f)
						continue; // wrong pointing direction

					runCutXMin = xf;

					bool isOkVal = std::isfinite(asLine.pi.x);
					BDTS_ASSERT(isOkVal, 0.0f);


					energy += integratePointEnergy(asLine, xi, xf, element.multEnergy);

				}
			}
		}

		return energy;

	}


	float EnergyForContactFrictionElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai, eEnergyComputePhase ePhase, float energySurface)
	{
		auto& element = o->elements[elementIndex];
		int firstSolid = GetSurfaceNeighIndexForContact(o, elementIndex);

		if (firstSolid == -1)
		{
			element.devEnergyPerArea = 0.0f;
			return 0.0f;// this contact has no surface
		}


		int sv0 = firstSolid;
		int sv1 = (firstSolid + 1) % 3;
		int sv2 = (firstSolid + 2) % 3;

#if 1 // temp disable friction
		float contactFriction = 0.0f;
		float minDistFriction = 0.0001f;
		auto& ematSurface = o->materials[o->elements[element.neighElem[firstSolid]].materialId];
		float fricMult = ematSurface.em_frictionSurface;
		float distNow = getElementSideDistance(o, rp, elementIndex, sv2, ri) - element.sideDistance[sv2];
		if (ePhase == eEnergyComputePhase::PREPERATION_COMPUTE_PHASE)
		{
			element.forceEstForFriction = estimateForceFromContactEnergy(energySurface, 1.0f, element.multEnergy);
		}

		if (element.forceEstForFriction > 0.00001)
		{

			contactFriction = fabsf(element.forceEstForFriction* distNow)*fricMult;// aka force * distance
		}
		else if (ePhase == eEnergyComputePhase::PREPERATION_COMPUTE_PHASE)
		{
			InitElementsSideDistances(o, elementIndex, ri);
		}


		if (ePhase == eEnergyComputePhase::LOSS_COMPUTE_PHASE && fabsf(distNow) > minDistFriction)
		{
#if BDTS_ENABLE_GLOBAL_ENERGY_TRACKING
			*o->energyFriction += contactFriction;
#endif

			element.sideDistance[sv2] = getElementSideDistance(o, rp, elementIndex, sv2, ri) - minDistFriction * (1 - 2 * signbit(distNow));
			float distNow2 = getElementSideDistance(o, rp, elementIndex, sv2, ri) - element.sideDistance[sv2];

			BDTS_ASSERT(fabsf(fabsf(distNow2) - minDistFriction) < 0.001f, 0.0f);
			contactFriction = fabsf(element.forceEstForFriction* distNow2)*fricMult;// force * distance
#if BDTS_ENABLE_GLOBAL_ENERGY_TRACKING
			*o->energyFriction -= contactFriction;
#endif
		}

		return contactFriction;
#else
		return 0.0f;
#endif
	}

	int GetSurfaceNeighIndexForContact(SimObj* o, int elementIndex)
	{

		const Element& element = o->elements[elementIndex];
		BDTS_ASSERT(element.neighElem[0] != -1, -1);
		BDTS_ASSERT(element.neighElem[1] != -1, -1);
		BDTS_ASSERT(element.neighElem[2] != -1, -1);
		int countNeighNonContact = 0;
		int firstSolid = -1;
		for (int i = 0; i < 3; i++)
		{
			if (element.neighElem[i] == -1) continue;

			countNeighNonContact += o->elements[element.neighElem[i]].elementType == Element::CONTACT_ELEMENT ? 0 : 1;
			if (countNeighNonContact > 0 && firstSolid == -1)
			{
				firstSolid = i;
			}
		}

		if (countNeighNonContact != 1)
		{
			return -1;
		}

		return firstSolid;
	}

	float EnergyForContactSurfaceElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai, eEnergyComputePhase ePhase)
	{
		float energy = 0.0f;
		auto& element = o->elements[elementIndex];
		int elementLayer = element.layerIndex;
		int firstSolid = GetSurfaceNeighIndexForContact(o, elementIndex);

		if (firstSolid == -1)
		{
			element.devEnergyPerArea = 0.0f;
			return energy;// this contact has no surface
		}

		int sv0 = firstSolid;
		int sv1 = (firstSolid + 1) % 3;
		int sv2 = (firstSolid + 2) % 3;

		float2 v0 = rcheck_replacePos(rp, element.vertindex[sv0], o->verts[element.vertindex[sv0]].r[ri]);
		float2 v1 = rcheck_replacePos(rp, element.vertindex[sv1], o->verts[element.vertindex[sv1]].r[ri]);
		float2 ss(v0 - v1); // super smooth
		float ssLen = Length(ss);
		float2 ssWindow(0.0f, ssLen);
		CononicalTrans2d conTrans(Normalize(ss.GetCross()), v1);

		int v2Index = element.vertindex[sv2];
		int currSurfVert = v2Index;

		// we step back twice to get to start 
		for (int i = 0; i < 2; i++)
		{
			BDTS_ASSERT(verts[currSurfVert].neighSurfLayer[elementLayer].neighSurfaceIn != -1, 0.0f);
			currSurfVert = o->verts[currSurfVert].neighSurfLayer[elementLayer].neighSurfaceIn;
		}

		for (int i = 0; i < 4; i++)
		{
			BDTS_ASSERT(verts[currSurfVert].neighSurfLayer[elementLayer].neighSurfaceOut != -1, 0.0f);
			int nextSurfIndex = o->verts[currSurfVert].neighSurfLayer[elementLayer].neighSurfaceOut;
			float2 vStart = rcheck_replacePos(rp, currSurfVert, o->verts[currSurfVert].r[ri]);
			float2 vEnd = rcheck_replacePos(rp, nextSurfIndex, o->verts[nextSurfIndex].r[ri]);
			Line2d testVec(vStart, vEnd);

			testVec = transLine(testVec, conTrans);
			auto windowLine = [&](Line2d &lineInOut)
			{
				float2 toWindow = float2(bdts_min(bdts_max(lineInOut.pi.x, ssWindow.x), ssWindow.y),
					bdts_min(bdts_max(lineInOut.pf.x, ssWindow.x), ssWindow.y));
				return toWindow;
			};


			auto addEnergy = [&](Line2d lineIn)
			{
				float2 winxixf = windowLine(lineIn);
				if ((winxixf.y - winxixf.x) > 0.00001f) // not flipped
				{
					float areaTotal = fabsf((winxixf.y - winxixf.x));
					energy += findEnergyIntegral(lineIn, winxixf, element.multEnergy);
				}
			};

			addEnergy(testVec);
			currSurfVert = nextSurfIndex;
		}

		return energy;
	}



	float EnergyForElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai, eEnergyComputePhase ePhase)
	{
		auto& element = o->elements[elementIndex];

		if (element.elementType == Element::CONTACT_ELEMENT)
		{
			return  EnergyForContactElement(o, rp, elementIndex, ri, ai, ePhase);
		}

		auto& emat = o->materials[element.materialId];

		float elementYoungsMod = emat.em_youngsMod;
		float elementPoissonRatio = emat.em_poissonRatio;

		float4 hvsa = rotFreeDeformForElement(o, rp, elementIndex, ri, ai);
		float matH = hvsa.x;
		float matV = hvsa.y;
		float matS = hvsa.z;
		float currArea = hvsa.w;
		//   [matH matS]
		//   [matS matV]    as deformation matrix without rotation

		float4 strainVec = { matH - 1.0f, matV - 1.0f, matS, matS };
		strainVec = strainVec - element.plasticStrainApplied; // remove the plastic strain 
		float4 hydroStrain = float4({ (strainVec.x + strainVec.y) * 0.5f, (strainVec.x + strainVec.y) * 0.5f, 0 , 0 });
		float hydroEnergyPerArea = UnitCellStrainEnergyRelation(hydroStrain, elementYoungsMod, elementPoissonRatio);
		float hydroEnergy = hydroEnergyPerArea * element.origArea;

		float4 devStrain = strainVec - hydroStrain; // remove the hydro (volume) strain
		element.strainVec = strainVec;
		element.strainDev = devStrain;
		float devEnergyPerArea = UnitCellStrainEnergyRelation(devStrain, elementYoungsMod, elementPoissonRatio);
		// plasticity 
		float4 maxStress = UnitCellFromStrainToStress(devStrain, elementYoungsMod, elementPoissonRatio);;// / (2.0f * G);
		float4 maxStrain = devStrain;
		float plasticEnergyPerArea = 0.0f;
		if (devEnergyPerArea > emat.em_maxElasticPerArea)
		{
			// elastic constraint to a constant
			// the shape of the constraint surface is basically an energy sphere for each strain/stress
			float sqVec = Dot(maxStress, maxStress);
			float G = elementYoungsMod / (2.0f* (1 + elementPoissonRatio));
			float t = sqrtf(emat.em_maxElasticPerArea * 2.0f * 2.0f * G / sqVec);
			t = bdts_min(1.0f, t);
			maxStress = t * maxStress; // stress is now constant (aka the forcec required to strain delta)
			maxStrain = UnitCellFromStressToStrain(maxStress, elementYoungsMod, elementPoissonRatio);
			devEnergyPerArea = UnitCellStrainEnergyRelation(maxStrain, elementYoungsMod, elementPoissonRatio);
			plasticEnergyPerArea = Dot(maxStress, devStrain - maxStrain);

			if (ePhase == eEnergyComputePhase::LOSS_COMPUTE_PHASE && plasticEnergyPerArea > 0.0f)
			{
				element.plasticStrainApplied = element.plasticStrainApplied + (devStrain - maxStrain);
#if BDTS_ENABLE_GLOBAL_ENERGY_TRACKING
				*o->energyPlastic += (plasticEnergyPerArea)*element.origArea;
#endif
			}
		}


		float devEnergy = devEnergyPerArea * element.origArea;
		float plasticEnergy = plasticEnergyPerArea * element.origArea;




		element.devEnergyPerArea = devEnergyPerArea;
		element.plasticEnergyPerArea = plasticEnergyPerArea;
		element.plasticEnergy = plasticEnergy;
		element.devEnergy = devEnergy;
		element.hydroEnergy = hydroEnergy;
		element.hydroEnergyPerArea = hydroEnergyPerArea;
		element.areaRatio = element.origArea / currArea;
		// we could use this for hydrostatic pressure?

		float compressiveEnergy = 0.0f;
		if (true)//element.origArea > currArea)
		{
			// this material is being compressed
			float compressionEnergyPerArea = abs(element.origArea - currArea) / element.origArea;
			element.compressionEnergyPerArea = compressionEnergyPerArea;
			compressiveEnergy = pow(compressionEnergyPerArea, 4)* sCompressiveMult * element.origArea;
			BDTS_ASSERT(compressiveEnergy >= 0.0f, 0.0f);
		}

		float sumEnergyOut = hydroEnergy + plasticEnergy + devEnergy + compressiveEnergy;


		//resitutive friciton
		float resitutitveMult = emat.em_frictionCof;
		float solidFrictionEnergy = 0.0f;
		float lossEnergy = 0.0f;
		solidFrictionEnergy = nced_evalEnergyState(&element.solidFrictionDissipator, ePhase, resitutitveMult, sumEnergyOut, lossEnergy);
#if BDTS_ENABLE_GLOBAL_ENERGY_TRACKING
		*o->energyPlastic += lossEnergy;
#endif
		sumEnergyOut += solidFrictionEnergy;



		return sumEnergyOut;
	}




	float EnergyForContactElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai, eEnergyComputePhase ePhase)
	{
		auto& element = o->elements[elementIndex];
		float energy = 0.0f;
		// surface to surface
		float energySurface = EnergyForContactSurfaceElement(o, rp, elementIndex, ri, ai, ePhase);

		// point to surface radial energy
		float energyPoint = EnergyForContactPointElement(o, rp, elementIndex, ri, ai, ePhase);

		// point to surface radial energy
		float energyFriction = EnergyForContactFrictionElement(o, rp, elementIndex, ri, ai, ePhase, energySurface);

		energy += energySurface;
		energy += energyPoint;

#if 1// resitutive friciton
		float resitutitveMult = 0.0f;
		{
			int sumWeight = 0;
			for (int i = 0; i < 3; i++)
			{
				if (o->elements[element.neighElem[i]].elementType != Element::CONTACT_ELEMENT)
				{
					int    neighMatIndex = o->elements[element.neighElem[i]].materialId;
					resitutitveMult += o->materials[neighMatIndex].em_frictionCof;
					sumWeight++;
				}
			}
			if (sumWeight > 1)
			{
				resitutitveMult = resitutitveMult / sumWeight; // rare contact requiring avg
			}
			if (resitutitveMult > 0.0001f)
			{

			}
			float lossEnergy = 0.0f;
			float solidFrictionEnergy = nced_evalEnergyState(&element.solidFrictionDissipator, ePhase, resitutitveMult, energy, lossEnergy);
#if BDTS_ENABLE_GLOBAL_ENERGY_TRACKING
			*o->energyPlastic += lossEnergy;
#endif

			energy += solidFrictionEnergy;
		}
#endif


		energy += energyFriction;
		element.devEnergyPerArea = energy;


		if (energy < 0.00001f)
		{
#if ADAPTIVE_ENERGY_FIX_ENABLED
			element.multEnergy = 1.0f;
#endif
		}
		return energy;
	}


	void InitElementsSideDistances(SimObj* o, int elementIndex, int ri)
	{
		Element& element = o->elements[elementIndex];
		for (int i = 0; i < 3; i++)
		{
			ReplacePos rp;
			rp.vertIndex = -1;// null replacement
			element.sideDistance[i] = getElementSideDistance(o, &rp, elementIndex, i, ri);
		}
	}


	float4 GetShapeSkewedMat(SimObj* o, int elementIndex, const float4& currShapedParams)
	{
		const Element& element = o->elements[elementIndex];
		float baseLength = currShapedParams.x;
		float barPerpCross = currShapedParams.y;
		float barPerpBase = currShapedParams.z;
		float currArea = (baseLength * barPerpCross)*0.5f;
		float scaleY = barPerpCross / element.barCrossLengthOrig;
		float scaleX = baseLength / element.baseLengthOrig;
		float scaleShear = (barPerpBase - element.barBaseLengthOrig*scaleX) / element.barCrossLengthOrig; // dx / tall-length oddly enough

		//  DEBUG verifies that this transformation will undo the deformation
		float  testBaseLength = element.baseLengthOrig * scaleX;
		float  testbarPerpCross = element.barCrossLengthOrig * scaleY;
		float  testbarPerpBase = element.barBaseLengthOrig * scaleX + element.barCrossLengthOrig  * scaleShear;

		// outputs
		// [ scaleX ScaleShear]
		// [    0      scaley ]
		return { scaleX, scaleY, scaleShear, currArea };
	}

	void ClearForce(SimObj* o, int ai)
	{
		for (int i = 0 ; i < o->vertsCount ; i++)
		{
			o->verts[i].a[ai] = float2();
		}
	}


	void ForceForElement(SimObj* o, int elementIndex, int ri, int vi, int ai)
	{
		auto& element = o->elements[elementIndex];
		static const float deltaSpatial = 0.0000331217f;
		int elementLayer = element.layerIndex;

		auto pdEnergy = [&](Vert& vertInOut, int vertIndex)
		{
			ReplacePos rp;
			rp.vertIndex = vertIndex;

			rp.pos = vertInOut.r[ri].v;
			auto& vert0 = rp.pos;
			float2 vertSaved = vertInOut.r[ri].v;

			vert0.x = (vertSaved.x + deltaSpatial);
			float energyPosX = EnergyForElement(o, &rp, elementIndex, ri, ai, eEnergyComputePhase::VARIATION_COMPUTE_PHASE);
			vert0 = vertSaved;
			vert0.x = (vertSaved.x - deltaSpatial);
			float energyNegX = EnergyForElement(o, &rp, elementIndex, ri, ai, eEnergyComputePhase::VARIATION_COMPUTE_PHASE);
			vert0 = vertSaved;

			vert0.y = (vertSaved.y + deltaSpatial);
			float energyPosY = EnergyForElement(o, &rp, elementIndex, ri, ai, eEnergyComputePhase::VARIATION_COMPUTE_PHASE);
			vert0 = vertSaved;
			vert0.y = (vertSaved.y - deltaSpatial);
			float energyNegY = EnergyForElement(o, &rp, elementIndex, ri, ai, eEnergyComputePhase::VARIATION_COMPUTE_PHASE);
			vert0 = vertSaved;



			vertInOut.a[ai] -= float2(0.5f* (energyPosX - energyNegX) / deltaSpatial, 0.5f*(energyPosY - energyNegY) / deltaSpatial);
		};


		tSetObj diffVertIndex;
		microset_init(&diffVertIndex);
		for (int i = 0; i < 3; i++)
		{
			pdEnergy(o->verts[element.vertindex[i]], element.vertindex[i]);
			microset_insert(&diffVertIndex, element.vertindex[i]);
		}

		if (element.elementType == Element::CONTACT_ELEMENT)
		{
			for (int i = 0; i < 3; i++)
			{
				auto& origVert = o->verts[element.vertindex[i]];

				if (origVert.neighSurfLayer[elementLayer].neighSurfaceIn != -1 &&
					microset_insert(&diffVertIndex, origVert.neighSurfLayer[elementLayer].neighSurfaceIn))
				{
					pdEnergy(o->verts[origVert.neighSurfLayer[elementLayer].neighSurfaceIn], origVert.neighSurfLayer[elementLayer].neighSurfaceIn);
				}

				if (origVert.neighSurfLayer[elementLayer].neighSurfaceOut != -1 &&
					microset_insert(&diffVertIndex, origVert.neighSurfLayer[elementLayer].neighSurfaceOut))
				{
					pdEnergy(o->verts[origVert.neighSurfLayer[elementLayer].neighSurfaceOut], origVert.neighSurfLayer[elementLayer].neighSurfaceOut);
				}
			}
		}
	}



#if BDTS_IS_CPP_COMPILE 
}
#endif