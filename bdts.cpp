
#include "bdts.h"
#include <Engine/SystemRendering2d.h>
#include <Engine/SystemInput.h>
#include <Engine/SystemRendering.h>
#include <Core/Utility.h>
#include <Core/Math/Random.h>
#include <algorithm> // clamp
#include <btds_math.h>
#include "btds_sim_fun.h"
#include "btds_sim_phys.h"
#define _USE_MATH_DEFINES

#include <math.h>


#define STRINGIFY(x) #x
#define AT  __FILE__  ":" STRINGIFY( __LINE__ );

#define BDTS_ASSERT( expression, rtnValue) \
if( !(expression)) { \
std::string assertFullName = AT ;\
if(!assertEnabled[ assertFullName ]) { assertEnabled[ assertFullName ] = true;     __debugbreak(); return rtnValue;   }}



#define ADAPTIVE_ENERGY_FIX_ENABLED 1

namespace bdts
{




	static 	Line2d transLine(const Line2d& line, const CononicalTrans2d& trans)
	{
		float2 pi = trans.transPos(line.pi);
		float2 pf = trans.transPos(line.pf);
		return Line2d(pi, pf);
	}



	TriangleSimulation::TriangleSimulation()
	{
		testIntegration();
	}

	void TriangleSimulation::createMeshGrid( int matIndex,  float2 bottom, float2 top, int dimx, int dimy)
	{
		GenerateBoxPoint(matIndex,dimx, dimy, bottom, top);
	}


	void TriangleSimulation::createMeshTube(bool bCentered, int  matIndex, float2 center, float radiusIn, float radiusOut, int segRad, int segScale, bool startNoAvec, bool geared, float gearRad, float rotOffset)
	{
		float rcurr = radiusIn;
		float rDelta = (radiusOut - radiusIn) / segScale;

		float thedaDelta = (M_PI*2.0f) / segRad;

		int centerVert = -1;
		if (bCentered)
		{
			centerVert = verts.size();
			verts.push_back({ center });
		}

		int vertNew = -1;
		int vertOld = -1;
		for (int i = 0; i < segScale+1; i++)
		{
			vertOld = vertNew;
			vertNew = (int)verts.size();
			for (int j = 0; j < segRad; j++)
			{
				float offsetForStable =  -i* thedaDelta/2.0f;
				float theda = (thedaDelta * j) + offsetForStable + rotOffset;
				float x = cos(theda) * rcurr;
				float y = sin(theda) * rcurr;

				if (geared && i == segScale)
				{
					if ((j % 4) < 2)
					{
						x = cos(theda) * gearRad;
						y = sin(theda) * gearRad;
					}
				}

				verts.push_back(center + float2( x,y ));
				if (startNoAvec && i == 0)
				{
					verts[vertNew + j].doNotAvect = true;
				}

			
			}

			rcurr += rDelta;

			if (vertOld != -1)
			{
				for (int j = 0; j < segRad; j++)
				{
					elements.push_back({ j + vertOld,  vertNew + ((j + 1) % segRad), vertNew + j , matIndex });
					elements.push_back({ j + vertOld,  vertOld + ((j + 1) % segRad), vertNew + ((j+1)% segRad),matIndex });
				}
			}
			else if (centerVert != -1)
			{
				for (int j = 0; j < segRad; j++)
				{
					elements.push_back({ centerVert,  vertNew + ((j + 1) % segRad), vertNew + j , matIndex });
				}
			}
		}
	}

	void TriangleSimulation::setSimulationBounds(float2 minCorner, float2 maxCorner)
	{
		float2 bottom = minCorner;
		float2 top = maxCorner;
		setLimits(bottom, top);
	}

	void TriangleSimulation::addBorder(int matIndex, int   dimX, int dimY, float thickness, int layer)
	{
		// create a border of triangles with the dims above
		std::vector<int> outerRing;
		std::vector<int> innerRing;
		float espilon = 0.001f;
		float2 outerMin = bottom - float2(espilon);
		float2 outerMax = top + float2(espilon);

		float2 outerStep = (outerMax - outerMin) / float2(dimX, dimY);

		float2 innerMin = outerMin + float2(thickness);
		float2 innerMax = outerMax - float2(thickness);

		float2 innerStep = (innerMax - innerMin) / float2(dimX - 1, dimY - 1);
		float2 currOuterPos = outerMin;
		float2 currInnerPos = innerMin;

		float2 currOuterStep = float2(outerStep.x, 0);
		float2 currInnerStep = float2(innerStep.x, 0);

		static const bool innerIsNoAvec = true;

		bool outerIsSurface = true;// flips;

		int prevOuterVert = verts.size();
		verts.push_back(currOuterPos);
		verts[prevOuterVert].doNotAvect = true;
		int prevInnerVert = verts.size();
		verts.push_back(currInnerPos);
		verts[prevInnerVert].doNotAvect = innerIsNoAvec;
		
		int startInnerVert = prevInnerVert;
		int startOuterVert = prevOuterVert;

		auto forwardBorder = [&]()
		{
			if (outerIsSurface)
			{
				currOuterPos += currOuterStep;
				int newOuterVert = verts.size();
				verts.push_back(currOuterPos);
				elements.push_back({ prevOuterVert, prevInnerVert, newOuterVert, matIndex , layer});
				prevOuterVert = newOuterVert;
				verts[prevOuterVert].doNotAvect = true;
			}
			else
			{
				currInnerPos += currInnerStep;
				int newInnerVert = verts.size();
				verts.push_back(currInnerPos);
				verts[newInnerVert].doNotAvect = innerIsNoAvec;
				elements.push_back({ prevInnerVert, newInnerVert, prevOuterVert, matIndex, layer});
				prevInnerVert = newInnerVert;
			}

			outerIsSurface = !outerIsSurface;
		};

		for (int i = 0; i < (dimX * 2) - 1; i++)
		{
			forwardBorder();
		}


		currOuterStep = float2(0, outerStep.y);
		currInnerStep = float2(0, innerStep.y);
		outerIsSurface = true;
		for (int i = 0; i < (dimY * 2) - 1; i++)
		{
			forwardBorder();
		}


		currOuterStep = float2(-outerStep.x, 0);
		currInnerStep = float2(-innerStep.x, 0);
		outerIsSurface = true;
		for (int i = 0; i < (dimX * 2) - 1; i++)
		{
			forwardBorder();
		}

		currOuterStep = float2(0, -outerStep.y);
		currInnerStep = float2(0, -innerStep.y);
		outerIsSurface = true;
		for (int i = 0; i < (dimY * 2) - 3; i++)
		{
			forwardBorder();
		}

		elements.push_back({ prevInnerVert, startInnerVert, prevOuterVert, matIndex , layer});
		elements.push_back({ prevOuterVert, startInnerVert, startOuterVert, matIndex , layer });

	}

	void TriangleSimulation::fillSimulation()
	{
		weldClose(0.0001f); // insures correct orientation as well
		for (int i = 0; i < sNumLayer; i++)
		{
			findSurface(i);
			while (connectOneSurface(i))
			{
				weldClose(0.0f);
			}
			conOpen.clear();
		}
		return;
	}

	void TriangleSimulation::deleteElement(int elementIndex)
	{
		for (auto&& each : elements)
		{
			for (int i = 0; i < 3; i++)
			{
				if (each.neighElem[i] == elementIndex)
				{
					each.neighElem[i] = -1;
				}
				else if (each.neighElem[i] > elementIndex)
				{
					each.neighElem[i] = each.neighElem[i] - 1;
				}
			}
		}

		elements.erase(elements.begin() + elementIndex);

		// find all the verts that have nothing pointing to them
		std::set<int> vertsToDelete;

		for (int vertIndex = 0; vertIndex < verts.size(); vertIndex++)
		{
			bool foundVert = false;
			for (auto&& each : elements)
			{
				if (each.vertindex[0] == vertIndex || each.vertindex[1] == vertIndex || each.vertindex[2] == vertIndex)
				{
					foundVert = true;
				}
			}

			if(!foundVert)
				vertsToDelete.insert(vertIndex);
		}

		if (vertsToDelete.size() == 0)
			return;

		std::vector<Vert> oldVerts = verts;// full copy very slow
		verts.clear();// bye bye
		int numLayedVerts = 0;
		std::map<int, int> remapVerts; // old to new
		for (int vertIndex = 0; vertIndex < oldVerts.size(); vertIndex++)
		{
			if (vertsToDelete.find(vertIndex) == vertsToDelete.end())
			{
				verts.push_back(oldVerts[vertIndex]);
				remapVerts[vertIndex] = numLayedVerts;
				numLayedVerts++;
			}
		}

		for (auto&& each : elements)
		{
			for (auto&&eachVert : each.vertindex)
			{
				eachVert = remapVerts[eachVert];// remap
			}
		}
	}




	void TriangleSimulation::GenerateBoxPoint(int matIndex, int countX, int countY, float2 low, float2 high)
	{
		float2 scaleStep = (high - low) / float2(countX, countY);

		std::vector<int> prevBox;
		prevBox.resize(countX + 1);
		for (int i = 0; i < countX + 1; i++)
		{
			prevBox[i] = verts.size();
			verts.push_back(low + float2(scaleStep.x * (float) i, 0) );
		}

	
		for (int j = 0; j < countY; j++)
		{
			std::vector<int> farBox;
			float farY = scaleStep.y * ( (float) (1 + j));
			if (true)
			{
				farBox.resize(countX + 1);
				for (int i = 0; i < countX + 1; i++)
				{
					farBox[i] = verts.size();
					verts.push_back(low + float2(scaleStep.x * ((float)i)  , farY));
				}
			}

			for (int i = 0; i < countX; i++)
			{
				// rows
				float2 midPointLocation = low + (scaleStep * float2(i, j)) + scaleStep / 2.0f;
#if 0// random midpoint location
				midPointLocation += (Random::GetVec2() - float2(0.5f, 0.5f)) * scaleStep * 0.5f;
#endif 
				verts.push_back(midPointLocation);
				int midPoint = verts.size() - 1;

				int lowXLowY = prevBox[i];
				int lowXhighY = farBox[i];
				int highXLowY = prevBox[i+1];
				int highXhighY = farBox[i + 1];
				elements.push_back({ midPoint, lowXLowY ,highXLowY , matIndex   });
				elements.push_back({ midPoint, lowXLowY , lowXhighY , matIndex });
				elements.push_back({ midPoint,  highXhighY , highXLowY, matIndex });
				elements.push_back({ midPoint , highXhighY ,lowXhighY , matIndex });
			}

			prevBox = farBox;
		}
	}

	void TriangleSimulation::assignNeighElements()
	{
		
		for (int i=0; i< elements.size(); i++)
		{
			auto& elemi = elements[i];
			for (int j = 0; j< elements.size(); j++)
			{
				auto& elemj = elements[j];
				if (i == j) continue;// no self
				if( elemi.layerIndex != elemj.layerIndex) continue;// no between layers
				// the logic here is that if all triangles have the same handedness 
				// one will have forward the other will have backward for adjoining
				for (int qi = 0; qi < 3; qi++)
				{
					int qi_forward = (qi + 1) % 3;
					for (int qj = 0; qj < 3; qj++)
					{
						int qj_forward = (qj + 1) % 3;


						if (elemi.vertindex[qi] == elemj.vertindex[qj_forward] &&
							elemi.vertindex[qi_forward] == elemj.vertindex[qj])
						{
							elemi.neighElem[qi] = j;
							elemj.neighElem[qj] = i;
						}
					}
				}


			}
		}

	}

	float getElementSideDistance(SimObj* o, ReplacePos* rp, int elementIndex, int triIndex, int ri)
	{
		Element& element = o->elements[elementIndex];
		float2 vert0 = rcheck_replacePos(rp, element.vertindex[triIndex], o->verts[element.vertindex[triIndex]].r[ri]);
		float2 vert1 = rcheck_replacePos(rp, element.vertindex[(triIndex + 1) % 3], o->verts[element.vertindex[(triIndex + 1) % 3]].r[ri]);
		float2 vert2 = rcheck_replacePos(rp, element.vertindex[(triIndex + 2) % 3], o->verts[element.vertindex[(triIndex + 2) % 3]].r[ri]);

		auto center = (vert1 + vert2) * 0.5f;
		auto ray = (vert1 - vert2);

		Line2d asLine(vert1, vert2);

		return asLine.projPointLineFloat(vert0);
	}



	Matrix22 TriangleSimulation::strainToDefMatrix(float4 strain)
	{
		Matrix22 mat(float4(strain.x, strain.z, strain.w, strain.y) +
			float4(1.0f, 0.0f, 0.0f, 1.0f));
		return mat;
	}

	SimObj TriangleSimulation::getSimObj()
	{
		SimObj o;
		o.elements = &elements.front();
		o.verts = &verts.front();
		o.vertsCount = verts.size();
		o.materials = &eMaterials.front();
#if BDTS_ENABLE_GLOBAL_ENERGY_TRACKING
		o.energyFriction = &energyFriction;
		o.energyOffset = &energyOffset;
		o.energyPlastic = &energyPlastic;
		o.simTime = &simTime;
#endif
		return o;
	}

	void TriangleSimulation::VertMassChangeForElement(Element& element, float density, float currArea)
	{
		if (element.elementType == Element::CONTACT_ELEMENT)
		{
			return;
		}
		for (int i = 0; i < 3; i++)
		{
			verts[element.vertindex[i]].mMassCount += (density * currArea) / 3.0f;
		}
	}

	void TriangleSimulation::InitElementsAsCurrent( int elementIndex, TriangleSimulation::tElementUnitAreaProp* initProps)
	{
		Element& element = elements[elementIndex];
		int ri = 0;
		float2 vert0 = verts[element.vertindex[0]].r[ri].v;// copies
		float2 vert1 = verts[element.vertindex[1]].r[ri].v;
		float2 vert2 = verts[element.vertindex[2]].r[ri].v;

		if (initProps)
		{
			//float4 strainVec = { matH - 1.0f, matV - 1.0f, matS, matS };
			Matrix22 whatAbout = rotRemoval(initProps->worldDefMat);
			initProps->worldDefMat = whatAbout * initProps->worldDefMat;
			vert0 = initProps->worldDefMat * vert0;
			vert1 = initProps->worldDefMat * vert1;
			vert2 = initProps->worldDefMat * vert2;
		}
		auto shapeVec = GetTriShapedParams(vert0, vert1, vert2);

		float barLength = shapeVec.x;
		float barPerpCross = shapeVec.y;
		float barPerpBase = shapeVec.z;

		if (false && initProps)
		{
			BDTS_ASSERT(fabsf(element.baseLengthOrig - barLength) < 0.001f,  );
			BDTS_ASSERT(fabsf(element.barCrossLengthOrig - barPerpCross) < 0.001f, );
			BDTS_ASSERT(fabsf(element.barBaseLengthOrig - barPerpBase) < 0.001f, );
		}

		element.baseLengthOrig = barLength;
		element.barCrossLengthOrig = barPerpCross;
		element.barBaseLengthOrig = barPerpBase;
		element.plasticStrainApplied = float4();
		float currArea = barLength * barPerpCross * 0.5f;
		element.origArea = currArea; // shouldnt matter if its regular or bent;

		float tempDensity = 1.0f;
		if (element.elementType == Element::CONTACT_ELEMENT)
		{
			tempDensity = 0.0f;
			SimObj o = getSimObj();
			InitElementsSideDistances(&o, elementIndex, ri);
		}
		
		if (initProps)
		{
			//tempDensity = 0.0f;
		}

		VertMassChangeForElement(element, tempDensity, currArea);
	}

	void TriangleSimulation::InitElements()
	{
		// lets go through all the elements and assign the correct values 
		// as our current positions are the relaxed ones
		assignNeighElements();

		for (int  elementIndex = 0 ; elementIndex < elements.size(); elementIndex++)
		{
 			InitElementsAsCurrent(elementIndex);
		}

		assignVertexSurfaces();
	}



	// these are all element properties NOT vertex
	//const float youngsMod_old = 10000.0101f;

	//const float maxElasticEnergy_old = 100000000.1f;
	//const float poissonRatio = 0.33333f;
	//const float G = youngsMod / (2.0f* (1 + poissonRatio));




	void TriangleSimulation::testInitElement(int elemIndex)
	{
		auto& elemi = elements[elemIndex];
		//auto initData = GetElmentUnitAreaProp(elemIndex);
		//InitElementsAsCurrent(elemi, &initData);
		adaptiveOne(elemi, elemIndex);
	}

	bool TriangleSimulation::refineCentered(int elemIndex)
	{
		auto& elemi = elements[elemIndex];

		//if (elemi.elementType != Element::CONTACT_ELEMENT )
		//	return false;
		SimObj o = getSimObj();
		auto repTri = toTriangle(&o, elemi.vertindex);
		float2 centerp = repTri.getCentroid();

		int newVertIndex = verts.size();
		verts.push_back(centerp);

		auto& newVert = verts[newVertIndex];
		newVert = verts[elemi.vertindex[0]];
		newVert.r[0].v = centerp;
		newVert.r[1].v = centerp;

		int n[3];
		int v[3];
		for (int i = 0; i < 3; i++)
		{
			v[i] = elemi.vertindex[i];
			n[i] = elemi.neighElem[i];
		}

		int e0Index = elemIndex;
		int e1Index = elements.size();
		elements.push_back({ -1,-1,-1 });
		int e2Index = elements.size();
		elements.push_back({ -1,-1,-1 });

		auto& e0 = elements[e0Index];
		auto& e1 = elements[e1Index];
		auto& e2 = elements[e2Index];


		e1.elementType = e0.elementType;
		e2.elementType = e0.elementType;

		e0.vertindex[0] = v[0];
		e0.vertindex[1] = v[1];
		e0.vertindex[2] = newVertIndex;
		e0.neighElem[0] = n[0];
		e0.neighElem[1] = e1Index;
		e0.neighElem[2] = e2Index;
	

		e1.vertindex[0] = v[1];
		e1.vertindex[1] = v[2];
		e1.vertindex[2] = newVertIndex;
		e1.neighElem[0] = n[1];
		e1.neighElem[1] = e2Index;
		e1.neighElem[2] = e0Index;
		e1.materialId = e0.materialId;

		e2.vertindex[0] = v[2];
		e2.vertindex[1] = v[0];
		e2.vertindex[2] = newVertIndex;
		e2.neighElem[0] = n[2];
		e2.neighElem[1] = e0Index;
		e2.neighElem[2] = e1Index;
		e2.materialId = e0.materialId;

		auto replaceSelfNeigh = [&](int indexOrig, int indexNew, int neighIndex)
		{
			if (neighIndex == -1) return;

			for (int i = 0; i < 3; i++)
			{
				if (elements[neighIndex].neighElem[i] == indexOrig)
				{
					elements[neighIndex].neighElem[i] = indexNew;
					return;
				}
			}
		};


		replaceSelfNeigh(elemIndex, e1Index, n[1]);
		replaceSelfNeigh(elemIndex, e2Index, n[2]);
	
		
		return true;
	}


	bool TriangleSimulation::refineBisplit(int elemIndexi, int elemIndexj, int neighVertIndexi, int neighVertIndexj)
	{
		// copies!!
		auto elemA = elements[elemIndexi];
		auto elemB = elements[elemIndexj];

		//if (elemA.elementType == Element::CONTACT_ELEMENT && elemB.elementType == Element::CONTACT_ELEMENT)
			//return false;
		SimObj o = getSimObj();
		auto elemAInitParam = GetElmentUnitAreaProp( elemIndexi);
		auto elemBInitParam = GetElmentUnitAreaProp( elemIndexj);

		// remove mass
		VertMassChangeForElement(elemA, -1.0f, elemA.origArea);
		VertMassChangeForElement(elemB, -1.0f, elemB.origArea);

		int oa0 = neighVertIndexi;
		int oa1 = (neighVertIndexi + 1) % 3;
		int oa2 = (neighVertIndexi + 2) % 3;

		int ob0 = neighVertIndexj;
		int ob1 = (neighVertIndexj + 1) % 3;
		int ob2 = (neighVertIndexj + 2) % 3;

		BDTS_ASSERT(elemA.vertindex[oa0] == elemB.vertindex[ob1], false);
		BDTS_ASSERT(elemA.vertindex[oa1] == elemB.vertindex[ob0], false);


		int xv = verts.size();
		float2 newSplitLocation = (verts[elemA.vertindex[oa0]].r[0].v + verts[elemB.vertindex[ob0]].r[0].v)* 0.5f;
		verts.push_back(newSplitLocation);

		verts[xv].v[0] = (verts[elemA.vertindex[oa0]].v[0] + verts[elemB.vertindex[ob0]].v[0])* 0.5f;


		int axIndex = elemIndexi;
		int ayIndex = elements.size();
		elements.push_back(elemA);
		int bxIndex = elemIndexj;
		int byIndex = elements.size();
		elements.push_back(elemB);

		auto& elemAX = elements[axIndex];
		auto& elemAY = elements[ayIndex];
		auto& elemBX = elements[bxIndex];
		auto& elemBY = elements[byIndex];

		// defined in "refine_split.png"
		elemAX.vertindex[oa1] = xv;
		elemAX.neighElem[oa0] = bxIndex;
		elemAX.neighElem[oa1] = ayIndex;

		elemAY.vertindex[oa0] = xv;
		elemAY.neighElem[oa0] = byIndex;
		elemAY.neighElem[oa2] = axIndex;


		elemBX.vertindex[ob0] = xv;
		elemBX.neighElem[ob0] = axIndex;
		elemBX.neighElem[ob2] = byIndex;

		elemBY.vertindex[ob1] = xv;
		elemBY.neighElem[ob0] = ayIndex;
		elemBY.neighElem[ob1] = bxIndex;

		auto replaceSelfNeigh = [&](int indexOrig, int indexNew, int neighIndex)
		{
			if (neighIndex == -1) return;

			for (int i = 0; i < 3; i++)
			{
				if (elements[neighIndex].neighElem[i] == indexOrig)
				{
					elements[neighIndex].neighElem[i] = indexNew;
					return;
				}
			}
		};

		replaceSelfNeigh(axIndex, ayIndex, elemA.neighElem[oa1]);
		replaceSelfNeigh(bxIndex, byIndex, elemB.neighElem[ob2]);

		InitElementsAsCurrent(axIndex, &elemAInitParam);
		InitElementsAsCurrent(ayIndex, &elemAInitParam);
		InitElementsAsCurrent(bxIndex, &elemBInitParam);
		InitElementsAsCurrent(byIndex, &elemBInitParam);
		return true;
	}




	void TriangleSimulation::adaptiveOne(Element& selElem, int elemIndex)
	{
		SimObj o = getSimObj();

		for (int selfNeighIndex = 0; selfNeighIndex < 3; selfNeighIndex++)
		{
			int firstNeighIndex = selElem.neighElem[selfNeighIndex];
			if (firstNeighIndex == -1)
				continue; // outer element

			auto& neighElem = elements[firstNeighIndex];

			int selfInNeighIndex;
			for (selfInNeighIndex = 0; selfInNeighIndex < 3; selfInNeighIndex++)
			{
				if (neighElem.neighElem[selfInNeighIndex] == elemIndex)
					break;
			}

			BDTS_ASSERT(selfInNeighIndex != 3, );

			if (selfInNeighIndex != 3)
			{
				if (tryAdaptive(&o, elemIndex, firstNeighIndex, selfNeighIndex, selfInNeighIndex))
				{
					
				}
			}
		}
	}
	void TriangleSimulation::adaptiveAll()
	{
#if 1// test for interlayer connection that is not valid
		for (auto&& each: elements)
		{
			for (int i = 0; i < 3; i++)
			{
				if (each.neighElem[i] != -1)
				{
					auto&& neighElem = elements[each.neighElem[i]];
					BDTS_ASSERT(each.layerIndex == neighElem.layerIndex, ); // dont point between layers
				}
			}
		}
#endif


		for (int elemIndex = 0; elemIndex < elements.size(); elemIndex++)
		{
			auto& selElem = elements[elemIndex];
			adaptiveOne(selElem, elemIndex);
		}

	}







	float4 rotFreeDeformForElement(SimObj* o, ReplacePos* rp, int elementIndex, int ri, int ai)
	{
		auto& element = o->elements[elementIndex];
		float2 vert0 = rcheck_replacePos(rp, element.vertindex[0], o->verts[element.vertindex[0]].r[ri]);
		float2 vert1 = rcheck_replacePos(rp, element.vertindex[1], o->verts[element.vertindex[1]].r[ri]);
		float2 vert2 = rcheck_replacePos(rp, element.vertindex[2], o->verts[element.vertindex[2]].r[ri]);
	
		float4 shapedParams = GetTriShapedParams(vert0, vert1, vert2);
		float4 matSkew = GetShapeSkewedMat(o, elementIndex, shapedParams);
		
		float scaleX = matSkew.x;
		float scaleY = matSkew.y;
		float strainShear = matSkew.z;
		float currArea = matSkew.w;
		// remove rotation																				  
		float asX = strainShear / (scaleX + scaleY);
		float phi_debug_only = atanf(asX);
		float cosAtan = 1 / sqrtf(asX*asX + 1.0f);  //cos(atan(x)) from wolfram
		float sinAtan = asX / sqrtf(asX*asX + 1.0f);  //sin(atan(x)) from wolfram
		float matH = cosAtan * scaleX;
		float matS = sinAtan * scaleX;
		float matV = sinAtan * strainShear + cosAtan * scaleY;

		return { matH, matV, matS, currArea };
	}


	




	TriangleSimulation::tElementUnitAreaProp TriangleSimulation::GetElmentUnitAreaProp( int elementIndex)
	{
		auto& element = elements[elementIndex];
		tElementUnitAreaProp out;
		float currArea = 0.0f;
		float2 vert0 = verts[element.vertindex[0]].r[0].v;
		float2 vert1 = verts[element.vertindex[1]].r[0].v;
		float2 vert2 = verts[element.vertindex[2]].r[0].v;
		for (int i = 0; i < 3; i++)
		{
			out.simVert[i] = verts[element.vertindex[i]].r[0].v;
		}
		
		SimObj o = getSimObj();
		float4 shapedParams = GetTriShapedParams(vert0, vert1, vert2);
		float4 matSkew = GetShapeSkewedMat(&o, elementIndex, shapedParams);
		out.matSkew = matSkew;
		Matrix22 defMat = Matrix22({matSkew.x, matSkew.z, 0.0f, matSkew.y});
		Matrix22 strainMat = defMat - Matrix22({ 1,0,0,1 });
		float2 triXDir = Normalize(vert1 - vert0);
		float2 triYDir = triXDir.GetCross() * shapedParams.w;// we might have mirrored y axis 

		Matrix22 rot = Matrix22({ triXDir.x, triYDir.x, 
								triXDir.y, triYDir.y });
		Matrix22 invRot = Invert(rot);
		Matrix22 invDef = Invert(defMat);
		out.invRot = invRot;
		out.rot = rot;
		out.invDef = invDef;
		out.worldDefMat = rot *  invDef * invRot; // go into element space then undo then rotate out again.
		Matrix22 rotRemovalMat = rotRemoval(out.worldDefMat);
		out.worldDefMat = rotRemovalMat * out.worldDefMat;

		out.expandedVert[0] = out.worldDefMat * vert0;
		out.expandedVert[1] = out.worldDefMat * vert1;
		out.expandedVert[2] = out.worldDefMat * vert2;

		return out;
	}


	
	

	void TriangleSimulation::DrawingOutput(int viewLayer, std::vector<float4>& posOut, std::vector<float4>& colOut, bool useUnique, bool useMaterial, float4 solidCol, float4 contactCol)
	{
		int element = 0;
		for (auto&& each : elements)
		{
			if(viewLayer != -1 && each.layerIndex != viewLayer) continue;

			for (int i = 0; i < 3; i++)
			{
				int v0 = each.vertindex[i];
				posOut.push_back({ verts[v0].r[0].v.x, verts[v0].r[1].v.y, 0.0f , 1.0f });
				float4 tempCol =   { Utility::rndFromIndex(element, 0)*0.5f + 0.5f, Utility::rndFromIndex(element, 1), Utility::rndFromIndex(element, 2), 0.5f };
				if (!useUnique)
				{
					tempCol = contactCol;
				}

			

				if (each.elementType == Element::SOLID_ELEMENT)
				{
					tempCol = tempCol+ float4({ 1.0,0,0.0f,0.5f });
					if (!useUnique)
					{
						tempCol = solidCol;
					}
					//
				}

				if (useMaterial )
				{
					if (each.materialId == -1)
					{
						tempCol = float4(0.1f, 0.1f, 0.1f, 0.1f);
					}
					else
					{
						tempCol = eMaterials[each.materialId].debugColor;
					}
				}

				//tempCol = tempCol * 0.6f;
				colOut.push_back(tempCol);
			}
			element++;
		}
	}


	void TriangleSimulation::setLimits(float2 bottomv, float2 topv)
	{
		bottom = bottomv;
		top = topv;
	}

	void  TriangleSimulation::findSurface(int forLayer)
	{
		conOpen.clear();
		for (int i = 0; i < elements.size(); i++)
		{
			auto& elemi = elements[i];
			if (elemi.layerIndex != forLayer) continue;; // only looking at one layer at a time

			for (int k = 0; k < 3; k++)
			{
				int vertA = elemi.vertindex[k];
				int vertB = elemi.vertindex[(k + 1) % 3];
				bool found = false;
				for (int j = 0; j < elements.size(); j++)
				{
					if (i == j) continue;

					auto& elemj = elements[j];

					bool hasVertA = elemj.vertindex[0] == vertA || elemj.vertindex[1] == vertA || elemj.vertindex[2] == vertA;
					bool hasVertB = elemj.vertindex[0] == vertB || elemj.vertindex[1] == vertB || elemj.vertindex[2] == vertB;
					if (hasVertB && hasVertA && elemi.layerIndex == elemj.layerIndex)
					{
						found = true;
						break;
					}
				}

				if (!found)
				{
					conOpen.push_back({ vertA, vertB, i });
				}
			}
		}

		
	}


	void TriangleSimulation::weldClose(float maxDist)
	{
		// let make sure these are all clockwise
		for (auto&& each : elements)
		{
			Line2d side = { verts[each.vertindex[0]].r[0].v , verts[each.vertindex[1]].r[0].v };
			if (side.isUpRelative(verts[each.vertindex[2]].r[0].v))
			{
				// wrong sided . switch
				int temp = each.vertindex[0];
				each.vertindex[0] = each.vertindex[1];
				each.vertindex[1] = temp;
			}
		}

		std::map<int, int> replacedToNewMap;
		for (int vertIndex = 0; vertIndex < verts.size(); vertIndex++)
		{
			if (replacedToNewMap.find(vertIndex) != replacedToNewMap.end())
			{
				continue;// its already replaced dont try to replace others
			}

			for (int otherIndex = 0; otherIndex < verts.size(); otherIndex++)
			{
				if (otherIndex == vertIndex) continue;// lol dont self replace

				if (Length(verts[vertIndex].r[0].v - verts[otherIndex].r[0].v) <= maxDist)
				{
					replacedToNewMap[otherIndex] = vertIndex;
				}
			}
		}



		if (replacedToNewMap.size() == 0) return; // nothing to weld


		std::vector<Vert> newVerts;
		std::map<int, int> oldToNewMapping;
		for (int vertIndex = 0; vertIndex < verts.size(); vertIndex++)
		{
			if (replacedToNewMap.find(vertIndex) != replacedToNewMap.end())
			{
				continue;
			}
			newVerts.push_back(verts[vertIndex]);
			oldToNewMapping[vertIndex] = newVerts.size() - 1;
		}

		// patch the elements
		for (auto&& each : elements)
		{
			for (int i = 0; i < 3; i++)
			{
				int needTrans = each.vertindex[i];
				if (replacedToNewMap.find(needTrans) != replacedToNewMap.end())
				{
					needTrans = replacedToNewMap[needTrans]; // translate to the replacemnt
				}

				int newTrans = oldToNewMapping[needTrans];
				each.vertindex[i] = newTrans;
			}
		}

		verts = newVerts; // force new verts;


	}

	float TriangleSimulation::CheckClosestSurface(float2* vertPoints, int* surfIndex)
	{
		auto oneOf = [&](int testIndex) {
			return surfIndex[0] == testIndex || surfIndex[1] == testIndex || surfIndex[2] == testIndex;
		};

		float closestSurface = 1000000.0f;
		for (auto&& each : conOpen)
		{
			Line2d triDir[3]={  { vertPoints[0], vertPoints[2] },
			{ vertPoints[2], vertPoints[1] }, // direction
			{ vertPoints[1], vertPoints[0] }}; // REVERSED TO THE ORIGINAL


			Line2d openLine = { verts[each.conVert[0]].r[0].v ,   verts[each.conVert[1]].r[0].v };
			float2 posOut;
			float2 tOut;
			bool lineCutAllowed = each.conVert[0] != surfIndex[2] && each.conVert[1] != surfIndex[2]; // neither can be endpoint
			bool aAllowed = each.conVert[0] != surfIndex[0] && each.conVert[1] != surfIndex[0]; // neither can be start
			bool bAllowed = each.conVert[0] != surfIndex[1] && each.conVert[1] != surfIndex[1]; // neither can be start

			if (lineCutAllowed && aAllowed &&  FiniteLineCutLine(openLine, triDir[0], posOut, tOut) )
			{
				closestSurface = 0.0f;
				break;
			}

			if (lineCutAllowed && bAllowed && FiniteLineCutLine(openLine, triDir[1], posOut, tOut) )
			{
				closestSurface = 0.0f;
				break;
			}


			for (int pi = 0; pi < 2; pi++)
			{
				if (!oneOf(each.conVert[pi]))
				{
					float maxDist = 0.0f;
					for (int k = 0; k < 3; k++)
					{
						maxDist = std::max(triDir[k].distToLine(verts[each.conVert[pi]].r[0].v), maxDist);
					}
					closestSurface = std::min(maxDist, closestSurface);
				}
			}	
		}

		return closestSurface;
	}

	bool TriangleSimulation::RemoveOpen(int index0, int index1)
	{
		for (auto iter = conOpen.begin(); iter != conOpen.end(); ++iter)
		{
			if (iter->conVert[0] == index0 && iter->conVert[1] == index1)
			{
				conOpen.erase(iter);
				return true;
			}
		}

		return false;
	}

	bool TriangleSimulation::joinSurface(int openIndex, int* vertIndex, float2& pTest, int forLayer)
	{
		elements.push_back({ vertIndex[1], vertIndex[0], vertIndex[2],-1, forLayer });
		int newElementIndex = elements.size() - 1;
		elements[newElementIndex].elementType = Element::CONTACT_ELEMENT;

		if (!RemoveOpen(vertIndex[2], vertIndex[0]))
		{
			conOpen.push_back({vertIndex[0], vertIndex[2], newElementIndex});
		}

		if (!RemoveOpen(vertIndex[1], vertIndex[2]))
		{
			conOpen.push_back({vertIndex[2], vertIndex[1], newElementIndex});
		}
		// remove the openIndex because we have closed it with the new element

		RemoveOpen(vertIndex[0], vertIndex[1]);
		return false;
	}


	bool TriangleSimulation::connectOneSurface(int forLayer)
	{
		float bestGood = 0.0f;
		float2 bestFarPoint;
		int  bestOpenIndex= 0;
		int bestVertIndex[3];

		for (int i = 0; i < conOpen.size(); i++)
		{
			auto& openi = conOpen[i];
			int surfIndex[3] = { openi.conVert[0] , openi.conVert[1] , -1 };
			float2 surfVert[3] = { verts[surfIndex[0]].r[0].v , verts[surfIndex[1]].r[0].v, float2() };
			Line2d openLine = { surfVert[0] , surfVert[1] };
			// search for a close acceptable  apex vertex

			if (openi.bestMatchVal >= 0.0f)
			{
				if (openi.bestMatchVal < bestGood)
					continue;
			}


			openi.bestMatchVal = 0.0f;
			for (int j = 0; j < conOpen.size(); j++)
			{
				if (i == j) continue;// dont connect to self

				auto& openj = conOpen[j];

				for (int k = 0; k < 2; k++)
				{
					surfIndex[2] = openj.conVert[k];
					surfVert[2] = verts[surfIndex[2]].r[0].v;
					// we need to make sure this new triangle doesnt pen another surface or point

					if (surfIndex[0] == surfIndex[1] || surfIndex[1] == surfIndex[2] || surfIndex[2] == surfIndex[0])
					{
						continue; // would be zero sized
					}


					if (openLine.distToLine(surfVert[2]) <= 0.0f)// should be outside orig triangle surface
					{
						continue;
					}

					Triangle2d tri = { surfVert[0], surfVert[1], surfVert[2] };
					float robustness = tri.robustness();

					float howGood = CheckClosestSurface(surfVert, surfIndex)* robustness;
					openi.bestMatchVal = std::max(openi.bestMatchVal, howGood);
					if (howGood > 0.0f && howGood > bestGood)
					{
						bestGood = howGood;
						bestFarPoint = surfVert[2];
						bestOpenIndex = i;
						for (int q = 0; q < 3; q++)
						{
							bestVertIndex[q] = surfIndex[q];
						}
					}
				}
			}
		}

		if (bestGood > 0.0f)
		{
			joinSurface(bestOpenIndex, bestVertIndex, bestFarPoint, forLayer);
			return true;
		}

		return false;
	}

	void TriangleSimulation::testIntegration()
	{
		return;



		return;


		{
			float energyAtContact = ContactEnergy(minContactDistance, 1.0f, 1.0f);
			float subContactDist = 0.0001f;
			float energyAtSunContact = ContactEnergy(minContactDistance - subContactDist, 1.0f, 1.0f);
			float forceAtContact = (energyAtSunContact - energyAtContact) / subContactDist;
			//assert(forceAtContact == 0);
			//assert(energyAtContact == 0);
			//assert(energyAtContact < energyAtSunContact);

			for (float x = minContactDistance; x > 0.01f; x -= 0.001f)
			{
				float forceForEnergy = (ContactEnergy(x, 1.0f, 3.3333f) - ContactEnergy(x + 0.0001f, 1.0f, 3.3333f)) / 0.0001f;
				float estForce = estimateForceFromContactEnergy(ContactEnergy(x, 1.0f, 3.3333f), 1.0f, 3.3333f);
				printf("temp");
			}

		}

		// line point
		{
			// line is expected to be +Y
			Line2d testLine(float2(-.3, minContactDistance * 0.7f), float2(1.0f, minContactDistance * 1.5005f));

			float avgY = (testLine.pi.y + testLine.pf.y) * 0.5f;
			auto toDistanceLine = [&](float x)
			{

				float a = testLine.getCompA();
				float b = testLine.getCompB();
	
				// t
				return   b / (cos(x) + a * sin(x));
			};

			// 0 to pi starting from x+ coord
			float2 start = testLine.pi;
			float2 stop = testLine.pf;
			testLine.getStartStopForDistance(minContactDistance, start, stop);

			float fromX = start.GetAngle() - (bdts_pi/2);
			float toX = stop.GetAngle() - (bdts_pi / 2);

			auto integrateDist = [&](float fromX, float toX)
			{
				double energyOut = 0.0f;
				int numSamples = 2000;
				float currX = fromX;
				float da = (toX - fromX) / (numSamples-1);
				for (int i = 0; i < numSamples; i++)
				{
					float dist = toDistanceLine(currX);
					energyOut += ContactEnergy(dist, 1.0f, 1.0f);
					currX += da;
				}

				return -(toX - fromX) * (energyOut / numSamples) ; // area * norm sample / pi
			};

			float computeEnergy = -integrateDist(fromX, toX);
			float intEnergy = (pureIntegratePointEnergyFunction(toX, testLine) - pureIntegratePointEnergyFunction(fromX, testLine)) ;

			//assert(computeEnergy == intEnergy);

		}


		// line line
		{
			float di = minContactDistance * 0.5f;
			float df = minContactDistance * 0.505f;
			Line2d testLine(float2(-1.f, di), float2(1.0f, df));
			float xi = 0.0f;
			float xf = testLine.length();

			auto integrateDist = [&]()
			{
				double sampledEnergy1 = 0;
				{
					double sum = 0.0f;
					int sampleCount = 1000;
					float runDist = di;
					float runInc = (df - di) / ((float)(sampleCount - 1));
					for (int i = 0; i < sampleCount; i++)
					{
						sum += ContactEnergy(runDist, 1.0f, 1.0f);
						runDist += runInc;
					}

					sampledEnergy1 = sum * ((xf - xi) / ((double)sampleCount));
				}
				return sampledEnergy1;
			};

			float computeEnergy = integrateDist();
			float intEnergy = (pureIntegrateLineEnergyFunction(xf, testLine) - pureIntegrateLineEnergyFunction(xi, testLine)) ;
			//assert(computeEnergy == intEnergy);
		}
	}

	void TriangleSimulation::debugRenderOut(std::vector<float4>& posOut, std::vector<float4>& colOut)
	{
#if 0
		for (auto&& vert : verts)
		{
			vert.vonCount = 0;
			vert.vonTotal = 0;
		}

		for (auto&& each : elements)
		{

			//each.tempVisEnergy = sqrt(each.compressionEnergyPerArea);
			each.tempVisEnergy = sqrtf(each.devEnergyPerArea);

			for (int i = 0; i < 3; i++)
			{
				verts[each.vertindex[i]].vonCount += 1.0f;
				verts[each.vertindex[i]].vonTotal += each.tempVisEnergy;
			}

		}

		int countElement = 0;
		for (auto&& each : elements)
		{
			//if (sCounter <= countElement)
			//	break;

			for (int i = 0; i < 3; i++)
			{
				posOut.push_back(float4(verts[each.vertindex[i]].r[0].x,
					verts[each.vertindex[i]].r[0].y
					, 1.0, 1.0));
			}

			for (int i = 0; i < 3; i++)
			{
				float smoothVon = verts[each.vertindex[i]].vonTotal / verts[each.vertindex[i]].vonCount;
				float4 singleColor = Utility::toFalseColor(smoothVon  *2.2);
				colOut.push_back({ 1,0,0,0.3f });// singleColor);// {1.f, 0.2f, 0.2f, 0.5f});
			}
			countElement++;
		}

#else
		SimObj o = getSimObj();
		ClearForce(&o, 0);
		ClearForce(&o, 1);
		for (int elementIndex = 0; elementIndex < elements.size(); elementIndex++)
		{
			ForceForElement(&o, elementIndex, 0, 0, 0);
		}

		for (auto&& each : elements)
		{

			for (int i = 0; i < 3; i++)
				posOut.push_back(float4(verts[each.vertindex[i]].r[0].v.x,
					verts[each.vertindex[i]].r[0].v.y
					, 1.0, 1.0));


			//each.tempVisEnergy = each.areaRatio;//  sqrt(each.hydroEnergyPerArea);//
			//each.tempVisEnergy = sqrt(Dot( each.plasticStrainApplied, each.plasticStrainApplied));
			each.tempVisEnergy = sqrtf(std::max(0.0f, each.devEnergyPerArea));
			//each.tempVisEnergy = sqrt(each.devEnergyPerArea);
			float4 singleColor = toFalseColor((each.tempVisEnergy  *0.5f));
			singleColor.w = 0.8f;
			colOut.push_back(singleColor);// {1.f, 0.2f, 0.2f, 0.5f});
			colOut.push_back(singleColor);// { 0.2f,1.f,0.2f,0.5f });
			colOut.push_back(singleColor); //{ 0.2f,0.2f,1.f,0.5f });
		}
#endif
	}

	bool TriangleSimulation::assignVertexSurfaces()
	{
		bool bWasChange = false;
		for (auto&& each : elements)
		{
			int elementLayer = each.layerIndex;
			if (each.elementType == Element::CONTACT_ELEMENT)
			{
				for (int i = 0 ; i <3; i++ )
				{
					int neighIndex = each.neighElem[i];
					if (neighIndex != -1 && elements[neighIndex].elementType == Element::SOLID_ELEMENT)
					{
						// my contact element has a solid surface running from neigh verts i to i + 1
						int vertIndexOut = each.vertindex[i];
						int vertIndexIn = each.vertindex[(i + 1) % 3];

						if (verts[vertIndexIn].neighSurfLayer[elementLayer].neighSurfaceIn != vertIndexOut
							|| verts[vertIndexOut].neighSurfLayer[elementLayer].neighSurfaceOut != vertIndexIn)
						{
							verts[vertIndexIn].neighSurfLayer[elementLayer].neighSurfaceIn = vertIndexOut;
							verts[vertIndexOut].neighSurfLayer[elementLayer].neighSurfaceOut = vertIndexIn;
							bWasChange = true;
						}
					}
				}
			}
		}
		return bWasChange;
	}


	
	double TriangleSimulation::energyTotal(int ri, int vi, int ai, bool includeField, bool includeKE, bool includePE)
	{
		SimObj o = getSimObj();
		ReplacePos rp;
		rp.vertIndex = -1;
		double energy = 0.0f;

		for (auto&& each : verts)
		{
			if (includeField)
			{
				energy += -Dot(each.r[ri].v , gravity) * each.mMassCount;
			}

			if (includeKE)
			{
				energy += Dot(each.v[vi], each.v[vi])  * each.mMassCount *0.5f;
			}
		}

		if(includePE)
		for (int i= 0 ; i < elements.size(); i++)
		{
			energy += EnergyForElement(&o, &rp, i, ri, ai,  eEnergyComputePhase::VARIATION_COMPUTE_PHASE);
		}

		return energy;
	}

	void TriangleSimulation::simAdvance(double timeForward)
	{
		static volatile bool bDoAvec = true;
		static bool bUseDampeningFriction = false;
		// called 
		SimObj o = getSimObj();
		ReplacePos rp;
		rp.vertIndex = -1;
		for (int microUpdate = 0; microUpdate < 1; microUpdate++)
		{
			simTime += deltaT * 1.0;

			ClearForce(&o, 0);
			ClearForce(&o, 1);

			for (int i = 0; i < elements.size(); i++)
			{
				EnergyForElement(&o, &rp, i, 0, 0, eEnergyComputePhase::PREPERATION_COMPUTE_PHASE); // apply phase
			}

			for (int i = 0; i < elements.size(); i++)
			{
				ForceForElement(&o,i, 0, 0, 0);
			}

			if (bDoAvec)
			{
				for (auto&& vert : verts)
				{
					if (!vert.doNotAvect)
					{
						vert.a[0] = vert.a[0] + gravity * vert.mMassCount;
						vert.r[1].v = vert.r[0].v + vert.v[0] * deltaT + vert.a[0] * deltaT * deltaT * 0.5f / vert.mMassCount;
					}
				}

				for (int i = 0; i < elements.size(); i++)
				{
					ForceForElement(&o, i, 1, 1, 1);
				}

				for (auto&& vert : verts)
				{
					if (!vert.doNotAvect)
					{
						vert.a[1] = vert.a[1] +  gravity * vert.mMassCount;
						float2 aDir = vert.a[1] + vert.a[0];
						aDir = aDir * 0.5f;
						vert.v[1] = vert.v[0] + aDir * deltaT  / vert.mMassCount;;
					}
				}


				float beforeE1 = energyTotal(1, 1, 1);
				for (int i = 0; i < elements.size(); i++)
				{
					EnergyForElement(&o,&rp, i, 1, 1, eEnergyComputePhase::LOSS_COMPUTE_PHASE); // apply phase
				}

				for (int i = 0; i < elements.size(); i++)
				{
					EnergyForElement(&o,&rp, i, 1, 1, eEnergyComputePhase::PREPERATION_COMPUTE_PHASE); // apply phase
				}

				float afterE1 = energyTotal(1, 1, 1);
				//energyFriction -= beforeE1 - afterE1;

				for (int i = 0; i < elements.size(); i++)
				{
					DiffusionForElement(&o, i, 1, 1, 1, deltaT); // viscosity
				}

				for (auto&& vert : verts)
				{
					if (!vert.doNotAvect)
					{
						vert.r[0] = vert.r[1];
						vert.v[0] = vert.v[1];
					}
				}

				if (bUseDampeningFriction)
				{
					for (auto&& element : elements)
					{
						if (element.elementType == Element::CONTACT_ELEMENT) continue;

						auto& vert0 = verts[element.vertindex[0]].v[0];
						auto& vert1 = verts[element.vertindex[1]].v[0];
						auto& vert2 = verts[element.vertindex[2]].v[0];

						float2 avg = (vert0 + vert1 + vert2) / 3.0f;

						vert0 = vert0 * 0.999f + avg * 0.001f;
						vert1 = vert1 * 0.999f + avg * 0.001f;
						vert2 = vert2 * 0.999f + avg * 0.001f;
					}
				}
			}
			adaptiveAll();
			BDTS_ASSERT(!assignVertexSurfaces(), );// right now we dont change surface topology
		}

	}

}