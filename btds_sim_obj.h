#ifndef BDTS_SIM_OBJ__INCLUDED
#define BDTS_SIM_OBJ__INCLUDED


#include "bdts_base.h"
#include "btds_math.h"
#include "btds_func.h"

#if BDTS_IS_CPP_COMPILE 
namespace bdts
{
#endif

	struct tSetObj
	{
		int ids[12];
		int count;
	};

	bool microset_contains(tSetObj* obj, int checkid);
	void microset_init(tSetObj* obj);
	bool microset_insert(tSetObj* obj, int checkid);

	enum class eEnergyComputePhase
	{
		VARIATION_COMPUTE_PHASE,
		PREPERATION_COMPUTE_PHASE,
		LOSS_COMPUTE_PHASE
	};

	struct NonConEnergyDissipator
	{
		float nominalEnergy = 0.0f;
		float prevEnergy = 0.0f;
		float frameEnergyVariation = 0.0f;
	};

	struct Element
	{
		enum ElementType
		{
			SOLID_ELEMENT,
			CONTACT_ELEMENT,
			ELEMENT_MAX,
		};

		Element(int a, int b, int c, int matIndex = -1, int layer = 0)
		{
			vertindex[0] = a;
			vertindex[1] = b;
			vertindex[2] = c;
			neighElem[0] = -1;
			neighElem[1] = -1;
			neighElem[2] = -1;
			stressX = 0;
			stressY = 0;
			stressShear = 0;
			devEnergyPerArea = 0.0f;
			plasticStrainApplied = float4();
			elementType = SOLID_ELEMENT;
			multEnergy = 1.0f;
			nominalDevEnergy = 0.0f;
			materialId = matIndex;
			layerIndex = layer;
		}
		int vertindex[3];
		int neighElem[3];// on surfaces
		float sideDistance[3];// for this vert where is it on the other verts (dynamic friction computations)
		float forceEstForFriction; // is it required that this force be stable ?? is this correct

		int layerIndex = 0;

		float baseLengthOrig; //  used to computer Exx
		float barCrossLengthOrig;// used to compute Eyy
		float barBaseLengthOrig; // used to compute Exy , Eyx
		float origArea;
		float stressX;
		float stressY;
		float stressShear;
		float multEnergy;
		ElementType elementType;
		float nominalDevEnergy; // we need to be able to tell how much is normal energy for restitutive
		NonConEnergyDissipator solidFrictionDissipator;
		int materialId;
		float4 plasticStrainApplied;

		// visualization stuff for now debugging
		float tempVisEnergy;
		float tempVisTestEnergy;
		float devEnergy;
		float devEnergyPerArea;
		float plasticEnergyPerArea;
		float plasticEnergy;
		float hydroEnergyPerArea;
		float hydroEnergy;
		float fullEnergy;
		float compressionEnergyPerArea;
		float areaRatio;
		float4 strainVec;
		float4 strainDev;
	};

	struct Pos
	{
		float2 v;
	};

	struct ReplacePos
	{
		int vertIndex;
		float2 pos;
	};



	float2 rcheck_replacePos(ReplacePos* rp, int vertIndex, Pos p);

	struct Vert
	{
		Vert(float2 initPos)
		{
			r[0].v = initPos;
			r[1].v = initPos;

			v[0] = float2();
			v[1] = float2();
			a[0] = float2();
			a[1] = float2();

			mMassCount = 0.0f;
			doNotAvect = false;
			for (auto&& each : neighSurfLayer)
			{
				each.neighSurfaceIn = -1;
				each.neighSurfaceOut = -1;
			}

		}

		Pos r[2];
		float2 v[2];
		float2 a[2];



		struct tNeighSurface
		{
			int neighSurfaceIn;
			int neighSurfaceOut;
		};

		tNeighSurface   neighSurfLayer[sNumLayer];

		float mMassCount;
		bool doNotAvect;

		// visual smoothing
		float vonTotal;
		float vonCount;
	};


	// material for vert (if required)
	struct VMat{};

	// material for an element (refered by index)
	struct EMat
	{
		EMat(float youngsMod, float maxElasticPerArea, float frictionCof, float frictionSurface, int colIndex);
		float em_youngsMod; // strength of spring
		float em_maxElasticPerArea; // when in energy does plastic deformation occure ( F * d = energy) so strength
		float em_frictionCof = 0.0f;
		float em_frictionSurface = 0.0f;// standard newton friction coefficent 
		float em_frictionResitRatio;
		float em_poissonRatio = 0.3333333333f;  // volume preservation through elasticity (youngs mod)
		float em_viscosity = 0.9f;
		// debug/visualize
		float4 debugColor;
		float em_maxPlasticPerArea;// NOT IMPL YET
	};

	struct SimObj
	{
		Vert * verts; 
		int vertsCount;
		EMat* materials;
		Element* elements;


		// global tracking (mostly for debugging),  BDTS_ENABLE_GLOBAL_ENERGY_TRACKING
		double* energyOffset;
		double* energyFriction;
		double* energyPlastic;
		double* simTime;

	};

	float nced_evalEnergyState(NonConEnergyDissipator * obj, eEnergyComputePhase ePhase, float resitutitveMult, float energyIn, float& lossOut);

	float nced_getCurrEnergy(NonConEnergyDissipator * obj, float resitutitveMult, float energyIn);

	void nced_startSample(NonConEnergyDissipator * obj, float resitutitveMult, float energyIn);

	float nced_endLoss(NonConEnergyDissipator * obj, float resitutitveMult, float energyIn);

#if BDTS_IS_CPP_COMPILE 
}
#endif

#endif