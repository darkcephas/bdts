//
// Copyright (c) Microsoft. All rights reserved.
// This code is licensed under the MIT License (MIT).
// THIS CODE IS PROVIDED *AS IS* WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING ANY
// IMPLIED WARRANTIES OF FITNESS FOR A PARTICULAR
// PURPOSE, MERCHANTABILITY, OR NON-INFRINGEMENT.
//
// Developed by Minigraph
//
// Author:  James Stanard 
//

#pragma once

#include "bdts_base.h"
#include <Engine/SystemsCore.h>
#include <vector>
#include <list>
#include <unordered_map>
#include <string>
#include <btds_math.h>
#include <btds_func.h>
#include <btds_sim_obj.h>
namespace bdts
{



	struct TraceArrow
	{
		float2 pi;
		float2 dir;
		float4 color;
	};

	struct TracePrint
	{
		float2 pos;
		std::string disp;
		float4 color;
	};
	

	class TriangleSimulation
	{
	public:



		// assymetric force integration issue hack
		// basically we need to lag energy in stored dissipators so that our force is non zero for friction 

		struct tElementUnitAreaProp
		{
			Matrix22 worldDefMat;// deformation in worldspace
			Matrix22 invRot;
			Matrix22 rot;
			Matrix22 invDef;
			float4 matSkew;
			float2 expandedVert[3];
			float2 simVert[3];
		};
		





		struct OpenSurface
		{
			OpenSurface(int verta, int vertb, int cone)
			{
				conVert[0] = verta;
				conVert[1] = vertb;
				conElement = cone;
			}

			bool hasVert(int index)
			{
				return conVert[0] == index || conVert[1] == index;
			}

			int conVert[2];
			int conElement;
			float bestMatchVal = -1.0f;
		};


		TriangleSimulation();


		// simulation
		SimObj getSimObj();
		void InitElementsAsCurrent(int elementIndex, TriangleSimulation::tElementUnitAreaProp* props = nullptr);
		void VertMassChangeForElement(Element & element, float density, float currArea);

		// timestepping simulation
		
		TriangleSimulation::tElementUnitAreaProp GetElmentUnitAreaProp(int elementIndex);
		void simAdvance(double timeForward);
		Matrix22 strainToDefMatrix(float4 strain);



		// sim mesh modification
		void InitElements();
		void adaptiveOne(Element & selElem, int elemIndex);
		void adaptiveAll();
		bool assignVertexSurfaces();
		bool refineBisplit(int elemIndexi, int elemIndexj, int neighVertIndexi, int neighVertIndexj);


		// mesh edit
		void GenerateBoxPoint(int matIndex, int countX, int countY, float2 low, float2 high);
		void addBorder(int matIndex, int dimX, int dimY, float thickness, int layer);
		void fillSimulation();
		void deleteElement(int elementIndex);


		// connection setup
		void findSurface(int forLayer);
		bool connectOneSurface(int forLayer);
		void assignNeighElements();
		bool joinSurface(int openIndex, int* vertIndex, float2& pTest, int forLayer);
		float CheckClosestSurface(float2* vertPoints, int* surfIndex);
		void weldClose(float maxDist);
		bool RemoveOpen(int index0, int index1);
		bool refineCentered(int elemIndex);

		// testing
		void testInitElement(int elemIndex);
		void testIntegration();

		// editor
		void createMeshGrid(int matIndex, float2 bottom, float2 top, int dimx, int dimy);
		void setLimits(float2 bottomv, float2 topv);
		void createMeshTube(bool bCentered, int matIndex,   float2 center, float radiusIn, float radiusOut, int segRad, int segScale,  bool startNoAvec = false, bool geared = false, float gearRad = 0.0f, float rotOffset = 0.0f);
		void setSimulationBounds(float2 minCorner, float2 maxCorner);

		// debug
		double energyTotal(int ri, int vi, int ai, bool includeField = true , bool includePE = true, bool includeKE = true);
		void DrawingOutput(int viewLayer, std::vector<float4>& posOut, std::vector<float4>& colOut, bool useUnique = true, bool useMaterial = false, float4 solidCol = float4(1.0f, 0.0f, 0.0f, 1.0f), float4 contactCol = float4(0.0f, 1.0f, 0.0f, 1.0f));
		void debugRenderOut(std::vector<float4>& posOut, std::vector<float4>& colOut);



		float2 bottom;
		float2 top;
		std::vector<OpenSurface> conOpen;
		std::vector<Element> elements;
		std::vector<Vert> verts;

		std::vector<EMat> eMaterials;

		// debugging
		double energyOffset = 0.0f;
		double energyFriction = 0.0f;
		double energyPlastic = 0.0f;
		double simTime = 0.0f;
		std::unordered_map<std::string, bool> assertEnabled;
		std::vector<TraceArrow> tracingArrows;
		std::vector<TracePrint> tracingPrints;
	};

	class FrameRecorder
	{
	public:

		void saveNewFrame(TriangleSimulation& inSim)
		{
			currFrame++;
			frames.resize(currFrame + 1);
			frames[currFrame] = inSim;
		}

		void loadFrame(TriangleSimulation& outSim)
		{
			outSim = frames[currFrame];
		}

		void advanceFrame() { setCurrFrame(currFrame + 1); }
		int getCurrFrame() { return currFrame; }
		int getLastFramint() { return frames.size()-1; }
		void setCurrFrame(int newFramint) {
			 currFrame = std::max(std::min(newFramint, getLastFramint()), 0);
		}
	private:
		int currFrame = -1;
		std::vector<TriangleSimulation> frames;
	};


}