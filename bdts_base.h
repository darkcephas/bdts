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

#ifndef BDTS_BASE__INCLUDED
#define BDTS_BASE__INCLUDED


#ifndef BDTS_IS_CPP_COMPILE 
#ifdef _WINDOWS
#define  BDTS_IS_CPP_COMPILE 1
#else
#define  BDTS_IS_CPP_COMPILE 0
#endif

#endif

#define BDTS_ENABLE_GLOBAL_ENERGY_TRACKING 1

#if BDTS_IS_CPP_COMPILE 
namespace bdts
{
#endif
	static const int sNumLayer = 2;
	static const float deltaT = 0.0001f *1.0f;
	static const float minContactDistance = 0.05f;
	static const float contactEnergyMult = 1.0f;
	static const float sCompressiveMult =  1000.0f;
#if BDTS_IS_CPP_COMPILE 
}
#endif

#endif //BDTS_BASE__INCLUDED