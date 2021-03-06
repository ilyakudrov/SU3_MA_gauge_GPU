/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
 *
 * Wraps the Philox4x32-10 counter based RNG (CBRNG).
 * The benefit of CBRNGs is the absence of large states which are costly to read and write in CUDA code.
 * For a detailed discussion see "Parallel random numbers: as easy as 1, 2, 3" (http://dl.acm.org/citation.cfm?doid=2063405)
 *
 * Any combination of the 4 * 32-bit counter and 2 * 32-bit key gives a uncorrelated random number.
 * As key we use the combination (thread id, seed), where seed can be chosen by the user.
 * The counter is (kernelCounter, globalCounter, *, * ) where
 *  - kernelCounter is a local variable of the kernel and is incremented by each call to generate new random numbers
 *  - globalCounter is a kernel parameter, incrementing has to be done on host side. (TAKE CARE TO DO IT!)
 *    This means that each kernel that uses random numbers has to have an counter parameter which has to be incremented on each kernel call.
 *  - the other two counters are set arbitrarily.
 *
 * Each call to philox4x32 calculates 4 32 bit values.
 *
 * We offer a rand() function that
 *  - returns a Real (float or double)
 *  - takes care to use already calculated numbers.
 * 	- and increments the kernelCounter.
 *
 * The static getNextCounter() function gives a global (runtime-wide) counter which can be given to the constructor.
 * I don't want to do this implicitly to avoid mixing of host and device variables.
 *
 *
 * TODO: get rid of the preprocessor statements: float/double has to be template argument (in all classes).
 * 		Because now there is no code possible that uses both single and double precision random numbers.
 * TODO: In CUDA5.0 a lot of stack-frame is used when Philox is involved (~ 200 bytes for the single precision SA kernel)
 */

#ifndef PHILOXWRAPPER_HXX_
#define PHILOXWRAPPER_HXX_

#include "../datatype/datatypes.h"
#include "../cuda/cuda_host_device.h"
#include "../../external/Random123/philox.h"
#include "../../external/Random123/u01.h"

#include <cuda.h>

#include <stdio.h>

class PhiloxWrapper
{
public:
	__device__ inline PhiloxWrapper( int tid, int seed, unsigned int globalCounter );
	__device__ inline virtual ~PhiloxWrapper();
	__device__ inline Real rand();

	static __host__ unsigned int getNextCounter();
	static __host__ unsigned int getCurrentCounter();

private:
	philox4x32_key_t k;
	philox4x32_ctr_t c;
	union
	{
		philox4x32_ctr_t res;
#ifdef DOUBLEPRECISION
		uint64_t i[2];
#else
		uint32_t i[4];
#endif
	} u;
	short localCounter;
//	int localCounter;
	static unsigned int globalCounter;
};


unsigned int PhiloxWrapper::globalCounter = 0;

__device__ PhiloxWrapper::PhiloxWrapper( int tid, int seed, unsigned int globalCounter )
{
	k[0] = tid;
	k[1] = seed;
	c[0] = 0; // kernelCounter
	c[1] = globalCounter;
	c[2] = 0x12345678;
	c[3] = 0xabcdef09;
	u.res = philox4x32(c,k);
	localCounter = 0;
}

__device__ PhiloxWrapper::~PhiloxWrapper()
{
}

__device__ Real PhiloxWrapper::rand()
{
#ifdef DOUBLEPRECISION
	if( localCounter == 0 )// we have to calculate two new doubles
	{
		localCounter = 1;
		c[0]++; // inkrement the kernel counter
		u.res = philox4x32(c,k);
		return u01_open_open_64_53( u.i[localCounter] );
	}
	else  // we have another double available.
	{
		localCounter--;
		return u01_open_open_64_53( u.i[localCounter] );
	}
#else
	if( localCounter == 0 ) // we have to calculate 4 new floats
	{
		localCounter = 3;
		c[0]++; // inkrement the kernel counter
		u.res = philox4x32(c,k);
		return u01_open_open_32_24( u.i[localCounter] );
	}
	else // we have another float available.
	{
		localCounter--;
		return u01_open_open_32_24( u.i[localCounter] );
	}
#endif
}

/**
 * Use this to get a global (static) counter to initialize the kernel-specific PhiloxWrapper.
 */
__host__ unsigned int PhiloxWrapper::getNextCounter()
{
	return globalCounter++;
}

__host__ unsigned int PhiloxWrapper::getCurrentCounter()
{
	return globalCounter;
}



#endif /* PHILOXWRAPPER_HXX_ */
