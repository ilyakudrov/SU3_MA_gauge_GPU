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
 * Stochastic Relaxation
 *
 * Flops: 32
 */

#ifndef SRUPDATE_HXX_
#define SRUPDATE_HXX_

#include "../../lattice/datatype/datatypes.h"
#include "../../lattice/rng/PhiloxWrapper.hxx"

class SrUpdate
{
public:
	__device__ inline SrUpdate();
	__device__ inline SrUpdate( float param, PhiloxWrapper *rng );
	__device__ inline void calculateUpdate( volatile Real (&shA)[4*NSB], short id );
	__device__ inline void setParameter( float param );
	__device__ inline float getParameter();
private:
	float srParameter;
	PhiloxWrapper *rng;
};

__device__ SrUpdate::SrUpdate()
{
}

__device__ SrUpdate::SrUpdate( float param, PhiloxWrapper *rng  ) : srParameter(param), rng(rng)
{
}

__device__ void SrUpdate::calculateUpdate( volatile Real (&shA)[4*NSB], short id )
{
#ifdef USE_DP_SRUPDATE
	double rand = rng->rand();
	double a0,a1,a2,a3,c;
	a0 = shA[id];
	a1 = shA[id+NSB];
	a2 = shA[id+2*NSB];
	a3 = shA[id+3*NSB];

	shA[id]    = (rand>=srParameter)*a0 + (rand<srParameter)*(a0*a0-a1*a1-a2*a2-a3*a3); // 12 flop
	shA[id+NSB] = (rand>=srParameter)*a1 + (rand<srParameter)*(2.0*a0*a1); // 7 flop
	shA[id+2*NSB] = (rand>=srParameter)*a2 + (rand<srParameter)*(2.0*a0*a2);
	shA[id+3*NSB] = (rand>=srParameter)*a3 + (rand<srParameter)*(2.0*a0*a3);

	c=rsqrt(shA[id]*shA[id]+shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]+shA[id+3*NSB]*shA[id+3*NSB]); // 9 flop
	shA[id]    *= c;
	shA[id+NSB] *= c;
	shA[id+2*NSB] *= c;
	shA[id+3*NSB] *= c;
	// 4 flop

	// sum: 32 flop
#else
	Real rand = rng->rand();
	Real a0,a1,a2,a3,c;
	a0 = shA[id];
	a1 = shA[id+NSB];
	a2 = shA[id+2*NSB];
	a3 = shA[id+3*NSB];
	
	shA[id]    = (rand>=srParameter)*a0 + (rand<srParameter)*(a0*a0-a1*a1-a2*a2-a3*a3);
	shA[id+NSB] = (rand>=srParameter)*a1 + (rand<srParameter)*(2.0*a0*a1);
	shA[id+2*NSB] = (rand>=srParameter)*a2 + (rand<srParameter)*(2.0*a0*a2);
	shA[id+3*NSB] = (rand>=srParameter)*a3 + (rand<srParameter)*(2.0*a0*a3);

	c=rsqrt(shA[id]*shA[id]+shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]+shA[id+3*NSB]*shA[id+3*NSB]);
	shA[id]    *= c;
	shA[id+NSB] *= c;
	shA[id+2*NSB] *= c;
	shA[id+3*NSB] *= c;
#endif
}

__device__ void SrUpdate::setParameter( float param )
{
	srParameter = param;
}

__device__ float SrUpdate::getParameter()
{
	return srParameter;
}
#endif /* SRUPDATE_HXX_ */
