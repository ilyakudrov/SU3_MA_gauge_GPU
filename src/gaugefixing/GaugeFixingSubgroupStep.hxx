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
 * Flops:
 * - summation of local update: 4*2*8 = 64 (Landau), 3*2*8 = 48 (Coulomb)
 * - calculating update (see update classes)
 * - apply transformation: 4*2*3*26 (mu*(up+down)*Nc*...) = 504
 */

#ifndef GAUGEFIXINGSUBGROUPSTEP_HXX_
#define GAUGEFIXINGSUBGROUPSTEP_HXX_

#include "../lattice/datatype/datatypes.h"
#include "../lattice/datatype/lattice_typedefs.h"
#include "../lattice/Quaternion.hxx"
#include "./GlobalConstants.h"
#include "../lattice/Matrix.hxx"

#define USEATOMIC
#ifdef DOUBLEPRECISION
// // this implementation of DP atomic add is slower than the explicit technique with syncthreads().
//__device__ double atomicAdd(double* address, double val)
//{
//    unsigned long long int* address_as_ull =
//                              (unsigned long long int*)address;
//    unsigned long long int old = *address_as_ull, assumed;
//    do {
//        assumed = old;
//        old = atomicCAS(address_as_ull, assumed,
//                        __double_as_longlong(val +
//                               __longlong_as_double(assumed)));
//    } while (assumed != old);
//    return __longlong_as_double(old);
//}

#undef USEATOMIC
#endif

template<class SUx, class Algorithm, GaugeType lc > class GaugeFixingSubgroupStep
{
public:
	__device__ inline GaugeFixingSubgroupStep( SUx* U, Algorithm algo, const short id, const short mu, const bool updown );
	__device__ inline void subgroup( const  int i, const  int j );
private:
	SUx* U;
	Algorithm algo;
	const short id;
	const short mu;
	const bool updown;
};

template<class SUx, class Algorithm, GaugeType lc > __device__ GaugeFixingSubgroupStep<SUx, Algorithm, lc>::GaugeFixingSubgroupStep( SUx* U, Algorithm algo, const short id, const short mu, const bool updown ): U(U), algo(algo), id(id), mu(mu), updown(updown)
{
}

template<class SUx, class Algorithm, GaugeType lc > __device__ void GaugeFixingSubgroupStep<SUx, Algorithm, lc>::subgroup( const int i, const int j )
{
	Quaternion<Real> q;
#ifdef USEATOMIC
	__shared__ Real shA[4*NSB];
#else
	__shared__ Real shA[4*NSB];
//	volatile __shared__ Real shA[4*NSB];
#endif
//	volatile __shared__ Real shA[4*NSB];

	__syncthreads();
	if( mu == 0 )
	{
		if( updown == 0 )
		{
			shA[id]	= 0;
			shA[id+NSB] = 0;
			shA[id+2*NSB] = 0;
			shA[id+3*NSB] = 0;
		}
	}
	__syncthreads();


	if( lc == LANDAU || lc == COULOMB )
	{

#ifdef USEATOMIC
		q = U->getSubgroupQuaternion( i, j );
		if( updown == 0 )
		{
			if( lc == LANDAU || mu != 0 )
			{
				atomicAdd( &shA[id], q[0] );
				atomicAdd( &shA[id+NSB], -q[1] );
				atomicAdd( &shA[id+2*NSB], -q[2] );
				atomicAdd( &shA[id+3*NSB], -q[3] );
			}
		}
		else
		{
			if( lc == LANDAU || mu != 0 )
			{
				atomicAdd( &shA[id], q[0] );
				atomicAdd( &shA[id+NSB], q[1] );
				atomicAdd( &shA[id+2*NSB], q[2] );
				atomicAdd( &shA[id+3*NSB], q[3] );
			}
		}


	#else
		if( lc == LANDAU )
		{
			if( updown == 0 && mu == 0 )
			{
			q = U->getSubgroupQuaternion( i, j );
				shA[id] += q[0];
				shA[id+NSB] -= q[1];
				shA[id+2*NSB] -= q[2];
				shA[id+3*NSB] -= q[3];
			}
			__syncthreads();
		}

		if( updown == 0 && mu == 1 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+NSB] -= q[1];
			shA[id+2*NSB] -= q[2];
			shA[id+3*NSB] -= q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 2 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+NSB] -= q[1];
			shA[id+2*NSB] -= q[2];
			shA[id+3*NSB] -= q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 3 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+NSB] -= q[1];
			shA[id+2*NSB] -= q[2];
			shA[id+3*NSB] -= q[3];
		}
		__syncthreads();


		if( lc == LANDAU )
		{
			if( updown == 1 && mu == 0 )
			{
			q = U->getSubgroupQuaternion( i, j );
				shA[id] += q[0];
				shA[id+NSB] += q[1];
				shA[id+2*NSB] += q[2];
				shA[id+3*NSB] += q[3];
			}
			__syncthreads();
		}

		if( updown == 1 && mu == 1 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+NSB] += q[1];
			shA[id+2*NSB] += q[2];
			shA[id+3*NSB] += q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 2 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+NSB] += q[1];
			shA[id+2*NSB] += q[2];
			shA[id+3*NSB] += q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 3 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+NSB] += q[1];
			shA[id+2*NSB] += q[2];
			shA[id+3*NSB] += q[3];
		}
#endif


		__syncthreads();
		if( mu == 0 )
		{
			if( updown == 0 )
			{
				algo.calculateUpdate( shA, id );
			}
		}
		__syncthreads();

		q[0] = shA[id];
		q[1] = shA[id+NSB];
		q[2] = shA[id+2*NSB];
		q[3] = shA[id+3*NSB];



		if( updown == 0 )
		{
			U->leftSubgroupMult( i, j, &q );
		}
		else
		{
			q.hermitian();
			U->rightSubgroupMult( i, j, &q );
		}
	}
	else if( lc == U1xU1 )
	{
		if( updown == 0 && mu == 0 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] -= q[1];
// 			shA[id+2*NSB] -= q[2];
			shA[id+3*NSB] -= q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 1 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] -= q[1];
// 			shA[id+2*NSB] -= q[2];
			shA[id+3*NSB] -= q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 2 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] -= q[1];
// 			shA[id+2*NSB] -= q[2];
			shA[id+3*NSB] -= q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 3 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] -= q[1];
// 			shA[id+2*NSB] -= q[2];
			shA[id+3*NSB] -= q[3];
		}
		__syncthreads();


		if( updown == 1 && mu == 0 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] += q[1];
// 			shA[id+2*NSB] += q[2];
			shA[id+3*NSB] += q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 1 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] += q[1];
// 			shA[id+2*NSB] += q[2];
			shA[id+3*NSB] += q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 2 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] += q[1];
// 			shA[id+2*NSB] += q[2];
			shA[id+3*NSB] += q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 3 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] += q[1];
// 			shA[id+2*NSB] += q[2];
			shA[id+3*NSB] += q[3];
		}
		__syncthreads();

				
		if( mu == 0 )
		{
			if( updown == 0 )
			{
				shA[id+NSB]=0.0;
				shA[id+2*NSB]=0.0;
				algo.calculateUpdate( shA, id );
			}
		}
		__syncthreads();

		q[0] = shA[id];
		q[1] = 0.0;//shA[id+NSB];
		q[2] = 0.0;//shA[id+2*NSB];
		q[3] = shA[id+3*NSB];


		Real norm = rsqrt( q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3] );

		q[0] *= norm;
		q[1] *= norm;
		q[2] *= norm;
		q[3] *= norm;



		if( updown == 0 )
		{
			U->leftSubgroupMult( i, j, &q );
		}
		else
		{
			q.hermitian();
			U->rightSubgroupMult( i, j, &q );
		}


	}
	else if( lc == MAG )
	{
		Matrix<Complex<Real>, 2 > subMat;

		if( updown == 0 && mu == 0 )
		{
			q.zero();
			subMat = U->getSubgroupMatrix( i, j );

			Complex<Real> c;

//			c = subMat.get( 0, 0 );
//			q[0] += c.x;
//			q[3] += c.y;
//			c = subMat.get( 1, 1 );
//			q[0] += c.x;
//			q[3] -= c.y;
//
//			c = subMat.get( 0, 1 );
//			q[2] += c.x;
//			q[1] += c.y;
//			c = subMat.get( 1, 0 );
//			q[2] -= c.x;
//			q[1] += c.y;
//
////			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
////			q*=.5;
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2]; //=D
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3]; //=E
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3]; //=F


			shA[id] += subMat.get(0,0).abs_squared() - subMat.get(0,1).abs_squared() - subMat.get(1,0).abs_squared() + subMat.get(1,1).abs_squared();
			c = subMat.get(0,0)*subMat.get(1,0).conj() - subMat.get(0,1)*subMat.get(1,1).conj();
			shA[id+NSB] += c.y;
			shA[id+2*NSB] += c.x;


		}
		__syncthreads();
		if( updown == 0 && mu == 1 )
		{
			subMat = U->getSubgroupMatrix( i, j );
			Complex<Real> c;
			shA[id] += subMat.get(0,0).abs_squared() - subMat.get(0,1).abs_squared() - subMat.get(1,0).abs_squared() + subMat.get(1,1).abs_squared();
			c = subMat.get(0,0)*subMat.get(1,0).conj() - subMat.get(0,1)*subMat.get(1,1).conj();
			shA[id+NSB] += c.y;
			shA[id+2*NSB] += c.x;




//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
////			q*=.5;
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 2 )
		{
			subMat = U->getSubgroupMatrix( i, j );
			Complex<Real> c;
			shA[id] += subMat.get(0,0).abs_squared() - subMat.get(0,1).abs_squared() - subMat.get(1,0).abs_squared() + subMat.get(1,1).abs_squared();
			c = subMat.get(0,0)*subMat.get(1,0).conj() - subMat.get(0,1)*subMat.get(1,1).conj();
			shA[id+NSB] += c.y;
			shA[id+2*NSB] += c.x;

//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
////			q*=.5;
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 3 )
		{
			subMat = U->getSubgroupMatrix( i, j );
			Complex<Real> c;
			shA[id] += subMat.get(0,0).abs_squared() - subMat.get(0,1).abs_squared() - subMat.get(1,0).abs_squared() + subMat.get(1,1).abs_squared();
			c = subMat.get(0,0)*subMat.get(1,0).conj() - subMat.get(0,1)*subMat.get(1,1).conj();
			shA[id+NSB] += c.y;
			shA[id+2*NSB] += c.x;


//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
////			q*=.5;
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();


		if( updown == 1 && mu == 0 )
		{

			subMat = U->getSubgroupMatrix( i, j );
			Complex<Real> c;
			shA[id] += subMat.get(0,0).abs_squared() - subMat.get(0,1).abs_squared() - subMat.get(1,0).abs_squared() + subMat.get(1,1).abs_squared();
			c = subMat.get(0,0).conj()*subMat.get(0,1) - subMat.get(1,0).conj()*subMat.get(1,1);
			shA[id+NSB] += c.y;
			shA[id+2*NSB] += c.x;


//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
////			q*=.5;
//			shA[id]    += q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 1 )
		{
			subMat = U->getSubgroupMatrix( i, j );
			Complex<Real> c;
			shA[id] += subMat.get(0,0).abs_squared() - subMat.get(0,1).abs_squared() - subMat.get(1,0).abs_squared() + subMat.get(1,1).abs_squared();
			c = subMat.get(0,0).conj()*subMat.get(0,1) - subMat.get(1,0).conj()*subMat.get(1,1);
			shA[id+NSB] += c.y;
			shA[id+2*NSB] += c.x;



//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
////			q*=.5;
//			shA[id]    += q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 2 )
		{
			subMat = U->getSubgroupMatrix( i, j );
			Complex<Real> c;
			shA[id] += subMat.get(0,0).abs_squared() - subMat.get(0,1).abs_squared() - subMat.get(1,0).abs_squared() + subMat.get(1,1).abs_squared();
			c = subMat.get(0,0).conj()*subMat.get(0,1) - subMat.get(1,0).conj()*subMat.get(1,1);
			shA[id+NSB] += c.y;
			shA[id+2*NSB] += c.x;


//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
////			q*=.5;
//			shA[id]    += q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 3 )
		{
			subMat = U->getSubgroupMatrix( i, j );
			Complex<Real> c;
			shA[id] += subMat.get(0,0).abs_squared() - subMat.get(0,1).abs_squared() - subMat.get(1,0).abs_squared() + subMat.get(1,1).abs_squared();
			c = subMat.get(0,0).conj()*subMat.get(0,1) - subMat.get(1,0).conj()*subMat.get(1,1);
			shA[id+NSB] += c.y;
			shA[id+2*NSB] += c.x;


//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
////			q*=.5;
//			shA[id]    += q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();


		if( mu == 0 )
		{
			if( updown == 0 )
			{
				//g0 -> D+sqrt(D*D+E*E+F*F):
//				shA[id] *= .5;
//				shA[id] -= 1./(Real)4.0;
				shA[id+NSB] *= 2.; // HV (2013-09-02) I think there should be a factor of 4 between E and D (and F and D) but factor of 2 finds better gff and Tuckers thesis says 2 is correct
				shA[id+2*NSB] *= 2.;
				shA[id+3*NSB] = 0;
				shA[id] = shA[id]+sqrt(shA[id]*shA[id]+shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]);


//				shA[id] *= .5;




				// make unitary
//				Real norm = rsqrt( shA[id]*shA[id]+shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB] );
//				if( norm > 1.0001 || norm < .999  ) printf("norm1 = %f\n", norm );
//				shA[id+NSB] *= norm;
//				shA[id+2*NSB] *= norm;
//				shA[id] *= norm;

				algo.calculateUpdate( shA, id );
			}
		}
		__syncthreads();

		q[0] = shA[id];
		q[1] = shA[id+NSB];
		q[2] = shA[id+2*NSB];
//		q[3] = 0; // should be zero here (because it was zero before) but it is not (if we use shA[id+3*NSB]! why?)
				// HV (2013-08-30): It is zero (only for RandomUpdate it is not! Therefore, don't set it to zero explicitly!
//		q[3] = shA[id+3*NSB];
		q[3] = 0;

		Real norm = rsqrt( q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3] );
//		if( norm > 1.0001 || norm < .999  ) printf("norm = %f\n", norm );
		if( q[3]*q[3] > 1e-12 ) printf("q[3]=%f\n", q[3] );



		q[0] *= norm;
		q[1] *= norm;
		q[2] *= norm;
		q[3] *= norm;

//		norm = rsqrt( q[0]*q[0]+q[1]*q[1]+q[2]*q[2] );

//		printf( "%f\n", norm );

		if( updown == 0 )
		{
			U->leftSubgroupMult( i, j, &q );
		}
		else
		{
			q.hermitian();
			U->rightSubgroupMult( i, j, &q );
		}
		//U->projectSU3withoutThirdRow();





//		if( updown == 0 && mu == 0 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3]; //=D
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3]; //=E
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3]; //=F
//		}
//		__syncthreads();
//		if( updown == 0 && mu == 1 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//		if( updown == 0 && mu == 2 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//		if( updown == 0 && mu == 3 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//
//
//		if( updown == 1 && mu == 0 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    += q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//		if( updown == 1 && mu == 1 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    += q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//		if( updown == 1 && mu == 2 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    += q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//		if( updown == 1 && mu == 3 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    += q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//
//
//		if( mu == 0 )
//		{
//			if( updown == 0 )
//			{
//				//g0 -> D+sqrt(D*D+E*E+F*F):
//				shA[id] /= 2.;
//				shA[id] -= 1./(Real)4.0;
//				shA[id] = shA[id]+sqrt(shA[id]*shA[id]+shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]);
////				shA[id+NSB] *= 2.;
////				shA[id+2*NSB] *= 2.;
//				algo.calculateUpdate( shA, id );
//			}
//		}
//		__syncthreads();
//
//		q[0] = shA[id];
//		q[1] = shA[id+NSB];
//		q[2] = shA[id+2*NSB];
//		q[3] = 0.0;
//
//		if( updown == 0 )
//		{
//			U->leftSubgroupMult( i, j, &q );
//		}
//		else
//		{
//			q.hermitian();
//			U->rightSubgroupMult( i, j, &q );
//		}


	}

}





#endif /* GAUGEFIXINGSUBGROUPSTEP_HXX_ */
