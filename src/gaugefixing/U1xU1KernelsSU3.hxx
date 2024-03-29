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
 */

#ifndef LANDAUKERNELSSU3_HXX_
#define LANDAUKERNELSSU3_HXX_

#include "GlobalConstants.h"
#include "kernel_launch_bounds.h"
#include "../lattice/datatype/datatypes.h"
#include "../lattice/datatype/lattice_typedefs.h"
#include "../lattice/rng/PhiloxWrapper.hxx"
#include "GaugeFixingSubgroupStep.hxx"
#include "algorithms/SaUpdate.hxx"
#include "algorithms/OrUpdate.hxx"
#include "algorithms/MicroUpdate.hxx"
#include "algorithms/RandomUpdate.hxx"
#include "../lattice/Matrix.hxx"
#include "../lattice/SU3.hxx"
#include "../lattice/Link.hxx"


// kernels as class members are not supported (even static): wrap the kernel calls and hide the kernels in namespace.
namespace U1xU1SU3
{
static const int Ndim = 4;
static const int Nc = 3;
__global__ void generateGaugeQualityPerSite( Real *U, double *dGff, double *dA );
__global__ void restoreThirdLine( Real* U, lat_index_t* nnt );
__global__ void randomTrafo( Real* U,lat_index_t* nnt, bool parity, int rngSeed, int rngCounter );
__global__ void orStep( Real* U, lat_index_t* nnt, bool parity, float orParameter );
__global__ void orStepSingleThread( Real* U, lat_index_t* nnt, bool parity, float orParameter );
__global__ void microStep( Real* U, lat_index_t* nnt, bool parity );
__global__ void saStep( Real* U, lat_index_t* nnt, bool parity, float temperature, int rngSeed, int rngCounter );
}

class U1xU1KernelsSU3
{
public:
//	__global__ static void heatbathStep( Real* UtDw, Real* Ut, Real* UtUp, lat_index_t* nnt, float beta, bool parity, int counter );

	// TODO remove static and make the init in constructor
	static void initCacheConfig()
	{
		cudaFuncSetCacheConfig( U1xU1SU3::generateGaugeQualityPerSite, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( U1xU1SU3::restoreThirdLine, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( U1xU1SU3::randomTrafo, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( U1xU1SU3::orStep, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( U1xU1SU3::microStep, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( U1xU1SU3::saStep, cudaFuncCachePreferL1 );
	}

	static void generateGaugeQualityPerSite( int a, int b, Real *U, double *dGff, double *dA )
	{
		U1xU1SU3::generateGaugeQualityPerSite<<<a,b>>>(U,dGff, dA);
	};
	static double getGaugeQualityPrefactorA()
	{
		return 1./(double)U1xU1SU3::Nc;
	};
	static double getGaugeQualityPrefactorGff()
	{
		return 1./(double)(U1xU1SU3::Nc*U1xU1SU3::Ndim);
	};

	static void restoreThirdLine( int a, int b, Real* U, lat_index_t* nnt )
	{
		U1xU1SU3::restoreThirdLine<<<a,b>>>(U,nnt);
	};
	static void randomTrafo( int a, int b, Real* U,lat_index_t* nnt, bool parity, int rngSeed, int rngCounter )
	{
		U1xU1SU3::randomTrafo<<<a,b>>>( U, nnt, parity, rngSeed, rngCounter );
	};
	static void orStep( int a, int b,  Real* U, lat_index_t* nnt, bool parity, float orParameter )
	{
		U1xU1SU3::orStep<<<a,b>>>( U, nnt, parity, orParameter );
	};
	static void microStep( int a, int b, Real* U, lat_index_t* nnt, bool parity )
	{
		U1xU1SU3::microStep<<<a,b>>>( U, nnt, parity );
	};
	static void saStep( int a, int b, Real* U, lat_index_t* nnt, bool parity, float temperature, int rngSeed, int rngCounter )
	{
		U1xU1SU3::saStep<<<a,b>>>( U, nnt, parity, temperature, rngSeed, rngCounter);
	};
private:
};





namespace U1xU1SU3
{

__global__ void generateGaugeQualityPerSite( Real *U, double *dGff, double *dA )
{
	typedef GpuPattern< SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;
	typedef Link<Gpu,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;

	SiteCoord<Ndim,FULL_SPLIT> s(DEVICE_CONSTANTS::SIZE);
	int site = blockIdx.x * blockDim.x + threadIdx.x;

	Matrix<Complex<Real>,Nc> locMatSum;
	SU3<Matrix<Complex<Real>,Nc> > Sum(locMatSum);

	Sum.zero();

	// TODO check if there is a faster way to compute DELTA
	for( int mu = 0; mu < 4; mu++ )
	{
		s.setLatticeIndex( site );

		Matrix<Complex<Real>,Nc> locMat;
		SU3<Matrix<Complex<Real>,Nc> > temp(locMat);

		TLink linkUp( U, s, mu );
		SU3<TLink> globUp( linkUp );

		temp.assignWithoutThirdLine( globUp );
//				temp.projectSU3withoutThirdRow();// TODO project here?
//				globUp.assignWithoutThirdLine( temp ); // TODO
		temp.reconstructThirdLine();
		Sum += temp;

		s.setNeighbour(mu,false);
		TLink linkDw( U, s, mu );
		SU3<TLink> globDw( linkDw );
		temp.assignWithoutThirdLine( globDw );
//				temp.projectSU3withoutThirdRow(); // TODO project here?
//				globDw.assignWithoutThirdLine( temp ); // TODO
		temp.reconstructThirdLine();
		Sum -= temp;
	}

	Sum -= Sum.trace()/Real(3.);

	Matrix<Complex<Real>,Nc> locMatSumHerm;
	SU3<Matrix<Complex<Real>,Nc> > SumHerm(locMatSumHerm);
	SumHerm = Sum;
	SumHerm.hermitian();

	Sum -= SumHerm;

	double prec = 0;

	for( int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			prec += Sum.get(i,j).abs_squared();
		}
	}
	dA[site] = prec;


	s.setLatticeIndex( site );
	double result = 0;

	Matrix<Complex<Real>,Nc> locTemp;
	SU3<Matrix<Complex<Real>,Nc> > temp(locTemp);
	for( int mu = 0; mu < 4; mu++ )
	{
		TLink linkUp( U, s, mu );
		SU3<TLink> globUp( linkUp );
		temp.assignWithoutThirdLine( globUp );  // TODO put this in the loop of dA to reuse globUp
		temp.reconstructThirdLine();
		result += temp.trace().x;
	}

	dGff[site] = result;

}

__global__ void restoreThirdLine( Real* U, lat_index_t* nnt )
{
	typedef GpuPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;
	typedef Link<Gpu,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;

//	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<Ndim,FULL_SPLIT> s(DEVICE_CONSTANTS::SIZE);
	s.nn = nnt;

	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	for( int mu = 0; mu < 4; mu++ )
	{
		TLink link( U, s, mu );
		SU3<TLink> glob( link );
		glob.projectSU3();
	}
}



template<class Algorithm> inline __device__ void apply( Real* U, lat_index_t* nn, bool parity, Algorithm algorithm  )
{
	typedef GpuPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(size);
	s.nn = nn;

	const bool updown = threadIdx.x / (4*NSB);
	const short mu = (threadIdx.x % (4*NSB)) / NSB;
	const short id = (threadIdx.x % (4*NSB)) % NSB;

	int site = blockIdx.x * blockDim.x/8 + id;
	if( parity == 1 ) site += s.getLatticeSize()/2;

	s.setLatticeIndex( site );
	if( updown==1 )
	{
		s.setNeighbour(mu,false);
	}

	Matrix<Complex<Real>,Nc> locMat;
	SU3<Matrix<Complex<Real>,Nc> > locU(locMat);

	TLinkIndex link( U, s, mu );

	SU3<TLinkIndex> globU( link );

	// make link local
	locU.assignWithoutThirdLine(globU);
	locU.reconstructThirdLine();

	GaugeFixingSubgroupStep<SU3<Matrix<Complex<Real>,Nc> >, Algorithm, U1xU1> subgroupStep( &locU, algorithm, id, mu, updown );

	// do the subgroup iteration
	SU3<Matrix<Complex<Real>,Nc> >::perSubgroup( subgroupStep );

	// copy link back
	globU.assignWithoutThirdLine(locU);
}

__global__ void __launch_bounds__(8*NSB,OR_MINBLOCKS) orStep( Real* U, lat_index_t* nnt, bool parity, float orParameter )
{
	OrUpdate overrelax( orParameter );
	apply( U, nnt, parity, overrelax );
}

__global__ void __launch_bounds__(8*NSB,MS_MINBLOCKS) microStep( Real* U, lat_index_t* nnt, bool parity )
{
	MicroUpdate micro;
	apply( U, nnt, parity, micro );
}

__global__ void __launch_bounds__(8*NSB,SA_MINBLOCKS) saStep( Real* U, lat_index_t* nnt, bool parity, float temperature, int rngSeed, int rngCounter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	SaUpdate sa( temperature, &rng );
	apply( U, nnt, parity, sa );
}

/**
 *  We do a lot of useless stuff here (gather a local functional value)
 *  but the random trafo is applied only once, so we don't care.
 */
__global__ void randomTrafo( Real* U, lat_index_t* nnt, bool parity, int rngSeed, int rngCounter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	RandomUpdate random( &rng );
	apply( U, nnt, parity, random );
}

}


#endif /* COULOMBKERNELSSU2_HXX_ */
