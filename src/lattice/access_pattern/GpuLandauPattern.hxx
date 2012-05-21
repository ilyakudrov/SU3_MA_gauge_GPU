/*
 * GpuTimeslicePattern.hxx
 *
 *  Created on: Apr 13, 2012
 *      Author: vogt
 */

#ifndef GPULANDAUPATTERN_HXX_
#define GPULANDAUPATTERN_HXX_

#include "../../util/cuda/cuda_host_device.h"
#include "../SiteCoord.hxx"
#include <assert.h>

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> class GpuLandauPattern
{
public:
	static const lat_group_dim_t Nc = T_Nc;
	static const lat_index_t Ndim = T_Ndim;
	CUDA_HOST_DEVICE static inline lat_index_t getSiteIndex( Site s );
	CUDA_HOST_DEVICE static inline lat_array_index_t getLinkIndex( Site s, lat_dim_t mu );
//	CUDA_HOST_DEVICE static inline lat_array_index_t getIndex( int linkIndex, int i, int j, int c );
	CUDA_HOST_DEVICE static inline lat_array_index_t getIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c );
	CUDA_HOST_DEVICE static inline lat_array_index_t getIndexByUnique( lat_array_index_t uniqueIndex, lat_coord_t size[T_Ndim] );
};

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_index_t GpuLandauPattern<Site, T_Ndim,T_Nc>::getSiteIndex( Site s )
{
	return  s.getLatticeIndex();
}

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t GpuLandauPattern<Site, T_Ndim,T_Nc>::getLinkIndex( Site s, lat_dim_t mu )
{
	int muSize = T_Nc*T_Nc*2*s.getLatticeSize();
	return s.getLatticeIndex()+mu*muSize;
}

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t GpuLandauPattern<Site, T_Ndim,T_Nc>::getIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c )
{
	return s.getLatticeIndex() + s.getLatticeSize()*( c + 2 * ( j + T_Nc *( i + T_Nc * mu ) ) );
}

/**
 * calculate the pattern index from unique index.
 */
template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t GpuLandauPattern<Site, T_Ndim, T_Nc>::getIndexByUnique( lat_array_index_t uniqueIndex, lat_coord_t size[T_Ndim] )
{

//	uniqueIndex /= site.getLatticeSize();
//	int mu = uniqueIndex /=


	bool c = uniqueIndex % 2;
	uniqueIndex /= 2;
	lat_group_dim_t j = uniqueIndex % T_Nc;
	uniqueIndex /= T_Nc;
	lat_group_dim_t i = uniqueIndex % T_Nc;
	uniqueIndex /= T_Nc;
	lat_dim_t mu = uniqueIndex % T_Ndim;
	uniqueIndex /= T_Ndim;
	lat_array_index_t latticeIndex = uniqueIndex;

//	printf( "c %d\n", c );
//	printf( "j %d\n", j );
//	printf( "i %d\n", i );
//	printf( "mu %d\n", mu );
//	printf( "latind %d\n", latticeIndex );


	Site s( size ); // parity must not be split for unique index
	s.setLatticeIndexFromNonParitySplitOrder( latticeIndex );

	return s.getLatticeIndex() + s.getLatticeSize()*( c + 2 * ( j + T_Nc *( i + T_Nc * mu ) ) );;
}

#endif /* GPUTIMESLICEPATTERN_HXX_ */