/**
 * Lattice Site using a Nd-dimensional array.
 *
 * This implementation uses a Nd-dimensional array to store the lattice coordinates of the current site.
 * We should use the "SiteIndex" implementation whenever we go for minimal register usage.
 *
 * Possible Optimizations:
 *  - We use on the fly neighbour calculation here. Check if a table is a more favorable choice.
 *  - See "SiteIndex".
 *
 *
 * @author Hannes Vogt (hannes@havogt.de) Universitaet Tuebingen - Institut fuer Theoretische Physik
 * @date 2012-04-16
 */

#ifndef SITECOORD_HXX_
#define SITECOORD_HXX_

#include "../util/log/Logger.hxx"
#include "../util/cuda/cuda_host_device.h"
#include "../util/datatype/lattice_typedefs.h"

/**
 * The template parameter "bool par" defines normal indexing (par==false) or parity split indexing (par==true).
 * Parity split indexing means that the first half of the indices are the even sites ((x+y+z+...)%2 == 0) and the second half the odd sites.
 * It has to be used for CUDA implementations to ensure coalesced memory reads for neighbouring (concerning its index) array elements.
 *
 * @param template Nd: Dimension of the lattice
 * @param template par
 */
template<lat_dim_t Nd, bool par> class SiteCoord
{
public:
	CUDA_HOST_DEVICE inline SiteCoord( const lat_coord_t size[Nd] );
	CUDA_HOST_DEVICE inline SiteCoord( const SiteCoord<Nd,par> &s);
	CUDA_HOST_DEVICE inline virtual ~SiteCoord();
	CUDA_HOST_DEVICE inline lat_coord_t& operator[](lat_dim_t i);
	CUDA_HOST_DEVICE inline lat_index_t getLatticeIndex();
	CUDA_HOST_DEVICE inline void setLatticeIndex( lat_index_t latticeIndex );
	CUDA_HOST_DEVICE inline void setLatticeIndexFromParitySplitOrder( lat_index_t latticeIndex );
	CUDA_HOST_DEVICE inline void setLatticeIndexFromNonParitySplitOrder( lat_index_t latticeIndex );
	CUDA_HOST_DEVICE inline lat_index_t getLatticeSize();
	CUDA_HOST_DEVICE inline lat_index_t getLatticeIndexTimeslice();
	CUDA_HOST_DEVICE inline lat_index_t getLatticeSizeTimeslice();
	CUDA_HOST_DEVICE inline void setNeighbour( lat_dim_t direction, lat_coord_t steps );

	lat_coord_t size[Nd];
	static const lat_dim_t Ndim = Nd;
	lat_coord_t site[Nd];
};



template <lat_dim_t Nd, bool par> SiteCoord<Nd, par>::SiteCoord( const lat_coord_t size[Nd] )
{
	for( int i = 0; i < Nd; i++ )
	{
		this->size[i] = size[i];
	}
}

template <lat_dim_t Nd, bool par> SiteCoord<Nd, par>::SiteCoord( const SiteCoord<Nd,par> &s )
{
	for( int i = 0; i < Nd; i++ )
	{
		this->size[i] = s.size[i];
		this->site[i] = s.site[i];
	}
}

template <lat_dim_t Nd, bool par> SiteCoord<Nd, par>::~SiteCoord()
{
	// TODO Auto-generated destructor stub
}

template<lat_dim_t Nd, bool par> lat_coord_t& SiteCoord<Nd,par>::operator[](lat_dim_t i)
{
	return site[i];
}

template<lat_dim_t Nd, bool par> lat_index_t SiteCoord<Nd, par>::getLatticeIndex() // TODO parity ordering
{
	if( par )
	{
		lat_index_t parity = 0;
		for(  lat_dim_t i = 0; i < Nd; i++ )
		{
			parity += site[i];
		}

		lat_index_t index = 0;
		for( lat_dim_t i = 0; i < Nd; i++ )
		{
			index *= size[i];
			index += site[i];
		}

		if( parity % 2 == 0 )
		{
			return index / 2;
		}
		else
		{
			return index / 2 + getLatticeSize()/2;
		}
	}
	else
	{
		lat_index_t index = 0;
		for( lat_dim_t i = 0; i < Nd; i++ )
		{
			index *= size[i];
			index += site[i];
		}
		return index;
	}
}

template<lat_dim_t Nd, bool par> void SiteCoord<Nd, par>::setLatticeIndex( lat_index_t latticeIndex )
{
	// TODO simplify this
	bool curParity;
	if( par )
	{
		if( latticeIndex >= getLatticeSize() / 2 ) // odd parity
		{
			latticeIndex -= getLatticeSize() / 2;
			latticeIndex *= 2;
//			latticeIndex++;
			curParity = 1;
		}
		else // even parity
		{
			latticeIndex *= 2;
			curParity = 0;
		}
	}

//	printf( "latticeIndex: %d\n" , latticeIndex );

	lat_index_t parity = 0;
	for( lat_dim_t i = Nd-1; i >= 0; i-- )
	{
		site[i] = latticeIndex % size[i];
		parity += site[i];
		latticeIndex /= size[i];
	}


	if( par && (parity % 2 != curParity) )
		site[Nd-1]++;
}


template<lat_dim_t Nd, bool par> void SiteCoord<Nd, par>::setLatticeIndexFromParitySplitOrder( lat_index_t latticeIndex )
{
	// TODO BUG! Compare to setLatticeIndex()!!!

	if( latticeIndex >= getLatticeSize() / 2 )
	{
		latticeIndex -= getLatticeSize() / 2;
		latticeIndex *= 2;
		latticeIndex++;
	}
	else
	{
		latticeIndex *= 2;
	}

	for( lat_dim_t i = Nd-1; i >= 0; i-- )
	{
		site[i] = latticeIndex % size[i];
		latticeIndex /= size[i];
	}
}

template<lat_dim_t Nd, bool par> void SiteCoord<Nd, par>::setLatticeIndexFromNonParitySplitOrder( lat_index_t latticeIndex )
{
	for( lat_dim_t i = Nd-1; i >= 0; i-- )
	{
		site[i] = latticeIndex % size[i];
		latticeIndex /= size[i];
	}
}


template<lat_dim_t Nd, bool par> lat_index_t SiteCoord<Nd, par>::getLatticeSize()
{
	int tmp = 1;
	for( lat_dim_t i = 0; i < Nd; i++ )
	{
		tmp *= size[i];
	}
	return tmp;
}

/**
 * Index within a timeslice
 */
template<lat_dim_t Nd, bool par> lat_index_t SiteCoord<Nd, par>::getLatticeIndexTimeslice()
{
	if( par )
	{
		lat_index_t parity = 0;
		for(  lat_dim_t i = 1; i < Nd; i++ )
		{
			parity += site[i];
		}

		lat_index_t index = 0;
		for( lat_dim_t i = 1; i < Nd; i++ )
		{
			index *= size[i];
			index += site[i];
		}

		if( parity % 2 == 0 )
		{
			return index / 2;
		}
		else
		{
			return index / 2 + getLatticeSizeTimeslice()/2;
		}
	}
	else
	{
		lat_index_t index = 0;
		for( lat_dim_t i = 1; i < Nd; i++ )
		{
			index *= size[i];
			index += site[i];
		}
		return index;
	}
}

template<lat_dim_t Nd, bool par> lat_index_t SiteCoord<Nd, par>::getLatticeSizeTimeslice()
{
	int tmp = 1;
	for( lat_dim_t i = 1; i < Nd; i++ )
	{
		tmp *= size[i];
	}
	return tmp;
}

/**
 * If you think this for-loop looks strange: This is a way the compiler can already index the array! -> The site[Ndim] array is not placed in (CUDA) memory
 */
template<lat_dim_t Nd, bool par> void SiteCoord<Nd, par>::setNeighbour( lat_dim_t direction, lat_coord_t steps )
{
	for( lat_dim_t i = 0; i < Nd; i++ )
	{
		if( direction == i )
		{
			site[i] += steps;

			if( site[i] >= size[i] )
			{
				site[i] -= size[i];
			}
			else
			{
				if( site[i] < 0 )
				{
					site[i] += size[i];
				}
			}
		}
	}
}


#endif /* SITECOORD_HXX_ */