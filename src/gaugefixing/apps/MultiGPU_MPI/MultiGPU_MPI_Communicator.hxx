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
 * compile with mpicc
 * 
 * This class takes care of the complete MPI communication:
 * 
 * -it offers methods to scatter and collect the gauge field from/to
 * the master host process to all process's devices.
 * -it applies an algorithm and takes care of the communication at 
 * the boundaries while hiding the time for the latter behind calculations
 * in the inner part of the domain.
 * -it offers methods setHot and projectSU3 that iterated over all
 * timeslices (no comm. necessary).
 * -it takes care of generating the gauge quality.
 * 
 */

#ifndef MULTIGPU_MPI_COMMUNICATOR_HXX_
#define MULTIGPU_MPI_COMMUNICATOR_HXX_

#include <stdio.h>
#include "../../../lattice/datatype/datatypes.h"
#include "../../../lattice/datatype/lattice_typedefs.h"
#include "../../GlobalConstants.h"
#include "./MultiGPU_MPI_LandauKernelsSU3.h"
#include "./MultiGPU_MPI_Reduce.h"
#include "./MultiGPU_MPI_AlgorithmOptions.h"

// MPI error handling macro
#define MPI_CHECK( call) \
    if((call) != MPI_SUCCESS) { \
        cerr << "MPI error calling \""#call"\"\n"; \
        MPI_Abort(MPI_COMM_WORLD, (-1) ); }
        

//TODO where to put these constants?
const int Ndim = 4;
const int Nc = 3;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;
const size_t timesliceSize = timesliceArraySize*sizeof(Real); 
        

template< class MultiGPU_MPI_GaugeKernels >
class MultiGPU_MPI_Communicator
{
public:
	// constructor
	MultiGPU_MPI_Communicator( int argc, char** argv );
	// destructor
	~MultiGPU_MPI_Communicator();
	// scatter the gauge field from 'master' to all other processes
	void scatterGaugeField( Real **dU, Real *U );
	// collect the gauge field from all processes to 'master'
	void collectGaugeField( Real **dU, Real *U );
	// get number of processes in MPI universe
	int getNumbProcs();
	// get rank
	int getRank();
	// get left neighbor rank
	int getLeft();
	// get right neighbor rank
	int getRight();
	// get min. timeslice the process takes care of
	int getMinTimeslice();
	// get max. timeslice the process takes care of
	int getMaxTimeslice();
	// is the process master?
	bool isMaster();
	// get numb. of timeslices the process takes care of
	int getNumbTimeslices();
	// get first slice of the six parts to hide the 6 parts of comm.
	int getStartPart( int beg );
	// get last slice of the six parts to hide the 6 parts of comm.
	int getEndPart( int end );
	// apply: applies algorithm, takes care of communication between devices
	void apply( Real** dU, lat_index_t** dNnt, bool evenodd, MultiGPU_MPI_AlgorithmOptions algoOptions );
	// projectSU3: iterates over the timeslices, no comm. needed.
	void projectSU3( Real** dU );
	// set dU to a random gauge field
	void setHot( Real** dU, MultiGPU_MPI_AlgorithmOptions algoOptions );
	// generate the gauge quality
	void generateGaugeQuality( Real** dU, lat_index_t** dNnt );
	// get the current value of the gauge functional
	double getCurrentGff();
	// get the current value of the gauge precision
	double getCurrentA();
	
private:
	// assign a device to the process
	void initDevice( const int device );
	int nprocs, rank, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	bool master;
	int lRank;
	int rRank;
	
	// MPI comm.
	MPI_Request request1, request2;
	MPI_Status  status;
	
	// useful variables
	int tmin;
	int tmax;
	int numbSlices;
	int startPart[6];
	int endPart[6];
	int theProcess[Nt];
	
	// device memory to collect the gauge fixing quality
	double *dGff;
	double *dA;
	double currentGff;
	double currentA;
	
	// halos
	Real* haloOut;
	Real* haloIn;
	Real* dHalo;
	
	// cudaStreams 
	cudaStream_t streamStd;
	cudaStream_t streamCpy;
};


template< class MultiGPU_MPI_GaugeKernels >
MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::MultiGPU_MPI_Communicator( int argc, char** argv )
{
	
	// initialize MPI communication
  MPI_CHECK( MPI_Init(&argc, &argv) );
  MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &nprocs) );
  MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &rank) );
  MPI_CHECK( MPI_Get_processor_name(processor_name, &namelen) );
	printf("Process %d on %s out of %d alive.\n", rank, processor_name, nprocs);
	
	// init. some variables
	lRank = ( rank - 1 + nprocs ) % nprocs;
	rRank = ( rank + 1 ) % nprocs;
	master = ( rank == 0 ? true : false );
	tmin = rank*Nt/nprocs;
	tmax = (rank+1)*Nt/nprocs;
	numbSlices = tmax-tmin;
	
	// set boarders of six parts to hide the comm. (apply,generateGaugeQuality)
	for( int l=0; l<6; l++ )
	{
		startPart[l] = tmin+1;
		endPart[l] = tmin+1;
	}
	for( int t=1; t<numbSlices; t++ )
	{
		for( int l=0; l<6; l++ )
		{
			if( l == (t-1)%6 )
			{
				endPart[l] += 1;
			}
			if( l > (t-1)%6 )
			{
				startPart[l] += 1;
				endPart[l] += 1;
			}
		}
	}
	
	// to which process belongs timeslice 't':
	for( int k=0; k<nprocs; k++ )
		for( int t = k*Nt/nprocs; t < (k+1)*Nt/nprocs; t++ )
			theProcess[t] = k;
	

	// print some information
	printf("Process %d: numbSlices %d\n", rank, numbSlices );
	for( int l=0; l<6; l++ )
	{
		printf("Process %d: startPart[%d] = %d, endPart[%d] = %d\n", rank, l, startPart[l], l, endPart[l] );
	}
	
	// init. the device
	initDevice( rank%4 );
	
	// init. cuda streams
	cudaStreamCreate( &streamStd );
	cudaStreamCreate( &streamCpy );
	
	// page-locked host memory for halo timeslices (two per thread)
 	if( nprocs > 1 ) cudaHostAlloc( &haloIn,  timesliceSize, 0 );
	if( nprocs > 1 ) cudaHostAlloc( &haloOut, timesliceSize, 0 );
	
	// device memory for halo timeslice (one per device)
	if( nprocs > 1 ) cudaMalloc( &dHalo, timesliceSize );
	
	// device memory for gauge quality
	cudaMalloc( &dGff, Nx*Ny*Nz*sizeof(double)/2 );
	cudaMalloc( &dA,   Nx*Ny*Nz*sizeof(double)/2 );
	
	// tell CUDA to prefer the L1 cache
	MultiGPU_MPI_GaugeKernels::initCacheConfig();
	
	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
}

template< class MultiGPU_MPI_GaugeKernels >
MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::~MultiGPU_MPI_Communicator()
{
	// free halo memory
 	cudaFreeHost( haloIn );
	cudaFreeHost( haloOut );
	cudaFree( dHalo );
	
	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
	MPI_CHECK( MPI_Finalize() );
}

template< class MultiGPU_MPI_GaugeKernels >
void MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::initDevice( const int device )
{
	cudaSetDevice(device);
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, device);
	printf("\nDevice %d: \"%s\"\n", device, deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);
}

template< class MultiGPU_MPI_GaugeKernels >
void MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::scatterGaugeField( Real **dU, Real *U )
{
	for( int t=0; t<Nt; t++ )
	{
		if( master && theProcess[t] == 0 )
		{
			cudaMemcpy( dU[t], &U[t*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
		}
		else if( theProcess[t] == rank )
		{
			MPI_CHECK( MPI_Recv( haloIn,  timesliceArraySize, MPI_Real, 0, 0, MPI_COMM_WORLD, &status ) );
			cudaMemcpy( dU[t], haloIn, timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
		}
		else if( master )
		{
			MPI_CHECK( MPI_Send( &U[t*timesliceArraySize], timesliceArraySize, MPI_Real, theProcess[t], 0, MPI_COMM_WORLD ) );
		}
		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
	}
}

template< class MultiGPU_MPI_GaugeKernels >
void MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::collectGaugeField( Real **dU, Real *U )
{
	// send back all timeslices to master
	for( int t=0; t<Nt; t++ )
	{
		if( master && theProcess[t] == 0 )
		{
			cudaMemcpy( &U[t*timesliceArraySize], dU[t], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
		}
		else if( theProcess[t] == rank )
		{
			cudaMemcpy( haloOut, dU[t], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
			MPI_CHECK( MPI_Send( haloOut,  timesliceArraySize, MPI_Real, 0, 0, MPI_COMM_WORLD ) );
		}
		else if( master )
		{
			MPI_CHECK( MPI_Recv( &U[t*timesliceArraySize], timesliceArraySize, MPI_Real, theProcess[t], 0, MPI_COMM_WORLD, &status ) );
		}
		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
	}
}

template< class MultiGPU_MPI_GaugeKernels >
int MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::getNumbProcs()
{
	return nprocs;
}

template< class MultiGPU_MPI_GaugeKernels >
int MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::getRank()
{
	return rank;
}

template< class MultiGPU_MPI_GaugeKernels >
int MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::getLeft()
{
	return lRank;
}

template< class MultiGPU_MPI_GaugeKernels >
int MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::getRight()
{
	return rRank;
}

template< class MultiGPU_MPI_GaugeKernels >
int MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::getMinTimeslice()
{
	return tmin;
}

template< class MultiGPU_MPI_GaugeKernels >
int MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::getMaxTimeslice()
{
	return tmax;
}

template< class MultiGPU_MPI_GaugeKernels >
bool MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::isMaster()
{
	return master;
}

template< class MultiGPU_MPI_GaugeKernels >
int MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::getNumbTimeslices()
{
	return numbSlices;
}

template< class MultiGPU_MPI_GaugeKernels >
int MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::getStartPart( int beg )
{
	return startPart[ beg ];
}

template< class MultiGPU_MPI_GaugeKernels >
int MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::getEndPart( int end )
{
	return endPart[ end ];
}

template< class MultiGPU_MPI_GaugeKernels >
inline void MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::apply( Real** dU, lat_index_t** dNnt, bool evenodd, MultiGPU_MPI_AlgorithmOptions algoOptions )
{
	static const int threadsPerBlock = NSB*8; // NSB sites are updated within a block (8 threads are needed per site)
	static const int numBlocks = Nx*Ny*Nz/2/NSB; // // half of the lattice sites (a parity) are updated in a kernel call
	
	// instantiate object of kernel wrapper class
	static MultiGPU_MPI_GaugeKernels kernelWrapper;

	
	if( nprocs > 1 )
	{
		int p_offset = evenodd?timesliceArraySize/2:0;
		
		// halo exchange forward step 1
		cudaMemcpyAsync( haloOut+p_offset, dU[tmax-1]+p_offset, timesliceSize/12, cudaMemcpyDeviceToHost, streamCpy );
		for( int t = startPart[2]; t < endPart[2]; t++ )
		{
			// call wrapper for one timelice
			kernelWrapper.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algoOptions );
		}
		cudaDeviceSynchronize(); // to ensure cudaMemcpyAsync finished
		
		// halo exchange forward step 2
		MPI_CHECK( MPI_Irecv( haloIn+p_offset,  timesliceArraySize/12, MPI_Real, lRank, 0, MPI_COMM_WORLD, &request2) );	
		MPI_CHECK( MPI_Isend( haloOut+p_offset, timesliceArraySize/12, MPI_Real, rRank, 0, MPI_COMM_WORLD, &request1) );
		for( int t = startPart[0]; t < endPart[0]; t++ )
		{
			kernelWrapper.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algoOptions );
		}
		MPI_CHECK( MPI_Wait( &request1, &status ) );
		MPI_CHECK( MPI_Wait( &request2, &status ) );
		
		// halo exchange forward step 3
		cudaMemcpyAsync( dHalo+p_offset, haloIn+p_offset, timesliceSize/12, cudaMemcpyHostToDevice, streamCpy );
		for( int t = startPart[3]; t < endPart[3]; t++ )
		{
			kernelWrapper.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algoOptions );		
		}
		
		// now call kernel wrapper for tmin with dU[t-1] replaced by dHalo
		kernelWrapper.applyOneTimeslice( numBlocks,threadsPerBlock, streamCpy, dU[tmin], dHalo, dNnt[rank], evenodd ^ (tmin%2) , algoOptions );
		
		// halo exchange back step 1
		cudaMemcpyAsync( haloOut+p_offset, dHalo+p_offset, timesliceSize/12, cudaMemcpyDeviceToHost, streamCpy );
		for( int t = startPart[4]; t < endPart[4]; t++ )
		{
			kernelWrapper.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algoOptions );
		}
		cudaDeviceSynchronize(); // to ensure cudaMemcpyAsync finished
		
		// halo exchange back step 2
		MPI_CHECK( MPI_Irecv( haloIn+p_offset,  timesliceArraySize/12, MPI_Real, rRank, 0, MPI_COMM_WORLD, &request2) );	
		MPI_CHECK( MPI_Isend( haloOut+p_offset, timesliceArraySize/12, MPI_Real, lRank, 0, MPI_COMM_WORLD, &request1) );
		for( int t = startPart[1]; t < endPart[1]; t++ )
		{
			kernelWrapper.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algoOptions );
		}
		MPI_CHECK( MPI_Wait( &request1, &status ) );
		MPI_CHECK( MPI_Wait( &request2, &status ) );
		
		// halo exchange back step 3
		cudaMemcpyAsync( dU[tmax-1]+p_offset, haloIn+p_offset, timesliceSize/12, cudaMemcpyHostToDevice, streamCpy );
		for( int t = startPart[5]; t < endPart[5]; t++ )
		{
			kernelWrapper.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algoOptions );
		}		
		cudaDeviceSynchronize(); // to ensure cudaMemcpyAsync finished
		
		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
	}
	else // nproc == 1
	{
		for( int t=tmin; t<tmax; t++ )
		{
			int tDw = ( t > 0 )?( t - 1 ):( Nt - 1 );
			kernelWrapper.applyOneTimeslice( numBlocks, threadsPerBlock, streamStd, dU[t], dU[tDw], dNnt[rank], evenodd ^ (t%2), algoOptions );
		}
		cudaDeviceSynchronize();
	} // end if nproc > 1
}


template< class MultiGPU_MPI_GaugeKernels >
inline void MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::projectSU3( Real** dU )
{
	int threadsPerBlock = 32;
	int numBlocks = Nx*Ny*Nz/32;
	
	// instantiate object of kernel wrapper class
	static MultiGPU_MPI_GaugeKernels kernelWrapper;

	for( int t=tmin; t<tmax; t++ )
	{
		kernelWrapper.projectSU3( numBlocks, threadsPerBlock, streamStd, dU[t] );
	}
	cudaDeviceSynchronize();
	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
}


template< class MultiGPU_MPI_GaugeKernels >
inline void MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::setHot( Real** dU, MultiGPU_MPI_AlgorithmOptions algoOptions )
{
	int threadsPerBlock = 32;
	int numBlocks = Nx*Ny*Nz/32;

	// instantiate object of kernel wrapper class
	static MultiGPU_MPI_GaugeKernels kernelWrapper;

	for( int t=tmin; t<tmax; t++ )
	{
		kernelWrapper.setHot( numBlocks, threadsPerBlock, streamStd, dU[t], algoOptions );
	}
	cudaDeviceSynchronize();
	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
}


template< class MultiGPU_MPI_GaugeKernels >
inline void MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::generateGaugeQuality( Real** dU, lat_index_t** dNnt )
{
	static const int threadsPerBlock = NSB; // NSB sites are updated within a block (8 threads are needed per site)
	static const int numBlocks = Nx*Ny*Nz/2/NSB; // // half of the lattice sites (a parity) are updated in a kernel call
	
	// instantiate object of kernel wrapper class
	static MultiGPU_MPI_GaugeKernels kernelWrapper;
	
	// reduce to collect dGff, dA
	static Reduce reduce( Nx*Ny*Nz/2 );
	
	double tempGff = 0.0;
	double tempA   = 0.0;
	
	if( nprocs > 1 )
	{
		for( int evenodd=0; evenodd<2; evenodd++ )
		{
			int p_offset = evenodd?timesliceArraySize/2:0;
			
			// halo exchange forward step 1
			cudaMemcpyAsync( haloOut+p_offset, dU[tmax-1]+p_offset, timesliceSize/12, cudaMemcpyDeviceToHost, streamCpy );
			for( int t = startPart[0]; t < endPart[1]; t++ )
			{
				// call wrapper for one timelice
				kernelWrapper.generateGaugeQualityPerSite( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , dGff, dA );
				tempGff += reduce.getReducedValue( streamStd, dGff );
				tempA   += reduce.getReducedValue( streamStd, dA );
			}
			cudaDeviceSynchronize(); // to ensure cudaMemcpyAsync finished
			
			// halo exchange forward step 2
			MPI_CHECK( MPI_Irecv( haloIn+p_offset,  timesliceArraySize/12, MPI_Real, lRank, 0, MPI_COMM_WORLD, &request2) );	
			MPI_CHECK( MPI_Isend( haloOut+p_offset, timesliceArraySize/12, MPI_Real, rRank, 0, MPI_COMM_WORLD, &request1) );
			for( int t = startPart[2]; t < endPart[3]; t++ )
			{
				kernelWrapper.generateGaugeQualityPerSite( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , dGff, dA );
				tempGff += reduce.getReducedValue( streamStd, dGff );
				tempA   += reduce.getReducedValue( streamStd, dA );
			}
			MPI_CHECK( MPI_Wait( &request1, &status ) );
			MPI_CHECK( MPI_Wait( &request2, &status ) );
			
			// halo exchange forward step 3
			cudaMemcpyAsync( dHalo+p_offset, haloIn+p_offset, timesliceSize/12, cudaMemcpyHostToDevice, streamCpy );
			for( int t = startPart[4]; t < endPart[5]; t++ )
			{
				kernelWrapper.generateGaugeQualityPerSite( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , dGff, dA );		
				tempGff += reduce.getReducedValue( streamStd, dGff );
				tempA   += reduce.getReducedValue( streamStd, dA );
			}
			
			// now call kernel wrapper for tmin with dU[t-1] replaced by dHalo
			kernelWrapper.generateGaugeQualityPerSite( numBlocks,threadsPerBlock, streamCpy, dU[tmin], dHalo, dNnt[rank], evenodd ^ (tmin%2) , dGff, dA );
			tempGff += reduce.getReducedValue( streamCpy, dGff );
			tempA   += reduce.getReducedValue( streamCpy, dA );
				
			cudaDeviceSynchronize();
			MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
		}
	}
	else // nproc == 1
	{
		for( int evenodd=0; evenodd<2; evenodd++ )
		{
			for( int t=tmin; t<tmax; t++ )
			{
				int tDw = ( t > tmin )?( t - 1 ):( tmax - 1 );
				kernelWrapper.generateGaugeQualityPerSite( numBlocks, threadsPerBlock, streamStd, dU[t], dU[tDw], dNnt[rank], evenodd ^ (t%2), dGff, dA );
				tempGff += reduce.getReducedValue( streamStd, dGff );
				tempA   += reduce.getReducedValue( streamStd, dA );
			}
			cudaDeviceSynchronize();
		}
	} // end if nproc > 1
	

	// collect dGff, dA from all devices
	if( nprocs > 1 )
	{
		MPI_CHECK( MPI_Allreduce( &tempGff, &currentGff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) );
		MPI_CHECK( MPI_Allreduce( &tempA,   &currentA,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) );
	}
	else
	{
		currentGff = tempGff;
		currentA = tempA;
	}
	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
}

template< class MultiGPU_MPI_GaugeKernels >
double MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::getCurrentGff()
{
	return currentGff/double(Nx*Ny*Nz*Nt)/(double)Ndim/(double)Nc;
}

template< class MultiGPU_MPI_GaugeKernels >
double MultiGPU_MPI_Communicator< MultiGPU_MPI_GaugeKernels >::getCurrentA()
{
	return currentA/double(Nx*Ny*Nz*Nt)/(double)Nc;
}



#endif /* MULTIGPU_MPI_COMMUNICATOR_HXX_ */
