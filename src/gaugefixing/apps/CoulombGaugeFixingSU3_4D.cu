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

#include <iostream>
#include <math.h>
#include <sstream>
#ifndef OSX
#include "malloc.h"
#endif
#include "../GlobalConstants.h"
#include "../GaugeFixingStats.hxx"
#include "../../lattice/access_pattern/StandardPattern.hxx"
#include "../../lattice/access_pattern/GpuPattern.hxx"
#include "../../lattice/access_pattern/GpuPatternTimeslice.hxx"
#include "../../lattice/SiteCoord.hxx"
#include "../../lattice/SiteIndex.hxx"
#include "../../lattice/LinkFile.hxx"
#include "../../util/timer/Chronotimer.h"
#include "../../lattice/filetypes/FileVogt.hxx"
#include "../../lattice/filetypes/FilePlain.hxx"
#include "../../lattice/filetypes/FileHeaderOnly.hxx"
#include "../../lattice/filetypes/filetype_typedefs.h"
#include "../../lattice/datatype/lattice_typedefs.h"
#include "../../lattice/rng/PhiloxWrapper.hxx"
#include "program_options/ProgramOptions.hxx"
#include "program_options/FileIterator.hxx"
#include "../CoulombKernelsSU3.hxx"
#include "../CommonKernelsSU3.hxx"

using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 3;

const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef StandardPattern<SiteCoord<Ndim,NO_SPLIT>,Ndim,Nc> Standard;
typedef GpuPattern< SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> GpuTimeslice;
typedef GpuPatternTimeslice<SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;

void readILDG_timeslice(SiteCoord<4,FULL_SPLIT> s, const char *file_name, const short SIZE[4], Real *U);
void writeILDG_timeslice(SiteCoord<4,FULL_SPLIT> s, const char *file_name, const char *output_name, const short SIZE[4], Real *U, int steps);

bool readQCDSTAG_timeslice(SiteCoord<4,FULL_SPLIT> s, const char *file_name, const short SIZE[4], Real *U);
bool writeQCDSTAG_timeslice(SiteCoord<4,FULL_SPLIT> s, const char *output_name, const short SIZE[4], Real *U);

int main(int argc, char* argv[])
{
	Chronotimer allTimer;
	allTimer.reset();
	allTimer.start();

	CoulombKernelsSU3::initCacheConfig();

	// read configuration from file or command line
	ProgramOptions options;
	int returncode = options.init( argc, argv );
	if( returncode != 0 ) return returncode;

	// Choose device and print device infos
	cudaDeviceProp deviceProp;
	int selectedDeviceNumber;
	if( options.getDeviceNumber() >= 0 )
	{
		cudaSetDevice( options.getDeviceNumber() );
		selectedDeviceNumber = options.getDeviceNumber();
	}
	else
	{
		cudaGetDevice( &selectedDeviceNumber );
	}
	cudaGetDeviceProperties(&deviceProp, selectedDeviceNumber );

	printf("\nDevice %d: \"%s\"\n", selectedDeviceNumber, deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);



	// SiteCoord is faster than SiteIndex when loading files
	SiteCoord<4,FULL_SPLIT> s(HOST_CONSTANTS::SIZE);


	// allocate Memory
	// host memory for configuration
	Real* U = (Real*)malloc( arraySize*sizeof(Real) );

	// device memory for timeslice t
	Real* dUtUp;
	cudaMalloc( &dUtUp, timesliceArraySize*sizeof(Real) );

	// device memory for timeslice t-1
	Real* dUtDw;
	cudaMalloc( &dUtDw, timesliceArraySize*sizeof(Real) );

	// host memory for the timeslice neighbour table (since we are working on timeslices, this is enough on device side)
	lat_index_t* nnt = (lat_index_t*)malloc( s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );

	// device memory for the timeslice neighbour table
	lat_index_t *dNnt;
	cudaMalloc( &dNnt, s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof( lat_index_t ) );

	// initialise the timeslice neighbour table
	SiteIndex<4,FULL_SPLIT> sTemp(HOST_CONSTANTS::SIZE_TIMESLICE);
	sTemp.calculateNeighbourTable( nnt );

	// copy neighbour table to device
	cudaMemcpy( dNnt, nnt, s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof( lat_index_t ), cudaMemcpyHostToDevice );



	// TODO maybe we should choose the filetype on compile time
	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfHeaderOnly( options.getReinterpret() );
	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfVogt( options.getReinterpret() );
	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfPlain( options.getReinterpret() );


	int threadsPerBlock = NSB*8; // NSB sites are updated within a block (8 threads are needed per site)
	int numBlocks = s.getLatticeSizeTimeslice()/2/NSB; // half of the lattice sites (a parity) are updated in a kernel call


	GaugeFixingStats<Ndim,Nc,CoulombKernelsSU3,AVERAGE> gaugeStats( dUtUp, HOST_CONSTANTS::SIZE_TIMESLICE );

	// timer to measure kernel times
	Chronotimer kernelTimer;
	kernelTimer.reset();
	kernelTimer.start();

	double orTotalKernelTime = 0; // sum up total kernel time for OR
	long orTotalStepnumber = 0;
	double saTotalKernelTime = 0;

	FileIterator fi( options );
	for( fi.reset(); fi.hasNext(); fi.next() )
	{
		bool loadOk;

		// ofstream output;
		// output.precision(17);
		// output.open("SA_test");
		// output << "temperature,functional,time"<<endl;

		if( !options.isSetHot() ) // load a file
		{

			switch(  options.getFType() )
			{
			case VOGT:
				loadOk = lfVogt.load( s, fi.getFilename(), U );
				break;
			case PLAIN:
				loadOk = lfPlain.load( s, fi.getFilename(), U );
				break;
			case HEADERONLY:
				loadOk = lfHeaderOnly.load( s, fi.getFilename(), U );
				break;
			case ILDG:
                loadOk = true;
                readILDG_timeslice(s, fi.getFilename().c_str(), HOST_CONSTANTS::SIZE, U);
                break;
            case QCDSTAG:
                loadOk = readQCDSTAG_timeslice(s, fi.getFilename().c_str(), HOST_CONSTANTS::SIZE, U);
                break;
			default:
				cout << "Filetype not set to a known value. Exiting...";
				exit(1);
			}

			if( !loadOk )
			{
				cout << "Error while loading. Trying next file." << endl;
				break;
			}
			else
			{
				cout << "File loaded." << endl;
			}
		}

		for( int t = 0; t < s.size[0]; t++ )
		{
			int tDw = (t > 0)?(t-1):(s.size[0]-1); // calculating t-1 (periodic boundaries)

			double bestGff = 0.0;
			for( int copy = 0; copy < options.getGaugeCopies(); copy++ )
			{
				if( !options.isSetHot() ) // if we want a hot random configuration we do not need to copy
				{
					// copying timeslice t ...
					cudaMemcpy( dUtUp, &U[t*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
					// ... and t-1 to device
					cudaMemcpy( dUtDw, &U[tDw*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
				}
				else // randomize the timeslices
				{
					CommonKernelsSU3::setHot( s.getLatticeSizeTimeslice()/32,32, dUtUp, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice(), options.getSeed(), PhiloxWrapper::getNextCounter() );
					CommonKernelsSU3::setHot( s.getLatticeSizeTimeslice()/32,32, dUtDw, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice(), options.getSeed(), PhiloxWrapper::getNextCounter() );
				}


				if( options.isRandomTrafo() )
				{
					CoulombKernelsSU3::randomTrafo(numBlocks,threadsPerBlock, dUtUp, dUtDw, dNnt, 0, options.getSeed(), PhiloxWrapper::getNextCounter() );
					CoulombKernelsSU3::randomTrafo(numBlocks,threadsPerBlock, dUtUp, dUtDw, dNnt, 1, options.getSeed(), PhiloxWrapper::getNextCounter() );
				}


				// calculate and print the gauge quality
				printf( "i:\t\tgff:\t\tdA:\n");
				gaugeStats.generateGaugeQuality();


				// SIMUALTED ANNEALING
				if( options.getSaSteps() > 0 ) printf( "SIMULATED ANNEALING\n" );
				float temperature = options.getSaMax();
				float tempStep = (options.getSaMax()-options.getSaMin())/(float)options.getSaSteps();

				kernelTimer.reset();
				kernelTimer.start();
				for( int i = 0; i < options.getSaSteps(); i++ )
				{
					CoulombKernelsSU3::saStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 0, temperature, options.getSeed(), PhiloxWrapper::getNextCounter() );
					CoulombKernelsSU3::saStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 1, temperature, options.getSeed(), PhiloxWrapper::getNextCounter() );

					for( int mic = 0; mic < options.getSaMicroupdates(); mic++ )
					{
						CoulombKernelsSU3::microStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 0 );
						CoulombKernelsSU3::microStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 1 );
					}

					if( i % options.getReproject() == 0 )
					{
						CommonKernelsSU3::projectSU3( s.getLatticeSizeTimeslice()/32, 32, dUtUp, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );
						CommonKernelsSU3::projectSU3( s.getLatticeSizeTimeslice()/32, 32, dUtDw, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );
					}

					if( i % options.getCheckPrecision() == 0 )
					{
						gaugeStats.generateGaugeQuality();
						printf( "%f\t\t%1.10f\t\t%e\n", temperature, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );
						// output<<temperature<<","<<gaugeStats.getCurrentGff()<<","<<t<<endl;
					}
					temperature -= tempStep;
				}
				cudaDeviceSynchronize();
				kernelTimer.stop();
				saTotalKernelTime += kernelTimer.getTime();


				// OVERRELAXATION
				if( options.getOrMaxIter() > 0 ) printf( "OVERRELAXATION\n" );
				kernelTimer.reset();
				kernelTimer.start();
				for( int i = 0; i < options.getOrMaxIter(); i++ )
				{
					CoulombKernelsSU3::orStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 0, options.getOrParameter() );
					CoulombKernelsSU3::orStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 1, options.getOrParameter() );

					if( i % options.getReproject() == 0 )
					{
						CommonKernelsSU3::projectSU3( s.getLatticeSizeTimeslice()/32, 32, dUtUp, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );
						CommonKernelsSU3::projectSU3( s.getLatticeSizeTimeslice()/32, 32, dUtDw, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );
					}

					if( i % options.getCheckPrecision() == 0 )
					{
						gaugeStats.generateGaugeQuality();
						printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );
						if( gaugeStats.getCurrentA() < options.getPrecision() ) break;
					}

					orTotalStepnumber++;
				}

				cudaDeviceSynchronize();
				kernelTimer.stop();
				orTotalKernelTime += kernelTimer.getTime();


				// reconstruct third line before copy back
				CommonKernelsSU3::projectSU3( s.getLatticeSizeTimeslice()/32, 32, dUtUp, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );
				CommonKernelsSU3::projectSU3( s.getLatticeSizeTimeslice()/32, 32, dUtDw, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );

				if( gaugeStats.getCurrentGff() > bestGff )
				{
					cout << "FOUND BETTER COPY" << endl;
					bestGff = gaugeStats.getCurrentGff();

					// copy back
					cudaMemcpy( &U[t*timesliceArraySize], dUtUp, timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
					cudaMemcpy( &U[tDw*timesliceArraySize], dUtDw, timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
				}
				else
				{
					cout << "NO BETTER COPY" << endl;
				}
			}
		}

		//saving file
		if( !options.isSetHot() )
		{
			cout << "saving " << fi.getOutputFilename() << " as " << options.getFType() << endl;

			switch( options.getFType() )
			{
			case VOGT:
				loadOk = lfVogt.save( s, fi.getOutputFilename(), U );
				break;
			case PLAIN:
				loadOk = lfPlain.save( s, fi.getOutputFilename(), U );
				break;
			case HEADERONLY:
				loadOk = lfHeaderOnly.save( s, fi.getOutputFilename(), U );
				break;
			case ILDG:
                loadOk = true;
                writeILDG_timeslice(s, fi.getFilename().c_str(), fi.getOutputFilename().c_str(), HOST_CONSTANTS::SIZE, U, options.getSaSteps());
                break;
            case QCDSTAG:
                loadOk = writeQCDSTAG_timeslice(s, fi.getOutputFilename().c_str(), HOST_CONSTANTS::SIZE, U);
                break;
			default:
				cout << "Filetype not set to a known value. Exiting";
				exit(1);
			}
		}
	}

	allTimer.stop();
	cout << "total time: " << allTimer.getTime() << " s" << endl;


	long hbFlops = 2252+86-16;
	long microFlops = 2252+14-16;
	cout << "Simulated Annealing (HB+Micro): " << (double)((long)(hbFlops+microFlops*options.getSaMicroupdates())*(long)s.getLatticeSize()*(long)options.getSaSteps()*(long)options.getGaugeCopies())/saTotalKernelTime/1.0e9 << " GFlops at "
					<< (double)((long)192*(long)s.getLatticeSize()*options.getSaSteps()*(options.getSaMicroupdates()+1)*(long)sizeof(Real))/saTotalKernelTime/1.0e9 << "GB/s memory throughput." << endl;


	long orFlops = 2252+22-16;
	cout << "Overrelaxation: " << (double)((long)orFlops*(long)s.getLatticeSize()*(long)orTotalStepnumber)/orTotalKernelTime/1.0e9/(double)s.size[0] << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(orTotalStepnumber)*(long)sizeof(Real))/orTotalKernelTime/1.0e9/(double)s.size[0] << "GB/s memory throughput." << endl;
}
