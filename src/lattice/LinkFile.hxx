/**
 * Class for loading and storing configurations.
 * TODO I named this class LinkFile because ConfigurationFile may be missleading but I don't like the naming.
 *      Suggestions for a more precise name are very welcome.
 *
 * This class serves as a middleman between different array patterns on file and memory side.
 * The pattern for saving configurations is usually such that matrices are stored in contiguous regions the memory pattern
 * for CUDA is completely different. This template class accepts as arguments the FilePattern, i.e. the array pattern on file side,
 * the MemoryPattern.
 * To allow a high flexibility you can define Header and Footer classes that are called before/after the read of the configuration.
 *
 * The intension of this class is to cover all possible File- and MemoryPatterns. This causes some computation overhead which can
 * be avoided by spezializations of this class.
 * Check the getUnique()/setUnique() function calls in load() and save(): Each Pattern needs to be able to provide a unique index
 * wich does not depend on the pattern. Calculating back and forth this index maybe costly.
 * TODO Check how much it costs actually.
 *
 * TODO:
 *  - Is this class really flexible? Do we have to move the loading of the configuration to an extra class?
 *  - Implement a switch for promotion to double / demotion to float
 *  - Think about partially loading when we want to use the reconstruction technique for SU3.
 *
 * @author Hannes Vogt (hannes@havogt.de) Universitaet Tuebingen - Institut fuer Theoretische Physik
 * @date 2012-04-9
 */

#ifndef LINKFILE_HXX_
#define LINKFILE_HXX_

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../util/datatype/datatypes.h"

template<class Header, class Footer, class FilePattern, class MemoryPattern, class TheSite> class LinkFile
{
public:
	LinkFile();
	virtual ~LinkFile();
	bool load( TheSite site, std::string filename, Real *U );
	bool save( std::string filename, Real U );
};

template <class Header, class Footer, class FilePattern, class MemoryPattern, class TheSite> LinkFile<Header, Footer, FilePattern, MemoryPattern, TheSite>::LinkFile()
{
}

template <class Header, class Footer, class FilePattern, class MemoryPattern, class TheSite> LinkFile<Header, Footer, FilePattern, MemoryPattern, TheSite>::~LinkFile()
{
}

/**
 *
 */
template <class Header, class Footer, class FilePattern, class MemoryPattern, class TheSite> bool LinkFile<Header, Footer, FilePattern, MemoryPattern, TheSite>::load( TheSite site, std::string filename, Real *U )
{
	// open file
	std::fstream file;
	file.open( filename.c_str(), std::ios::in | std::ios::binary);

	if( !file )
	{
		util::Logger::log( util::ERROR, "Can't open file");
		util::Logger::log( util::ERROR, filename.c_str() );
		return false;
	}

	// load header

	Header header;
	if( !header.load( &file ) )
	{
		util::Logger::log( util::ERROR, "Can't read header");
	}

	// load config

	// calculate number of reals in configuration:
	int configSize = 1;
	// calculate lattice size
	for( int i = 0; i < site.Ndim; i++ )
	{
		configSize *= site.size[i];
	}

	std::cout << "lattice size: " << configSize << std::endl;

	// * number of directions mu * Nc^2 * 2 (complex number)
	configSize *= FilePattern::Ndim * FilePattern::Nc * FilePattern::Nc * 2;

	std::cout << "Ndim: " << FilePattern::Ndim << std::endl;

	std::cout << "total size: " << configSize << std::endl;

	Real temp;
	for( int i = 0; i < configSize; i++ )
	{
		// load number
		file.read( (char*)&temp, sizeof(Real) );
//		std::cout << temp << std::endl;
//		if( i < 18 ) std::cout << temp << std::endl;
		int index = MemoryPattern::getIndexByUnique( FilePattern::getUniqueIndex(i, site.size) , site.size );
		U[index] = temp;
//		break;
	}


	// load footer
//	if( !footer.load( file ) )
//	{
//		util::Logger::log( util::ERROR, "Can't read footer");
//	}


	file.close();
	return true;
}


#endif /* LINKFILE_HXX_ */