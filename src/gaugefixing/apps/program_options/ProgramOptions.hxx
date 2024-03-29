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

#ifndef PROGRAMOPTIONS_HXX_
#define PROGRAMOPTIONS_HXX_

#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>

#include "../../../lattice/filetypes/filetype_typedefs.h"

#include <fstream>
#include <string>
#include <iostream>

class ProgramOptions
{
public:
	ProgramOptions();
	int init( int argc, char* argv[] );

	int getDeviceNumber() const {
		return deviceNumber;
	}

	std::string getFBasename() const {
		return fBasename;
	}

	std::string getFEnding() const {
		return fEnding;
	}

	int getFNumberformat() const {
		return fNumberformat;
	}

	int getFStartnumber() const {
		return fStartnumber;
	}

	int getFStepnumber() const {
		return fStepnumber;
	}

	FileType getFType() const {
		return fType;
	}

	std::string getFOutputAppendix() const {
		return fOutputAppendix;
	}

	int getGaugeCopies() const {
		return gaugeCopies;
	}

	bool isSetHot() const {
		return setHot;
	}

	int getNconf() const {
		return nconf;
	}

	bool isRandomTrafo() const {
		return randomTrafo;
	}

	int getCheckPrecision() const {
		return checkPrecision;
	}

	int getOrMaxIter() const {
		return orMaxIter;
	}

	float getOrParameter() const {
		return orParameter;
	}

	int getSrMaxIter() const {
		return srMaxIter;
	}

	float getSrParameter() const {
		return srParameter;
	}

	float getPrecision() const {
		return precision;
	}

	ReinterpretReal getReinterpret() const {
		return reinterpret;
	}

	int getReproject() const {
		return reproject;
	}

	float getSaMax() const {
		return saMax;
	}

	int getSaMicroupdates() const {
		return saMicroupdates;
	}

	float getSaMin() const {
		return saMin;
	}

	int getSaSteps() const {
		return saSteps;
	}

	long getSeed() const {
		return seed;
	}

	std::string get_output_SA_functional() const {
		return output_SA_functional;
	}

	std::string getOutConfPath() const {
        return outputConf;
    }

	std::string getOutputEnding() const {
		return outputEnding;
	}

	bool getSaveEach() const {
		return saveEach;
	}

	bool getDoSA() const {
		return doSA;
	}

private:
	boost::program_options::variables_map options_vm;
	boost::program_options::options_description options_desc;

//	int argc;
//	char* argv[];

	// variables
	std::string configFile;
	std::string output_SA_functional;
	std::string outputConf;
	std::string outputEnding;
	bool saveEach;
	bool doSA;

	int deviceNumber;

	FileType fType;
	std::string fBasename;
	std::string fEnding;
	int fNumberformat;
	int fStartnumber;
	int fStepnumber;
	int nconf;
	std::string fOutputAppendix;

	ReinterpretReal reinterpret;

	bool setHot;

	long seed;

	int gaugeCopies;
	bool randomTrafo;
	int reproject;

	int saSteps;
	float saMin;
	float saMax;
	int saMicroupdates;

	int orMaxIter;
	float orParameter;

	int srMaxIter;
	float srParameter;


	float precision;
	int checkPrecision;



};

ProgramOptions::ProgramOptions() : options_desc("Allowed options")//, argc(argc), argv(*argv)
{
}

int ProgramOptions::init( int argc, char* argv[] )
{
	options_desc.add_options()
			("help", "produce help message")

			("config-file", boost::program_options::value<std::string>(&configFile), "config file (command line arguments overwrite config file settings)")
			("output_SA_functional", boost::program_options::value<std::string>(&output_SA_functional), "output for temperature-functional data (part before numbering starts)")
			("output_conf", boost::program_options::value<std::string>(&outputConf), "path for output configuration (part before numbering starts)")
			("output_ending", boost::program_options::value<std::string>(&outputEnding)->default_value(""), "file ending to append to output_conf (default: "")")
			("save_each", boost::program_options::value<bool>(&saveEach)->default_value(false), "true - save each gauge copy, false - not (default: false)")
			("doSA", boost::program_options::value<bool>(&doSA)->default_value(true), "true - do simulated annealing, false - don't do (default: true)")

			("devicenumber,D", boost::program_options::value<int>(&deviceNumber)->default_value(-1), "number of the CUDA device (or -1 for auto selection)")

			("ftype", boost::program_options::value<FileType>(&fType), "type of configuration (PLAIN, HEADERONLY, VOGT, ILDG, QCDSTAG)")
			("fbasename", boost::program_options::value<std::string>(&fBasename), "file basename (part before numbering starts)")
			("fending", boost::program_options::value<std::string>(&fEnding)->default_value(".vogt"), "file ending to append to basename (default: .vogt)")
			("fnumberformat", boost::program_options::value<int>(&fNumberformat)->default_value(1), "number format for file index: 1 = (0,1,2,...,10,11), 2 = (00,01,...), 3 = (000,001,...),...")
			("fstartnumber", boost::program_options::value<int>(&fStartnumber)->default_value(0), "file index number to start from (startnumber, ..., startnumber+nconf-1")
			("fstepnumber", boost::program_options::value<int>(&fStepnumber)->default_value(1), "load every <fstepnumber>-th file")
			("nconf,m", boost::program_options::value<int>(&nconf)->default_value(1), "how many files to gaugefix")
			("fappendix", boost::program_options::value<std::string>(&fOutputAppendix)->default_value("gaugefixed_"), "appendix to be inserted beween (input-)filename and number")

			("reinterpret", boost::program_options::value<ReinterpretReal>(&reinterpret)->default_value(STANDARD), "reinterpret Real datatype (STANDARD = do nothing, FLOAT = read input as float and cast to Real, DOUBLE = ...)")

			("hotgaugefield", boost::program_options::value<bool>(&setHot)->default_value(false), "don't load gauge field; fill with random SU(3).")

			("seed", boost::program_options::value<long>(&seed)->default_value(1), "RNG seed")

			("gaugecopies", boost::program_options::value<int>(&gaugeCopies)->default_value(1), "Number of gauge copies")
			("randomtrafo", boost::program_options::value<bool>(&randomTrafo)->default_value(true), "do a random trafo before each gf run" )
			("reproject", boost::program_options::value<int>(&reproject)->default_value(100), "reproject every arg-th step")

			("sasteps", boost::program_options::value<int>(&saSteps)->default_value(1000), "number of SA steps")
			("samin", boost::program_options::value<float>(&saMin)->default_value(.01), "min. SA temperature")
			("samax", boost::program_options::value<float>(&saMax)->default_value(.4), "max. SA temperature")
			("microupdates", boost::program_options::value<int>(&saMicroupdates)->default_value(3), "number of microcanoncial updates at each SA temperature")

			("ormaxiter", boost::program_options::value<int>(&orMaxIter)->default_value(1000), "Max. number of OR iterations")
			("orparameter", boost::program_options::value<float>(&orParameter)->default_value(1.7), "OR parameter")

			("srmaxiter", boost::program_options::value<int>(&srMaxIter)->default_value(1000), "Max. number of SR iterations")
			("srparameter", boost::program_options::value<float>(&srParameter)->default_value(1.7), "SR parameter")

			("precision", boost::program_options::value<float>(&precision)->default_value(1E-7), "OR precision (dmuAmu)")
			("checkprecision", boost::program_options::value<int>(&checkPrecision)->default_value(100), "how often to check the gauge precision")
			;

	boost::program_options::positional_options_description options_p;
	options_p.add("config-file", -1);

	boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
			options(options_desc).positional(options_p).run(), options_vm);
	boost::program_options::notify(options_vm);

	std::ifstream cfg( configFile.c_str() );
	boost::program_options::store(boost::program_options::parse_config_file( cfg, options_desc), options_vm);
	boost::program_options::notify(options_vm);

	if (options_vm.count("help")) {
		std::cout << "Usage: " << argv[0] << " [options] [config-file]" << std::endl;
		std::cout << options_desc << "\n";
		return 1;
	}
	return 0;
}


#endif /* PROGRAMOPTIONS_HXX_ */
