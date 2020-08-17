//
//  main.cpp
//  LatticePoly
//
//  Created by mtortora on 21/10/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <memory>
#include <fstream>
#include <sys/stat.h>

#include "SimFactory.hpp"
#include "InputParser.hpp"


int main(int argc, const char** argv)
{
	if ( argc != 2 )
	{
		std::cout << "Syntax is " << argv[0] << " input_file" << std::endl;
		
		return 1;
	}
	
	try
	{
		// Parse input file
		InputParser parser(argv[1]);
		
		parser.ParseVars();
		
		// Create output directory if necessary
		if ( mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1 )
		{
			if ( errno != EEXIST )
				throw std::runtime_error("main: Could not create directory " + outputDir);
		}
		
		std::cout << "Writing output to directory " << outputDir << std::endl;

		// Initialise simulation
		SimFactory::CheckInputOpt();

		std::unique_ptr<IMCSim> sim(SimFactory::GetSimulationInstance());
		
		sim->Init();
		
		// Run
		for ( int i = sim->Ninit; i < sim->Nfinal; ++i )
		{
			if ( i >= Nrelax )
				sim->DumpVTK(i);
			
			for ( int j = 0; j < Ninter; ++j )
				sim->Run();
			
			sim->PrintStats();
		}
		
		sim->DumpVTK(sim->Nfinal);
	}
	
	catch ( std::exception& e )
	{
		std::cout << e.what() << std::endl;
		
		return 1;
	}
	
	return 0;
}
