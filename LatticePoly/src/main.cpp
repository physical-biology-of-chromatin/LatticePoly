//
//  main.cpp
//  LatticePoly
//
//  Created by mtortora on 21/10/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <memory>

#include "SimFactory.hpp"
#include "InputParser.hpp"


int main(int argc, const char** argv)
{

	if ( argc != 2 )
	{
		std::cerr << "Syntax is " << argv[0] << " input_file" << std::endl;
		
		return 1;
	}
	
	try
	{

		// Parse input file
		InputParser parser(argv[1]);
		
		parser.ParseVars();
		
		// Initialise simulation
		std::unique_ptr<IMCSim> sim(SimFactory::GetSimulationInstance());
		
		sim->Init();

		// Run
		for ( int frame = sim->Ninit; frame < sim->Nfinal; ++frame )
		{


			// Print VTK frames every Ninter-th MC cycle beyond Nrelax
			if ( frame >= Nrelax )
			{
				sim->DumpVTK(frame);

			}

			for ( int i = 0; i < Ninter; ++i )
				sim->Run(frame);


			sim->PrintStats();


		}
		
		sim->DumpVTK(sim->Nfinal);

	}
	
	catch ( std::exception& e )
	{
		std::cerr << e.what() << std::endl;
		
		return -1;
	}
	
	return 0;
}
