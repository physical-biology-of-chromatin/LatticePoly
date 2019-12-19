//
//  main.cpp
//  LatticePoly
//
//  Created by mtortora on 21/10/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include "MCSim.hpp"


int main(int, const char**)
{
	try
	{
		MCSim<MCLiqLattice, MCHeteroPoly> sim;
		
		sim.Init();
		
		for ( int i = 1; i < Nmeas; i++ )
		{
			for ( int j = 0; j < Ninter; j++ )
			{
				sim.Run();
			}
						
			std::cout << "Performed " << sim.step << " out of "
			          << (Nmeas-1)*Ninter << " MC steps" << std::endl;
			
			sim.DumpVTK(i);
		}
	}
	
	catch ( std::exception& e )
	{
		std::cout << e.what() << std::endl;
		
		exit(1);
	}
	
	return 0;
}
