//
//  main.cpp
//  LatticePoly
//
//  Created by mtortora on 21/10/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <memory>
#include <fstream>
#include <dirent.h>
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
		
		std::string dataDir = __DATA_PATH;
		
		outputDir = dataDir + "/" + outputDir;
		
		// Create output folder if necessary
		if ( mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1 )
		{
			if( errno != EEXIST )
			{
				throw std::runtime_error("Could not create folder " + outputDir);
			}
		}
		
		std::cout << "Writing output to folder " << outputDir << std::endl;
		
		// Clear output folder
		DIR* folder = opendir(outputDir.c_str());
		
		struct dirent* nextFile;
		char filePath[256];

		while ( (nextFile = readdir(folder)) != NULL )
		{
			sprintf(filePath, "%s/%s", outputDir.c_str(), nextFile->d_name);
			remove(filePath);
		}
		
		closedir(folder);
		
		// Copy input file to output folder
		std::string x;
		std::ifstream inFile(argv[1]);
		
		std::ofstream outFile(outputDir + "/input.cfg");
		
		while ( std::getline(inFile, x) )
		{
			outFile << x << std::endl;
		}
		
		inFile.close();
		outFile.close();
		
		// Initialise simulation
		SimFactory factory;
		
		std::unique_ptr<IMCSim> sim(factory.GetSimulationInstance());
		
		sim->Init();
		
		// Run
		for ( int i = 1; i < Nmeas; i++ )
		{
			for ( int j = 0; j < Ninter; j++ )
			{
				sim->Run();
			}
						
			std::cout << "Performed " << sim->step << " out of " << (Nmeas-1)*Ninter << " MC steps" << std::endl;
			
			sim->DumpVTK(i);
		}
	}
	
	catch ( std::exception& e )
	{
		std::cout << e.what() << std::endl;
		
		return 1;
	}
	
	return 0;
}
