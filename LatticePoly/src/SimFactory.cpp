//
//  SimFactory.cpp
//  LatticePoly
//
//  Created by mtortora on 22/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <set>
#include <fstream>
#include <sys/stat.h>

#include "SimFactory.hpp"


IMCSim* SimFactory::GetSimulationInstance()
{
	SimFactory::CheckInputOpts();
	SimFactory::CreateOutputDir();

	if ( (latticeType == "MCLattice") && ( polyType == "MCPoly" ) )
		return new MCSim<MCLattice, MCPoly>;
	
	else if ( (latticeType == "MCLattice") && ( polyType == "MCHeteroPoly" ) )
		return new MCSim<MCLattice, MCHeteroPoly>;
	
	else if ( (latticeType == "MCLattice") && ( polyType == "MCLivingPoly" ) )
		return new MCSim<MCLattice, MCLivingPoly>;
	
	else if ( (latticeType == "MCLattice") && ( polyType == "MCReplicPoly" ) )
		return new MCSim<MCLattice, MCReplicPoly>;
	
	else if ( (latticeType == "MCLiqLattice") && ( polyType == "MCHeteroPoly" ) )
		return new MCSim<MCLiqLattice, MCHeteroPoly>;
	
	else if ( (latticeType == "MCLiqLattice") && ( polyType == "MCLivingPoly" ) )
		return new MCSim<MCLiqLattice, MCLivingPoly>;
	
	else if ( (latticeType == "MCLiqLattice") && ( polyType == "MCReplicPoly" ) )
		return new MCSim<MCLiqLattice, MCReplicPoly>;
	
	else if ( (latticeType == "MCLattice") && ( polyType == "MCHeteroPoly_looped" ) )
		return new MCSim<MCLattice, MCHeteroPoly_looped>;
	
	else
		throw std::runtime_error("SimFactory: Unsupported combination of polyType and latticeType");
	
	return nullptr;
}

void SimFactory::CheckInputOpts()
{
	std::set<std::string> polyTypes;
	std::set<std::string> latticeTypes;
	
	polyTypes.insert("MCPoly");
	
	polyTypes.insert("MCHeteroPoly");
	polyTypes.insert("MCHeteroPoly_looped");

	polyTypes.insert("MCLivingPoly");
	polyTypes.insert("MCReplicPoly");

	latticeTypes.insert("MCLattice");
	latticeTypes.insert("MCLiqLattice");

	auto polyFind = polyTypes.find(polyType);
	auto latticeFind = latticeTypes.find(latticeType);
	
	if ( polyFind == polyTypes.end() )
		throw std::runtime_error("SimFactory: Invalid polyType");
	
	if ( latticeFind == latticeTypes.end() )
		throw std::runtime_error("SimFactory: Invalid latticeType");
}

void SimFactory::CreateOutputDir()
{
	int status = 0;
	size_t pos = 0;
		
	while ( (status == 0) && (pos != std::string::npos) )
	{
		pos = outputDir.find("/", pos+1);
		std::string path = outputDir.substr(0, pos);

		if ( (status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)) == -1 )
		{
			if ( errno == EEXIST )
				status = 0;
			else
				throw std::runtime_error("SimFactory: Could not create directory " + outputDir + " (error code " + std::to_string(errno) + ")");
		}
	}
}
