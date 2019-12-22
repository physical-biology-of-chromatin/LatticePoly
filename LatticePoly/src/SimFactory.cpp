//
//  SimFactory.cpp
//  LatticePoly
//
//  Created by mtortora on 22/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <set>

#include "SimFactory.hpp"


SimFactory::SimFactory()
{
	std::set<std::string> polyTypes;
	std::set<std::string> latticeTypes;
	
	polyTypes.insert("MCPoly");
	polyTypes.insert("MCHeteroPoly");
	
	latticeTypes.insert("MCLattice");
	latticeTypes.insert("MCLiqLattice");
	
	std::set<std::string>::iterator polyFind = polyTypes.find(polyType);
	std::set<std::string>::iterator latticeFind = latticeTypes.find(latticeType);
	
	if ( polyFind == polyTypes.end() ) throw std::runtime_error("SimFactory: invalid polyType");
	if ( latticeFind == latticeTypes.end() ) throw std::runtime_error("SimFactory: invalid latticeType");
}

IMCSim* SimFactory::GetSimulationInstance()
{
	if ( latticeType == "MCLattice" )
	{
		if ( polyType == "MCPoly" ) return new MCSim<MCLattice, MCPoly>;
		if ( polyType == "MCHeteroPoly" ) return new MCSim<MCLattice, MCHeteroPoly>;
	}
	
	if ( latticeType == "MCLiqLattice" )
	{
		if ( polyType == "MCPoly" ) return new MCSim<MCLiqLattice, MCPoly>;
		if ( polyType == "MCHeteroPoly" ) return new MCSim<MCLiqLattice, MCHeteroPoly>;
	}

	return NULL;
}
