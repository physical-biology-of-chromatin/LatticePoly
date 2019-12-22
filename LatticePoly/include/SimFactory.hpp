//
//  SimFactory.hpp
//  LatticePoly
//
//  Created by mtortora on 22/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef SimFactory_hpp
#define SimFactory_hpp

#include <set>

#include "MCSim.hpp"


class SimFactory
{
public:
	SimFactory();
	
	IMCSim* GetSimulationInstance();
};

#endif /* SimFactory_hpp */
