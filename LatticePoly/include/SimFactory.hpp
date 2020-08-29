//
//  SimFactory.hpp
//  LatticePoly
//
//  Created by mtortora on 22/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef SimFactory_hpp
#define SimFactory_hpp

#include "MCSim.hpp"


class SimFactory
{
public:
	static IMCSim* GetSimulationInstance();
	
private:
	static void CheckInputOpts();	
	static void CreateOutputDir();
};


#endif /* SimFactory_hpp */
