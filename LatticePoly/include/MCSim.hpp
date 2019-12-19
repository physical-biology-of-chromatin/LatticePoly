//
//  MCSim.hpp
//  LatticePoly
//
//  Created by mtortora on 29/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCSIM_HPP_
#define MCSIM_HPP_

#include "MCHeteroPoly.hpp"


template<class lattice, class polymer>
class MCSim
{
public:
	MCSim();
	~MCSim();

	void Init();
	void InitRNG();
	void Run();
	void Update(int);
	void DumpVTK(int);
	
	bool MetropolisMove(double);
	bool ArrheniusMove(double, double);

	int step;
	
private:
	lattice* lat;
	polymer* pol;
	
	std::mt19937_64 rngEngine;
	std::uniform_real_distribution<double> rngDistrib{0., 1.};
};


#endif /* MCSim_hpp */
