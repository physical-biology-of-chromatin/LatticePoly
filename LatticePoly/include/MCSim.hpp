//
//  MCSim.hpp
//  LatticePoly
//
//  Created by mtortora on 29/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCSIM_HPP_
#define MCSIM_HPP_

#include <chrono>

#include "MCUpdater.hpp"



class IMCSim
{
public:
	virtual ~IMCSim() {};

	virtual void Init() = 0;
	virtual void Run() = 0;
	virtual void DumpVTK(int) = 0;
	virtual void PrintStats() = 0;

	int Ninit;
	int Nfinal;
};


template<class lattice, class polymer>
class MCSim: public IMCSim
{
public:
	MCSim();
	~MCSim();
	
	void Init();
	void Run();
	void DumpVTK(int);
	void PrintStats();
	
private:
	void InitRNG();
	void InitSimRange();

	double acceptAveLiq;
	double acceptAvePoly;
	
	int NliqMoves;
	
	unsigned long long cycle;
	unsigned long long acceptCountLiq;
	unsigned long long acceptCountPoly;
	
	lattice* lat;
	polymer* pol;
	
    std::chrono::high_resolution_clock::time_point tStart;
    std::chrono::high_resolution_clock::time_point tCycle;
};


#endif /* MCSim_hpp */
