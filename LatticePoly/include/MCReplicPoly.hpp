//
//  MCReplicPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#ifndef MCReplicPoly_hpp
#define MCReplicPoly_hpp

#include "MCHeteroPoly.hpp"


class MCReplicPoly: public MCHeteroPoly
{
public:
	MCReplicPoly(MCLattice*);
	
	void Init(int);
	void TrialMove(double*);
	
protected:
	void Replicate(MCTad*);
	void ReplicateTADs(MCTad*);
	void ReplicateBonds(MCTad*);
	
	void UnsetFork(MCTad*);
	void Update();

	virtual std::vector<double3> GetPBCConf();
	
	int Nfork;
	int Norigin;

	std::vector<MCTad*> activeForks;
	std::vector<MCTad*> inactiveOrigins;
	
	double3 centerMass2;

private:
	void BuildPBCPair(std::vector<MCTad*>&, std::vector<double3>&, MCTad*, MCTad*);
};


#endif /* MCReplicPoly_hpp */
