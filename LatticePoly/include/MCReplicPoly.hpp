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
	void ForkextraTrialMove(double*);
	double GetEffectiveEnergy() const;
	void AcceptMove();
	void OriginMove();
	void ForkMove();
	

protected:
	void Replicate(MCTad*);
	void ReplicateTADs(MCTad*);
	void ReplicateBonds(MCTad*);
	int ReplTable[3][Ntot];
	void UpdateReplTable(MCTad*);


	
	void UnsetFork(MCTad*);
	void Update();

	virtual std::vector<double3> GetPBCConf();
	
	double Jint = 0.0;
	int Nfork;
	int MCsteps;
	int MCrepl;
	std::vector<int> origins;

	
	
	
private:
	void BuildPBCPair(std::vector<MCTad*>&, std::vector<double3>&, MCTad*, MCTad*);
};


#endif /* MCReplicPoly_hpp */
