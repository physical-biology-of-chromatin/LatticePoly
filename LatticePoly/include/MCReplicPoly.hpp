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
#include "MCLiqLattice.hpp"



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
	void MoveChoesin(MCTad*);

	int ReplTable[3][Ntot];
	void UpdateReplTable(MCTad*);


	
	void UnsetFork(MCTad*);
	void Update();

	virtual std::vector<double3> GetPBCConf();
	
	double Jint = 0.0;
	int Nfork;
	int MCsteps;
	int MCrepl;
	bool neigh;
	
	std::vector<int> origins;
	std::vector<double> mrt;
	std::vector<double> weights;
	std::vector<int> CAR;

	std::vector<int> anchor1;
	std::vector<int> anchor2;



	

	
	
	
private:
	void BuildPBCPair(std::vector<MCTad*>&, std::vector<double3>&, MCTad*, MCTad*);
};


#endif /* MCReplicPoly_hpp */
