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
	double GetCouplingForkEnergy(const std::vector<int>) const;
	std::vector<std::vector<int>> binded_particles;


	
	void AcceptMove();
	void OriginMove(MCTad*);
	void ForkMove();
	std::vector<int> dangling_ends;
	std::vector<MCTad*> cohesive_CARs;



protected:
	void Replicate(MCTad*);
	void TurnCohesive();
	void Find_cohesive_CAR();


	void ReplicateTADs(MCTad*);
	void ReplicateBonds(MCTad*);
	int ReplTable[3][Ntot];
	void UpdateReplTable(MCTad*);


	
	void UnsetFork(MCTad*);
	void Update();
	
	virtual vtkSmartPointer<vtkPolyData> GetVTKData();
	virtual void SetVTKData(const vtkSmartPointer<vtkPolyData>);
	
	double Jint = 0.0;
	int Nfork;
	int MCsteps;
	int MCrepl;
	

	int lattice_neigh1[55];
	int lattice_neigh2[55];
	std::vector<int> origins;
	std::vector<double> mrt;
	std::vector<double> weights;
	

	std::vector<int> anchor1;
	std::vector<int> anchor2;



	

	
	
	
private:
	void BuildPBCPair(std::vector<MCTad*>&, std::vector<double3>&, MCTad*, MCTad*);
	int Norigin;

	std::vector<MCTad*> inactiveOrigins;	
};


#endif /* MCReplicPoly_hpp */
