//
//  MCReplicPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#ifndef MCReplicPoly_hpp
#define MCReplicPoly_hpp

#include "MCLivingPoly.hpp"

#include "MCLiqLattice.hpp"



class MCReplicPoly: public MCLivingPoly
{
public:
	MCReplicPoly(MCLattice*);
	
	void Init(int,int,int[3]);
	void TrialMove(double*);
	void ForkextraTrialMove(double*);
	double GetEffectiveEnergy();
	double GetCouplingForkEnergy(const std::vector<int>) const;
	std::vector<std::vector<int>> binded_particles;


	std::vector<int> active_cars;
	void AcceptMove();
	void OriginMove_explicit(const int[Ntot]);
	void OriginMove_implicit();

	void ForkMove();
	std::vector<int> dangling_ends;
	std::vector<MCTad*> cohesive_CARs;
	void LoadExtruders();
	void unLoadExtruders();

	std::vector<int> Spin_pos_toDelete;
	std::vector<int> Spin_pos_toCreate;
	
	std::vector<MCTad*> activeForks;
	std::vector<MCTad*> activeOrigins;
	std::vector<MCTad*> active_extruders;

	int NbindedForks;
	int NbindedCohesin;
	std::vector<int> MergedForkPos;
	int total_activated_cars;
	int individual_Nchain;
	int individual_Ndf;



protected:
	void Replicate(MCTad*);
	void TurnCohesive(MCTad*);
	void Find_cohesive_CAR();

	std::vector<double> PODLS;
	std::vector<double> ChIP;



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
	std::vector<int> loaded_mcms;

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
