//
//  MCPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCPoly_hpp
#define MCPoly_hpp

#include "MCTadUpdater.hpp"


class MCPoly
{
public:
	MCPoly(MCLattice*);
	virtual ~MCPoly();

	void Init(int);
	void GenerateRandom(int);
	
	void ToVTK(int);
	void FromVTK(int);

	void TrialMove(double*);
	void AcceptMove();
	void  OriginMove(MCTad*);
	void ForkMove();
	void GenerateCAR();
	
	int Ntad;
	int Nbond;
	std::vector<MCTad*> activeForks;
	std::vector<MCTad*> activeOrigins;

	std::vector<MCTad*> CAR;
	std::vector<MCTad*> interCAR;
	std::vector<int>  nreptation74;
	std::vector<int>  nreptation73;
	std::vector<int>  nreptation64;
	std::vector<int>  nreptation44;
	std::vector<int>  nreptation14;







	MCLattice* lat;
	MCTadUpdater* tadUpdater;

		
protected:
	MCTad* tadTrial;
	
	std::vector<MCTad> tadConf;
	std::vector<MCBond> tadTopo;

	void CreateBond(MCBond&);
	void FixPBCPair(std::vector<double3>&, int, int);
	
	virtual std::vector<double3> GetPBCConf();
	double3 GetPBCCenterMass(std::vector<double3>::iterator, std::vector<double3>::iterator);
	
	double3 centerMass;
};


#endif /* MCPoly_hpp */
