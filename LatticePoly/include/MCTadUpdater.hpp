//
//  MCTadUpdater.hpp
//  LatticePoly
//
//  Created by mtortora on 01/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#ifndef MCTadUpdater_hpp
#define MCTadUpdater_hpp

#include "MCTad.hpp"


class MCTadUpdater
{
public:
	MCTadUpdater(MCLattice*);
	
	void TrialMove(const MCTad*, double*);
	void AcceptMove(MCTad*) const;

	void TrialMoveLeftEnd(const MCTad*, double*);
	void TrialMoveRightEnd(const MCTad*, double*);
	void TrialMoveLinear(const MCTad*, double*);
	void TrialMoveFork(const MCTad*, double*);
	
	void TrialMoveTopo(const MCTad*, double*);
	void TrialSwapTopo(const MCTad*,const MCTad*, double*);
	void AcceptMoveTopo(MCTad*,MCTad*) const;
	
	bool legal;
	bool legalTopo1;
	bool legalTopo2;
	
	int vo;
	int vn;
	
	int vin;
	int vjn;


private:
	int dn1;
	int dn2;
	int dn3;
	
	int din1;
	int din2;
	int djn1;
	int djn2;

	MCLattice* lat;
};


#endif /* MCTadUpdater_hpp */
