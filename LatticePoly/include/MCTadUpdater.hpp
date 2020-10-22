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
	void AcceptMovePos(MCTad*) const;

	void TrialMoveLeftEnd(const MCTad*, double*);
	void TrialMoveRightEnd(const MCTad*, double*);
	void TrialMoveLinear(const MCTad*, double*);
	//void TrialMoveFork(const MCTad*, double*);

	
	bool legal;
	
	int vo;
	int vn;
	
private:
	int nv1;
	int nv2;
	
	MCLattice* lat;
};


#endif /* MCTadUpdater_hpp */
