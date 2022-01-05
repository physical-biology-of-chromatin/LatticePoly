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
	void TrialReptationMove(const MCTad*, int) ;


	void TrialMoveLeftEnd(const MCTad*, double*);
	void TrialMoveRightEnd(const MCTad*, double*);
	void TrialMoveLinear(const MCTad*, double*);
	void TrialMoveFork(const MCTad*, double*);
	void SaveSpecialMonomers(const MCTad* tad);


	bool legal;
	int reptation_step;
	int rept_dir;
	
	std::vector<std::vector<int>> reptation_values;
	std::vector<std::vector<int>>  reptating_choesins;

private:
	
	MCLattice* lat;
};


#endif /* MCTadUpdater_hpp */
