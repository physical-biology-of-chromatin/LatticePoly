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


	void TrialMoveLeftEnd(const MCTad*);
	void TrialMoveRightEnd(const MCTad*);
	void TrialMoveLinear(const MCTad*);
	void TrialMoveFork(const MCTad*);
	void CheckForkLegal(const MCTad*, int , int , int , int, int);



	
	bool legal;

	
	std::vector<std::vector<int>> reptation_values;
	std::vector<const MCTad*>  reptating_tads;


private:
	
	MCLattice* lat;
};


#endif /* MCTadUpdater_hpp */
