//
//  MCLivingPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 27/04/2021.
//  Copyright © 2021 ENS Lyon. All rights reserved.
//

#ifndef MCLivingPoly_hpp
#define MCLivingPoly_hpp

#include "MCHeteroPoly_looped.hpp"


class MCLivingPoly: public MCHeteroPoly_looped
{
public:
	MCLivingPoly(MCLattice*);
	
	void Init(int,int,int[3]);
	void TrialMove(double*);
	
	void ToVTK(int,std::string);
	
private:
	void UpdateFromFile(int);
	
	std::vector<std::vector<int> > colorData;
};


#endif /* MCLivingPoly_hpp */
