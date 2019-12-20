//
//  globals.hpp
//  LatticePoly
//
//  Created by mtortora on 29/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_


// Chain parameters
#define Nchain 9827
#define Kint 1.2

// Box/MC parameters
#define L 27
#define Nmeas 200
#define Ninter 1000

// Number of liquid moves per TAD MC move
#define NliqMC 1

// Set InitDrop to 1 for initial drops, 0 for uniform liquid distribution
#define InitDrop 1

#define Ndrop 4 // Droplet number (if InitDrop == 1)
#define R 5. // Size of drops (must have have L >= 2*(R+1))

// Uniform starting density (if InitDrop == 0)
#define Ldens 0.07

// Liquid bleach radius (measured from the box center)
#define Rbleach 2.

// Set to 1 for Arrhenius dynamics, 0 for Metropolis
#define Arrhenius 1

#define Trel 0 // Relaxation time to turn on attractive interactions
#define Tbleach 49999 // Spin bleach time

// Liquid-liquid, liquid-polymer and polymer-polymer interaction strengths
#define Jll 0.7
#define Jlp 0.5
#define Jpp 0.3

// Heterochromatin domain dimensions
#define Ndom 5 // Number of domains
#define Nloc 350 // Number of heterochromatic loci per domain (must be <= Nchain)


#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define __DATA_PATH TOSTRING(__DPATH__)


#define L2 SQR(L)
#define L3 CUB(L)
#define Ntot (4*L3)


#endif /* globals_hpp */
