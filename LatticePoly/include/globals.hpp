//
//  globals.hpp
//  LatticePoly
//
//  Created by mtortora on 29/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//
//

#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_

#include <array>
#include <string>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCubeSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>


// Box linear dimension
#define L 40

// Runtime global parameters
extern std::string outputDir;

extern std::string domainPath;
extern std::string colorPath;

extern std::string latticeType;
extern std::string polyType;

extern int Nrelax;
extern int Nmeas;
extern int Ninter;
extern int NG1;
extern int Nchain;

extern int NliqMC;
extern int Ndrop;

extern int Qcg;

extern bool InitDrop;
extern bool RestartFromFile;
extern bool RablConf;

extern bool neigh;


extern double Kint;
extern double krouse;

extern double R;
extern double Ldens;

extern double Jll;
extern double Jlp;
extern double Jpp;
extern double Jpair;
extern double stop_replication;


extern double Jf;
extern double Jf_sister;
extern int enhancement_sister;
extern int enhancement_fork;

extern int Centromere;




extern double inactiveRatio;
extern double propRate;

extern int propagationMode;

extern double originRate;
extern double replicRate;
extern int Ndf;
extern int Ncohesins;


// Custom macros, compile-time constants & typedefs
#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))

#define L2 SQR(L)
#define L3 CUB(L)
#define Ntot (4*L3)

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define __DATA_PATH TOSTRING(__DPATH__)

typedef std::array<double, 3> double3;


#endif /* globals_hpp */
