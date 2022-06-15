//
//  MCLattice.cpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <fstream>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCubeSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include "MCLattice.hpp"


MCLattice::MCLattice()
{
	ReadInputArrays();

	double sum = 0.;
	
	for ( int v1 = 1; v1 < 13; ++v1 )
	{
		for ( int v2 = 1; v2 < 13; ++v2 )
			sum += exp(-cTheta[v1][v2]);
	}
	
	sum = -log(sum / ((double) 12*12));
	
	for ( int v = 0; v < 13; ++v )
	{
		cTheta[0][v] = sum;
		cTheta[v][0] = sum;
	}
	
	opp[0] = 0;
	
	for ( int i = 0; i < 6; ++i )
	{
		opp[2*i+1]   = 2*(i+1);
		opp[2*(i+1)] = 2*i+1;
	}
}

void MCLattice::ReadInputArrays()
{
	std::string dataDir = __DATA_PATH;

	std::string cosPath(dataDir + "/costhet.in");
	std::string xyzPath(dataDir + "/nbxyz.in");
	std::string nnPath (dataDir + "/nbnn.in");
	
	std::ifstream cosFile(cosPath);
	std::ifstream xyzFile(xyzPath);
	std::ifstream nnFile(nnPath);
	
	if ( !cosFile.good() )
		throw std::runtime_error("MCLattice: Couldn't open file " + cosPath);
	if ( !xyzFile.good() )
		throw std::runtime_error("MCLattice: Couldn't open file " + xyzPath);
	if ( !nnFile.good() )
		throw std::runtime_error("MCLattice: Couldn't open file " + nnPath);
	
	for ( int v1 = 0; v1 < 13; ++v1 )
	{
		for ( int v2 = 0; v2 < 13; ++v2 )
		{
			for ( int v3 = 0; v3 < 13; ++v3 )
				nnFile >> nbNN[v3][v1][v2];
			
			cosFile >> cTheta[v1][v2];
			
			cTheta[v1][v2] *= Kint;
		}
		
		for ( int i = 0; i < 3; ++i )
			xyzFile >> nbXYZ[i][v1];
	}
	
	cosFile.close();
	xyzFile.close();
	nnFile.close();
}

void MCLattice::Init(int)
{
	for ( int vi = 0; vi < Ntot; ++vi )
	{
		int iz = vi/(2*L2);
		int iy = vi/L - 2*L*iz;
		int ix = vi - L*(2*L*iz + iy);
		
		int mod = iz % 2;
		
		double x = ix + 0.5*(1-((iy+1+mod) % 2));
		double y = iy*0.5;
		double z = iz*0.5;
		
		xyzTable[0][vi] = x;
		xyzTable[1][vi] = y;
		xyzTable[2][vi] = z;

		bitTable[0][vi] = 0;
		
		if (Rconfinement > 0)
		{
			double c = (L-0.5)/2;
			double d2 = SQR(xyzTable[0][vi]-c)+SQR(xyzTable[1][vi]-c)+SQR(xyzTable[2][vi]-c);
			if ( d2 > SQR(Rconfinement) ) 
				bitTable[0][vi] = -100;
		}		

		for ( int v = 0; v < 12; ++v )
		{
			double xp = x + nbXYZ[0][v+1];
			double yp = y + nbXYZ[1][v+1];
			double zp = z + nbXYZ[2][v+1];
			
			if ( xp >= L ) xp -= L;
			if ( xp < 0 )  xp += L;
			
			if ( yp >= L ) yp -= L;
			if ( yp < 0 )  yp += L;
			
			if ( zp >= L ) zp -= L;
			if ( zp < 0 )  zp += L;
			
			int ixp = (int) 1*xp;
			int iyp = (int) 2*yp;
			int izp = (int) 4*zp;
			
			bitTable[v+1][vi] = ixp + iyp*L + izp*L2;
		}
	}
	
	if ( RestartFromFile )
		BoxFromVTK();
	else
		BoxToVTK();
}

void MCLattice::BoxToVTK()
{
	std::string path = outputDir + "/box.vtp";
	
	auto cubeSource = vtkSmartPointer<vtkCubeSource>::New();
	
	cubeSource->SetCenter((L-0.5)/2., (L-0.5)/2., (L-0.5)/2.);
	
	cubeSource->SetXLength(L+0.5);
	cubeSource->SetYLength(L+0.5);
	cubeSource->SetZLength(L+0.5);
	
	cubeSource->Update();
	
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	
	writer->SetFileName(path.c_str());
	writer->SetInputConnection(cubeSource->GetOutputPort());
	
	writer->Write();
}

void MCLattice::BoxFromVTK()
{
	std::string path = outputDir + "/box.vtp";
	
	auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

	reader->SetFileName(path.c_str());
	reader->Update();
	
	vtkPolyData* polyData = reader->GetOutput();
	
	int Lx = (int) polyData->GetBounds()[1];
	int Ly = (int) polyData->GetBounds()[3];
	int Lz = (int) polyData->GetBounds()[5];
	
	if ( (Lx != L) || (Ly != L) || (Lz != L) )
		throw std::runtime_error("MCLattice: Found box file with incompatible dimensions");
}
