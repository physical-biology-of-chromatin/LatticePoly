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
#include <vtkXMLPolyDataWriter.h>

#include "MCLattice.hpp"


MCLattice::MCLattice() {InitConstArrays();}

void MCLattice::InitConstArrays()
{
	ReadInputArrays();

	double sum = 0.;
	
	for ( int i = 1; i < 13; i++ )
	{
		for ( int j = 1; j < 13; j++ )
		{
			sum += exp(-cTheta[i][j]);
		}
	}
	
	sum = -log(sum / ((double) 12*12));
	
	for ( int i = 0; i < 13; i++ )
	{
		cTheta[0][i] = sum;
		cTheta[i][0] = sum;
	}
	
	opp[0] = 0;
	
	for ( int i = 0; i < 6; i++ )
	{
		opp[2*i+1]   = 2*(i+1);
		opp[2*(i+1)] = 2*i+1;
	}
}

void MCLattice::ReadInputArrays()
{
	std::string dataDir = __DATA_PATH;

	std::string cosPath(dataDir + "/costhet.out");
	std::string xyzPath(dataDir + "/voisxyz.out");
	std::string nnPath (dataDir + "/voisnn.out");
	
	std::ifstream cosFile(cosPath);
	std::ifstream xyzFile(xyzPath);
	std::ifstream nnFile(nnPath);
	
	if ( !cosFile.good() ) throw std::runtime_error("MCLattice: Couldn't open file " + cosPath);
	if ( !xyzFile.good() ) throw std::runtime_error("MCLattice: Couldn't open file " + xyzPath);
	if ( !nnFile.good()  ) throw std::runtime_error("MCLattice: Couldn't open file " + nnPath);
	
	for ( int i = 0; i < 13; i++ )
	{
		for ( int j = 0; j < 13; j++ )
		{
			for ( int k = 0; k < 13; k++ )
			{
				nnFile >> nbNN[k][i][j];
				nbNN[k][i][j] -= 1;
			}
			
			cosFile >> cTheta[i][j];
			cTheta[i][j] *= Kint;
		}
		
		for ( int j = 0; j < 3; j++ )
		{
			xyzFile >> nbXYZ[j][i];
		}
	}
	
	cosFile.close();
	xyzFile.close();
	nnFile.close();
}

void MCLattice::Init(std::mt19937_64&)
{
	for ( int i = 0; i < Ntot; i++ )
	{
		int iz = i/(2*L2);
		int iy = i/L - 2*L*iz;
		int ix = i - L*(2*L*iz + iy);
		
		int mod = iz % 2;
		
		double x = ix + 0.5*(1-((iy+1+mod) % 2));
		double y = iy*0.5;
		double z = iz*0.5;
		
		xyzTable[0][i] = x;
		xyzTable[1][i] = y;
		xyzTable[2][i] = z;

		bitTable[0][i] = 0;
		
		for ( int j = 0; j < 12; j++ )
		{
			double xp = x + nbXYZ[0][j+1];
			double yp = y + nbXYZ[1][j+1];
			double zp = z + nbXYZ[2][j+1];
			
			if ( xp >= L ) xp = xp - L;
			if ( xp < 0 )  xp = xp + L;
			
			if ( yp >= L ) yp = yp - L;
			if ( yp < 0 )  yp = yp + L;
			
			if ( zp >= L ) zp = zp - L;
			if ( zp < 0 )  zp = zp + L;
			
			int ixp = (int) xp;
			int iyp = (int) 2*yp;
			int izp = (int) 4*zp;
			
			bitTable[j+1][i] = ixp + iyp*L + izp*L2;
		}
	}
	
	ToVTK(0);
}

void MCLattice::ToVTK(int)
{
	std::string filename = outputDir + "/box.vtp";
	
	vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();
	
	cubeSource->SetCenter((L-0.5)/2., (L-0.5)/2., (L-0.5)/2.);
	
	cubeSource->SetXLength(L+0.5);
	cubeSource->SetYLength(L+0.5);
	cubeSource->SetZLength(L+0.5);
	
	cubeSource->Update();
	
	// Write file
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	
	writer->SetFileName(filename.c_str());
	writer->SetInputConnection(cubeSource->GetOutputPort());
	
	writer->Write();
}
