//
//  MCCGLattice.cpp
//  LatticePoly
//
//  Created by mtortora on 29/01/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCubeSource.h>
#include <vtkXMLPolyDataWriter.h>

#include "MCCGLattice.hpp"


void MCCGLattice::Init(std::mt19937_64& rngEngine)
{
	MCLattice::Init(rngEngine);
	
	nLiq = 0;
	legal = true;

	for ( int i = 0; i < Ntot; i++ )
	{
		spinTable[i] = 0;
		spinIdTable[i] = -1;
	}
	
	if ( InitDrop )
	{
		int centers[3][Ndrop];
		
		int r = std::floor(R)+1; // Set to 1 to allow initial droplets to cross PBCs
		
		for ( int i = 0; i < Ndrop; i++ )
		{
			centers[0][i] = (rngEngine() % (L-2*r+1)) + r;
			centers[1][i] = (rngEngine() % (L-2*r+1)) + r;
			centers[2][i] = (rngEngine() % (L-2*r+1)) + r;
		}
		
		for ( int i = 0; i < Ntot; i++ )
		{
			for ( int j = 0; j < Ndrop; j++ )
			{
				double dx = xyzTable[0][i] - centers[0][j];
				double dy = xyzTable[1][i] - centers[1][j];
				double dz = xyzTable[2][i] - centers[2][j];

				if      ( SQR(dx-L) + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx-L) + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx-L) + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[i] = Qcg;

				else if ( SQR(dx-L) + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx-L) + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx-L) + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[i] = Qcg;
				
				else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[i] = Qcg;
				
				else if ( SQR(dx)   + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx)   + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx)   + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[i] = Qcg;

				else if ( SQR(dx)   + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx)   + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx)   + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[i] = Qcg;
				
				else if ( SQR(dx)   + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx)   + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx)   + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[i] = Qcg;
				
				else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[i] = Qcg;

				else if ( SQR(dx+L) + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx+L) + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx+L) + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[i] = Qcg;
				
				else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[i] = Qcg;
				else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[i] = Qcg;
				
				if ( spinTable[i] > 0 )
				{
					nLiq++;
					break;
				}
			}
		}
	}
	
	else
	{
		int nLiqTot = std::floor(Ntot*Ldens);
		
		while ( nLiq < nLiqTot )
		{
			int idx = rngEngine() % Ntot;
			
			if ( spinTable[idx] < Qcg )
			{
				spinTable[idx]++;
				nLiq++;
			}
		}
	}

	for ( int i = 0; i < nLiq; i++ )
	{
		disp initDisp;
		
		initDisp.dx = 0.;
		initDisp.dy = 0.;
		initDisp.dz = 0.;

		spinDisp.push_back(initDisp);
	}
	
	std::cout << "Set up lattice with fixed liquid density " << nLiq / ((double) Ntot) << std::endl;
}

void MCCGLattice::TrialMove(std::mt19937_64& rngEngine, double* dE)
{
	int v = rngEngine() % 12;

	idx1 = rngEngine() % Ntot;
	idx2 = bitTable[v+1][idx1];
	
	legal = (spinTable[idx1] > 0) & (spinTable[idx2] < Qcg);
	*dE = GetSpinEnergy();
}

void MCCGLattice::AcceptMove()
{
	spinTable[idx1] -= 1;
	spinTable[idx2] += 1;
}

double MCCGLattice::GetSpinEnergy() const
{
	double dE = (spinTable[idx2] - spinTable[idx1] + 1.) / (double) Qcg;

	if ( legal )
	{
		for ( int i = 0; i < 12; i++ )
		{
			if ( bitTable[i+1][idx1] != idx2 )
				dE += spinTable[bitTable[i+1][idx1]];
		
			if ( bitTable[i+1][idx2] != idx1 )
				dE -= spinTable[bitTable[i+1][idx2]];
		}
	}
	
	dE *= Jll;
	
	return dE;
}

double MCCGLattice::GetSpinDensity(int idx) const
{
	double aveLiq = 0.;

	for ( int j = 0; j < 12; j++ )
		aveLiq += spinTable[idx] * spinTable[bitTable[j+1][idx]] / (12. * SQR(Qcg));
		
	return aveLiq;
}

double MCCGLattice::GetCouplingEnergy(const int tadHetTable[Ntot]) const
{
	double E1 = 0.;
	double E2 = 0.;
	
	for ( int i = 0; i < 13; i++ )
	{
		int v1 = (i == 0) ? idx1 : bitTable[i][idx1];
		int v2 = (i == 0) ? idx2 : bitTable[i][idx2];
		
		E1 -= tadHetTable[v1];
		E2 -= tadHetTable[v2];
	}
	
	return Jlp * (E2-E1);
}

void MCCGLattice::ToVTK(int idx)
{
	char buf[128];
	sprintf(buf, "%05d", idx);
	
	std::string filename = outputDir + "/liq" + buf + ".vtp";
	
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	
	vtkSmartPointer<vtkFloatArray> liqDensity = vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkFloatArray> liqDisplacement = vtkSmartPointer<vtkFloatArray>::New();
	
	liqDensity->SetName("Density");
	liqDensity->SetNumberOfComponents(1);
	
	liqDisplacement->SetName("Displacement");
	liqDisplacement->SetNumberOfComponents(3);
	
	int id = -1;
	
	for ( int i = 0; i < Ntot; i++ )
	{
		if ( spinTable[i] > 0 )
		{
			double aveDensity = GetSpinDensity(i);

			double x = xyzTable[0][i];
			double y = xyzTable[1][i];
			double z = xyzTable[2][i];
			
			id = (spinIdTable[i] == -1) ? id + 1 : spinIdTable[i];

			double dx = (spinIdTable[i] == -1) ? 0. : spinDisp[id].dx;
			double dy = (spinIdTable[i] == -1) ? 0. : spinDisp[id].dy;;
			double dz = (spinIdTable[i] == -1) ? 0. : spinDisp[id].dz;;
					
			points->InsertPoint(id, x, y, z);
		
			liqDisplacement->InsertTuple3(id, dx, dy, dz);
			liqDensity->InsertValue(id, aveDensity);
		}
	}
	
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	polydata->SetPoints(points);
	
	polydata->GetPointData()->AddArray(liqDensity);
	polydata->GetPointData()->AddArray(liqDisplacement);

	writer->SetFileName(filename.c_str());
	writer->SetInputData(polydata);
	
	writer->Write();
}
