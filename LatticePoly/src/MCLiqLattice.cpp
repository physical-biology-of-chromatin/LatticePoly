//
//  MCLiqLattice.cpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCubeSource.h>
#include <vtkXMLPolyDataWriter.h>

#include "MCLiqLattice.hpp"


void MCLiqLattice::Init(std::mt19937_64& rngEngine)
{
	MCLattice::Init(rngEngine);
	
	for ( int i = 0; i < Ntot; i++ )
		spinTable[i] = 0;
	
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

				if      ( SQR(dx-L) + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx-L) + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx-L) + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;

				else if ( SQR(dx-L) + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx-L) + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx-L) + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
				
				else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
				
				else if ( SQR(dx)   + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx)   + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx)   + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;

				else if ( SQR(dx)   + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx)   + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx)   + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
				
				else if ( SQR(dx)   + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx)   + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx)   + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
				
				else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;

				else if ( SQR(dx+L) + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx+L) + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx+L) + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
				
				else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
				else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
				
				if ( spinTable[i] == 1 )
				{
					spinConf.push_back(i);
					
					break;
				}
			}
		}
		
		nLiq = (int) spinConf.size();
	}
	
	else
	{
		nLiq = std::floor(Ntot*Ldens);
		
		while ( spinConf.size() < (size_t) nLiq )
		{
			int idx = rngEngine() % Ntot;
			
			if ( spinTable[idx] == 0 )
			{
				spinTable[idx] = 1;
				spinConf.push_back(idx);
			}
		}
		
		std::sort(spinConf.begin(), spinConf.end());
	}
	
	int idx = 0;
	
	for ( int i = 0; i < Ntot; i++ )
	{
		spinIdTable[i] = (spinTable[i] == 1) ? idx++ : -1;
		spinTypeTable[i] = 0;
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

void MCLiqLattice::BleachSpins()
{
	for ( int i = 0; i < Ntot; i++ )
	{
		double dx = xyzTable[0][i] - (L-0.5)/2.;
		double dy = xyzTable[1][i] - (L-0.5)/2.;
		
		if ( spinTable[i] == 1 )
			spinTypeTable[i] = (SQR(dx)+SQR(dy) < SQR(Rbleach));
	}
}

void MCLiqLattice::TrialMoveSpin(std::mt19937_64& rngEngine, double* dE)
{
	int v1 = rngEngine() % nLiq;
	int v2 = rngEngine() % 12;

	idxSpin1 = spinConf[v1];
	idxSpin2 = bitTable[v2+1][idxSpin1];
	
	*dE = GetSpinEnergy();
}

void MCLiqLattice::AcceptMoveSpin()
{
	int spin1 = spinTable[idxSpin1];
	int spin2 = spinTable[idxSpin2];
	
	if ( (spin1 == 1) || (spin2 == 1) )
	{
		DisplaceSpins();

		int type1 = spinTypeTable[idxSpin1];
		int type2 = spinTypeTable[idxSpin2];

		spinTable[idxSpin1] = spin2;
		spinTable[idxSpin2] = spin1;
		
		spinTypeTable[idxSpin1] = type2;
		spinTypeTable[idxSpin2] = type1;
				
		if ( (spin1 == 1) && (spin2 == 0) )
		{
			int id1 = spinIdTable[idxSpin1];
		
			spinIdTable[idxSpin1] = -1;
			spinIdTable[idxSpin2] = id1;
			
			spinConf[id1] = idxSpin2;
		}
		
		else if ( (spin1 == 0) && (spin2 == 1) )
		{
			int id2 = spinIdTable[idxSpin2];
			
			spinIdTable[idxSpin1] = id2;
			spinIdTable[idxSpin2] = -1;
				
			spinConf[id2] = idxSpin1;
		}
		
		else
		{
			int id1 = spinIdTable[idxSpin1];
			int id2 = spinIdTable[idxSpin2];
			
			spinIdTable[idxSpin1] = id2;
			spinIdTable[idxSpin2] = id1;
			
			spinConf[id1] = idxSpin2;
			spinConf[id2] = idxSpin1;
		}
	}
}

void MCLiqLattice::DisplaceSpins()
{
	int spin1 = spinTable[idxSpin1];
	int spin2 = spinTable[idxSpin2];
	
	double dx = xyzTable[0][idxSpin2] - xyzTable[0][idxSpin1];
	double dy = xyzTable[1][idxSpin2] - xyzTable[1][idxSpin1];
	double dz = xyzTable[2][idxSpin2] - xyzTable[2][idxSpin1];
	
	if ( std::abs(dx) > L/2. ) dx -= std::copysign(L, dx);
	if ( std::abs(dy) > L/2. ) dy -= std::copysign(L, dy);
	if ( std::abs(dz) > L/2. ) dz -= std::copysign(L, dz);
	
	if ( spin1 == 1 )
	{
		int id1 = spinIdTable[idxSpin1];

		spinDisp[id1].dx += dx;
		spinDisp[id1].dy += dy;
		spinDisp[id1].dz += dz;
	}
	
	if ( spin2 == 1 )
	{
		int id2 = spinIdTable[idxSpin2];

		spinDisp[id2].dx -= dx;
		spinDisp[id2].dy -= dy;
		spinDisp[id2].dz -= dz;
	}
}

double MCLiqLattice::GetSpinEnergy() const
{
	double dE = 0.;
	
	if ( spinTable[idxSpin1] != spinTable[idxSpin2] )
	{
		for ( int i = 0; i < 12; i++ )
		{
			if ( bitTable[i+1][idxSpin1] != idxSpin2 )
				dE += spinTable[bitTable[i+1][idxSpin1]];
		
			if ( bitTable[i+1][idxSpin2] != idxSpin1 )
				dE -= spinTable[bitTable[i+1][idxSpin2]];
		}
	
		dE *= Jll * (spinTable[idxSpin1]-spinTable[idxSpin2]);
	}
	
	return dE;
}

double MCLiqLattice::GetBindingEnergy(const int tadTable[Ntot]) const
{
	double E1 = 0.;
	
	for ( int i = 0; i < 13; i++ )
	{
		int v1 = (i == 0) ? idxSpin1 : bitTable[i][idxSpin1];
		E1 -= tadTable[v1];
		
		if ( spinTable[idxSpin2] == 1 )
		{
			int v2 = (i == 0) ? idxSpin2 : bitTable[i][idxSpin2];
			E1 -= tadTable[v2];
		}
	}
	
	return Jlp * E1;
}

double MCLiqLattice::GetCouplingEnergy(const int tadTable[Ntot]) const
{
	double E1 = 0.;
	double E2 = 0.;
	
	if ( spinTable[idxSpin2] == 0 )
	{
		for ( int i = 0; i < 13; i++ )
		{
			int v1 = (i == 0) ? idxSpin1 : bitTable[i][idxSpin1];
			int v2 = (i == 0) ? idxSpin2 : bitTable[i][idxSpin2];
		
			E1 -= tadTable[v1];
			E2 -= tadTable[v2];
		}
	}
	
	return Jlp * (E2-E1);
}

void MCLiqLattice::ToVTK(int idx)
{
	char buf[256];
	sprintf(buf, "%04d", idx);
	
	std::string filename = outputDir + "/liq" + buf + ".vtp";
	
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	
	vtkSmartPointer<vtkIntArray> frapTypes = vtkSmartPointer<vtkIntArray>::New();
	vtkSmartPointer<vtkFloatArray> liqDensity = vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkFloatArray> liqDisplacement = vtkSmartPointer<vtkFloatArray>::New();

	frapTypes->SetName("FRAP type");
	frapTypes->SetNumberOfComponents(1);
	
	liqDensity->SetName("Density");
	liqDensity->SetNumberOfComponents(1);
	
	liqDisplacement->SetName("Displacement");
	liqDisplacement->SetNumberOfComponents(3);
	
	for ( int i = 0; i < nLiq; i++ )
	{
		int v = spinConf[i];
		
		double x = xyzTable[0][v];
		double y = xyzTable[1][v];
		double z = xyzTable[2][v];
		
		double dx = spinDisp[i].dx;
		double dy = spinDisp[i].dy;
		double dz = spinDisp[i].dz;

		double aveLiq = 0.;
		int frapType = spinTypeTable[v];

		for ( int j = 0; j < 12; j++ )
			aveLiq += spinTable[bitTable[j+1][v]] / 12.;
		
		points->InsertNextPoint(x, y, z);
		
		frapTypes->InsertNextValue(frapType);
		liqDensity->InsertNextValue(aveLiq);
		
		liqDisplacement->InsertNextTuple3(dx, dy, dz);
	}
	
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	polydata->SetPoints(points);
	
	polydata->GetPointData()->AddArray(frapTypes);
	polydata->GetPointData()->AddArray(liqDensity);
	polydata->GetPointData()->AddArray(liqDisplacement);

	writer->SetFileName(filename.c_str());
	writer->SetInputData(polydata);
	
	writer->Write();
}
