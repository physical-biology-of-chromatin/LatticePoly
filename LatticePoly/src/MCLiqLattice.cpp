//
//  MCLiqLattice.cpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

// #include <vtkLine.h>
// #include <vtkPointData.h>
// #include <vtkFloatArray.h>
// #include <vtkCubeSource.h>
// #include <vtkXMLPolyDataReader.h>
// #include <vtkXMLPolyDataWriter.h>

#include "MCLiqLattice.hpp"


void MCLiqLattice::Init(int Ninit)
{
	MCLattice::Init(Ninit);

	nLiq = 0;
	
	for ( int vi = 0; vi < Ntot; ++vi )
	{
		spinTable[vi] = 0;
		spinNeighborhood[vi] = 0;
		lookupTable[vi] = -1;
	}
	
	// if ( RestartFromFile ) // FromVTK(Ninit);

	// 	FromHDF5(Ninit);

	// else
	// {
	if ( InitDrop )
		GenerateDroplets();
	else
		GenerateRandom();
	
	spinConf.resize(nLiq);
	spinDisp.resize(nLiq);

	std::fill(spinDisp.begin(), spinDisp.end(), (double3) {0., 0., 0.});

	int ctr = 0;
	
	for ( int vi = 0; vi < Ntot; ++vi )
	{
		if ( spinTable[vi] > 0 )
		{
			lookupTable[vi] = ctr;
			spinConf[ctr] = vi;
			
			++ctr;
		}
	}
	// }
	
	std::cout << "Set up lattice with fixed liquid density " << nLiq / ((double) Ntot) << std::endl;
}

void MCLiqLattice::GenerateDroplets()
{
	std::vector<double3> centers(Ndrop);
	
	int r = std::floor(R) + 1; // Set to 1 to allow initial droplets to cross PBCs
	
	for ( int i = 0; i < Ndrop; ++i )
	{
		centers[i][0] = 12; // (rngEngine() % (L-2*r+1)) + r;
		centers[i][1] = 12; // (rngEngine() % (L-2*r+1)) + r;
		centers[i][2] = 12; // (rngEngine() % (L-2*r+1)) + r;
	}
	
	for ( int vi = 0; vi < Ntot; ++vi )
	{
		for ( int i = 0; i < Ndrop; ++i )
		{
			double dx = xyzTable[0][vi] - centers[i][0];
			double dy = xyzTable[1][vi] - centers[i][1];
			double dz = xyzTable[2][vi] - centers[i][2];

			if      ( SQR(dx-L) + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;

			else if ( SQR(dx-L) + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			else if ( SQR(dx)   + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;

			else if ( SQR(dx)   + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			else if ( SQR(dx)   + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;

			else if ( SQR(dx+L) + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			if ( spinTable[vi] > 0 )
			{

				for (int v = 0; v<13; ++v)
				{
					
				int pos = (v == 0) ? vi : bitTable[v][vi];

				spinNeighborhood[pos] += 1;
			
				}
				++nLiq;
				break;
			}
		}
	}
}

void MCLiqLattice::GenerateRandom()
{
	while ( nLiq < std::floor(Ntot*Ldens) )
	{
		int vi = rngEngine() % Ntot;
		
		if ( spinTable[vi] == 0 )
		{
			++spinTable[vi];
			++nLiq;

			for (int v = 0; v<13; ++v)
			{
				
			int pos = (v == 0) ? vi : bitTable[v][vi];

			spinNeighborhood[pos] += 1;
		
			}
		}
	}
}

void MCLiqLattice::TrialMove(double* dE)
{
	int n = rngEngine() % nLiq;
	int v = rngEngine() % 12;

	v1 = spinConf[n];
	v2 = bitTable[v+1][v1];
	
	*dE = GetSpinEnergy();	
}

void MCLiqLattice::AcceptMove()
{
	DisplaceSpins();
					
	if ( spinTable[v2] == 0 )
	{
		int id1 = lookupTable[v1];
	
		lookupTable[v1] = -1;
		lookupTable[v2] = id1;
		
		spinConf[id1] = v2;

		for (int v = 0; v < 13; ++v)
		{
			int vi1 = (v == 0) ? v1 : bitTable[v][v1];
			int vi2 = (v == 0) ? v2 : bitTable[v][v2];
			
			spinNeighborhood[vi1] -= 1;
			spinNeighborhood[vi2] += 1;

		}
	}
	
	else
	{
		int id1 = lookupTable[v1];
		int id2 = lookupTable[v2];
		
		lookupTable[v1] = id2;
		lookupTable[v2] = id1;
		
		spinConf[id1] = v2;
		spinConf[id2] = v1;
	}
	
	spinTable[v1] = spinTable[v2];
	spinTable[v2] = 1;
}

void MCLiqLattice::DisplaceSpins()
{
	for ( int i = 0; i < 3; ++i )
	{
		double disp = xyzTable[i][v2] - xyzTable[i][v1];

		if ( std::abs(disp) > L/2. )
			disp -= std::copysign(L, disp);

		int id1 = lookupTable[v1];
		spinDisp[id1][i] += disp;
		
		if ( spinTable[v2] == 1 )
		{
			int id2 = lookupTable[v2];
			spinDisp[id2][i] -= disp;
		}
	}
}

double MCLiqLattice::GetSpinEnergy() const
{
	if ( Jll > 0. )
	{
		if ( spinTable[v2] == 0 )
		{
			double dE = 0.;
			
			for ( int v = 1; v < 13; ++v )
			{
				if ( bitTable[v][v1] == v2)
				{
					for ( int i = 0; i<7; ++i)
					{
						if (spinTable[bitTable[enNN[v][i]][v1]] == 1.)
						{							
							dE += ((Jll_Valency < spinNeighborhood[bitTable[enNN[v][i]][v1]]-1) ? 0. : 1.);
						}
					}
				}
				if (bitTable[v][v2] == v1)
				{
					for ( int i = 0; i<7; ++i)
					{
						if (spinTable[bitTable[enNN[v][i]][v2]] == 1.)
						{
							dE -= ((Jll_Valency < spinNeighborhood[bitTable[enNN[v][i]][v2]]) ? 0. : 1.);
						}
					}
				}
			}
		return Jll / 2 * (((Jll_Valency < spinNeighborhood[v1]-1) ? Jll_Valency : (spinNeighborhood[v1]-1)) - ((Jll_Valency < spinNeighborhood[v2]-1) ? Jll_Valency : (spinNeighborhood[v2]-1)) + dE );
		}
	}
	
	return 0.;
}

double MCLiqLattice::GetCouplingEnergy(const double hetTable[Ntot], const double hetNeighborhood[Ntot]) const
{
	if ( Jlp > 0. )
	{
		if ( spinTable[v2] == 0 )
		{
			double dN = 0.;
			for ( int v = 1; v < 13; ++v )
			{
				if ( bitTable[v][v1] == v2)
				{
					for ( int i = 0; i<7; ++i)
					{
						if (hetTable[bitTable[enNN[v][i]][v1]] != 0.)
							dN += ((Jpl_Valency < spinNeighborhood[bitTable[enNN[v][i]][v1]]) ? 0. : hetTable[bitTable[enNN[v][i]][v1]]);
					}
				}
				if (bitTable[v][v2] == v1)
				{
					for ( int i = 0; i<7; ++i)
					{
						if (hetTable[bitTable[enNN[v][i]][v2]] != 0.)
							dN -= ((Jpl_Valency < spinNeighborhood[bitTable[enNN[v][i]][v2]]+1) ? 0. : hetTable[bitTable[enNN[v][i]][v2]]);
					}
				}
			}

			return Jlp / 2 * (((Jlp_Valency < hetNeighborhood[v1]) ? Jlp_Valency : hetNeighborhood[v1]) - ((Jlp_Valency < hetNeighborhood[v2]) ? Jlp_Valency : hetNeighborhood[v2]) + dN) + EV * (bitTable[0][v2]-bitTable[0][v1]);
		}
	}
	
	return 0.;
}

double MCLiqLattice::GetCouplingEnergyPainter(const double hetTable[Ntot], const double hetNeighborhood[Ntot],const double painterTable[Ntot], const double painterNeighborhood[Ntot]) const
{

	double dE = MCLiqLattice::GetCouplingEnergy(hetTable, hetNeighborhood);
	
	if ( ( Jlpp > 0. ) )
	{
		if ( spinTable[v2] == 0 )
		{
			double dN = 0.;
			for ( int v = 1; v < 13; ++v )
			{
				if ( bitTable[v][v1] == v2)
				{
					for ( int i = 0; i<7; ++i)
					{
						if (painterTable[bitTable[enNN[v][i]][v1]] != 0.)
							dN += ((Jppl_Valency < spinNeighborhood[bitTable[enNN[v][i]][v1]]) ? 0. : painterTable[bitTable[enNN[v][i]][v1]]);
					}
				}
				if (bitTable[v][v2] == v1)
				{
					for ( int i = 0; i<7; ++i)
					{
						if (painterTable[bitTable[enNN[v][i]][v2]] != 0.)
							dN -= ((Jppl_Valency < spinNeighborhood[bitTable[enNN[v][i]][v2]]+1) ? 0. : painterTable[bitTable[enNN[v][i]][v2]]);
					}
				}
			}	
			dE += Jlpp / 2 * (((Jlpp_Valency < painterNeighborhood[v1]) ? Jlpp_Valency : painterNeighborhood[v1]) - ((Jlpp_Valency < painterNeighborhood[v2]) ? Jlpp_Valency : painterNeighborhood[v2]) + dN);
		}
	}	
    
	
	return dE;
}

void MCLiqLattice::ToHDF5(int frame)
{	
	std::clock_t c_start = std::clock();
	
	using namespace H5;
	
	double data_position[nLiq][3];
	double data_density[nLiq];
	double data_displacement[nLiq][3];

	for ( int i = 0; i < nLiq; ++i )
	{
		int vi = spinConf[i];
		double aveDensity = 0.;

		for ( int v = 0; v < 12; ++v )
			aveDensity += spinTable[bitTable[v+1][vi]] / 12.;
				
		data_position[i][0] = xyzTable[0][vi];
		data_position[i][1] = xyzTable[1][vi];
		data_position[i][2] = xyzTable[2][vi];

		data_displacement[i][0] = spinDisp[i][0];
		data_displacement[i][1] = spinDisp[i][1];
		data_displacement[i][2] = spinDisp[i][2];
		
		data_density[i] = aveDensity;
	}

	const H5std_string FILE_PATH( H5filePath );
	const int RANK = 2;

	// char Frame_group_Name[32];
	// sprintf(Frame_group_Name, "/Liq/Frame%05d", frame);
	// const H5std_string FRAME_GROUP_NAME( Frame_group_Name );

	H5File file(FILE_PATH, H5F_ACC_RDWR);
	Group liq_group(file.openGroup("/Liq"));
	// Group *frame_group = new Group(file.openGroup(FRAME_GROUP_NAME));

	char Position_dataset_Name[32];
	sprintf(Position_dataset_Name, "%05d_position_dataset", frame);
	const H5std_string POSITION_DATASET_NAME(Position_dataset_Name);

	hsize_t dims_position[RANK];
	dims_position[0] = nLiq;
	dims_position[1] = 3;
	
	DataSpace *dataspace_position = new DataSpace( RANK, dims_position);
	
	DataSet *dataset_position = new DataSet(liq_group.createDataSet( POSITION_DATASET_NAME, PredType::NATIVE_DOUBLE, *dataspace_position ));
	
	dataset_position->write(data_position , PredType::NATIVE_DOUBLE);

	delete dataset_position;
	delete dataspace_position;

	char Density_dataset_Name[32];
	sprintf(Density_dataset_Name, "%05d_density_dataset", frame);
	const H5std_string DENSITY_DATASET_NAME(Density_dataset_Name);

	hsize_t dims_density[1];
	dims_density[0] = nLiq;
	
	DataSpace *dataspace_density = new DataSpace( 1, dims_density);
	
	DataSet *dataset_density = new DataSet(liq_group.createDataSet( DENSITY_DATASET_NAME, PredType::NATIVE_DOUBLE, *dataspace_density ));
	
	dataset_density->write(data_density , PredType::NATIVE_DOUBLE);

	delete dataset_density;
	delete dataspace_density;

	
	char Displacement_dataset_Name[32];
	sprintf(Displacement_dataset_Name, "%05d_displacement_dataset", frame);
	const H5std_string DISPLACEMENT_DATASET_NAME(Displacement_dataset_Name);

	hsize_t dims_displacement[RANK];
	dims_displacement[0] = nLiq;
	dims_displacement[1] = 3;
	
	DataSpace *dataspace_displacement = new DataSpace( RANK, dims_displacement);
	
	DataSet *dataset_displacement = new DataSet(liq_group.createDataSet( DISPLACEMENT_DATASET_NAME, PredType::NATIVE_DOUBLE, *dataspace_displacement ));
	
	dataset_displacement->write(data_displacement , PredType::NATIVE_DOUBLE);

	delete dataset_displacement;
	delete dataspace_displacement;
	liq_group.close();
	// frame_group->close();
	file.close();

	std::clock_t c_end = std::clock();

	double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	std::cout << "CPU time used for MCLiqLattice::ToHDF5: " 
          	<< time_elapsed_ms 
          	<< " ms\n";
}

// void MCLiqLattice::ToVTK(int frame)
// {
// 	std::clock_t c_start = std::clock();

// 	char fileName[32];
// 	sprintf(fileName, "liq%05d.vtp", frame);
	
// 	std::string path = outputDir + "/" + fileName;
	
// 	auto points = vtkSmartPointer<vtkPoints>::New();
// 	auto liqDensity = vtkSmartPointer<vtkFloatArray>::New();
// 	auto liqDisplacement = vtkSmartPointer<vtkFloatArray>::New();
	
// 	liqDensity->SetName("Density");
// 	liqDensity->SetNumberOfComponents(1);
	
// 	liqDisplacement->SetName("Displacement");
// 	liqDisplacement->SetNumberOfComponents(3);
		
// 	for ( int i = 0; i < nLiq; ++i )
// 	{
// 		int vi = spinConf[i];
// 		double aveDensity = 0.;

// 		for ( int v = 0; v < 12; ++v )
// 			aveDensity += spinTable[bitTable[v+1][vi]] / 12.;
		
// 		double x = xyzTable[0][vi];
// 		double y = xyzTable[1][vi];
// 		double z = xyzTable[2][vi];
		
// 		double dx = spinDisp[i][0];
// 		double dy = spinDisp[i][1];
// 		double dz = spinDisp[i][2];
				
// 		points->InsertNextPoint(x, y, z);
	
// 		liqDensity->InsertNextValue(aveDensity);
// 		liqDisplacement->InsertNextTuple3(dx, dy, dz);
// 	}
	
// 	auto polyData = vtkSmartPointer<vtkPolyData>::New();
// 	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

// 	polyData->SetPoints(points);
	
// 	polyData->GetPointData()->AddArray(liqDensity);
// 	polyData->GetPointData()->AddArray(liqDisplacement);

// 	writer->SetFileName(path.c_str());
// 	writer->SetInputData(polyData);
	
// 	writer->Write();
	
// 	std::clock_t c_end = std::clock();

// 	double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
// 	std::cout << "CPU time used for MCLiqLattice::ToVTK: " 
//           	<< time_elapsed_ms 
//           	<< " ms\n";
// }

// void MCLiqLattice::FromVTK(int frame)
// {
// 	char fileName[32];
// 	sprintf(fileName, "liq%05d.vtp", frame);
	
// 	std::string path = outputDir + "/" + fileName;
	
// 	auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

// 	reader->SetFileName(path.c_str());
// 	reader->Update();
	
// 	vtkPolyData* polyData = reader->GetOutput();
// 	vtkDataArray* dispData = polyData->GetPointData()->GetArray("Displacement");

// 	nLiq = (int) polyData->GetNumberOfPoints();
	
// 	spinConf.reserve(nLiq);
// 	spinDisp.reserve(nLiq);

// 	if ( (InitDrop == 0) && (nLiq != std::floor(Ntot*Ldens)) )
// 		throw std::runtime_error("MCLiqLattice: Found liquid configuration file with incompatible dimension " + std::to_string(nLiq));
	
// 	std::cout << "Starting from liquid configuration file " << path << std::endl;
	
// 	for ( int i = 0; i < nLiq; ++i )
// 	{
// 		double3 initDisp;
// 		double point[3];
		
// 		polyData->GetPoint(i, point);
		
// 		for ( int j = 0; j < 3; ++j )
// 			initDisp[j] = dispData->GetComponent(i, j);

// 		int ixp = (int) 1*point[0];
// 		int iyp = (int) 2*point[1];
// 		int izp = (int) 4*point[2];
		
// 		int vi = ixp + iyp*L + izp*L2;
		
// 		lookupTable[vi] = i;
		
// 		spinConf.push_back(vi);
// 		spinDisp.push_back(initDisp);
		
// 		++spinTable[vi];
// 	}
// }
