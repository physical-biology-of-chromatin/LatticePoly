//
//  MCPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <iterator>
#include <algorithm>

#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCubeSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include "MCPoly.hpp"


MCPoly::MCPoly(MCLattice* _lat): lat(_lat)
{
	tadUpdater = new MCTadUpdater(lat);
}

MCPoly::~MCPoly()
{
	delete tadUpdater;
}

void MCPoly::Init(int Ninit)
{
	tadConf.reserve(2*Ntot);
	tadTopo.reserve(2*Ntot);
	
	std::fill(centerMass.begin(), centerMass.end(), 0.);

	if ( RestartFromFile )
		FromVTK(Ninit);
	else
		GenerateRandom(L/2);

	for ( auto bond = tadTopo.begin(); bond != tadTopo.end(); ++bond )
		CreateBond(*bond);
	
	/*
	CAR.push_back(&tadConf.at(580));
	CAR.back()->choesin_binding_site = &tadConf.at(583);
	CAR.push_back(&tadConf.at(584));
	CAR.back()->choesin_binding_site = &tadConf.at(589);
	CAR.push_back(&tadConf.at(590));
	CAR.back()->choesin_binding_site = &tadConf.at(599);

*/
	//CAR.back()->choesin_binding_site = &tadConf.at(125);
	//GenerateCAR();

	

	
	std::cout << "Running with initial polymer density " << Ntad / ((double) Ntot) << std::endl;
	std::cout << "Using " << Ntad << " TADs, including main chain of length " << Nchain << std::endl;
}

void MCPoly::CreateBond(MCBond& bond)
{
	MCTad* tad1 = &tadConf[bond.id1];
	MCTad* tad2 = &tadConf[bond.id2];
	
	int id1 = ( !bond.set && (tad1->links == 2) ) ? 2 : 1;
	int id2 = ( !bond.set && (tad2->links == 2) ) ? 2 : 0;

	tad1->neighbors[id1] = tad2;
	tad2->neighbors[id2] = tad1;

	tad1->bonds[id1] = &bond;
	tad2->bonds[id2] = &bond;
	
	if ( !bond.set || (tad1->links == 0) )
		++tad1->links;
	if ( !bond.set || (tad2->links == 0) )
		++tad2->links;
	
	bond.set = true;
}

void MCPoly::GenerateRandom(int lim)
{
	Ntad = Nchain;
	Nbond = Nchain-1;
	
	tadConf.resize(Ntad);
	tadTopo.resize(Nbond);
	
	for ( int b = 0; b < Nbond; ++b )
	{
		tadTopo[b].id1 = b;
		tadTopo[b].id2 = b+1;
	}
	
	int turn1[7];
	
	turn1[0] = 12;
	turn1[1] = 12;
	turn1[2] = 1;
	turn1[3] = 1;
	turn1[4] = 11;
	turn1[5] = 11;
	turn1[6] = 2;

	int turn2[7];
	
	turn2[0] = 12;
	turn2[1] = 1;
	turn2[2] = 1;
	turn2[3] = 11;
	turn2[4] = 11;
	turn2[5] = 2;
	turn2[6] = 2;
	
	int vi = 2*CUB(L) + SQR(L) + L/2; // Set to lat->rngEngine() % Ntot for random chromosome placement
	
	tadConf[0].pos = vi;
	lat->bitTable[0][vi] = 1;
	
	int ni = 1;
	
	for ( int i = 0; i < lim; ++i )
	{
		for ( int j = 0; j < 7; ++j )
		{
			int turn = ((i % 2) == 0) ? turn1[j] : turn2[j];
			
			tadTopo[ni-1].dir = turn;
			tadConf[ni].pos = lat->bitTable[turn][tadConf[ni-1].pos];
			
			lat->bitTable[0][tadConf[ni].pos] = 1;
			
			++ni;
		}
		
		tadTopo[ni-1].dir = 10;
		tadConf[ni].pos = lat->bitTable[10][tadConf[ni-1].pos];
		
		lat->bitTable[0][tadConf[ni].pos] = 1;
		
		++ni;
	}
	
	--ni;
	
	while ( ni < Nbond )
	{
		int t = lat->rngEngine() % ni;
		int iv = lat->rngEngine() % lat->nbNN[0][0][tadTopo[t].dir];
		
		int nd1 = lat->nbNN[2*iv+1][0][tadTopo[t].dir];
		int nd2 = lat->nbNN[2*(iv+1)][0][tadTopo[t].dir];
		
		int en2 = tadConf[t].pos;
		int v1 = (nd1 == 0) ? en2 : lat->bitTable[nd1][en2];
		
		int b = lat->bitTable[0][v1];
					
		if ( b == 0 )
		{
			for ( int i = ni+1; i > t+1; --i )
			{
				tadConf[i].pos = tadConf[i-1].pos;
				tadTopo[i].dir = tadTopo[i-1].dir;
			}
			
			tadConf[t+1].pos = v1;
			
			tadTopo[t].dir = nd1;
			tadTopo[t+1].dir = nd2;

			lat->bitTable[0][v1] = 1;
			
			++ni;
		}
	}
	
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		for ( int i = 0; i < 3; ++i )
			centerMass[i] += lat->xyzTable[i][tad->pos] / Ntad;
	}
}
void MCPoly::GenerateCAR()
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::vector<double>  weights = { 0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  14.27626036,  14.27626036,   0.        ,
										   0.        ,   6.67564303,   6.67564303,   0.        ,
										   0.        ,   7.49082708,   7.49082708,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  50.83845243,  50.83845243,   0.        ,
										   0.        ,  13.34882586,  13.34882586,   0.        ,
										   0.        ,   5.38868585,   5.38868585,   0.        ,
										   0.        ,   0.        ,   0.        ,   3.62944389,
										   3.62944389,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  48.94537977,
										   48.94537977,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  20.18919962,
										   20.18919962,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   6.86066502,
										   6.86066502,   0.        ,   0.        ,  38.14000157,
										   38.14000157,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  11.38773356,  11.38773356,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  82.74387584,  82.74387584,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  14.83973734,  14.83973734,   0.        ,
										   0.        ,   0.        ,   0.        ,  32.64170696,
										   32.64170696,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   2.35807053,
										   2.35807053,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  68.23137377,
										   68.23137377,   0.        ,   0.        ,   0.        ,
										   0.        ,   5.9424787 ,   5.9424787 ,   0.        ,
										   0.        ,   1.77550621,   1.77550621,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  18.38090779,
										   18.38090779,   0.        ,   0.        ,   0.        ,
										   0.        ,  47.21823968,  47.21823968,   0.        ,
										   0.        ,   6.0830042 ,   6.0830042 ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  31.79419245,  31.79419245,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   2.70838661,
										   2.70838661,   0.        ,   0.        ,   0.        ,
										   0.        ,   7.43685078,   7.43685078,   0.        ,
										   0.        ,  10.84705691,  10.84705691,   0.        ,
										   0.        ,   5.26489901,   5.26489901,   0.        ,
										   0.        ,  57.06845937,  57.06845937,   0.        ,
										   0.        ,   4.13590402,   4.13590402,   0.        ,
										   0.        ,   0.        ,   0.        ,  26.62566264,
										   26.62566264,   0.        ,   0.        ,   0.        ,
										   0.        ,  22.59782564,  22.59782564,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   3.2438692 ,
										   3.2438692 ,   0.        ,   0.        ,   0.        ,
										   0.        ,  52.43152972,  52.43152972,   0.        ,
										   0.        ,   0.        ,   0.        ,  10.62196365,
										   10.62196365,   0.        ,   0.        ,   6.08010273,
										   6.08010273,   0.        ,   0.        ,  27.5173382 ,
										   27.5173382 ,   0.        ,   0.        ,   0.        ,
										   0.        ,   4.31336569,   4.31336569,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  36.42149416,
										   36.42149416,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   4.59219405,
										   4.59219405,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  21.25370714,
										   21.25370714,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  59.13374107,
										   59.13374107,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  12.09382979,  12.09382979,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  40.67129366,  40.67129366,   0.        ,
										   0.        ,   0.        ,   0.        ,  25.75381863,
										   25.75381863,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   4.31862223,
										   4.31862223,   0.        ,   0.        ,  42.81995705,
										   42.81995705,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   7.18842353,
										   7.18842353,   0.        ,   0.        ,   0.        ,
										   0.        ,  30.77592557,  30.77592557,   0.        ,
										   0.        ,  20.06957348,  20.06957348,   0.        ,
										   0.        ,  12.21239846,  12.21239846,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  16.28098852,
										   16.28098852,   0.        ,   0.        ,  14.69736163,
										   14.69736163,   0.        ,   0.        ,   0.        ,
										   0.        ,  29.86362233,  29.86362233,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  36.81557448,  36.81557448,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   9.51038517,   9.51038517,   0.        ,
										   0.        ,   0.        ,   0.        ,   3.53905236,
										   3.53905236,   0.        ,   0.        ,  64.45170046,
										   64.45170046,   0.        ,   0.        ,  37.36763813,
										   37.36763813,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        , 123.2460735 ,
										   123.2460735 ,   0.        ,   0.        ,   0.        ,
										   0.        ,  81.19332237,  81.19332237,   0.        ,
										   0.        ,   0.        ,   0.        ,  95.99293587,
										   95.99293587,   0.        ,   0.        ,  73.48476231,
										   73.48476231,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  23.00618088,  23.00618088,   0.        ,
										   0.        ,   7.60203323,   7.60203323,   0.        ,
										   0.        ,  55.62943702,  55.62943702,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  16.41008893,
										   16.41008893,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  73.99881282,  73.99881282,   0.        ,
										   0.        ,  18.65211009,  18.65211009,   0.        ,
										   0.        ,   0.        ,   0.        ,  13.94028737,
										   13.94028737,   0.        ,   0.        ,   0.        ,
										   0.        ,   6.02414586,   6.02414586,   0.        ,
										   0.        ,   5.10896347,   5.10896347,   0.        ,
										   0.        ,  39.60084123,  39.60084123,   0.        ,
										   0.        ,   0.        ,   0.        ,   4.96229659,
										   4.96229659,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  17.4096751 ,
										   17.4096751 ,   0.        ,   0.        ,   0.        ,
										   0.        ,  68.72975746,  68.72975746,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   7.40871035,
										   7.40871035,   0.        ,   0.        ,   0.        ,
										   0.        ,  33.24269127,  33.24269127,   0.        ,
										   0.        ,  12.13854765,  12.13854765,   0.        ,
										   0.        ,   0.        ,   0.        ,  21.79969402,
										   21.79969402,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   3.90731672,   3.90731672,   0.        ,
										   0.        ,   5.47353201,   5.47353201,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  13.79992119,  13.79992119,   0.        ,
										   0.        ,  64.20658772,  64.20658772,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  26.0420695 ,  26.0420695 ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  24.48564135,  24.48564135,   0.        ,
										   0.        ,  13.26902743,  13.26902743,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  20.7656081 ,
										   20.7656081 ,   0.        ,   0.        ,  34.27498647,
										   34.27498647,   0.        ,   0.        ,   0.        ,
										   0.        ,  21.84644346,  21.84644346,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  77.92069811,
										   77.92069811,   0.        ,   0.        ,   0.        ,
										   0.        ,   2.95056537,   2.95056537,   0.        ,
										   0.        ,   0.        ,   0.        ,  11.8125279 ,
										   11.8125279 ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  20.92211171,
										   20.92211171,   0.        ,   0.        ,  30.82386768,
										   30.82386768,   0.        ,   0.        ,  35.93547695,
										   35.93547695,   0.        ,   0.        ,   0.        ,
										   0.        ,  18.90583297,  18.90583297,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  28.28064854,  28.28064854,   0.        ,
										   0.        ,   7.69223731,   7.69223731,   0.        ,
										   0.        ,   0.        ,   0.        ,  16.19796187,
										   16.19796187,   0.        ,   0.        ,   9.99566822,
										   9.99566822,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  30.88520044,
										   30.88520044,   0.        ,   0.        ,   0.        ,
										   0.        ,  15.73928145,  15.73928145,   0.        ,
										   0.        ,  69.85321963,  69.85321963,   0.        ,
										   0.        ,   3.17062835,   3.17062835,   0.        ,
										   0.        ,   0.        ,   0.        ,   4.93342523,
										   4.93342523,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  28.40941272,  28.40941272,   0.        ,
										   0.        ,   4.60117337,   4.60117337,   0.        ,
										   0.        ,  30.11536822,  30.11536822,   0.        ,
										   0.        ,   0.        ,   0.        ,   4.67043308,
										   4.67043308,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  17.27665656,
										   17.27665656,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  16.53537024,
										   16.53537024,   0.        ,   0.        ,   0.        ,
										   0.        ,  19.27773877,  19.27773877,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  52.78562262,  52.78562262,   0.        ,
										   0.        ,  15.95216076,  15.95216076,   0.        ,
										   0.        ,   0.        ,   0.        ,  52.22702069,
										   52.22702069,   0.        ,   0.        ,   6.87252377,
										   6.87252377,   0.        ,   0.        ,   0.        ,
										   0.        ,  39.42432575,  39.42432575,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   8.01235673,   8.01235673,   0.        ,
										   0.        ,   0.        ,   0.        ,   5.73486297,
										   5.73486297,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,  11.30434909,  11.30434909,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  23.23952551,
										   23.23952551,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  49.93764114,
										   49.93764114,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,   3.93019451,
										   3.93019451,   0.        ,   0.        ,   0.        ,
										   0.        ,   1.95903063,   1.95903063,   0.        ,
										   0.        ,   0.        ,   0.        ,  38.99492352,
										   38.99492352,   0.        ,   0.        ,   0.        ,
										   0.        ,   0.        ,   0.        ,  15.49085337,
										   15.49085337,   0.        ,   0.        ,   0.        ,
										   0.        ,  49.6599029 ,  49.6599029 ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.        ,
										   0.        ,   5.46991939,   5.46991939,   0.        ,
										   0.        ,   0.        ,   0.        ,  22.23939083,
										   22.23939083,   0.        ,   0.        ,   0.        ,
										   0.        ,  17.6775166 ,  17.6775166 ,   0.        ,
										   0.        ,   0.        ,   0.        ,   0.    };
	
	std::discrete_distribution<> d(weights.begin(), weights.end());
	
	std::vector<int> CAR_sample;
	
	while((int) CAR_sample.size() < 160)
	{
		
		int sampled_car=d(gen);
		if(std::find(CAR_sample.begin(),CAR_sample.end(),sampled_car) == CAR_sample.end())
			CAR_sample.push_back(sampled_car);
		
	}

	std::sort (CAR_sample.begin(), CAR_sample.end());
	auto CAR_sample_copy = CAR_sample;
	while(CAR_sample.size()>0)
	{
		
		int t = (lat->rngEngine() % (int) CAR_sample.size());

		double rnd = lat->rngDistrib(lat->rngEngine);
		int i = rnd > 0.5 ? 1 : -1;
		if( CAR_sample.at(t)==CAR_sample_copy.at(0))
			i=1;
		if( CAR_sample.at(t)==CAR_sample_copy.back())
			i=-1;
		
		int t2=(int)std::distance(CAR_sample_copy.begin(), std::find(CAR_sample_copy.begin(), CAR_sample_copy.end(),CAR_sample.at(t)));
		
		
		if(std::find(CAR_sample.begin(),CAR_sample.end(),CAR_sample_copy.at(t2+i)) == CAR_sample.end())
		{
			if(CAR_sample.at(t)==CAR_sample_copy.back() or CAR_sample.at(t)==CAR_sample_copy.at(0))
				CAR_sample.erase(CAR_sample.begin()+t);
			else if(std::find(CAR_sample.begin(),CAR_sample.end(),CAR_sample_copy.at(t2-i)) == CAR_sample.end() or CAR_sample_copy.at(t2-i)==CAR_sample_copy.at(t2)-i )
				CAR_sample.erase(CAR_sample.begin()+t);
			
		}
		else if(CAR_sample_copy.at(t2+i)!=CAR_sample_copy.at(t2)+i)
		{
			
			CAR.push_back(&tadConf.at(CAR_sample_copy.at(t2)));
			CAR.back()->choesin_binding_site = &tadConf.at(CAR_sample_copy.at(t2+i));
			CAR.back()->choesin_binding_site->choesin_binding_site=CAR.back();
			CAR_sample.erase(CAR_sample.begin()+t);
			int next_anchor=(int)std::distance(CAR_sample.begin(), std::find(CAR_sample.begin(), CAR_sample.end(),CAR_sample_copy.at(t2+i)));
			CAR_sample.erase(CAR_sample.begin()+ next_anchor);
			
			std::cout <<"choesin "<< CAR_sample_copy.at(t2) <<" and "<< CAR_sample_copy.at(t2+i) << std::endl;

		}
	}
	for ( int t = 0; t < (int) CAR.size(); ++t )
		std::cout << CAR.at(t) << std::endl;

}
void MCPoly::TrialMove(double* dE)
{

	double rnd = lat->rngDistrib(lat->rngEngine);
	if(rnd < (double) Ntad/(Ntad+0*activeForks.size())){
		int t = lat->rngEngine() % Ntad;
		tadTrial = &tadConf[t];
		tadUpdater->TrialMove(tadTrial, dE);
		*dE = tadUpdater->legal ? *dE : 0.;
		if(t==74)
			nreptation74.push_back((int) tadUpdater->reptation_values.size());
		if(t==73)
			nreptation73.push_back((int) tadUpdater->reptation_values.size());
		if(t==64)
			nreptation64.push_back((int) tadUpdater->reptation_values.size());
		if(t==44)
			nreptation44.push_back((int) tadUpdater->reptation_values.size());
		if(t==14)
			nreptation14.push_back((int) tadUpdater->reptation_values.size());



			
	}else{
		int forkID = lat->rngEngine() % (int) activeForks.size();
		tadTrial = activeForks[forkID];
		tadUpdater->TrialMove(tadTrial, dE);
		*dE = tadUpdater->legal ? *dE : 0.;
	}
	
}

void MCPoly::AcceptMove()
{
	
	tadUpdater->AcceptMove(tadTrial);
	

	for ( int t = 0; t < (int) tadUpdater->reptation_values.size(); ++t )
	{
		--lat->bitTable[0][tadUpdater->reptation_values[t][0]];
		++lat->bitTable[0][tadUpdater->reptation_values[t][1]];
	}
	
	
	
	for ( int i = 0; i < (int) CAR.size(); ++i )
	{
		bool co_lococalized=false;
		if(CAR.at(i)->pos == CAR.at(i)->choesin_binding_site->pos)
			co_lococalized=true;
		for ( int v = 0; v < 12; ++v )
			if(CAR.at(i)->pos == lat->bitTable[v+1][CAR.at(i)->choesin_binding_site->pos])
				co_lococalized=true;
		
		
		if (co_lococalized==true)
		{
			for ( int t = 0; t < (int) tadConf.size(); ++t )
			{
				if(&tadConf.at(t)==CAR.at(i))
					std::cout <<"BINDED "<< t <<"with "<<std::endl;
				
			}
			for ( int t = 0; t < (int) tadConf.size(); ++t )
			{
				if(&tadConf.at(t)==CAR.at(i)->choesin_binding_site)
					std::cout << t <<std::endl;
				
			}
			
			CAR.at(i)->isChoesin=true;
			CAR.at(i)->choesin_binding_site->isChoesin=true;
			CAR.at(i)->neighbors[2] = CAR.at(i)->choesin_binding_site;
			CAR.at(i)->neighbors[2]->neighbors[2] = CAR.at(i);
			
		}
	}
	
	for ( int i = 0; i < (int) interCAR.size(); ++i )
	{
		if(interCAR.at(i)->SisterID!=-1)
		{
			interCAR.at(i)->choesin_binding_site=&tadConf.at(interCAR.at(i)->SisterID);
			tadConf.at(interCAR.at(i)->SisterID).choesin_binding_site=interCAR.at(i);
		}
	}
	
	for ( int i = 0; i < (int) interCAR.size(); ++i )
	{
		if(interCAR.at(i)->SisterID!=-1)
		{
			
			bool co_lococalized=false;
			if(interCAR.at(i)->pos == interCAR.at(i)->choesin_binding_site->pos)
				co_lococalized=true;
			for ( int v = 0; v < 12; ++v )
				if(interCAR.at(i)->pos == lat->bitTable[v+1][interCAR.at(i)->choesin_binding_site->pos])
					co_lococalized=true;
			
			
			if (co_lococalized==true)
			{
				interCAR.at(i)->isChoesin=true;
				interCAR.at(i)->choesin_binding_site->isChoesin=true;
				interCAR.at(i)->neighbors[2] = interCAR.at(i)->choesin_binding_site;
				interCAR.at(i)->neighbors[2]->neighbors[2] = interCAR.at(i);
			}
		}
	}
	
	for ( int i = 0; i < (int) CAR.size(); ++i )
		if(CAR.at(i)->isChoesin)
			CAR.erase(CAR.begin()+i);
	for ( int i = 0; i < (int) interCAR.size(); ++i )
		if(interCAR.at(i)->isChoesin)
			interCAR.erase(interCAR.begin()+i);
	
}

void MCPoly::ToVTK(int frame)
{
	char fileName[32];
	sprintf(fileName, "poly%05d.vtp", frame);
	
	std::string path = outputDir + "/" + fileName;
	
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto lines = vtkSmartPointer<vtkCellArray>::New();
	
	auto types = vtkSmartPointer<vtkIntArray>::New();
	auto forks = vtkSmartPointer<vtkIntArray>::New();
	auto status = vtkSmartPointer<vtkIntArray>::New();
	auto SisterIDs = vtkSmartPointer<vtkIntArray>::New();
	auto choesins = vtkSmartPointer<vtkIntArray>::New();


	
	types->SetName("TAD type");
	types->SetNumberOfComponents(1);
	
	forks->SetName("Fork type");
	forks->SetNumberOfComponents(1);
	
	status->SetName("Replication status");
	status->SetNumberOfComponents(1);
	
	SisterIDs->SetName("SisterID");
	SisterIDs->SetNumberOfComponents(1);
	
	choesins->SetName("Choesin");
	choesins->SetNumberOfComponents(1);
	
	std::vector<double3> conf = GetPBCConf();

	for ( int t = 0; t < Ntad; ++t )
	{
		int type = tadConf[t].type;
		int state = tadConf[t].status;
		int fork = tadConf[t].isFork() ? (tadConf[t].isLeftFork() ? -1 : 1) : 0;
		int sisterID = tadConf[t].SisterID;
		int choesin = tadConf[t].isChoesin;

		
		points->InsertNextPoint(conf[t][0], conf[t][1], conf[t][2]);
		
		types->InsertNextValue(type);
		forks->InsertNextValue(fork);
		status->InsertNextValue(state);
		SisterIDs->InsertNextValue(sisterID);
		choesins->InsertNextValue(choesin);

	}
	
	for ( auto bond = tadTopo.begin(); bond != tadTopo.end(); ++bond )
	{
		auto line = vtkSmartPointer<vtkLine>::New();
		
		line->GetPointIds()->SetId(0, bond->id1);
		line->GetPointIds()->SetId(1, bond->id2);
	
		lines->InsertNextCell(line);
	}
	
	auto polyData = vtkSmartPointer<vtkPolyData>::New();
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	polyData->SetPoints(points);
	polyData->SetLines(lines);
	
	polyData->GetPointData()->AddArray(types);
	polyData->GetPointData()->AddArray(forks);
	polyData->GetPointData()->AddArray(status);
	polyData->GetPointData()->AddArray(SisterIDs);
	polyData->GetPointData()->AddArray(choesins);


	
	writer->SetFileName(path.c_str());
	writer->SetInputData(polyData);
	
 	writer->Write();
}

void MCPoly::FromVTK(int frame)
{
	char fileName[32];
	sprintf(fileName, "poly%05d.vtp", frame);
	
	std::string path = outputDir + "/" + fileName;

	std::cout << "Starting from polymer configuration file " << path << std::endl;

	auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

	reader->SetFileName(path.c_str());
	reader->Update();
	
	vtkPolyData* polyData = reader->GetOutput();
	
	vtkCellArray* lineData = polyData->GetLines();
	vtkDataArray* typeData = polyData->GetPointData()->GetArray("TAD type");
	vtkDataArray* statusData = polyData->GetPointData()->GetArray("Replication status");


	Ntad = (int) polyData->GetNumberOfPoints();
	Nbond = (int) polyData->GetNumberOfLines();
	
	tadConf.resize(Ntad);
	tadTopo.resize(Nbond);
			
	for ( int t = 0; t < Ntad; ++t )
	{
		double point[3];
		
		polyData->GetPoint(t, point);
		
		tadConf[t].type = (int) typeData->GetComponent(t, 0);
		tadConf[t].status = (int) statusData->GetComponent(t, 0);

		
		for ( int i = 0; i < 3; ++i )
		{
			centerMass[i] += point[i] / ((double) Ntad);

			while ( point[i] >= L ) point[i] -= L;
			while ( point[i] < 0 )  point[i] += L;
		}

		int ixp = (int) 1*point[0];
		int iyp = (int) 2*point[1];
		int izp = (int) 4*point[2];
		
		tadConf[t].pos = ixp + iyp*L + izp*L2;
		
		++lat->bitTable[0][tadConf[t].pos];
	}
	
	for ( int b = 0; b < Nbond; ++b )
	{
		auto cell = vtkSmartPointer<vtkIdList>::New();
		
		lineData->GetCellAtId(b, cell);

		int t1 = (int) cell->GetId(0);
		int t2 = (int) cell->GetId(1);

		tadTopo[b].id1 = t1;
		tadTopo[b].id2 = t2;
		
		if ( tadConf[t1].pos == tadConf[t2].pos )
			tadTopo[b].dir = 0;
		
		else
		{
			for ( int v = 0; v < 12; ++v )
			{
				if ( lat->bitTable[v+1][tadConf[t1].pos] == tadConf[t2].pos )
				{
					tadTopo[b].dir = v+1;
					break;
				}
			}
		}
	}
	
	auto lastBond = std::find_if(tadTopo.begin(), tadTopo.end(), [](const MCBond& b){return b.id2 != b.id1+1;});
	int length = (int) std::distance(tadTopo.begin(), lastBond) + 1;
	
	if ( length != Nchain )
		throw std::runtime_error("MCPoly: Found incompatible main chain dimension " + std::to_string(length));
}

std::vector<double3> MCPoly::GetPBCConf()
{
	std::vector<double3> conf(Ntad);

	for ( int t = 0; t < Ntad; ++t )
	{
		for ( int i = 0; i < 3; ++i )
			conf[t][i] = lat->xyzTable[i][tadConf[t].pos];
	}
	
	for ( auto bond = tadTopo.begin(); bond != tadTopo.end(); ++bond )
		FixPBCPair(conf, bond->id1, bond->id2);
	
	centerMass = GetPBCCenterMass(conf.begin(), conf.end());
	
	return conf;
}

void MCPoly::FixPBCPair(std::vector<double3>& conf, int id1, int id2)
{
	double3* pos1 = &conf[id1];
	double3* pos2 = &conf[id2];

	for ( int i = 0; i < 3; ++i )
	{
		double deltaTad = (*pos2)[i] - (*pos1)[i];
		
		while ( std::abs(deltaTad) > L/2. )
		{
			double pbcShift = std::copysign(L, deltaTad);

			(*pos2)[i] -= pbcShift;
			deltaTad -= pbcShift;
		}
	}
}

double3 MCPoly::GetPBCCenterMass(std::vector<double3>::iterator end1, std::vector<double3>::iterator end2)
{
	double3 center = {0., 0., 0.};

	for ( int i = 0; i < 3; ++i )
	{
		for ( auto tadPos = end1; tadPos != end2; ++tadPos )
			center[i] += (*tadPos)[i] / ((double) std::distance(end1, end2));

		double deltacenterMass = center[i] - centerMass[i];
		
		// Translate chain center of mass into the same box as the previous conformation
		while ( std::abs(deltacenterMass) > L/2. )
		{
			double pbcShift = std::copysign(L, deltacenterMass);
			
			for ( auto tadPos = end1; tadPos != end2; ++tadPos )
				(*tadPos)[i] -= pbcShift;
			
			deltacenterMass -= pbcShift;
			center[i] -= pbcShift;
		}
	}
	
	return center;
}
void MCPoly::OriginMove(MCTad*)
{}
void MCPoly::ForkMove()
{}
