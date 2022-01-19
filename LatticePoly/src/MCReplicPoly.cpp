//
//  MCReplicPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include "MCReplicPoly.hpp"

#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <random>

MCReplicPoly::MCReplicPoly(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCReplicPoly::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);


	
	MCsteps=0;
	MCrepl=0;
	activeForks.reserve(Nchain);
	neigh=false;

	// Locate existing forks
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		if ( tad->isFork() )
			activeForks.push_back(&(*tad));
	}
	//chr4
	/*
	origins={0,7,12,17,36,40,68,76,80,99,110,126,170,185,188,203,205,253,263,281,326,348,355,370,381,387,
		404,444,454,503,511,515,562,576,598,602,644,676,687,703,719,723,731,737,754,780,813,817,
		826,846,888,904,927,932,963,992,1021,1042,1046,1082,1103,1108,1123,1158,1163,1169,1189,1199,1202,1204,1219};
	
 mrt={1000,1000,1000,0.9208379238843918,0.7043074965476991,0.7438885271549225,0.7962751537561417,0.8440044820308685,0.8253782838582993,0.4947620034217834,0.642608255147934,0.8812578544020653,0.5261937975883484,0.6554138362407684,0.6670551002025604,0.7322472929954529,0.6728757321834564,0.3876601457595825,0.20954644680023196,0.6647264659404755,0.34575146436691284,0.2060534954071045,0.3015140891075134,0.16298043727874756,0.32828861474990845,0.4039586782455444,0.4051220417022705,0.07683342695236206,0.2630963921546936,0.7741559892892838,0.7334106266498566,0.7415598928928375,0.6938305497169495,0.8661236464977264,0.6821883320808411,0.6880099475383759,0.7229337096214294,0.9208379238843918,0.850989431142807,1000,0.22118771076202395,0.1548311710357666,0.0,0.09313195943832396,0.637950986623764,0.9289871901273729,0.3725259304046631,0.3597203493118286,0.6123398542404175,0.551804929971695,0.7881258875131607,0.6682194173336029,0.0861470103263855,0.18509960174560547,0.8672879636287689,0.6949939131736755,0.7660067230463028,0.5506406128406525,0.7834695726633072,0.4761348366737366,0.8416768163442612,0.7415598928928375,0.7334106266498566,0.7462162077426909,0.6821883320808411,0.5250294804573059,0.8125727027654648,0.850989431142807,0.8195576518774033,0.8428401648998259,0.8800935298204422};
	*/
	
	//chr 7
	origins={1,6,14,20,26,51,55,89,91,94,130,133,137,149,163,185,192,213,215,220,228,230,255,270,282,311,336,365,385,388,407,431,442,454,459,460,465,473,485,497,501,523,527,572,594,612,622,636,64,667,677,690,706,710,733,776,782,799,801,804,829,850,859,866,871};
	
	weights={0.0008441 , 0.00113861, 0.00615295, 0.0011208 , 0.00126901,0.01631902, 0.00450673, 0.00759752, 0.00806823, 0.00269195,0.0199772 , 0.01153813, 0.002372  , 0.00223142, 0.03786356,0.00247823, 0.00170537, 0.00148846, 0.00186694, 0.00142549,0.04605837, 0.04407694, 0.00312641, 0.00132944, 0.00448002,0.0357918 , 0.07323744, 0.00183768, 0.05248039, 0.08644023,0.06210324, 0.00209721, 0.00130972, 0.02321365, 0.01614855,0.01131486, 0.00320274, 0.00343364, 0.00220153, 0.00144076,0.00138033, 0.00346672, 0.04149439, 0.04672691, 0.00159914,0.00221997, 0.06408213, 0.00273266, 0.00133516, 0.08104678,0.00311496, 0.00147574, 0.01425362, 0.07897057, 0.00378222,0.00140005, 0.00499844, 0.01212906, 0.01337517, 0.00358758,0.00097005, 0.00169265, 0.00099994, 0.00379431, 0.00139114};
	
		
	
	//origins={20,40,60,80,100,120,140,160,180};
	//origins={10,30,50,70,90,110,130,150,170,190,20,40,60,80,100,120,140,160,180};
	


	/*
	for (int i = 1; i < (int)origins.size()-1; ++i)
	tadConf[origins[i]].isCAR=true;

	for (int i = 1; i < (int) Nchain-1; ++i)
		origins.push_back(i);

	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<> d({0.1294,  0.1327,  0.1359,  0.1433,  0.1474,  0.1536,  0.179 ,
		0.2045,  0.2494,  0.3211,  0.4307,  0.7594,  0.9349,  1.0427,
		0.9673,  0.7989,  0.6131,  0.451 ,  0.3311,  0.2043,  0.1762,
		0.1622,  0.1572,  0.1626,  0.1733,  0.187 ,  0.1995,  0.2096,
		0.2026,  0.1867,  0.1523,  0.1408,  0.1337,  0.1298,  0.1287,
		0.1322,  0.1362,  0.1408,  0.1451,  0.1428,  0.1381,  0.1324,
		0.1285,  0.1263,  0.1297,  0.1395,  0.2108,  0.3129,  0.5246,
		0.9052,  1.4703,  2.5655,  2.6783,  2.3125,  1.1384,  0.7085,
		0.4488,  0.3079,  0.2366,  0.1857,  0.1812,  0.185 ,  0.2007,
		0.208 ,  0.2099,  0.2026,  0.1883,  0.1559,  0.1432,  0.1344,
		0.128 ,  0.129 ,  0.132 ,  0.1363,  0.1417,  0.1516,  0.1532,
		0.1517,  0.1434,  0.1398,  0.1379,  0.1383,  0.1411,  0.1641,
		0.192 ,  0.2382,  0.4542,  0.6638,  0.9241,  1.1944,  1.3815,
		1.2684,  1.0351,  0.7836,  0.4232,  0.3232,  0.2636,  0.2253,
		0.2036,  0.1858,  0.1842,  0.1852,  0.1899,  0.191 ,  0.1901,
		0.1874,  0.1834,  0.1732,  0.1689,  0.1662,  0.1625,  0.162 ,
		0.163 ,  0.1638,  0.1653,  0.1691,  0.1721,  0.1745,  0.1784,
		0.181 ,  0.1838,  0.189 ,  0.1985,  0.2472,  0.3065,  0.4193,
		0.9499,  1.4275,  2.0423,  2.6719,  3.1406,  3.0067,  2.4702,
		1.8139,  0.8627,  0.608 ,  0.4601,  0.3729,  0.3243,  0.2834,
		0.2751,  0.267 ,  0.248 ,  0.238 ,  0.2291,  0.225 ,  0.228 ,
		0.2602,  0.2966,  0.3508,  0.4867,  0.5403,  0.5508,  0.5229,
		0.468 ,  0.3836,  0.3727,  0.3997,  0.6853,  1.0632,  1.7711,
		2.8392,  4.1044,  5.9525,  5.9932,  5.2951,  2.6755,  1.5836,
		0.9079,  0.5529,  0.3819,  0.2692,  0.2668,  0.2879,  0.3844,
		0.4469,  0.51  ,  0.5534,  0.5647,  0.5098,  0.4711,  0.4377,
		0.3968,  0.3921,  0.3903,  0.3896,  0.3849,  0.3677,  0.3524,
		0.3302,  0.2896,  0.2767,  0.2681,  0.2615,  0.2559,  0.2419,
		0.2315,  0.2204,  0.2049,  0.2008,  0.1986,  0.1983,  0.2006,
		0.2075,  0.213 ,  0.2171,  0.2211,  0.2204,  0.2161,  0.21  ,
		0.2072,  0.2113,  0.219 ,  0.234 ,  0.2753,  0.2935,  0.3058,
		0.3044,  0.289 ,  0.2388,  0.2241,  0.2199,  0.2742,  0.3814,
		0.6194,  1.1629,  2.1903,  5.6969,  7.2408,  7.9619,  6.9293,
		5.2026,  3.3755,  1.868 ,  1.0359,  0.4093,  0.3168,  0.2903,
		0.322 ,  0.3665,  0.4256,  0.4734,  0.5013,  0.4485,  0.394 ,
		0.3425,  0.2887,  0.2887,  0.3041,  0.3398,  0.3924,  0.5284,
		0.5823,  0.6157,  0.5611,  0.4915,  0.4232,  0.3682,  0.3274,
		0.2709,  0.2526,  0.2347,  0.2099,  0.2017,  0.1997,  0.1999,
		0.2024,  0.2078,  0.2102,  0.2112,  0.209 ,  0.2054,  0.2021,
		0.1979,  0.1968,  0.2074,  0.2229,  0.2489,  0.3407,  0.4133,
		0.4998,  0.6065,  0.7043,  0.8492,  0.8561,  0.8219,  0.6675,
		0.5759,  0.4936,  0.4155,  0.3487,  0.2605,  0.2355,  0.2197,
		0.2121,  0.217 ,  0.2295,  0.2469,  0.2701,  0.3182,  0.3323,
		0.3358,  0.3201,  0.3151,  0.3182,  0.3415,  0.3984,  0.7649,
		1.1918,  1.9349,  4.4039,  5.6268,  6.3535,  6.2523,  5.2894,
		2.3278,  1.2543,  0.6457,  0.2704,  0.2254,  0.21  ,  0.2211,
		0.2488,  0.3689,  0.4453,  0.49  ,  0.4246,  0.3644,  0.315 ,
		0.2913,  0.3028,  0.5661,  1.0867,  2.2798,  7.5163, 10.0658,
		11.5136, 11.9793, 11.5885,  7.7993,  4.9719,  2.7204,  0.8267,
		0.5598,  0.4428,  0.4238,  0.4554,  0.7042,  0.8991,  1.0523,
		1.0816,  0.9393,  0.7576,  0.5878,  0.4633,  0.3102,  0.273 ,
		0.2497,  0.2346,  0.2376,  0.2412,  0.2516,  0.2556,  0.2708,
		0.2809,  0.2889,  0.275 ,  0.2559,  0.2374,  0.225 ,  0.222 ,
		0.262 ,  0.328 ,  0.4406,  0.7892,  0.9295,  0.9112,  0.8034,
		0.6421,  0.4374,  0.4266,  0.4676,  1.2007,  2.4871,  4.8634,
		8.2504, 11.1229, 13.5455, 13.5892, 13.0005,  8.349 ,  4.833 ,
		2.3422,  1.0389,  0.5106,  0.2   ,  0.1644,  0.1524,  0.1598,
		0.1815,  0.2327,  0.3444,  0.5959,  2.1752,  3.9134,  6.0175,
		9.2096,  9.7632,  9.4783,  7.9791,  5.7335,  1.6992,  0.8439,
		0.4846,  0.2894,  0.2952,  0.351 ,  0.4744,  0.6525,  1.0821,
		1.1514,  1.0629,  0.668 ,  0.5219,  0.4139,  0.3461,  0.3123,
		0.3088,  0.3162,  0.3283,  0.3393,  0.3297,  0.3026,  0.2746,
		0.2507,  0.2164,  0.2058,  0.2004,  0.2009,  0.2001,  0.2006,
		0.2027,  0.2059,  0.2215,  0.2303,  0.2348,  0.2405,  0.2497,
		0.2697,  0.3085,  0.3789,  0.7667,  1.2103,  1.8809,  3.6494,
		4.3559,  4.659 ,  4.6156,  4.1635,  2.5387,  1.7788,  1.2601,
		0.6957,  0.5733,  0.5198,  0.5035,  0.5167,  0.5938,  0.6337,
		0.6547,  0.6501,  0.6189,  0.5788,  0.5398,  0.5148,  0.4973,
		0.4994,  0.5082,  0.519 ,  0.5107,  0.4931,  0.4704,  0.438 ,
		0.3797,  0.3606,  0.3461,  0.3288,  0.3248,  0.3158,  0.3069,
		0.2965,  0.288 ,  0.2869,  0.283 ,  0.2682,  0.2549,  0.2403,
		0.2265,  0.2153,  0.2074,  0.2103,  0.217 ,  0.2394,  0.2511,
		0.2597,  0.2597,  0.2543,  0.2238,  0.2066,  0.1945,  0.1859,
		0.1882,  0.1973,  0.2115,  0.2286,  0.2517,  0.2517,  0.2395,
		0.206 ,  0.1966,  0.1961,  0.2089,  0.2475,  0.545 ,  0.9853,
		1.8676,  5.0205,  6.5233,  7.3451,  7.3492,  6.4516,  3.1467,
		1.8009,  0.9995,  0.4074,  0.32  ,  0.2839,  0.2767,  0.2903,
		0.3402,  0.3588,  0.3657,  0.3324,  0.3037,  0.2754,  0.2538,
		0.2429,  0.2421,  0.2449,  0.2471,  0.2307,  0.212 ,  0.1931,
		0.1764,  0.1643,  0.1557,  0.1602,  0.1704,  0.2093,  0.2343,
		0.254 ,  0.2595,  0.2529,  0.2385,  0.2407,  0.2555,  0.3925,
		0.5929,  1.0004,  1.7875,  3.0599,  6.211 ,  7.3459,  7.8112,
		5.9555,  4.06  ,  2.3488,  1.2135,  0.6358,  0.2546,  0.2045,
		0.1868,  0.1993,  0.222 ,  0.2595,  0.3114,  0.373 ,  0.4551,
		0.461 ,  0.4374,  0.3589,  0.319 ,  0.2884,  0.2665,  0.2514,
		0.2375,  0.2341,  0.2313,  0.2203,  0.2151,  0.2079,  0.1996,
		0.1899,  0.1784,  0.1769,  0.181 ,  0.2108,  0.2426,  0.2874,
		0.3378,  0.379 ,  0.3848,  0.349 ,  0.3145,  0.2859,  0.3001,
		0.3586,  0.507 ,  0.8405,  2.9644,  5.0976,  7.3555, 10.0743,
		10.0913,  9.0767,  7.1289,  4.6431,  1.2919,  0.6905,  0.4291,
		0.2734,  0.2667,  0.2853,  0.3169,  0.36  ,  0.4269,  0.4296,
		0.4069,  0.3353,  0.3151,  0.3143,  0.3284,  0.3664,  0.5334,
		0.6701,  0.8124,  0.9877,  0.9636,  0.8806,  0.8016,  0.7583,
		0.8111,  0.9249,  1.1102,  1.601 ,  1.7714,  1.7763,  1.6426,
		1.4111,  0.9989,  0.9136,  0.9148,  1.5958,  2.6463,  4.5494,
		7.2016,  9.8157, 12.7413, 13.0402, 12.6204,  8.7857,  5.6579,
		3.0135,  1.5533,  0.8608,  0.455 ,  0.4405,  0.4897,  0.7492,
		0.9149,  1.0353,  1.0391,  0.9594,  0.6083,  0.4582,  0.3609,
		0.2622,  0.2448,  0.2368,  0.2337,  0.232 ,  0.2354,  0.2452,
		0.2588,  0.321 ,  0.3821,  0.4794,  0.6053,  0.7537,  0.9414,
		0.909 ,  0.829 ,  0.712 ,  0.76  ,  0.9273,  1.3601,  2.2408,
		6.1193,  8.5693, 10.6511, 12.4149, 12.0923, 11.0556,  8.8431,
		6.0272,  1.6712,  0.8583,  0.516 ,  0.3232,  0.317 ,  0.3497,
		0.4055,  0.5011,  0.6953,  0.738 ,  0.7494,  0.6485,  0.5826,
		0.5313,  0.4954,  0.4886,  0.5235,  0.5589,  0.5946,  0.6224,
		0.6075,  0.575 ,  0.533 ,  0.4956,  0.4268,  0.3931,  0.3634,
		0.31  ,  0.2878,  0.2684,  0.256 ,  0.246 ,  0.2327,  0.2256,
		0.2162,  0.1954,  0.1859,  0.1776,  0.1713,  0.1667,  0.1609,
		0.1593,  0.1592,  0.1584,  0.1577,  0.1551,  0.1524,  0.1502,
		0.1504,  0.153 ,  0.1585,  0.1724,  0.1775,  0.1808,  0.1798,
		0.1766,  0.1622,  0.1571,  0.1562,  0.1695,  0.1879,  0.2201,
		0.2696,  0.3456,  0.5506,  0.6691,  0.7672,  0.7858,  0.7026,
		0.5713,  0.4369,  0.3245,  0.1915,  0.1585,  0.139 ,  0.1247,
		0.1242,  0.1275,  0.1364,  0.1558,  0.2737,  0.415 ,  0.6627,
		1.489 ,  1.9068,  2.1524,  2.1027,  1.7434,  0.8649,  0.564 ,
		0.3736,  0.2202,  0.1946,  0.1846,  0.1848,  0.1906,  0.2137,
		0.2235,  0.2298,  0.2181,  0.2066,  0.1944,  0.1835,  0.1764,
		0.1708,  0.17  ,  0.1707,  0.1701,  0.1684,  0.1662,  0.1631,
		0.16  ,  0.1551,  0.1537,  0.1525,  0.1537,  0.1551,  0.1567,
		0.1581,  0.1615,  0.168 ,  0.1712,  0.1756,  0.1784,  0.1786,
		0.1797,  0.1818,  0.1844,  0.1918,  0.1968,  0.2031,  0.2186,
		0.2298,  0.2419,  0.2551,  0.2661,  0.2753,  0.268 ,  0.2534,
		0.2099,  0.1883,  0.1706,  0.1593,  0.153 ,  0.1572,  0.1696,
		0.1928,  0.288 ,  0.3617,  0.4575,  0.5417,  0.5965,  0.5387,
		0.4554,  0.3706,  0.25  ,  0.2187});
	origins={};
	for(int n=0; n<int(Nchain*1.25)/5; ++n)
	{
		int origin=d(gen);
		origins.push_back(origin);

	}*/

	
	origins={50};
	int i=5;
	while(i+5<Nchain)
	{
		interCAR.push_back(&tadConf.at(i));
		i=i+1;
	}


	
}


void MCReplicPoly::TrialMove(double* dE)
{
	MCHeteroPoly::TrialMove(dE);
	

}

void MCReplicPoly::OriginMove()
{
	
	/*
	if(Ntad>=int(.95*Nchain+Nchain))
	{
		
		std::ostringstream streamObj;
		streamObj << originRate;
		std::string strObj = streamObj.str();

		std::ofstream outfile(outputDir+"repltime"+std::to_string(Ndf)+ "_" + strObj+".res", std::ios_base::app | std::ios_base::out);

		outfile << MCsteps << std::endl;


		
		
		exit(0);
		
	}*/

	if ( origins.size() > 0 and MCsteps> (Nrelax)*Ninter )
	{

		auto originsCopy =origins;
		//auto weightsCopy =weights;

		std::vector<int> indexes; //create a indexes vector
		indexes.reserve(originsCopy.size());
		for (int i = 0; i < (int)originsCopy.size(); ++i)
			indexes.push_back(i); //populate
		std::random_shuffle(indexes.begin(), indexes.end()); // randomize
		
		for ( int i=0 ; i < (int)indexes.size(); i++) //for every element in indexes
		{
			MCTad* origin = &tadConf[originsCopy[indexes[i]]]; //select origin taf
			
			//delete passivated origin
			if(origin->status!=0)
			{
				std::vector<int>::iterator itr = std::find(origins.begin(), origins.end(), originsCopy[indexes[i]]);
				
				origins.erase(origins.begin() + std::distance(origins.begin(), itr));
				//weights.erase(weights.begin()+ std::distance(origins.begin(), itr));
				
			}
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			if ( rndReplic < (Ndf- int(double(Nfork)/2 + 0.5))*originRate and origin->status==0)
			{

				Replicate(origin);

				
				std::vector<int>::iterator itr = std::find(origins.begin(), origins.end(), originsCopy[indexes[i]]);

				origins.erase(origins.begin()+std::distance(origins.begin(), itr));
				//weights.erase(weights.begin()+ std::distance(origins.begin(), itr));
				
			}
		}
	}
	MCsteps+=1;
}
void MCReplicPoly::ForkMove()
{
	if ( Nfork > 0 )
	{
		auto activeForksCopy =activeForks;
		for ( int i=0 ; i < (int)activeForksCopy.size(); i++)
		{
			MCTad* fork = activeForks[i];
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			if ( rndReplic < replicRate)
			{
				Replicate(fork);
			}
		}
	}
}
void MCReplicPoly::Replicate(MCTad* tad)
{
	if (tad->isRightEnd() || tad->isLeftEnd())
		return;
	
	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];
		
	double rnd = lat->rngDistrib(lat->rngEngine);
	
	//origin replication
	if ( !tad->isFork() )
	{
		// Can't replicate tad if it's already adjacent to a fork
		if ( nb1->isFork() || nb2->isFork() )
			return;
		
		else
		{
			// Replicate extremities at half the normal rate
			if ( nb1->isLeftEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}

			else
			{
				activeForks.push_back(nb1);
				UpdateReplTable(nb1);//increase energy around new fork
			}
			
			if ( nb2->isRightEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}
			
			else
			{
				activeForks.push_back(nb2);
				UpdateReplTable(nb2);//increase energy around new fork
			}
		}
	}
	
	else //fork displacement
	{
		// Replicating left fork means displacing it to its left neighbor, so we need to check if it's already a fork or chain end
		if ( tad->isLeftFork() )
		{
			if ( nb1->isLeftFork() || nb1->isRightEnd() )
				// Probably should never happen, do nothing
				return;
			
			if ( nb1->isRightFork() || nb1->isLeftEnd() )
			{
				// Merge forks/replicate extremities at half the normal rate
				if ( rnd < 0.5 )
					return;
			}
		
			else
			{
				activeForks.push_back(nb1);
				UpdateReplTable(nb1);//increase energy around new fork
			}
		}
		
		// Same for right forks
		else if ( tad->isRightFork() )
		{
			if ( nb2->isRightFork() || nb2->isLeftEnd() )
				return;
			
			if  ( nb2->isLeftFork() || nb2->isRightEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}
		
			else
			{
				activeForks.push_back(nb2);
				UpdateReplTable(nb2);//increase energy around new fork

			}
		}
		
		// Delete old forks
		auto fork = std::find(activeForks.begin(), activeForks.end(), tad);
		activeForks.erase(fork);
		UpdateReplTable(tad);

		
		if ( nb1->isFork() || nb2->isFork() )
		{
			auto fork2 = std::find(activeForks.begin(), activeForks.end(), nb1->isFork() ? nb1 : nb2);
			activeForks.erase(fork2);
			if(nb1->isFork()){
				UpdateReplTable(nb1);//since it's a fork about to be replicated, energy decreases
			}
			else{
				UpdateReplTable(nb2);//since it's a fork about to be replicated, energy decreases
			}
			
		}
	}
	
	ReplicateTADs(tad);
	ReplicateBonds(tad);

	Update();

	
}

void MCReplicPoly::ReplicateTADs(MCTad* tad)
{
	
	MCTad tadReplic;
	
	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];

	// Replicate left end/fork, if applicable
	if ( nb1->isLeftEnd() || nb1->isRightFork() )
	{

		tadReplic = *nb1;
		tadConf.push_back(tadReplic);
		nb1->SisterID= (int) tadConf.size()-1;
		tadConf.back().SisterID = (int) std::distance(tadConf.data(), nb1);
		
		if(nb1->isCAR)
		{

			tadConf.back().isCAR=true;
			tadConf.back().isChoesin=true;
			nb1->isChoesin=true;
			nb1->choesin_binding_site = &tadConf.back();
			tadConf.back().choesin_binding_site=nb1;
		}


	}
	
	// Replicate TAD
	tadReplic = *tad;
	tadConf.push_back(tadReplic);
	tad->SisterID= (int) tadConf.size()-1;
	tadConf.back().SisterID = (int) std::distance(tadConf.data(), tad);
	
	if(tad->isCAR)
	{

		tadConf.back().isCAR=true;
		tadConf.back().isChoesin=true;
		tad->isChoesin=true;
		tad->choesin_binding_site = &tadConf.back();
		tadConf.back().choesin_binding_site=tad;
	}

	
	// Same for right end/fork
	if ( nb2->isRightEnd() || nb2->isLeftFork() )
	{
		tadReplic = *nb2;
		tadConf.push_back(tadReplic);
		nb2->SisterID= (int) tadConf.size()-1;
		tadConf.back().SisterID = (int) std::distance(tadConf.data(), nb2);

		
		if(nb2->isCAR)
		{

			tadConf.back().isCAR=true;
			tadConf.back().isChoesin=true;
			nb2->isChoesin=true;
			nb2->choesin_binding_site = &tadConf.back();
			tadConf.back().choesin_binding_site=nb2;
		}

	}
}

void MCReplicPoly::ReplicateBonds(MCTad* tad)
{
	MCBond bondReplic1;
	MCBond bondReplic2;

	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];
	
	// Create/modify bonds between relevant neighbors and replicated tad
	MCBond* bond1 = tad->isRightFork() ? tad->bonds[2] : tad->bonds[0];
	MCBond* bond2 = tad->isLeftFork() ? tad->bonds[2] : tad->bonds[1];
	
	bondReplic1.id1 = (nb1->isLeftEnd() || nb1->isRightFork()) ? Ntad : bond1->id1;
	bondReplic1.id2 = (nb1->isLeftEnd() || nb1->isRightFork()) ? Ntad+1 : Ntad;
				
	bondReplic1.dir = bond1->dir;

	bondReplic2.id1 = (nb1->isLeftEnd() || nb1->isRightFork()) ? Ntad+1 : Ntad;
	bondReplic2.id2 = (nb2->isRightEnd() || nb2->isLeftFork()) ? Ntad+1 : bond2->id2;
		
	bondReplic2.dir = bond2->dir;
	
	// For right forks, update bond1 to link left (replicated) neighbor to new tad
	if ( tad->isRightFork() )
	{
		bond1->id2 = bondReplic1.id2;
		
		CreateBond(*bond1);
		UnsetFork(tad);
		
		// Merge forks if necessary
		if ( nb2->isLeftFork() )
		{
			MCBond* bond3 = nb2->bonds[2];
			bond3->id1 = bondReplic2.id2;
			
			CreateBond(*bond3);
			UnsetFork(nb2);
		}
	}
	
	else
		tadTopo.push_back(bondReplic1);
			
	// Same for left forks
	if ( tad->isLeftFork() )
	{
		bond2->id1 = bondReplic2.id1;
		
		CreateBond(*bond2);
		UnsetFork(tad);
		
		if ( nb1->isRightFork() )
		{
			MCBond* bond3 = nb1->bonds[2];
			bond3->id2 = bondReplic1.id1;
			
			CreateBond(*bond3);
			UnsetFork(nb1);
		}
	}
	
	else
		tadTopo.push_back(bondReplic2);
}

void MCReplicPoly::UnsetFork(MCTad* tad)
{
	if ( tad->isFork() )
	{
		tad->bonds[2] = nullptr;
		tad->neighbors[2] = nullptr;
		
		--tad->links;
	}
}

void MCReplicPoly::Update()
{
	// Update bonds
	if ( (int) tadTopo.size() > Nbond )
	{
		for ( auto bond = tadTopo.begin()+Nbond; bond != tadTopo.end(); ++bond )
			CreateBond(*bond);
		
		Nbond = (int) tadTopo.size();
	}
	
	// Update tads
	if ( (int) tadConf.size() > Ntad )
	{
		for ( auto tad = tadConf.begin()+Ntad; tad != tadConf.end(); ++tad )
		{
			if ( tad->type == 1 )
			{
				for ( int v = 0; v < 13; ++v )
				{
					int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
					
					++hetTable[vi];
				}
			}
			
			++lat->bitTable[0][tad->pos];
			//Update the choesin status: if in the original chain is CAR so it is in the new chain.
			if(tadConf[tad->SisterID].isCAR)
			{
				tad->isCAR=true;
				tad->isChoesin=true;
				tadConf[tad->SisterID].isChoesin=true;
				tad->choesin_binding_site=&tadConf[tad->SisterID];
				tadConf[tad->SisterID].choesin_binding_site = &tadConf[tadConf[tad->SisterID].SisterID];
				
			}
		}
		
		Ntad = (int) tadConf.size();
	}
	
	Nfork = (int) activeForks.size();
}
void MCReplicPoly::MoveChoesin(MCTad* tad)
{

	
}
double MCReplicPoly::GetEffectiveEnergy() const
{
	if ( Jf > 0.  )
	{
		if (tadTrial->isFork() and neigh==true)
		{
			
			std::vector<int> indexes; //create a indexes vector
			indexes.reserve(activeForks.size());
			for (int i = 0; i < (int)activeForks.size(); ++i)
				indexes.push_back(i); //populate
			std::random_shuffle(indexes.begin(), indexes.end()); // randomize
			
			double Jf1=0.0;
			double Jf2=0.0;
			
			for ( int i = 0; i < (int) indexes.size(); ++i )
			{
				int forkpos = activeForks[indexes[i]]->pos;
				if(forkpos!=tadUpdater->vo)
				{
					for ( int dir = 0; dir < 3; ++dir )
						if(abs(lat->xyzTable[dir][forkpos]-lat->xyzTable[dir][tadUpdater->vo])>6)
							break;
					
					bool neigh1=false;
					bool neigh2=false;
					for ( int v = 0; v < 13; ++v )
					{
						int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
						int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];
						
						for ( int v1 = 0; v1 < 13; ++v1)
						{
							int vi1 = (v1 == 0) ? vo: lat->bitTable[v1][vo];
							int vi2 = (v1 == 0) ? vn : lat->bitTable[v1][vn];
							
							if(forkpos==vi2 )
								neigh2=true;

							if(forkpos==vi1)
								neigh1=true;
						}
						if(neigh1==neigh2 and neigh2==true)
							break;
					}
					if(neigh2==true)
						Jf2+=Jf;
					if(neigh1==true)
						Jf1+=Jf;
				}
			}
			return 	MCHeteroPoly::GetEffectiveEnergy() -Jf2+Jf1;
		}
		if (tadTrial->isFork() and neigh==false)
		{
			
			std::vector<int> indexes; //create a indexes vector
			indexes.reserve(activeForks.size());
			for (int i = 0; i < (int)activeForks.size(); ++i)
				indexes.push_back(i); //populate
			std::random_shuffle(indexes.begin(), indexes.end()); // randomize
			
			double Jf1=0.0;
			double Jf2=0.0;
			
			for ( int i = 0; i < (int) indexes.size(); ++i )
			{
				int forkpos = activeForks[indexes[i]]->pos;
				if(forkpos!=tadUpdater->vo)
				{
					for ( int dir = 0; dir < 3; ++dir )
						if(abs(lat->xyzTable[dir][forkpos]-lat->xyzTable[dir][tadUpdater->vo])>6)
							break;
					
					bool neigh1=false;
					bool neigh2=false;
					for ( int v = 0; v < 13; ++v )
					{
						int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
						int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];

							
						if(forkpos==vn )
							neigh2=true;
							
						if(forkpos==vo)
							neigh1=true;
					}
					if(neigh1==neigh2 and neigh2==true)
						break;
				if(neigh2==true)
					Jf2+=Jf;
				if(neigh1==true)
					Jf1+=Jf;
				}
			}
			return 	MCHeteroPoly::GetEffectiveEnergy() -Jf2+Jf1;
		}
	}
	if ( Jpair > 0.  )
	{
		if (tadTrial->isChoesin == true and neigh==true)
		{


			double Jbott1=0.0;
			double Jbott2=0.0;
			
			for ( int v = 0; v < 13; ++v )
			{
				int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
				int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];


				
				for ( int v1 = 0; v1 < 13; ++v1)
				{
					int vi1 = (v1 == 0) ? vo: lat->bitTable[v1][vo];
					int vi2 = (v1 == 0) ? vn : lat->bitTable[v1][vn];
				
					if(tadTrial->choesin_binding_site->pos==vi2 ){
						Jbott2=Jpair;
					}
					if(tadTrial->choesin_binding_site->pos==vi1){
						Jbott1=Jpair;
					}
				}
				if(Jbott1==Jbott2 and Jbott2==Jpair)
					break;
			}

			return 	MCHeteroPoly::GetEffectiveEnergy() -Jbott2+Jbott1;
		}
		if (tadTrial->isChoesin == true and neigh==false)
		{
			
			
			double Jbott1=0.0;
			double Jbott2=0.0;
			
			for ( int v = 0; v < 13; ++v )
			{
				int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
				int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];

				if(tadTrial->choesin_binding_site->pos==vn )
					Jbott2=Jpair;
				
				if(tadTrial->choesin_binding_site->pos==vo)
					Jbott1=Jpair;

			
				if(Jbott1==Jbott2 and Jbott2==Jpair)
					break;
			}
			return 	MCHeteroPoly::GetEffectiveEnergy() -Jbott2+Jbott1;
		}
		
	}
	
	return 	MCHeteroPoly::GetEffectiveEnergy();
}

void MCReplicPoly::AcceptMove()
{
	MCHeteroPoly::AcceptMove();

	if ( tadTrial->isFork()) //increase energy at fork site
	{
		std::vector<int> updatedposition1;
		std::vector<int> updatedposition2;

		
		for ( int v = 0; v < 13; ++v )
		{
			int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
			int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];
			
			for ( int v1 = 0; v1 < 13; ++v1)
			{
				int vi1 = (v1 == 0) ? vo: lat->bitTable[v1][vo];
				int vi2 = (v1 == 0) ? vn : lat->bitTable[v1][vn];
				if ( std::find(updatedposition1.begin(), updatedposition1.end(), vi1) != updatedposition1.end() )
				{
					--ReplTable[0][vi1];
					updatedposition1.push_back(vi1);
				}
				if ( std::find(updatedposition2.begin(), updatedposition2.end(), vi2) != updatedposition2.end() )
				{
					++ReplTable[0][vi2];
					updatedposition2.push_back(vi2);
					
				}
			}
		}
	}
	
}

void MCReplicPoly::UpdateReplTable(MCTad* tad)
{
/*
	if(tad->isFork())
	{
		std::vector<int> updatedposition1;
		
		
		for ( int v = 0; v < 13; ++v )
		{
			int vo =(v == 0) ?  tad->pos : lat->bitTable[v][tad->pos];
			
			for ( int v1 = 0; v1 < 13; ++v1)
			{
				int vi1 = (v1 == 0) ? vo: lat->bitTable[v1][vo];
				if ( std::find(updatedposition1.begin(), updatedposition1.end(), vi1) != updatedposition1.end() )
				{
					--ReplTable[0][vi1];
					updatedposition1.push_back(vi1);
				}
			}
		}
	}
	else
	{
		std::vector<int> updatedposition1;

		for ( int v = 0; v < 13; ++v )
		{
			int vo =(v == 0) ?  tad->pos : lat->bitTable[v][tad->pos];
			
			for ( int v1 = 0; v1 < 13; ++v1)
			{
				int vi1 = (v1 == 0) ? vo: lat->bitTable[v1][vo];
				if ( std::find(updatedposition1.begin(), updatedposition1.end(), vi1) != updatedposition1.end() )
				{
					++ReplTable[0][vi1];
					updatedposition1.push_back(vi1);
				}
			}
		}
	}*/
}
std::vector<double3> MCReplicPoly::GetPBCConf()
{
	std::vector<MCTad*> leftEnds;
	std::vector<MCTad*> builtTads;
	
	std::vector<double3> conf(Ntad);
	
	builtTads.reserve(Ntad);
	
	for ( int t = 0; t < Ntad; ++t )
	{
		for ( int i = 0; i < 3; ++i )
			conf[t][i] = lat->xyzTable[i][tadConf[t].pos];
		
		if ( tadConf[t].isLeftEnd() )
			leftEnds.push_back(&tadConf[t]);
	}
	
	// Grow chains recursively, starting from their respective left extremities
	auto leftEnd = leftEnds.begin();
	
	while ( (int) builtTads.size() < Ntad )
	{
		MCTad *tad1, *tad2;
		tad1 = *leftEnd;
		
		bool builtTad1 = (std::find(builtTads.begin(), builtTads.end(), tad1) != builtTads.end());

		if ( !builtTad1 )
		{
			builtTads.push_back(tad1);

			// Traverse main branch
			while ( (tad2 = tad1->neighbors[1]) )
			{
				bool builtTad2 = (std::find(builtTads.begin(), builtTads.end(), tad2) != builtTads.end());

				if ( !builtTad2 )
				{
					BuildPBCPair(builtTads, conf, tad1, tad2);
					
					// Traverse side branches
					if ( tad2->isFork() )
					{
						MCTad *tad3, *tad4;
						tad3 = tad2->neighbors[2];
						
						BuildPBCPair(builtTads, conf, tad2, tad3);
					
						while ( (tad4 = (tad2->isLeftFork() ? tad3->neighbors[1] : tad3->neighbors[0])) )
						{
							BuildPBCPair(builtTads, conf, tad3, tad4);
						
							if ( tad4->isFork() )
								break;
							
							tad3 = tad4;
						}
					}
				}
				
				tad1 = tad2;
			}
		}
		
		++leftEnd;
	}
	
	int chainNum = Ntad / Nchain;
	int chainLength = (chainNum == 1) ? Ntad : Nchain;
	
	std::vector<double3> centers(chainNum);
	
	for ( int c = 0; c < chainNum; ++c )
	{
		auto end1 = conf.begin() + c*chainLength;
		auto end2 = conf.begin() + (c+1)*chainLength;

		centers[c] = GetPBCCenterMass(end1, end2);
	}
	
	centerMass = {0., 0., 0.};
	
	for ( int c = 0; c < chainNum; ++c )
	{
		for ( int i = 0; i < 3; ++i )
			centerMass[i] += centers[c][i] / ((double) chainNum);
	}
	
	return conf;
}

void MCReplicPoly::BuildPBCPair(std::vector<MCTad*>& builtTads, std::vector<double3>& conf, MCTad* tad1, MCTad* tad2)
{
	int id1 = (int) std::distance(tadConf.data(), tad1);
	int id2 = (int) std::distance(tadConf.data(), tad2);
	
	FixPBCPair(conf, id1, id2);
	
	builtTads.push_back(tad2);
}
