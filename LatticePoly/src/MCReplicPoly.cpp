//
//  MCReplicPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <random>

#include "MCReplicPoly.hpp"


MCReplicPoly::MCReplicPoly(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCReplicPoly::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);

	std::vector<int> lattice_neigh_load1={0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1, 1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4, 4,  4,  4,  5,  5,  5,  6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9, 10, 10, 11, 12};
	std::vector<int> lattice_neigh_load2={0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  1,  3,  5,  7, 9, 11, 12,  2,  4,  6,  8, 10, 11, 12,  3,  5,  8, 10, 11,  4,  6, 7,  9, 12,  5, 10, 12,  6,  9, 11,  7,  9, 12,  8, 10, 11,  9, 11, 10, 12, 11, 12};
	
	for(int n=0; n< 55 ; ++n)
	{
		lattice_neigh1[n]=lattice_neigh_load1[n];
		lattice_neigh2[n]=lattice_neigh_load2[n];

	}
	
	for ( int vi = 0; vi < Ntot; ++vi )
		ReplTable[0][vi] = 0;

	
	activeForks.reserve(Nchain);
	binded_particles.reserve(Ndf);
	for (int i = 0; i < (int) Ndf; ++i)
		binded_particles.push_back({});

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
	
	
	//chr 7
	origins={1,6,14,20,26,51,55,89,91,94,130,133,137,149,163,185,192,213,215,220,228,230,255,270,282,311,336,365,385,388,407,431,442,454,459,460,465,473,485,497,501,523,527,572,594,612,622,636,64,667,677,690,706,710,733,776,782,799,801,804,829,850,859,866,871};
	
	weights={0.0008441 , 0.00113861, 0.00615295, 0.0011208 , 0.00126901,0.01631902, 0.00450673, 0.00759752, 0.00806823, 0.00269195,0.0199772 , 0.01153813, 0.002372  , 0.00223142, 0.03786356,0.00247823, 0.00170537, 0.00148846, 0.00186694, 0.00142549,0.04605837, 0.04407694, 0.00312641, 0.00132944, 0.00448002,0.0357918 , 0.07323744, 0.00183768, 0.05248039, 0.08644023,0.06210324, 0.00209721, 0.00130972, 0.02321365, 0.01614855,0.01131486, 0.00320274, 0.00343364, 0.00220153, 0.00144076,0.00138033, 0.00346672, 0.04149439, 0.04672691, 0.00159914,0.00221997, 0.06408213, 0.00273266, 0.00133516, 0.08104678,0.00311496, 0.00147574, 0.01425362, 0.07897057, 0.00378222,0.00140005, 0.00499844, 0.01212906, 0.01337517, 0.00358758,0.00097005, 0.00169265, 0.00099994, 0.00379431, 0.00139114};
	
		
	
	//origins={20,40,60,80,100,120,140,160,180};
	//origins={10,30,50,70,90,110,130,150,170,190,20,40,60,80,100,120,140,160,180};
	

	
	
	for (int i = 1; i < (int)origins.size()-1; ++i)
	tadConf[origins[i]].isCAR=true;

	for (int i = 1; i < (int) Nchain-1; ++i)
		origins.push_back(i);
*/

	//chr 7 weights

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
	/*std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<> d({0.3587,  0.3272,  0.2957,  0.2477,  0.2318,  0.2217,  0.2171,
		0.2244,  0.2369,  0.2533,  0.2724,  0.3082,  0.3149,  0.3113,
		0.2709,  0.2445,  0.2207,  0.2007,  0.187 ,  0.1744,  0.1742,
		0.1766,  0.1854,  0.1873,  0.1856,  0.18  ,  0.1708,  0.1528,
		0.1481,  0.1468,  0.162 ,  0.1866,  0.2367,  0.3273,  0.4905,
		1.0969,  1.4782,  1.8098,  1.7697,  1.4424,  1.0602,  0.75  ,
		0.5186,  0.2902,  0.2431,  0.2159,  0.1962,  0.1932,  0.1918,
		0.1912,  0.1898,  0.1842,  0.1805,  0.1754,  0.169 ,  0.1697,
		0.172 ,  0.1767,  0.1822,  0.1915,  0.1955,  0.1999,  0.2107,
		0.2201,  0.2329,  0.248 ,  0.2648,  0.2983,  0.3104,  0.3218,
		0.3271,  0.324 ,  0.3207,  0.313 ,  0.3067,  0.2813,  0.2608,
		0.2402,  0.2018,  0.1882,  0.1804,  0.1768,  0.1763,  0.1862,
		0.1954,  0.2056,  0.2181,  0.2151,  0.2049,  0.1926,  0.1831,
		0.1795,  0.1952,  0.2367,  0.5373,  0.956 ,  1.7927,  3.061 ,
		4.5548,  6.0988,  5.5689,  4.2191,  1.5021,  0.8344,  0.4817,
		0.3315,  0.2579,  0.2408,  0.2701,  0.3208,  0.4819,  0.5531,
		0.5826,  0.5491,  0.473 ,  0.3075,  0.2551,  0.2248,  0.2082,
		0.2167,  0.2341,  0.2626,  0.2946,  0.3605,  0.3757,  0.3769,
		0.3449,  0.328 ,  0.3111,  0.295 ,  0.2831,  0.259 ,  0.2458,
		0.2354,  0.2258,  0.2245,  0.2246,  0.2232,  0.2213,  0.2149,
		0.214 ,  0.2147,  0.2185,  0.2238,  0.2303,  0.2374,  0.2421,
		0.2413,  0.2337,  0.2226,  0.1999,  0.1915,  0.1876,  0.1862,
		0.1888,  0.2053,  0.2174,  0.2276,  0.2329,  0.2274,  0.219 ,
		0.2117,  0.209 ,  0.2287,  0.2649,  0.3394,  0.7104,  1.1035,
		1.6721,  2.3911,  3.097 ,  3.6801,  3.2439,  2.534 ,  1.1359,
		0.7182,  0.4884,  0.359 ,  0.294 ,  0.2554,  0.2653,  0.2889,
		0.3809,  0.4498,  0.5226,  0.5988,  0.6644,  0.7358,  0.7476,
		0.7348,  0.6281,  0.5392,  0.4449,  0.3616,  0.2986,  0.2308,
		0.2215,  0.2264,  0.291 ,  0.3692,  0.4887,  0.6753,  0.9083,
		1.3558,  1.4435,  1.3837,  0.9766,  0.7495,  0.5624,  0.4241,
		0.3312,  0.2396,  0.2184,  0.2068,  0.1984,  0.2001,  0.2022,
		0.2044,  0.2058,  0.2058,  0.2052,  0.2048,  0.2087,  0.2122,
		0.2171,  0.2234,  0.2303,  0.2443,  0.2479,  0.2511,  0.2538,
		0.2524,  0.2517,  0.2475,  0.2458,  0.2482,  0.2535,  0.2604,
		0.2811,  0.2942,  0.3046,  0.3076,  0.306 ,  0.2798,  0.2672,
		0.2637,  0.2857,  0.3212,  0.378 ,  0.4546,  0.5439,  0.6808,
		0.6672,  0.6001,  0.442 ,  0.3984,  0.3954,  0.4427,  0.6017,
		1.8759,  3.677 ,  6.6408, 12.7123, 14.1904, 14.4653, 13.4995,
		11.0302,  4.299 ,  2.1865,  1.122 ,  0.4485,  0.379 ,  0.3697,
		0.4112,  0.4947,  0.7587,  0.8833,  0.9733,  0.9464,  0.8527,
		0.7808,  0.7205,  0.6937,  0.6862,  0.6883,  0.6827,  0.5966,
		0.5243,  0.4437,  0.3697,  0.3121,  0.2469,  0.231 ,  0.2228,
		0.2168,  0.2161,  0.2131,  0.2127,  0.2119,  0.2282,  0.2448,
		0.2622,  0.299 ,  0.3097,  0.3084,  0.2983,  0.281 ,  0.2406,
		0.2264,  0.2208,  0.2298,  0.2532,  0.2873,  0.3293,  0.3742,
		0.4336,  0.4382,  0.4261,  0.4174,  0.4423,  0.5193,  0.6768,
		1.015 ,  2.5986,  4.0373,  5.668 ,  7.8165,  7.556 ,  6.3183,
		4.7074,  3.0562,  1.1152,  0.7459,  0.5566,  0.4267,  0.4278,
		0.445 ,  0.4571,  0.4672,  0.4849,  0.5098,  0.5687,  0.9959,
		1.5855,  2.9075,  5.1431,  8.3491, 13.5536, 14.3958, 13.8208,
		8.9612,  5.6627,  3.0496,  1.5991,  0.9312,  0.4901,  0.4327,
		0.4042,  0.3946,  0.381 ,  0.3728,  0.3713,  0.4001,  0.5727,
		0.8041,  1.31  ,  4.1623,  7.1333, 10.4011, 13.3627, 15.0636,
		14.56  , 12.004 ,  7.8623,  1.7247,  0.7822,  0.4257,  0.2958,
		0.2609,  0.3401,  0.4832,  0.7496,  1.773 ,  2.5614,  3.2611,
		3.6818,  3.7479,  3.0914,  2.4968,  1.9276,  1.0149,  0.7211,
		0.5164,  0.3957,  0.3212,  0.2606,  0.2655,  0.2986,  0.5256,
		0.8042,  1.3129,  2.1483,  3.2053,  5.2559,  5.8811,  5.9342,
		4.305 ,  3.019 ,  1.8983,  1.1041,  0.6655,  0.3161,  0.2548,
		0.2252,  0.2125,  0.2205,  0.2378,  0.2609,  0.2844,  0.3175,
		0.3177,  0.3085,  0.2677,  0.2457,  0.2289,  0.2189,  0.2243,
		0.2761,  0.3404,  0.4507,  0.8456,  1.0931,  1.2745,  1.2946,
		1.1888,  0.8771,  0.7736,  0.7695,  1.3462,  2.4076,  4.6937,
		8.6706, 13.1162, 18.2555, 18.9798, 18.7957, 14.8512, 10.2126,
		5.4164,  2.4309,  1.1177,  0.4201,  0.3651,  0.3752,  0.5838,
		0.7951,  1.0523,  1.2952,  1.4205,  1.2852,  1.098 ,  0.9139,
		0.677 ,  0.6302,  0.6147,  0.6223,  0.6249,  0.6108,  0.5844,
		0.5432,  0.4578,  0.4296,  0.4109,  0.4019,  0.3979,  0.3958,
		0.3872,  0.3728,  0.3402,  0.3319,  0.3216,  0.3114,  0.2951,
		0.2569,  0.2377,  0.2223,  0.204 ,  0.199 ,  0.1957,  0.1921,
		0.188 ,  0.1738,  0.1658,  0.1582,  0.1503,  0.1503,  0.1528,
		0.1599,  0.1719,  0.2111,  0.2356,  0.2565,  0.2802,  0.2774,
		0.273 ,  0.2694,  0.276 ,  0.3303,  0.3936,  0.4968,  0.8675,
		1.1416,  1.4126,  1.5806,  1.5715,  1.0877,  0.7855,  0.5555,
		0.3092,  0.2553,  0.2298,  0.2186,  0.2225,  0.2541,  0.2765,
		0.2982,  0.316 ,  0.3031,  0.2785,  0.2508,  0.2235,  0.1845,
		0.1742,  0.1683,  0.1657,  0.1685,  0.1731,  0.1796,  0.1853,
		0.1926,  0.1936,  0.1912,  0.1853,  0.1831,  0.1815,  0.182 ,
		0.1838,  0.1932,  0.198 ,  0.2001,  0.1897,  0.1769,  0.1658,
		0.1575,  0.1528,  0.1582,  0.1756,  0.2136,  0.436 ,  0.6977,
		1.1009,  1.6302,  2.148 ,  2.3693,  1.9569,  1.3871,  0.5654,
		0.374 ,  0.2688,  0.2133,  0.1861,  0.1728,  0.1777,  0.1884,
		0.2214,  0.2382,  0.2458,  0.2424,  0.2267,  0.1845,  0.1672,
		0.1559,  0.1478,  0.1486,  0.1514,  0.1555,  0.16  ,  0.1659,
		0.1647,  0.1599,  0.147 ,  0.142 ,  0.1393,  0.1396,  0.1424,
		0.1701,  0.2085,  0.2788,  0.6286,  0.9693,  1.3778,  1.7869,
		2.0748,  1.8313,  1.4302,  1.0447,  0.522 ,  0.3884,  0.3098,
		0.2649,  0.2404,  0.2176,  0.2086,  0.2011,  0.1879,  0.1803,
		0.1747,  0.17  ,  0.1675,  0.1668,  0.1689,  0.1701,  0.1669,
		0.1624,  0.1568,  0.1522,  0.1493,  0.1483,  0.1507,  0.1548,
		0.1645,  0.1669,  0.1656,  0.1605,  0.1521,  0.1376,  0.1342,
		0.1338,  0.1507,  0.1784,  0.2359,  0.3489,  0.5695,  1.4755,
		2.0492,  2.5011,  2.1926,  1.6394,  1.0964,  0.6958,  0.4479,
		0.2415,  0.2057,  0.1886,  0.1822,  0.1856,  0.1925,  0.1969,
		0.2022,  0.2   ,  0.1921,  0.1843,  0.1697,  0.1669,  0.1668,
		0.1697,  0.1744,  0.1853,  0.1911,  0.1954,  0.1998,  0.2027,
		0.2064,  0.2097,  0.2154,  0.2287,  0.2362,  0.2431,  0.2557,
		0.2639,  0.2731,  0.2841,  0.2981,  0.3237,  0.3282,  0.3292,
		0.3167,  0.3085,  0.3043,  0.3017,  0.3063,  0.3409,  0.3707,
		0.4085,  0.4978,  0.526 ,  0.5418,  0.5418,  0.523 ,  0.454 ,
		0.424 ,  0.4042,  0.4162,  0.4553,  0.5194,  0.6128,  0.7277,
		0.9545,  0.98  ,  0.9248,  0.6373,  0.5042,  0.4218,  0.3746,
		0.3641,  0.4606,  0.5837,  0.7854,  1.3432,  1.5573,  1.6119,
		1.5329,  1.408 ,  1.2853,  1.4462,  1.9165,  5.1476,  8.529 ,
		12.3523, 15.4793, 17.3727, 18.3919, 17.5626, 15.585 ,  7.4076,
		3.6998,  1.7706,  0.8926,  0.544 ,  0.3757,  0.3951,  0.4591,
		0.7248,  0.8901,  1.0173,  1.0819,  1.042 ,  0.8286,  0.743 ,
		0.6543,  0.601 ,  0.6166,  0.6259,  0.6528,  0.6514,  0.6142,
		0.5751,  0.5291,  0.4557,  0.4272,  0.4002,  0.3805,  0.3603,
		0.3179,  0.2979,  0.2819,  0.2627,  0.2566,  0.2495,  0.2413,
		0.2317,  0.212 ,  0.2052,  0.2015,  0.199 ,  0.1992,  0.2008,
		0.2015,  0.2024,  0.1978,  0.1945,  0.192 ,  0.193 ,  0.1999,
		0.2124,  0.2332,  0.2624,  0.3499,  0.3964,  0.431 ,  0.4094,
		0.3628,  0.3107,  0.26  ,  0.2227,  0.1823,  0.1765,  0.1792,
		0.2072,  0.2335,  0.2619,  0.289 ,  0.3064,  0.3006,  0.2877,
		0.2823,  0.3097,  0.3697,  0.5041,  0.7791,  1.3082,  3.7717,
		5.5501,  7.0528,  7.6548,  6.5813,  5.0148,  3.4117,  2.1269,
		0.8955,  0.6679,  0.5636,  0.5462,  0.5874,  0.6233,  0.6485,
		0.6346,  0.5216,  0.4581,  0.396 ,  0.3073,  0.2848,  0.2685,
		0.2567,  0.2444,  0.2115,  0.1939,  0.1823,  0.1738,  0.1798,
		0.2005,  0.2433,  0.3339,  0.8299,  1.376 ,  2.1323,  3.5387,
		3.6823,  3.1886,  2.4463,  1.6442,  0.6722,  0.4578,  0.3447,
		0.2634,  0.2607,  0.2683,  0.2854,  0.3054,  0.3397,  0.3416,
		0.3256,  0.2695,  0.2442,  0.2224,  0.2088,  0.2009,  0.1984,
		0.2018,  0.2047,  0.2161,  0.2226,  0.2295,  0.2361,  0.2424,
		0.2511,  0.2522,  0.2546,  0.2595,  0.2649,  0.2736,  0.2876,
		0.3084,  0.3708,  0.4079,  0.4435,  0.4954,  0.4994,  0.4877,
		0.4663,  0.4357,  0.3729,  0.348 ,  0.3272,  0.301 ,  0.2916,
		0.2862,  0.2841,  0.2815,  0.2788,  0.2789,  0.2799,  0.2902,
		0.301 ,  0.3136,  0.3287,  0.3508,  0.4143,  0.4533,  0.4956,
		0.5698,  0.5984,  0.622 ,  0.653 ,  0.7053,  0.8391,  0.9508,
		1.1036,  1.5531,  1.815 ,  2.0682,  2.2416,  2.401 ,  2.8082,
		3.2665,  3.9858,  6.9501,  9.1788, 11.4789, 13.4404, 14.9848,
		15.997 , 15.4711, 14.2979,  9.2541,  6.1513,  3.7121,  2.1395,
		1.2835,  0.6279,  0.5314,  0.4893,  0.5279,  0.5869,  0.6571,
		0.7171,  0.7645,  0.7964,  0.7805,  0.7517,  0.6966,  0.6737,
		0.6535,  0.6264,  0.594 ,  0.5136,  0.47  ,  0.4342,  0.3937,
		0.3887,  0.3844,  0.3809,  0.3811,  0.3765,  0.3691,  0.3585,
		0.3355,  0.3245,  0.3134,  0.2982,  0.2741,  0.2195,  0.1965,
		0.1774,  0.159 ,  0.156 ,  0.1557,  0.1584,  0.1641,  0.1718,
		0.1717,  0.1655,  0.1473,  0.1399,  0.1356,  0.134 ,  0.1385,
		0.1768,  0.2269,  0.3305,  0.8194,  1.2613,  1.7643,  2.1535,
		2.1829,  1.4327,  0.9516,  0.5957,  0.2648,  0.2015,  0.171 ,
		0.1558,  0.152 ,  0.16  ,  0.1707,  0.1845,  0.2136,  0.2195,
		0.2141,  0.1978,  0.1822,  0.1598,  0.1578,  0.1617,  0.1957,
		0.234 ,  0.2918,  0.3768,  0.4816,  0.6488,  0.6557,  0.5934,
		0.3988,  0.3142,  0.251 ,  0.2083,  0.1816,  0.157 ,  0.1516,
		0.1494,  0.1457,  0.1426,  0.1414,  0.141 ,  0.1442,  0.1665,
		0.2011,  0.2718,  0.6904,  1.1976,  1.9623,  2.8221,  3.5654,
		3.4737,  2.7073,  1.8849,  0.7586,  0.506 ,  0.37  ,  0.2945,
		0.2514,  0.2187,  0.2124,  0.2103,  0.21  ,  0.2074,  0.2039,
		0.1997,  0.1936,  0.181 ,  0.175 ,  0.1687,  0.1574,  0.1528,
		0.1508,  0.1506,  0.1523,  0.1594,  0.166 ,  0.173 ,  0.1832,
		0.1822,  0.1766,  0.1706,  0.1663,  0.169 ,  0.1857,  0.2257,
		0.4973,  0.8673,  1.5694,  2.6029,  3.8463,  5.3662,  4.9939,
		3.9404,  1.573 ,  0.9129,  0.5321,  0.3559,  0.2618,  0.2017,
		0.1975,  0.2003,  0.2185,  0.227 ,  0.2338,  0.2337,  0.23  ,
		0.2161,  0.2096,  0.2057,  0.2072,  0.2117,  0.2173,  0.2277,
		0.2391,  0.2601,  0.2687,  0.273 ,  0.2596,  0.2435,  0.2235,
		0.2043,  0.1887,  0.1718,  0.1746,  0.1869,  0.2787,  0.3872,
		0.5759,  0.859 ,  1.2553,  1.9808,  2.0017,  1.7694,  0.9184,
		0.5923,  0.3944,  0.2824,  0.2203,  0.1719,  0.1668,  0.169 ,
		0.186 ,  0.196 ,  0.205 ,  0.2069,  0.2027,  0.1856,  0.1781,
		0.1728,  0.1776,  0.1846,  0.1936,  0.203 ,  0.2089,  0.202 ,
		0.1907,  0.1798,  0.1668,  0.1667,  0.1723,  0.1823,  0.1975,
		0.2388,  0.2563,  0.2627,  0.2325,  0.2101,  0.1919,  0.1789,
		0.1742,  0.2081,  0.2688,  0.4059,  1.222 ,  2.0829,  3.1794,
		4.1752,  4.7046,  3.362 ,  2.2426,  1.2946,  0.4363,  0.2981,
		0.23  ,  0.203 ,  0.1932,  0.2096,  0.2296,  0.2481,  0.277 ,
		0.278 ,  0.2688,  0.2564,  0.2417,  0.2269,  0.2258,  0.2268,
		0.2214,  0.2122,  0.198 ,  0.1848,  0.1754,  0.1685,  0.1738,
		0.1895,  0.265 ,  0.3349,  0.4394,  0.5641,  0.6928,  0.804 ,
		0.7617,  0.6676,  0.4285,  0.3265,  0.2546,  0.2085,  0.1812,
		0.1597,  0.1591,  0.1635,  0.1922,  0.2158,  0.2458,  0.2793,
		0.311 ,  0.3414,  0.3346,  0.3145,  0.2585,  0.2361});*/
	origins={};
	for(int n=0; n<int(Nchain*1.25)/5; ++n)
	{

		int origin=d(gen);
		origins.push_back(origin);

	}

	origins={50};

	for (int i=0 ; i < (int) origins.size();++i)
		std::cout <<origins[i]<<  std::endl;

		


	for (int i=0 ; i < (int) origins.size();++i)
		activeOrigins.push_back( &tadConf[origins[i]]);

	
	if ( !RestartFromFile )
	{
		std::ifstream domainFile(CARpath);
		
		if ( !domainFile.good() )
		throw std::runtime_error("MCReplicPoly: Couldn't open file " + CARpath);
		
		std::string line;
		
		while ( std::getline(domainFile, line) )
		{
			std::istringstream ss(line);
			
			int d1;
			
			if ( ss >> d1 )
			{
				if ( (d1 >= 0)  && (d1 <= Nchain)  )
				{
					tadConf[d1].isCAR = true;
					
					
				}
				
				else
					throw std::runtime_error("MCReplicPoly: Found inconsistent domain boundaries '" + line + "' in file " + CARpath);
			}

		}
		
		domainFile.close();
	}

}


void MCReplicPoly::TrialMove(double* dE)
{

	MCHeteroPoly::TrialMove(dE);
	/*for (int i=0 ; i < Ntot;++i)
		if(ReplTable[0][i]>2)
			std::cout <<ReplTable[0][i]<<  std::endl;*/

			
}

void  MCReplicPoly::OriginMove(MCTad* origin_tad)
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
	if(latticeType!="MCLiqLattice")
	{
		if ( (int) activeOrigins.size() > 0 )
		{
			auto originsCopy =activeOrigins;
			std::shuffle (originsCopy.begin(), originsCopy.end(), lat->rngEngine);


			for ( int i=0 ; i < (int)originsCopy.size(); i++) //for every element in indexes
			{
				
				MCTad* origin = originsCopy[i]; //select origin taf
				double rndReplic = lat->rngDistrib(lat->rngEngine);

				int Nocc = activeForks.size() % 2 == 0 ? int(activeForks.size()) : int(activeForks.size())+ 1;
				
				// -1 since origin firing implicate 2 new monomer in the system
				if ( rndReplic < double(2*Ndf- Nocc) * originRate and origin->status==0 and  Ntad < Nchain -1 + int(Nchain * stop_replication))
				{

					Replicate(origin);

				}
			}
		}
	}
	else
	{
			Replicate(origin_tad);
	}
}
void MCReplicPoly::ForkMove()
{
	/*int tot=0;
	for ( int v = 0; v < Ntot ; ++v )
		tot=tot+ReplTable[0][v];
	if(tot!=int(activeForks.size())*55)
		std::cout <<tot<<  std::endl;
	*/
	if ( activeForks.size() > 0   )
	{
		auto activeForksCopy =activeForks;
		for ( int i=0 ; i < (int)activeForksCopy.size(); i++)
		{
			MCTad* fork = activeForks[i];
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			if ( fork->status==0 and rndReplic < replicRate and Ntad < Nchain + int(Nchain*stop_replication) )
				Replicate(fork);
			
			
				
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
				nb1->binding_site = nb2;

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
				nb2->binding_site = nb1;


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
			
			if ( nb1->isRightFork() || nb1->isLeftEnd() ) // SPECIAL CASE
			{
				// Merge forks/replicate extremities at half the normal rate
				if ( rnd < 0.5 )
					return;
				
				//MERGING

			}
		
			else
			{

				activeForks.push_back(nb1); // STANDARD CASE
				UpdateReplTable(nb1);//increase energy around new fork
				nb1->binding_site=tad->binding_site;
				tad->binding_site->binding_site=nb1;

			}
		}
		
		// Same for right forks
		else if ( tad->isRightFork() )
		{
			if ( nb2->isRightFork() || nb2->isLeftEnd() )
				return;
			
			if  ( nb2->isLeftFork() || nb2->isRightEnd() ) // SPECIAL CASE
			{
				if ( rnd < 0.5 )
					return;
				//MERGING

				
			}
		
			else
			{

				activeForks.push_back(nb2);
				UpdateReplTable(nb2);//increase energy around new fork
				nb2->binding_site=tad->binding_site;
				tad->binding_site->binding_site=nb2;

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
			//MODIFY
			//tadConf.back().isChoesin=true;
			//nb1->isChoesin=true;
			//nb1->binding_site = &tadConf.back();
			//tadConf.back().binding_site=nb1;
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
		//MODIFY
		//tadConf.back().isChoesin=true;
		//tad->isChoesin=true;
		//tad->choesin_binding_site = &tadConf.back();
		//tadConf.back().choesin_binding_site=tad;
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
			//MODIFY

			//tadConf.back().isChoesin=true;
			//nb2->isCohesin=true;
			//nb2->binding_site = &tadConf.back();
			//tadConf.back().binding_site=nb2;
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
		
		SetBond(*bond1);
		UnsetFork(tad);
		
		// Merge forks if necessary
		if ( nb2->isLeftFork() )
		{
			MCBond* bond3 = nb2->bonds[2];
			bond3->id1 = bondReplic2.id2;
			
			SetBond(*bond3);
			UnsetFork(nb2);
		}
	}
	
	else
		tadTopo.push_back(bondReplic1);
			
	// Same for left forks
	if ( tad->isLeftFork() )
	{
		bond2->id1 = bondReplic2.id1;
		
		SetBond(*bond2);
		UnsetFork(tad);
		
		if ( nb1->isRightFork() )
		{
			MCBond* bond3 = nb1->bonds[2];
			bond3->id2 = bondReplic1.id1;
			
			SetBond(*bond3);
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
			SetBond(*bond);
		
		Nbond = (int) tadTopo.size();
	}
	
	// Update tads
	if ( (int) tadConf.size() > Ntad )
	{
		for ( auto tad = tadConf.begin()+Ntad; tad != tadConf.end(); ++tad )
		{
			
			if(tadConf.at(tad->SisterID).isCentromere)
				tad->isCentromere=true;
			
			if ( tad->type == 1 )
			{
				for ( int v = 0; v < 13; ++v )
				{
					int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
					
					++hetTable[vi];
				}
			}
			
			++lat->bitTable[0][tad->pos];
		}
		
		Ntad = (int) tadConf.size();
	}
	
	// Update origins
	inactiveOrigins.erase(std::remove_if(inactiveOrigins.begin(), inactiveOrigins.end(), [](const MCTad* tad){return tad->status != 0;}),
						  inactiveOrigins.end());
	activeOrigins.erase(std::remove_if(activeOrigins.begin(), activeOrigins.end(), [](const MCTad* tad){return tad->status != 0;}),
						  activeOrigins.end());

	
	// Update fork/origin counters
	Nfork = (int) activeForks.size();
	Norigin = (int) inactiveOrigins.size();
	
	//check how many forks are binded to their sister 
	NbindedForks = 0;
	for (int i=0; i < (int) activeForks.size();++i)
	{
		if (activeForks.at(i)->binding_site->isFork())
		{
			int pos_binded=activeForks.at(i)->binding_site->pos;
			if ( Jf_sister > 0.  and neigh==1)
			{
				
				for ( int v = 0; v < 55 ; ++v )
				{
					
					int vo =(lattice_neigh1[v] == 0) ? activeForks.at(i)->pos : lat->bitTable[lattice_neigh1[v]][activeForks.at(i)->pos];
					int v1 = (lattice_neigh2[v] == 0) ? vo: lat->bitTable[lattice_neigh2[v]][vo];
					if(v1==pos_binded)
					{
						++NbindedForks;
						break;
					}
				}
			}
			
			if (Jf_sister > 0. and neigh==0 )
			{

				for ( int v = 0; v < 13 ; ++v )
				{
					
					
					int vo =(v == 0) ?  activeForks.at(i)->pos : lat->bitTable[v][activeForks.at(i)->pos];
					if(vo==pos_binded)
					{
						++NbindedForks;
						break;
					}
				}
			}
		}

	}
	
	if(Jlp>0 and latticeType == "MCLiqLattice")
	{

		for ( int i = 0; i < (int) binded_particles.size(); ++i )
			binded_particles.at(i).clear();

		for ( int j = 0; j < (int) activeForks.size(); ++j )
		{
			if(activeForks.at(j)->binding_particle!=-1)
				binded_particles.at(activeForks.at(j)->binding_particle).push_back(j);
		}
		
		for ( int i = 0; i < (int) binded_particles.size(); ++i )
			if(binded_particles.at(i).size()==2)
				std::cout << "binded_particles at  " <<i<<"HAS SIZE "<< binded_particles.at(i).size()  <<"and bindings "<< activeForks.at(binded_particles.at(i).at(0))->binding_particle<<"and "<< activeForks.at(binded_particles.at(i).at(1))->binding_particle<<std::endl;
			else
				std::cout << "binded_particles at  " <<i<<"HAS SIZE "<< binded_particles.at(i).size()<<std::endl;
				

	}
}

double MCReplicPoly::GetEffectiveEnergy() const
{
	if (tadTrial->isFork() )
	{
		
		double Etot = 0.;

		if ( Jf > 0.  )
			Etot=Etot+Jf*(ReplTable[0][tadUpdater->vo]-ReplTable[0][tadUpdater->vn]);


		if ( Jf_sister > 0.  and tadTrial->binding_site->isFork())
		{

			double Jsister_replisome1=0.0;
			double Jsister_replisome2=0.0;

			double old_dist=0.0;
			double new_dist=0.0;
			for ( int dir = 0; dir < 3; ++dir )
			{
				double distance=lat->xyzTable[dir][tadUpdater->vo]-lat->xyzTable[dir][tadTrial->binding_site->pos];
				while ( std::abs(distance) > L/2. )
				{
					double pbcShift = std::copysign(L, distance);
					distance -= pbcShift;
				}
					
				old_dist=old_dist+SQR(distance);
					
				double distance1=lat->xyzTable[dir][tadUpdater->vn]-lat->xyzTable[dir][tadTrial->binding_site->pos];
				while ( std::abs(distance1) > L/2. )
				{
					double pbcShift = std::copysign(L, distance1);
					distance1 -= pbcShift;
				}
				new_dist=new_dist+SQR(distance1);
			}
			double thr_distance = (neigh==1) ? 2 : 0.5;

			Jsister_replisome1= old_dist<=thr_distance ? 1 : old_dist/2;
			Jsister_replisome2= new_dist<=thr_distance ? 1 : new_dist/2;
			
			Etot=Etot-Jf_sister*(Jsister_replisome1-Jsister_replisome2);
			

	
		}

		return MCHeteroPoly::GetEffectiveEnergy()+Etot;
	}

	if (tadTrial->isCohesin )
	{
		
		double Etot = 0.;
		
		
		
		if ( Jpair > 0.  )
		{
			
			double Jpair_anchors1=0.0;
			double Jpair_anchors2=0.0;
			
			double old_dist=0.0;
			double new_dist=0.0;
			for ( int dir = 0; dir < 3; ++dir )
			{
				double distance=lat->xyzTable[dir][tadUpdater->vo]-lat->xyzTable[dir][tadTrial->binding_site->pos];
				while ( std::abs(distance) > L/2. )
				{
					double pbcShift = std::copysign(L, distance);
					distance -= pbcShift;
				}
				
				old_dist=old_dist+SQR(distance);
				
				double distance1=lat->xyzTable[dir][tadUpdater->vn]-lat->xyzTable[dir][tadTrial->binding_site->pos];
				while ( std::abs(distance1) > L/2. )
				{
					double pbcShift = std::copysign(L, distance1);
					distance1 -= pbcShift;
				}
				new_dist=new_dist+SQR(distance1);
			}
			
			double thr_distance = (neigh==1) ? 2 : 0.5;
			
			Jpair_anchors1= old_dist<=thr_distance ? 1 : old_dist/2;
			Jpair_anchors2= new_dist<=thr_distance ? 1 : new_dist/2;
			
			Etot=Etot-Jf_sister*(Jpair_anchors1-Jpair_anchors2);
			
			
			
		}
		return MCHeteroPoly::GetEffectiveEnergy()+Etot;
	}
	
	return 	MCHeteroPoly::GetEffectiveEnergy();
}

void MCReplicPoly::TurnCohesive()
{
	if(!tadTrial->isCohesin and tadTrial->status!=0 and activeForks.size()>0 and std::find(cohesive_CARs.begin(),cohesive_CARs.end(),tadTrial) == cohesive_CARs.end())
	{
		//std::cout <<  "start turnCohesive func"<<  std::endl;

		if(ReplTable[0][tadUpdater->vn]>0)//should also act as weight
		{
			double rnd = lat->rngDistrib(lat->rngEngine);
			if(rnd<keco1*ReplTable[0][tadUpdater->vn])
			{
				cohesive_CARs.push_back(tadTrial);
				//std::cout <<  "car turned cohesive"<<  std::endl;

			}
		}
	}
}

void MCReplicPoly::Find_cohesive_CAR()
{
	if(cohesive_CARs.size()>1 and std::find(cohesive_CARs.begin(),cohesive_CARs.end(),tadTrial) != cohesive_CARs.end())
	{
		//std::cout <<  "start Find_cohesive_CAR func"<<  std::endl;

		auto cohesive_CARs_copy=cohesive_CARs;
		std::shuffle (cohesive_CARs_copy.begin(), cohesive_CARs_copy.end(), lat->rngEngine);
		
		for ( int i = 0; i < (int) cohesive_CARs_copy.size(); ++i )
		{
			if(cohesive_CARs_copy.at(i)!=tadTrial and !cohesive_CARs_copy.at(i)->isCohesin)
			{
				for ( int v = 0; v < 13 ; ++v )
				{
					int pos =(v == 0) ?  tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
					if(cohesive_CARs_copy.at(i)->pos==pos)
					{
						cohesive_CARs_copy.at(i)->isCohesin=true;
						tadTrial->isCohesin=true;
						tadTrial->binding_site=cohesive_CARs_copy.at(i);
						cohesive_CARs_copy.at(i)->binding_site=tadTrial;
						cohesive_CARs.erase(std::remove_if(cohesive_CARs.begin(), cohesive_CARs.end(), [](const MCTad* tad){return tad->isCohesin;}), cohesive_CARs.end());
						NbindedCohesin+=2;
						//std::cout <<  "car found partner "<<  std::endl;

					}
				}
			}
		}
	}
}

void MCReplicPoly::AcceptMove()
{
	MCHeteroPoly::AcceptMove();
	
	if ( tadTrial->isFork()) //increase energy at fork site
	{
			if( neigh==true)
			{
			for ( int v = 0; v < 55 ; ++v )
			{
				int vo =(lattice_neigh1[v] == 0) ? tadUpdater->vo: lat->bitTable[lattice_neigh1[v]][tadUpdater->vo];
				int vn =(lattice_neigh1[v] == 0) ? tadUpdater->vn: lat->bitTable[lattice_neigh1[v]][tadUpdater->vn];
		
				int vi1 = (lattice_neigh2[v] == 0) ? vo: lat->bitTable[lattice_neigh2[v]][vo];
				int vi2 = (lattice_neigh2[v] == 0) ? vn: lat->bitTable[lattice_neigh2[v]][vn];

				--ReplTable[0][vi1];
				++ReplTable[0][vi2];
			}
		}
		else
		{
			for ( int v = 0; v < 13 ; ++v )
			{
				int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
				int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];
				
				--ReplTable[0][vo];
				++ReplTable[0][vn];
			}
		}
	}
	//if CAR is replicated can be turned into cohesive
	if(Jpair>0. and tadTrial->isCAR)
	{
		//if CAR is replicated can be turned into cohesive
		TurnCohesive();
		
		//if CAR check if monomer is in turned cohesive and uncoupled:
		Find_cohesive_CAR();
	}
}

double MCReplicPoly::GetCouplingForkEnergy(const std::vector<int> spinConf) const
{
	if ( Jlp > 0. )
	{
		if ( tadTrial->isFork() )
		{
			double Jlp1=0.0;
			double Jlp2=0.0;
			
			for ( int v = 0; v < 13; ++v )
			{
				int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
				int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];
				
				if(spinConf[tadTrial->binding_particle]==vn )
					Jlp2=Jlp;
				
				if(spinConf[tadTrial->binding_particle]==vo)
					Jlp1=Jlp;
				
				
				if(Jlp1==Jlp2 and Jlp2==Jlp)
					break;
			}
			
			return Jlp1-Jlp2;
		}
	}
	
	return 0.;
}
void MCReplicPoly::UpdateReplTable(MCTad* tad)
{
	
	if(tad->isFork())
	{
		if(neigh==true)
		{
			for ( int v = 0; v < 55 ; ++v )
			{
				
				int vo =(lattice_neigh1[v] == 0) ?  tad->pos : lat->bitTable[lattice_neigh1[v]][tad->pos];
				int vi1 = (lattice_neigh2[v] == 0) ? vo: lat->bitTable[lattice_neigh2[v]][vo];

				--ReplTable[0][vi1];
			}
		}
		else
		{
			for ( int v = 0; v < 13 ; ++v )
			{
				int vo =(v == 0) ?  tad->pos : lat->bitTable[v][tad->pos];
				--ReplTable[0][vo];
			}
		}
	}
	else
	{
		if(neigh==true)
		{
			for ( int v = 0; v < 55 ; ++v )
			{
				
				
				int vo =(lattice_neigh1[v] == 0) ?  tad->pos : lat->bitTable[lattice_neigh1[v]][tad->pos];
				int vi1 = (lattice_neigh2[v] == 0) ? vo: lat->bitTable[lattice_neigh2[v]][vo];

				
				++ReplTable[0][vi1];


			}
		}
		else{
			for ( int v = 0; v < 13 ; ++v )
			{
				int vo =(v == 0) ?  tad->pos : lat->bitTable[v][tad->pos];
				++ReplTable[0][vo];
			}
		}
	}
}


vtkSmartPointer<vtkPolyData> MCReplicPoly::GetVTKData()
{
	vtkSmartPointer<vtkPolyData> polyData = MCHeteroPoly::GetVTKData();
	
	auto forks = vtkSmartPointer<vtkIntArray>::New();
	auto status = vtkSmartPointer<vtkIntArray>::New();
	auto sisterID = vtkSmartPointer<vtkIntArray>::New();
	auto cohesin = vtkSmartPointer<vtkIntArray>::New();
	auto cars = vtkSmartPointer<vtkIntArray>::New();



	
	forks->SetName("Fork type");
	forks->SetNumberOfComponents(1);
	
	status->SetName("Replication status");
	status->SetNumberOfComponents(1);
	
	sisterID->SetName("Sister ID");
	sisterID->SetNumberOfComponents(1);
	
	cohesin->SetName("Cohesin");
	cohesin->SetNumberOfComponents(1);
	
	cars->SetName("CAR");
	cars->SetNumberOfComponents(1);


	

	
	for ( int t = 0; t < Ntad; ++t )
	{
		int fork = tadConf[t].isFork() ? (tadConf[t].isLeftFork() ? -1 : 1) : 0;
		
		forks->InsertNextValue(fork);
		status->InsertNextValue(tadConf[t].status);
		sisterID->InsertNextValue(tadConf[t].SisterID);
		//MODIFY THE OUTPUT
		cohesin->InsertNextValue(tadConf[t].isCohesin);
		cars->InsertNextValue(tadConf[t].isCAR);
	}
	
	polyData->GetPointData()->AddArray(forks);
	polyData->GetPointData()->AddArray(status);
	polyData->GetPointData()->AddArray(sisterID);
	polyData->GetPointData()->AddArray(cohesin);
	polyData->GetPointData()->AddArray(cars);

	


	
	return polyData;
}

void MCReplicPoly::SetVTKData(const vtkSmartPointer<vtkPolyData> polyData)
{
	MCHeteroPoly::SetVTKData(polyData);

	vtkDataArray* status = polyData->GetPointData()->GetArray("Replication status");
	vtkDataArray* sisterID = polyData->GetPointData()->GetArray("Sister ID");



	for ( int t = 0; t < Ntad; ++t )
	{
		tadConf[t].status = (int) status->GetComponent(t, 0);
		tadConf[t].SisterID = (int) sisterID->GetComponent(t, 0);

	}
}
