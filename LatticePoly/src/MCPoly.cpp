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
	
	GenerateCAR();
	
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
	std::discrete_distribution<> d({0.00148905, 0.00091369, 0.00134444, 0.00129954, 0.00094887,
		0.0008361 , 0.0008967 , 0.00130561, 0.00134687, 0.00112967,
		0.00146032, 0.00256203, 0.00167084, 0.00140109, 0.0009309 ,
		0.00050712, 0.00042954, 0.00040754, 0.00040194, 0.00059222,
		0.00213679, 0.00339143, 0.00104173, 0.00080058, 0.00051873,
		0.00067132, 0.00177964, 0.00100932, 0.00050621, 0.00077024,
		0.0020318 , 0.00148459, 0.00056911, 0.00043075, 0.00024015,
		0.00041917, 0.00105505, 0.00244438, 0.00140925, 0.0012595 ,
		0.00134421, 0.0006799 , 0.0016666 , 0.00474803, 0.00093349,
		0.00046755, 0.00076747, 0.00098083, 0.00193637, 0.00109016,
		0.00107562, 0.00112664, 0.00122856, 0.00196752, 0.0010149 ,
		0.00070908, 0.00091573, 0.0006187 , 0.00049926, 0.00088729,
		0.00130882, 0.00110027, 0.00039171, 0.00089913, 0.0006431 ,
		0.00112859, 0.00156477, 0.00189972, 0.00119064, 0.0004543 ,
		0.00031346, 0.00079566, 0.00116723, 0.00304486, 0.00150582,
		0.00072873, 0.00234246, 0.00105929, 0.00040538, 0.00067465,
		0.00166147, 0.00155557, 0.00105276, 0.00063703, 0.00035808,
		0.00034157, 0.00042026, 0.00058584, 0.00073896, 0.00099116,
		0.00087908, 0.00097514, 0.00057699, 0.00069328, 0.00171329,
		0.00377608, 0.00097435, 0.00043606, 0.00080812, 0.00113162,
		0.00045502, 0.00040194, 0.00065915, 0.00253599, 0.00174304,
		0.00062123, 0.00043743, 0.00070166, 0.00064165, 0.00078446,
		0.00133536, 0.00075561, 0.0009783 , 0.00091384, 0.00152559,
		0.00500221, 0.00136628, 0.00082234, 0.00076722, 0.00164172,
		0.00086661, 0.0005345 , 0.00041559, 0.00043788, 0.0003577 ,
		0.0006255 , 0.0012654 , 0.00126395, 0.00364564, 0.00097335,
		0.00053908, 0.00071221, 0.00110153, 0.00065358, 0.00087023,
		0.00156884, 0.00087845, 0.00032511, 0.00039435, 0.00055592,
		0.00045755, 0.0004049 , 0.00081245, 0.00122128, 0.00284359,
		0.00165507, 0.00024729, 0.00056977, 0.00104594, 0.00186889,
		0.00198073, 0.00060733, 0.00028749, 0.00046564, 0.00107895,
		0.00286538, 0.00071469, 0.00042469, 0.00030031, 0.00028123,
		0.00041641, 0.00086178, 0.00071982, 0.00270043, 0.00371541,
		0.00066737, 0.00040004, 0.00039314, 0.00024268, 0.00028692,
		0.00035437, 0.00056977, 0.00061593, 0.00050242, 0.00159132,
		0.00220125, 0.00066547, 0.0003585 , 0.00039193, 0.0003978 ,
		0.00058164, 0.00062629, 0.00060933, 0.00052534, 0.00091005,
		0.00169306, 0.00152729, 0.00183951, 0.00242497, 0.00128438,
		0.00048978, 0.00026114, 0.00032841, 0.00033962, 0.00046594,
		0.00111118, 0.00127224, 0.00167448, 0.00060525, 0.00066105,
		0.00069902, 0.00070073, 0.00147245, 0.00316514, 0.00160016,
		0.0012231 , 0.00035492, 0.00033308, 0.0004564 , 0.00097209,
		0.00375545, 0.00142028, 0.00064398, 0.00048475, 0.00061112,
		0.00151737, 0.00142637, 0.00045238, 0.00053158, 0.00088982,
		0.00081841, 0.00041647, 0.00060472, 0.00165641, 0.00331914,
		0.00094631, 0.00045988, 0.00093917, 0.00070845, 0.00030467,
		0.00052833, 0.00078681, 0.00076911, 0.00172513, 0.00265114,
		0.00191236, 0.00087845, 0.00100408, 0.00228586, 0.00084574,
		0.00037557, 0.00040735, 0.00070529, 0.00054735, 0.00095373,
		0.00230418, 0.00041084, 0.00273317, 0.00077492, 0.00055592,
		0.00063901, 0.00142321, 0.00122395, 0.00111544, 0.00070377,
		0.00047226, 0.00042882, 0.00037457, 0.00207551, 0.00081841,
		0.00046387, 0.00032736, 0.00059586, 0.00086586, 0.00035431,
		0.00068187, 0.00170178, 0.00113228, 0.00184982, 0.00175816,
		0.00165009, 0.00076626, 0.00046776, 0.00028024, 0.00048206,
		0.00084392, 0.00045755, 0.00047481, 0.00040889, 0.0005345 ,
		0.00090762, 0.00114451, 0.00076343, 0.00156346, 0.00468294,
		0.00162797, 0.00091598, 0.00158057, 0.00097416, 0.0006278 ,
		0.00034819, 0.00071563, 0.00159925, 0.00092126, 0.00072743,
		0.00041458, 0.00119691, 0.00280294, 0.00468067, 0.001939  ,
		0.00129998, 0.00136705, 0.00162355, 0.00145475, 0.00086454,
		0.00097072, 0.00088856, 0.00085886, 0.00073941, 0.00125826,
		0.0011048 , 0.00085443, 0.00071946, 0.00048293, 0.00046123,
		0.00033747, 0.0007239 , 0.00080261, 0.00068633, 0.00103666,
		0.00557069, 0.00290077, 0.00079629, 0.00084748, 0.00087006,
		0.00087304, 0.00066484, 0.00123614, 0.00138206, 0.00112997,
		0.00252811, 0.00278959, 0.00154581, 0.00080058, 0.00127533,
		0.00099599, 0.00080453, 0.00064165, 0.0007432 , 0.00097893,
		0.00203623, 0.00063324, 0.00093975, 0.00130164, 0.00250505,
		0.00158466, 0.00170191, 0.00259236, 0.00329012, 0.00253357,
		0.00145   , 0.00255634, 0.00661421, 0.00282721, 0.00191489,
		0.00339181, 0.00534553, 0.00371966, 0.00231394, 0.00654868,
		0.00533772, 0.00318091, 0.001336  , 0.00204518, 0.00348365,
		0.00908407, 0.00379246, 0.00231809, 0.00086454, 0.00088037,
		0.00163682, 0.00134396, 0.0006383 , 0.00097767, 0.00094645,
		0.00146745, 0.0020002 , 0.00073688, 0.00131957, 0.0014039 ,
		0.00177838, 0.00361834, 0.00491045, 0.00187634, 0.00097198,
		0.00057085, 0.00080387, 0.0011938 , 0.00076664, 0.00066673,
		0.00138276, 0.00047904, 0.00051701, 0.00175639, 0.00616555,
		0.00163394, 0.00106172, 0.00072876, 0.00073056, 0.00060912,
		0.00049855, 0.00049357, 0.00057291, 0.00058122, 0.00049087,
		0.00038551, 0.00060002, 0.00312765, 0.00287675, 0.00067978,
		0.00120614, 0.00180031, 0.00075711, 0.0008083 , 0.00064824,
		0.00092719, 0.00040332, 0.00068739, 0.00220292, 0.00085082,
		0.00061302, 0.00363593, 0.00406487, 0.00131891, 0.00147187,
		0.00200453, 0.00143761, 0.00151371, 0.00213011, 0.00132484,
		0.00160838, 0.0011325 , 0.00254496, 0.00095009, 0.00050742,
		0.00056878, 0.00079566, 0.00145146, 0.00199262, 0.00090689,
		0.00040269, 0.00043986, 0.00055045, 0.00151917, 0.00314148,
		0.00073878, 0.00073598, 0.00087971, 0.00055856, 0.00044606,
		0.00049363, 0.00085399, 0.00056908, 0.00094038, 0.00150726,
		0.00072804, 0.00039435, 0.0005058 , 0.00065879, 0.0008164 ,
		0.00057347, 0.00056602, 0.00065159, 0.00056483, 0.001069  ,
		0.00181061, 0.00185588, 0.00351694, 0.0022827 , 0.0008611 ,
		0.00069968, 0.00052012, 0.00040056, 0.00029609, 0.00036534,
		0.00051437, 0.0005909 , 0.00068921, 0.00117945, 0.0006206 ,
		0.00049327, 0.00110912, 0.00112555, 0.00090337, 0.001008  ,
		0.00381955, 0.00208868, 0.0005674 , 0.00036907, 0.00056713,
		0.00200453, 0.00181124, 0.0012393 , 0.00094948, 0.00047247,
		0.00041316, 0.00088161, 0.00042758, 0.00100764, 0.00372936,
		0.00213739, 0.000779  , 0.00041647, 0.00030461, 0.00030335,
		0.00035067, 0.00048536, 0.00060101, 0.00072172, 0.00066275,
		0.00081601, 0.00200574, 0.00391926, 0.0019026 , 0.00109142,
		0.00094948, 0.00247532, 0.00106235, 0.00065409, 0.00059704,
		0.00082318, 0.00071224, 0.00032334, 0.00120455, 0.00224781,
		0.00059417, 0.00033914, 0.00050621, 0.00250894, 0.00313905,
		0.00092323, 0.00064019, 0.0005513 , 0.00045916, 0.0011857 ,
		0.0006187 , 0.00059514, 0.00075109, 0.00324836, 0.00197872,
		0.00054932, 0.00112724, 0.00210132, 0.0009583 , 0.0007108 ,
		0.00071818, 0.00082554, 0.00056625, 0.0003716 , 0.0004502 ,
		0.00043561, 0.00104379, 0.00241405, 0.00184497, 0.00036402,
		0.00034759, 0.00047525, 0.00157804, 0.00337869, 0.00073246,
		0.00063571, 0.00053718, 0.00041963, 0.00036655, 0.00047444,
		0.00047452, 0.00054258, 0.00066421, 0.00057151, 0.00151431,
		0.00582125, 0.00189593, 0.00063083, 0.00070034, 0.00082852,
		0.00074123, 0.00070377, 0.00064007, 0.00136886, 0.00079398,
		0.00052163, 0.00041018, 0.00030928, 0.00035264, 0.00045709,
		0.0006364 , 0.00214303, 0.00328817, 0.00190351, 0.00062945,
		0.00045565, 0.00056246, 0.00067594, 0.00055637, 0.00055487,
		0.00047626, 0.00042279, 0.00046513, 0.00264823, 0.00097198,
		0.00038677, 0.00037781, 0.0010326 , 0.0008652 , 0.00027235,
		0.00032247, 0.00027433, 0.00035264, 0.00038056, 0.00058971,
		0.00190685, 0.00265499, 0.00081245, 0.00268221, 0.00047777,
		0.0003082 , 0.00033299, 0.00037987, 0.00054141, 0.0005258 ,
		0.00071937, 0.00066041, 0.00075647, 0.00159827, 0.00357102,
		0.00188365, 0.00081019, 0.00043572, 0.00163808, 0.0025614 ,
		0.00061221, 0.00041923, 0.0006915 , 0.00282246, 0.00125465,
		0.00068847, 0.00159258, 0.00088666, 0.00049264, 0.00037574,
		0.00036264, 0.00028742, 0.00050368, 0.00039172, 0.00062975,
		0.0007517 , 0.00046881, 0.00061609, 0.00202232, 0.0006541 ,
		0.00163687, 0.0005546 , 0.00051443, 0.00060791, 0.00125953,
		0.00088578, 0.00068193, 0.00071814, 0.00064651, 0.00130061,
		0.00379628, 0.00116347, 0.00068375, 0.00190321, 0.00186195,
		0.00076458, 0.00097261, 0.00055108, 0.00030137, 0.00033065,
		0.00065725, 0.00169572, 0.00391016, 0.00060417, 0.00143497,
		0.00094368, 0.00080023, 0.00052419, 0.00029757, 0.00072459,
		0.00129253, 0.00198883, 0.00184158, 0.00061744, 0.00047349,
		0.00056815, 0.0004413 , 0.00036006, 0.00049908, 0.00047019,
		0.00065289, 0.00278716, 0.00254054, 0.00055499, 0.00036138,
		0.00046152, 0.00072873, 0.00120005, 0.00060139, 0.00057131,
		0.00072459, 0.00131009, 0.00295303, 0.00131891, 0.00064189,
		0.00061842, 0.00049475, 0.00095997, 0.00042848, 0.0004803 ,
		0.00057454, 0.00077859, 0.00078507, 0.00071425, 0.00060354,
		0.00070804, 0.00087845, 0.00093659, 0.00184031, 0.00507077,
		0.00186889, 0.00053402, 0.00059343, 0.00063703, 0.00066863,
		0.00035643, 0.00054273, 0.00052454, 0.00103828, 0.00106415,
		0.00073661, 0.00046655, 0.00034642, 0.00037792, 0.00049908,
		0.0007531 , 0.00139672, 0.00338597, 0.000929  , 0.00042535,
		0.0005269 , 0.00049639, 0.001139  , 0.00270526, 0.00051569,
		0.00105793, 0.00128659, 0.00172315, 0.00258162, 0.00197844,
		0.00223325, 0.00057818, 0.00045832, 0.00048657, 0.00049989,
		0.00044537, 0.00095744, 0.00185801, 0.00122856, 0.00078673,
		0.00036336, 0.00050242, 0.00046425, 0.00088414, 0.00067892,
		0.00149083, 0.00231366, 0.00092016, 0.00052492, 0.0007153 ,
		0.00098127, 0.00111511, 0.00130313, 0.00138448, 0.00119823,
		0.00287089, 0.00138719, 0.00093722, 0.00095618, 0.00083547,
		0.00073107, 0.00035006, 0.00058607, 0.00056726, 0.00068981,
		0.00087971, 0.00059027, 0.00060606, 0.00101622, 0.00247861,
		0.00104898, 0.00058032, 0.00060945, 0.00053528, 0.00111025,
		0.00294491, 0.00147488, 0.00059021, 0.00042785, 0.00100039,
		0.00115394, 0.00039498, 0.00032034, 0.00034885, 0.00042331,
		0.00069163, 0.00095699, 0.00107567, 0.00137902, 0.0009448 ,
		0.00149652, 0.00094417, 0.00098864, 0.00142574, 0.00117153,
		0.00061337, 0.00055085, 0.00204562, 0.00305714, 0.00185285,
		0.00117685, 0.0006043 , 0.00043047, 0.00043158, 0.00035188,
		0.00044777, 0.00049921, 0.00069081, 0.00089245, 0.00093553,
		0.00233093, 0.00078673, 0.00065346, 0.00217463, 0.00656385,
		0.00157802, 0.00083421, 0.0006453 , 0.00044744, 0.00040952,
		0.00044112, 0.0005091 , 0.00051316, 0.00045141, 0.00035896,
		0.00057257, 0.00086248, 0.00041641, 0.00141264, 0.00227693,
		0.00034354, 0.00044997, 0.0004822 , 0.00122161, 0.0035783 ,
		0.00201059, 0.00066869, 0.0004348 , 0.00049926, 0.00050052,
		0.0007107 , 0.0017127 , 0.002573  , 0.00093027, 0.00050712,
		0.00076201, 0.0005019 , 0.00036546, 0.00037124, 0.00056119,
		0.00114198});
	std::vector<int> CAR_sample;

	for(int n=0; n< int(Nchain*1.25/5); ++n)
		CAR_sample.push_back(d(gen));
	std::sort (CAR_sample.begin(), CAR_sample.end());
	CAR_sample.erase( unique( CAR_sample.begin(), CAR_sample.end() ), CAR_sample.end() );
	double rnd = lat->rngDistrib(lat->rngEngine);

	int i=0;
	if(rnd>0.5)
		i=1;

	while(i+1< (int)CAR_sample.size())
	{
		CAR.push_back(&tadConf.at(CAR_sample.at(i)));
		CAR.back()->choesin_binding_site = &tadConf.at(CAR_sample.at(i+1));
		CAR.back()->choesin_binding_site->choesin_binding_site=CAR.back();
		i=i+2;
	}


}
void MCPoly::TrialMove(double* dE)
{
	
	double rnd = lat->rngDistrib(lat->rngEngine);
	if(rnd < (double) Ntad/(Ntad+0*activeForks.size())){
		int t = lat->rngEngine() % Ntad;
		tadTrial = &tadConf[t];
		tadUpdater->TrialMove(tadTrial, dE);
		*dE = tadUpdater->legal ? *dE : 0.;
	}else{
		int forkID = lat->rngEngine() % (int) activeForks.size();
		tadTrial = activeForks[forkID];
		tadUpdater->TrialMove(tadTrial, dE);
		*dE = tadUpdater->legal ? *dE : 0.;
	}
	
	
	for ( int i = 0; i < (int) CAR.size(); ++i )
	{
		bool co_lococalized=false;
		if(CAR.at(i)->pos == CAR.at(i)->choesin_binding_site->pos)
			co_lococalized=true;
		for ( int v = 0; v < 12; ++v )
			if(CAR.at(i)->pos == lat->bitTable[v+1][CAR.at(i)->choesin_binding_site->pos])
				co_lococalized=true;

				
			
		CAR.at(i)->isChoesin=true;
		CAR.at(i)->choesin_binding_site->isChoesin=true;
	}
	for ( int i = 0; i < (int) CAR.size(); ++i )
		if(CAR.at(i)->isChoesin)
			CAR.erase(CAR.begin()+i);
	
}

void MCPoly::AcceptMove()
{
	tadUpdater->AcceptMove(tadTrial);
	
	--lat->bitTable[0][tadUpdater->vo];
	++lat->bitTable[0][tadUpdater->vn];
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

	
	types->SetName("TAD type");
	types->SetNumberOfComponents(1);
	
	forks->SetName("Fork type");
	forks->SetNumberOfComponents(1);
	
	status->SetName("Replication status");
	status->SetNumberOfComponents(1);
	
	SisterIDs->SetName("SisterID");
	SisterIDs->SetNumberOfComponents(1);
	
	std::vector<double3> conf = GetPBCConf();

	for ( int t = 0; t < Ntad; ++t )
	{
		int type = tadConf[t].type;
		int state = tadConf[t].status;
		int fork = tadConf[t].isFork() ? (tadConf[t].isLeftFork() ? -1 : 1) : 0;
		int sisterID = tadConf[t].SisterID;
		
		points->InsertNextPoint(conf[t][0], conf[t][1], conf[t][2]);
		
		types->InsertNextValue(type);
		forks->InsertNextValue(fork);
		status->InsertNextValue(state);
		SisterIDs->InsertNextValue(sisterID);
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
void MCPoly::OriginMove()
{}
void MCPoly::ForkMove()
{}
