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
	
	//CAR.push_back(&tadConf.at(75));
	//CAR.push_back(&tadConf.at(125));

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
	std::vector<double>  weights = {23.10356699,   4.41884159,   2.38859219,   6.29402075,
		13.18449657,  15.29337363,   2.90914774,   7.21021356,
		7.0547994 ,   6.33119205,   3.04476467,   1.19892222,
		2.97140717,   9.34155845,   4.24487656,   1.59057799,
		1.74052746,   1.88541666,   1.62192617,   9.79761547,
		67.5953341 ,   9.50748353,   2.25669612,   2.78268028,
		11.49572611,  14.93843798,   5.06990178,   1.96351625,
		0.86504861,   6.28251296,   4.75214829,   1.0589514 ,
		2.06814799,   1.87217021,   2.24079808,   4.47946345,
		1.98502456,   2.20816552,   1.98326658,   1.84587102,
		2.21135838,  31.00639736,  69.45616748,  11.2807858 ,
		2.58744168,   1.57214809,   1.90166627,   1.74484711,
		2.47409013,  12.33311118,  25.33533709,  13.78025112,
		3.70300509,   1.64023361,   1.19797335,   1.64867204,
		1.8834736 ,   5.82499593,   8.12079751,   4.33547764,
		2.03481059,   6.95905037,  53.35738615,  12.97181476,
		2.81038774,   1.6239302 ,   1.53380077,   1.12369803,
		1.51429197,   1.84887006,   1.59791418,   7.46521999,
		13.10178808,   1.10340642,   1.41654262,   2.37203825,
		1.20976709,   2.2969785 ,   2.0611629 ,   1.4491716 ,
		4.95012192, 109.61269855,  36.91472138,   3.35898862,
		1.64567551,   1.58053711,   2.16582483,   2.07404451,
		20.43913156,   7.47647476,   2.26987113,   1.064884  ,
		2.62612978,  26.94061108,  40.76906626,   3.04250325,
		1.37077222,   1.96224609,   1.05217486,   0.85208247,
		1.49029467,   1.93515665,   1.66774449,   2.78555277,
		2.31459771,   1.40277609,   1.53223473,   1.42444508,
		1.88413743,   2.49683188,   4.82180403,  89.62487603,
		39.02733715,   3.90228816,   2.10105446,   2.34412691,
		7.40913871,   2.29079459,   1.79160551,   1.33792368,
		1.59832865,   1.91811254,   1.7293222 ,   1.30552643,
		1.48218169,   1.45128982,   0.95837191,   1.30652652,
		1.89201768,   4.58400826,  23.68025456,   8.94519305,
		1.59540779,   0.99481496,   1.5059105 ,   2.85740614,
		33.6885306 ,  58.21112831,   3.78370737,   5.31537054,
		7.66839977,   2.67264219,   2.07199955,   1.28042033,
		0.80107931,   1.33028842,   3.30453239,  18.87600673,
		19.7022825 ,  42.11893029,  26.93440205,   2.53976087,
		1.38150326,   2.9727632 ,   1.9303435 ,   2.12807672,
		1.44244486,   0.81985591,   1.97050433,   3.12841191,
		1.59721344,   1.60423961,   1.9516258 ,   3.86803429,
		8.73154367,   5.71926826,   6.96047055,   5.86855787,
		14.20358488,   3.90608859,   3.21477786,   3.98041087,
		6.52526921,   3.03424381,   2.69967152,   4.69624016,
		75.3111374 ,  32.62631973,   3.00698784,   1.8623005 ,
		1.3957412 ,   5.38321827,   4.88218731,   1.1320326 ,
		2.12567699,  26.78869894,  33.82368618,   2.93216849,
		2.19300317,   1.68281518,   0.99996183,   2.64678757,
		25.81977242,  18.81774635,   2.18486142,   1.28778844,
		1.07665306,   1.53432449,   2.2456184 ,   1.35015984,
		1.58402737,   2.40729549,   3.72932122,   2.45371857,
		3.02750431,   3.02927526,   5.42195606,  27.85522284,
		49.28496651,  55.64616019,   7.57960098,   3.12656627,
		2.4295946 ,  12.78001525,  13.76014902,   2.62033459,
		1.87492925,   2.10423545,   4.13318818,   7.39990586,
		6.75121705,   2.36783506,  38.76237931,   7.83851625,
		1.46116776,   2.09858702,   1.62364396,   1.22300041,
		3.15716083,   5.09151991,   1.91813621,   1.18068015,
		2.04478228,   2.01553824,   2.2324936 ,   2.6424037 ,
		3.22701055,  11.43661966,  48.22936991,   7.38433815,
		18.45746544,  36.18899682,   3.56369098,   1.73158063,
		1.10867631,   0.70219362,   3.70972164,   5.19071635,
		1.66874247,   4.37141323,   3.02064166,   1.77085053,
		2.09953903,   2.5187698 ,   5.79801528,  28.10751898,
		12.15166738,   4.20791143,   1.96416624,   1.42851746,
		2.51780789,  10.54850661,  67.87222495,  48.85494838,
		3.45460669,   2.04852566,   1.55600638,   1.06932659,
		1.43544679,   1.3682233 ,   1.20248578,   1.83370863,
		3.05964829,  15.15116722,   3.94673832,   2.08484047,
		2.71516927,   2.38557247,   2.42587288,  21.66963347,
		51.30386558,   5.75059178,   1.602666  ,   2.38325563,
		2.07924695,   6.72619982,   8.09920057,  36.66576899,
		29.49229921,  10.69770372,   4.08292082,   3.62048659,
		1.55195038,   0.90123146,   0.75946475,   1.19563459,
		1.79125289,   1.98125446,   3.2235712 ,   5.12470146,
		2.46438668,   2.1070747 ,   5.6465685 ,  56.1229338 ,
		36.326588  ,   3.21478436,   1.48095699,   1.09194169,
		1.61089276,   2.73471289,   6.20292149,   8.06973153,
		2.2191198 ,   2.63563743,   4.09425507,  35.83190271,
		39.9513923 ,   5.12438855,   1.44464676,  11.65599453,
		24.55928789,  14.03539732,   5.42281317,   4.25433302,
		13.66043062,  10.78590727,   3.85818881,   1.701393  ,
		1.49189844,   3.33682435,   3.17429512,   2.81685601,
		20.41907966,   7.33462419,   7.04612425,  20.79452842,
		2.81036146,   2.35651194,   6.01897309,  20.19982926,
		17.22724199,   1.9248664 ,   1.07397293,   1.58694177,
		24.24285932,  33.81921036,   3.02588565,   5.55488753,
		9.27697034,   7.53470818,   1.95878675,  21.04772187,
		50.41778056,  16.01466242,   3.14458099,   2.88973331,
		2.30410034,   2.56152943,   1.52835802,   1.59915527,
		1.49956087,   1.27799614,   1.47050057,   2.22134349,
		2.32529344,  12.16202616,   6.29777261,   3.39961281,
		2.0473913 ,   2.15982607,   3.38692751,   3.67854401,
		3.14941662,   3.65320865,  40.20993228,  83.7919956 ,
		7.49569913,   3.3427634 ,   8.78682956,  50.48251858,
		26.79834292,  35.89671345,   9.5928346 ,  11.96330111,
		6.98440304,   6.10109246,   2.1118364 ,   1.80311725,
		2.61489663,   6.16734412,  17.52349444, 187.77803734,
		74.89207117,  34.04575469,  31.72292095,  21.64728359,
		46.43965341, 112.4620328 ,  64.482078  ,  21.13522118,
		27.53687926,  38.67559961,  47.24516716, 140.68573647,
		35.36730211,  16.99069077,  98.08130821,  44.7144162 ,
		8.43144204,   4.78364356,   3.73018791,   3.75332923,
		4.17401997,   4.02350893,   4.74333278,   7.98421094,
		24.949209  ,  20.95658361,   6.97936544,   7.57034421,
		6.94886603,   8.18877665,   4.59423187,   6.2081044 ,
		36.53550906,  73.10835395,  25.4503526 ,  39.73276011,
		15.48324748,   5.72034716,   5.17509854,   4.87920854,
		9.61261366,  17.62969199,  10.99030487,  21.20451329,
		2.4791904 ,   1.86448642,   1.86351084,   3.00317614,
		3.14635534,   4.87930313,  20.96876229,  32.59651764,
		6.24197013,  99.98951011,  15.42618066,   3.73553857,
		6.16760059,  24.04882974,   6.97663377,   5.74719186,
		3.45708515,   4.99230624,  11.74731191,  15.78250876,
		6.68455553,   4.28628942,   3.40422199,   3.05878779,
		2.49926018,   7.56284874,   4.40393657,   1.93856717,
		6.16424096,   2.89581207,   3.99239322,   3.09885239,
		25.88304096,  49.00061872,   4.39556626,   1.17956505,
		1.57589009,   1.77674867,   3.07217925,   6.0277165 ,
		2.3829771 ,   1.1484242 ,   1.99976932,   2.27839534,
		2.60465923,   2.38406038,   4.68951737,  22.16473407,
		5.35222847,   4.24081091,  11.82574273,   2.72681003,
		8.46499759,  94.66152681,  16.84565521,   3.30925742,
		1.76992436,   1.81443154,   1.79083741,   2.00371461,
		3.04120627,   8.3342599 ,   7.98252644,   6.76316724,
		2.7536484 ,   1.08339968,   1.22438932,   3.10044059,
		23.45227761,  41.17936448,   4.51725665,   2.9041668 ,
		14.03746601,   9.9697004 ,   3.23280941,   1.9486528 ,
		1.38602354,   7.9185926 ,  24.57937227,  18.86596957,
		11.10199778,   1.59903803,   3.35568678,   2.08573119,
		1.22895225,   1.50401261,   1.39851016,   1.19188046,
		1.74158084,   1.15583331,   1.91348176,   3.14099753,
		3.32903429,   4.56112991,   1.74904379,   2.58935235,
		5.00704452,   5.87212694,   1.64169325,   2.23654758,
		3.46358132,   2.30519998,   2.09502239,   4.10475303,
		17.73570532,   4.48155376,   1.80519058,   2.44908726,
		5.47322288,  87.93169564,  73.30466968,   4.9436169 ,
		2.43737168,   1.48963637,   1.42857946,   1.54461776,
		1.63718231,   2.1572039 ,   2.41816549,   3.40909564,
		2.23191442,  31.63546922,  29.63870943,   4.48899043,
		10.3858782 ,   4.27655745,   1.59239499,   9.16325669,
		31.47540172,  16.25707481,   7.32851492,   3.67929387,
		6.36875197,  17.42759345,  10.63906203,   2.35645179,
		1.30803702,   1.5262628 ,   1.10815785,   1.2876066 ,
		1.28024051,   1.67384095,   3.2237877 ,  26.04348415,
		20.16047665,   5.6536341 ,  44.35713263,   5.36688596,
		2.74745708,   2.24801715,   2.31339028,   6.94362584,
		27.53602962,  11.73162362,   1.52240853,   2.0252401 ,
		1.75011365,   1.69775446,   2.23638279,   2.47529406,
		3.78618657,  23.51345484, 105.3077878 ,   7.76561762,
		4.66596344,   4.7570494 ,   2.24685481,   3.0662971 ,
		3.63306898,   0.99405503,   1.58684962,   2.60760565,
		2.40607371,   2.65238933,  14.87015418,   7.76815174,
		2.80207689,   1.59692907,   2.0382143 ,   2.20956869,
		2.95788633,   2.31946855,   8.55164236,  29.26457006,
		1.78497907,   3.97292153,   9.96248431,  43.32327692,
		5.58569734,   3.80851767,  35.67623101,  36.21409588,
		2.7132504 ,   2.64040728,   1.7543141 ,   8.58464081,
		24.18630616,   5.62464282,   2.30543397,   7.57832903,
		2.31231343,   2.04349109,   1.24362376,   1.3909571 ,
		1.74182657,   1.85715732,   3.65736705,  28.03277475,
		37.7590417 ,   9.9444713 ,   2.11855473,   1.68453583,
		1.60088695,   9.74732882,   3.77865442,   4.34296871,
		5.66510058,  20.99838602,  21.62335434,   6.81135189,
		4.33905019,   3.66865799,  12.06716548,   1.31777428,
		2.52788496,   4.1804385 ,   3.39839753,   2.05336835,
		2.23488853,   4.63198382,  42.38502289,  14.23647153,
		4.55772472,   1.55551734,   2.16516601,   3.58214363,
		19.57288568,  10.2745162 ,   3.13567955,  18.90887812,
		103.81138229,  18.77208513,   2.69736102,   2.88461151,
		2.49156725,   3.67291799,   2.11399333,   1.92838134,
		2.10409221,   2.51745391,   5.77767987,   3.69342624,
		5.06507795,   4.01567182,   2.54439345,   1.68214458,
		1.52001493,   1.32592274,   1.62759189,  24.72760367,
		37.9748755 ,   4.36726446,   2.7670946 ,   2.18133223,
		4.3120735 ,   4.85501715,   1.94154399,   4.67071017,
		4.16403169,  38.9080908 ,  22.05442775,   3.08528685,
		1.5014462 ,   2.02357453,   2.06714179,   5.88158999,
		4.52572138,   1.66790231,   1.56323224,   1.08938279,
		1.61808266,   5.50012495,   6.88109222,  23.65898934,
		18.2278431 ,   2.86338891,   5.56881344,   2.77020393,
		1.0632503 ,   2.15154856,   3.73592198,   8.51407732,
		20.50732619,   6.3777504 ,  19.79601417,   3.01227858,
		1.51462111,   2.18773713,   2.39031663,   2.0871887 ,
		1.45497491,  22.52489986,  16.26177274,   1.49218824,
		1.82676676,   2.12152201,   2.08285477,   2.08069331,
		3.98274676,   2.85906069,   1.50647997,   4.97282464,
		72.55893647,  28.31026919,  11.53436037,   9.69627097,
		19.73039743,   9.67527331,   2.46602795,   1.29731529,
		1.82755336,   3.22966256,  26.39041754,  69.23372307,
		5.10647427,   4.30211646,   3.06094854,   8.66705071,
		4.62685564,   1.44413421,   1.419146  ,   1.03111892,
		14.92925193,  53.02478825,  48.62482962,   5.05246351,
		4.05097184,   2.53810547,   2.841587  ,   4.64639275,
		10.04050518,   4.01551696,   2.1195186 ,   1.77704202,
		1.64737657,   2.66647369,   7.04006454,   3.56868635,
		2.53836624,   1.63412431,   1.99854857,   1.60450363,
		5.93132258,  11.05363122,  10.45745233,  10.32226126,
		14.57257655,   5.63828371,   1.95483618,   1.29693228,
		2.14999953,   8.22611366,  18.16044354,   6.11653371,
		3.82976375,  21.41911036,  32.03920994,  12.2333136 ,
		21.67000744,   5.62857409,   2.13503586,   1.66981702,
		1.85455015,   1.88026578,  16.5855346 ,  67.57546092,
		6.612948  ,   2.17117698,   2.41967232,   3.53165093,
		1.86486015,   1.12598913,   2.01210148,   4.78871774,
		2.77113909,   0.97371343,   1.21855723,   1.44012813,
		1.69542547,   2.17921847,   1.76214016,   1.99831518,
		2.45408768,  13.96798577,  36.41162306,  41.50025717,
		7.07888348,   2.31798446,   1.40108199,   1.82471505,
		1.35488517,   2.59470823,  19.52930605,   4.7119617 ,
		1.0139345 ,   0.74047809,   1.72180118,   2.27839545,
		4.47338937,  65.29769845,  34.63897855,   2.39012451,
		4.19461788,   2.08109974,   1.40150021,   2.73075782,
		4.10400466,   6.41141422,   4.23247819,   1.63390485,
		2.58421553,   8.2589466 ,  29.18520421,  11.5173331 ,
		2.5588713 ,   3.60129745,   2.71838306,   7.24616705,
		20.00752781,  13.89805044,   3.1603168 ,   6.23976368,
		2.8829051 ,   1.58106828,   1.92071182 };
	
	std::discrete_distribution<> d(weights.begin(), weights.end());
	
	std::vector<int> CAR_sample;
	
	while((int) CAR_sample.size() < 400)
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

		
		
		if(std::find(CAR_sample.begin(),CAR_sample.end(),CAR_sample_copy.at(t2+i)) == CAR_sample.end() or abs(CAR_sample_copy.at(t2)-CAR_sample_copy.at(t2+i)) <= 1 )
		{
			
			if(CAR_sample.at(t)==CAR_sample_copy.back() or CAR_sample.at(t)==CAR_sample_copy.at(0))
				CAR_sample.erase(CAR_sample.begin()+t);
			else if(std::find(CAR_sample.begin(),CAR_sample.end(),CAR_sample_copy.at(t2+(i*-1))) == CAR_sample.end() )
				CAR_sample.erase(CAR_sample.begin()+t);
			else if (abs(CAR_sample_copy.at(t2)-CAR_sample_copy.at(t2-i)) <= 1)
				CAR_sample.erase(CAR_sample.begin()+t);
			
		}else
		{
			if(abs(CAR_sample_copy.at(t2)-CAR_sample_copy.at(t2+i)) > 1)
			   {
				CAR.push_back(&tadConf.at(CAR_sample_copy.at(t2)));
				CAR.back()->choesin_binding_site = &tadConf.at(CAR_sample_copy.at(t2+i));
				CAR.back()->choesin_binding_site->choesin_binding_site=CAR.back();
				CAR_sample.erase(CAR_sample.begin()+t);
				int next_anchor=(int)std::distance(CAR_sample.begin(), std::find(CAR_sample.begin(), CAR_sample.end(),CAR_sample_copy.at(t2+i)));
				CAR_sample.erase(CAR_sample.begin()+ next_anchor);
				std::cout<< "choesin is at " << CAR_sample_copy.at(t2) << std::endl;
				std::cout<< "with anchor at" << CAR_sample_copy.at(t2+i) << std::endl;
			   }
		}
	}
	std::cout<< "Found Nloop= " << 2*CAR.size() << std::endl;

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
		
		
		if (co_lococalized==true)
		{

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

void MCPoly::AcceptMove()
{
	
	tadUpdater->AcceptMove(tadTrial);
	

	for ( int t = 0; t < (int) tadUpdater->reptation_values.size(); ++t )
	{
		--lat->bitTable[0][tadUpdater->reptation_values[t][0]];
		++lat->bitTable[0][tadUpdater->reptation_values[t][1]];
	}

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
