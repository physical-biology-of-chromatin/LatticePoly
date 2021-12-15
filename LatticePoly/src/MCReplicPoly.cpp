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


MCReplicPoly::MCReplicPoly(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCReplicPoly::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);


	
	MCsteps=0;
	MCrepl=0;
	activeForks.reserve(Nchain);
	neigh=true;

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
	

	//mrt={0.8,0.6,0.4,0.19999999999999996,0.0,0.19999999999999996,0.3999999999999999,0.6000000000000001,0.8};
	//mrt={0.9,0.8,0.7,0.6,0.5,0.4,0.30000000000000004,0.19999999999999996,0.09999999999999998,0.0,0.10000000000000009,0.19999999999999996,0.30000000000000004,0.3999999999999999,0.5,0.600000000000001,0.7,0.8,0.8999999999999999};
	//CAR={5,15,25,35,45,65,75,85,95,105,125,135,155,175,185,195,10,30,50,70,90,110,130,150,170,190,20,40,60,80,100,120,140,160,180};
	
	origins={};
	/*
	for (int i = 1; i < (int)origins.size()-1; ++i)
		tadConf[origins[i]].isCAR=true;
*/
	for (int i = 1; i < (int) Nchain-1; ++i)
		origins.push_back(i);
	
	
	weights={0.00013853, 0.00014187, 0.0001496 , 0.00015388, 0.00016035,
		0.00018686, 0.00021349, 0.00026036, 0.00033521, 0.00044962,
		0.00079277, 0.00097598, 0.00108851, 0.0010098 , 0.000834  ,
		0.00064004, 0.00047082, 0.00034565, 0.00021328, 0.00018394,
		0.00016933, 0.00016411, 0.00016974, 0.00018091, 0.00019522,
		0.00020827, 0.00021881, 0.0002115 , 0.0001949 , 0.00015899,
		0.00014699, 0.00013957, 0.0001355 , 0.00013435, 0.00013801,
		0.00014218, 0.00014699, 0.00015148, 0.00014907, 0.00014417,
		0.00013822, 0.00013415, 0.00013185, 0.0001354 , 0.00014563,
		0.00022006, 0.00032665, 0.00054765, 0.00094497, 0.0015349 ,
		0.00267822, 0.00279598, 0.00241411, 0.00118842, 0.00073963,
		0.00046852, 0.00032143, 0.000247  , 0.00019386, 0.00018916,
		0.00019313, 0.00020952, 0.00021714, 0.00021912, 0.0002115 ,
		0.00019657, 0.00016275, 0.00014949, 0.00014031, 0.00013362,
		0.00013467, 0.0001378 , 0.00014229, 0.00014793, 0.00015826,
		0.00015993, 0.00015837, 0.0001497 , 0.00014594, 0.00014396,
		0.00014438, 0.0001473 , 0.00017131, 0.00020044, 0.00024867,
		0.00047416, 0.00069297, 0.0009647 , 0.00124688, 0.0014422 ,
		0.00132413, 0.00108058, 0.00081803, 0.00044179, 0.0003374 ,
		0.00027518, 0.0002352 , 0.00021255, 0.00019396, 0.00019229,
		0.00019334, 0.00019824, 0.00019939, 0.00019845, 0.00019563,
		0.00019146, 0.00018081, 0.00017632, 0.0001735 , 0.00016964,
		0.00016912, 0.00017016, 0.000171  , 0.00017256, 0.00017653,
		0.00017966, 0.00018217, 0.00018624, 0.00018895, 0.00019188,
		0.0001973 , 0.00020722, 0.00025806, 0.00031997, 0.00043772,
		0.00099164, 0.00149022, 0.00213203, 0.0027893 , 0.00327859,
		0.00313881, 0.00257873, 0.0018936 , 0.0009006 , 0.00063471,
		0.00048032, 0.00038928, 0.00033855, 0.00029585, 0.00028719,
		0.00027873, 0.0002589 , 0.00024846, 0.00023917, 0.00023489,
		0.00023802, 0.00027163, 0.00030963, 0.00036621, 0.00050808,
		0.00056404, 0.000575  , 0.00054587, 0.00048856, 0.00040045,
		0.00038908, 0.00041726, 0.00071541, 0.00110991, 0.00184892,
		0.00296395, 0.00428474, 0.00621404, 0.00625653, 0.00552775,
		0.00279305, 0.00165318, 0.00094779, 0.00057719, 0.00039868,
		0.00028103, 0.00027852, 0.00030055, 0.00040129, 0.00046654,
		0.00053241, 0.00057771, 0.00058951, 0.0005322 , 0.0004918 ,
		0.00045693, 0.00041423, 0.00040933, 0.00040745, 0.00040672,
		0.00040181, 0.00038386, 0.00036788, 0.00034471, 0.00030232,
		0.00028886, 0.00027988, 0.00027299, 0.00026714, 0.00025253,
		0.00024167, 0.00023008, 0.0002139 , 0.00020962, 0.00020733,
		0.00020701, 0.00020941, 0.00021662, 0.00022236, 0.00022664,
		0.00023081, 0.00023008, 0.00022559, 0.00021923, 0.0002163 ,
		0.00022058, 0.00022862, 0.00024428, 0.0002874 , 0.0003064 ,
		0.00031924, 0.00031777, 0.0003017 , 0.00024929, 0.00023395,
		0.00022956, 0.00028625, 0.00039816, 0.00064661, 0.00121399,
		0.00228654, 0.00594721, 0.00755894, 0.00831172, 0.00723376,
		0.00543119, 0.00352381, 0.00195007, 0.00108141, 0.00042728,
		0.00033072, 0.00030306, 0.00033615, 0.0003826 , 0.0004443 ,
		0.0004942 , 0.00052333, 0.00046821, 0.00041131, 0.00035755,
		0.00030138, 0.00030138, 0.00031746, 0.00035473, 0.00040964,
		0.00055162, 0.00060788, 0.00064275, 0.00058575, 0.0005131 ,
		0.00044179, 0.00038438, 0.00034179, 0.0002828 , 0.0002637 ,
		0.00024501, 0.00021912, 0.00021056, 0.00020847, 0.00020868,
		0.00021129, 0.00021693, 0.00021944, 0.00022048, 0.00021818,
		0.00021442, 0.00021098, 0.0002066 , 0.00020545, 0.00021651,
		0.00023269, 0.00025984, 0.00035567, 0.00043146, 0.00052176,
		0.00063315, 0.00073525, 0.00088651, 0.00089371, 0.00085801,
		0.00069683, 0.0006012 , 0.00051529, 0.00043376, 0.00036402,
		0.00027195, 0.00024585, 0.00022935, 0.00022142, 0.00022653,
		0.00023958, 0.00025775, 0.00028197, 0.00033218, 0.0003469 ,
		0.00035055, 0.00033416, 0.00032894, 0.00033218, 0.0003565 ,
		0.0004159 , 0.00079851, 0.00124416, 0.00201991, 0.0045974 ,
		0.00587403, 0.00663266, 0.00652701, 0.0055218 , 0.00243008,
		0.00130941, 0.00067407, 0.00028228, 0.0002353 , 0.00021923,
		0.00023081, 0.00025973, 0.00038511, 0.00046487, 0.00051153,
		0.00044326, 0.00038041, 0.00032884, 0.0003041 , 0.0003161 ,
		0.00059097, 0.00113445, 0.00237997, 0.00784655, 0.01050806,
		0.01201948, 0.01250564, 0.01209767, 0.00814198, 0.00519035,
		0.00283993, 0.00086302, 0.0005844 , 0.00046226, 0.00044242,
		0.00047541, 0.00073514, 0.0009386 , 0.00109854, 0.00112912,
		0.00098057, 0.00079089, 0.00061363, 0.00048366, 0.00032383,
		0.00028499, 0.00026067, 0.00024491, 0.00024804, 0.0002518 ,
		0.00026265, 0.00026683, 0.0002827 , 0.00029324, 0.00030159,
		0.00028708, 0.00026714, 0.00024783, 0.00023489, 0.00023175,
		0.00027351, 0.00034241, 0.00045996, 0.00082388, 0.00097034,
		0.00095124, 0.0008387 , 0.00067031, 0.00045662, 0.00044534,
		0.00048815, 0.00125346, 0.00259638, 0.00507708, 0.0086129 ,
		0.01161161, 0.01414065, 0.01418627, 0.01357171, 0.00871583,
		0.00504535, 0.00244511, 0.00108455, 0.00053303, 0.00020879,
		0.00017162, 0.0001591 , 0.00016682, 0.00018947, 0.00024292,
		0.00035953, 0.00062208, 0.00227077, 0.00408534, 0.00628189,
		0.00961425, 0.01019217, 0.00989475, 0.00832968, 0.00598541,
		0.00177386, 0.00088098, 0.00050589, 0.00030212, 0.00030817,
		0.00036642, 0.00049524, 0.00068117, 0.00112964, 0.00120199,
		0.0011096 , 0.00069735, 0.00054483, 0.00043209, 0.00036131,
		0.00032602, 0.00032237, 0.00033009, 0.00034272, 0.00035421,
		0.00034419, 0.0003159 , 0.00028667, 0.00026172, 0.00022591,
		0.00021484, 0.00020921, 0.00020973, 0.00020889, 0.00020941,
		0.00021161, 0.00021495, 0.00023123, 0.00024042, 0.00024512,
		0.00025107, 0.00026067, 0.00028155, 0.00032205, 0.00039555,
		0.00080039, 0.00126348, 0.00196354, 0.00380974, 0.00454729,
		0.0048637 , 0.0048184 , 0.00434643, 0.00265024, 0.00185696,
		0.00131547, 0.00072627, 0.00059849, 0.00054264, 0.00052562,
		0.0005394 , 0.00061989, 0.00066154, 0.00068347, 0.00067866,
		0.00064609, 0.00060423, 0.00056352, 0.00053742, 0.00051915,
		0.00052134, 0.00053053, 0.0005418 , 0.00053314, 0.00051477,
		0.00049107, 0.00045724, 0.00039638, 0.00037644, 0.00036131,
		0.00034325, 0.00033907, 0.00032968, 0.00032038, 0.00030953,
		0.00030065, 0.00029951, 0.00029543, 0.00027998, 0.0002661 ,
		0.00025086, 0.00023645, 0.00022476, 0.00021651, 0.00021954,
		0.00022653, 0.00024992, 0.00026213, 0.00027111, 0.00027111,
		0.00026547, 0.00023363, 0.00021568, 0.00020305, 0.00019407,
		0.00019647, 0.00020597, 0.00022079, 0.00023864, 0.00026276,
		0.00026276, 0.00025002, 0.00021505, 0.00020524, 0.00020472,
		0.00021808, 0.00025837, 0.00056895, 0.00102859, 0.00194966,
		0.00524109, 0.00680992, 0.00766782, 0.0076721 , 0.00673507,
		0.00328496, 0.00188003, 0.00104342, 0.0004253 , 0.00033406,
		0.00029637, 0.00028886, 0.00030306, 0.00035515, 0.00037456,
		0.00038177, 0.000347  , 0.00031704, 0.0002875 , 0.00026495,
		0.00025357, 0.00025274, 0.00025566, 0.00025796, 0.00024084,
		0.00022131, 0.00020158, 0.00018415, 0.00017152, 0.00016254,
		0.00016724, 0.00017789, 0.0002185 , 0.00024459, 0.00026516,
		0.0002709 , 0.00026401, 0.00024898, 0.00025128, 0.00026673,
		0.00040975, 0.00061895, 0.00104435, 0.00186604, 0.00319434,
		0.00648389, 0.00766866, 0.0081544 , 0.00621717, 0.00423839,
		0.002452  , 0.00126682, 0.00066374, 0.00026579, 0.00021349,
		0.00019501, 0.00020806, 0.00023175, 0.0002709 , 0.00032508,
		0.00038939, 0.0004751 , 0.00048126, 0.00045662, 0.00037467,
		0.00033302, 0.00030107, 0.00027821, 0.00026245, 0.00024794,
		0.00024439, 0.00024146, 0.00022998, 0.00022455, 0.00021703,
		0.00020837, 0.00019824, 0.00018624, 0.00018467, 0.00018895,
		0.00022006, 0.00025326, 0.00030003, 0.00035264, 0.00039565,
		0.00040171, 0.00036433, 0.00032832, 0.00029846, 0.00031329,
		0.00037436, 0.00052928, 0.00087743, 0.00309465, 0.00532157,
		0.00767868, 0.01051694, 0.01053468, 0.00947551, 0.00744212,
		0.00484711, 0.00134866, 0.00072084, 0.00044795, 0.00028541,
		0.00027842, 0.00029784, 0.00033082, 0.00037582, 0.00044566,
		0.00044848, 0.00042478, 0.00035003, 0.00032894, 0.00032811,
		0.00034283, 0.0003825 , 0.00055684, 0.00069954, 0.00084809,
		0.0010311 , 0.00100594, 0.00091929, 0.00083682, 0.00079162,
		0.00084674, 0.00096554, 0.00115898, 0.00167134, 0.00184923,
		0.00185435, 0.00171477, 0.0014731 , 0.00104279, 0.00095374,
		0.00095499, 0.00166592, 0.00276257, 0.00474929, 0.00751802,
		0.01024698, 0.01330112, 0.01361315, 0.01317491, 0.00917172,
		0.00590649, 0.00314591, 0.00162155, 0.00089862, 0.00047499,
		0.00045985, 0.00051122, 0.00078212, 0.0009551 , 0.00108079,
		0.00108476, 0.00100155, 0.00063503, 0.00047833, 0.00037676,
		0.00027372, 0.00025556, 0.0002472 , 0.00024397, 0.00024219,
		0.00024574, 0.00025597, 0.00027017, 0.0003351 , 0.00039889,
		0.00050046, 0.0006319 , 0.00078682, 0.00098276, 0.00094894,
		0.00086542, 0.00074328, 0.00079339, 0.00096804, 0.00141986,
		0.00233925, 0.00638817, 0.00894581, 0.01111908, 0.01296038,
		0.0126236 , 0.01154135, 0.00923164, 0.00629202, 0.00174463,
		0.00089601, 0.00053867, 0.0003374 , 0.00033093, 0.00036506,
		0.00042332, 0.00052312, 0.00072585, 0.00077043, 0.00078233,
		0.00067699, 0.0006082 , 0.00055464, 0.00051717, 0.00051007,
		0.0005465 , 0.00058346, 0.00062073, 0.00064975, 0.00063419,
		0.00060026, 0.00055642, 0.00051738, 0.00044555, 0.00041037,
		0.00037937, 0.00032362, 0.00030045, 0.00028019, 0.00026725,
		0.00025681, 0.00024292, 0.00023551, 0.0002257 , 0.00020399,
		0.00019407, 0.0001854 , 0.00017883, 0.00017402, 0.00016797,
		0.0001663 , 0.00016619, 0.00016536, 0.00016463, 0.00016191,
		0.0001591 , 0.0001568 , 0.00015701, 0.00015972, 0.00016546,
		0.00017997, 0.0001853 , 0.00018874, 0.0001877 , 0.00018436,
		0.00016933, 0.000164  , 0.00016306, 0.00017695, 0.00019616,
		0.00022977, 0.00028145, 0.00036078, 0.00057479, 0.0006985 ,
		0.00080091, 0.00082033, 0.00073347, 0.0005964 , 0.0004561 ,
		0.00033876, 0.00019991, 0.00016546, 0.00014511, 0.00013018,
		0.00012966, 0.0001331 , 0.00014239, 0.00016265, 0.00028573,
		0.00043323, 0.00069182, 0.00155442, 0.00199058, 0.00224697,
		0.00219509, 0.00182   , 0.0009029 , 0.00058878, 0.00039001,
		0.00022988, 0.00020315, 0.00019271, 0.00019292, 0.00019897,
		0.00022309, 0.00023332, 0.0002399 , 0.00022768, 0.00021568,
		0.00020294, 0.00019156, 0.00018415, 0.0001783 , 0.00017747,
		0.0001782 , 0.00017757, 0.0001758 , 0.0001735 , 0.00017027,
		0.00016703, 0.00016191, 0.00016045, 0.0001592 , 0.00016045,
		0.00016191, 0.00016358, 0.00016505, 0.0001686 , 0.00017538,
		0.00017872, 0.00018332, 0.00018624, 0.00018645, 0.0001876 ,
		0.00018979, 0.0001925 , 0.00020023, 0.00020545, 0.00021202,
		0.0002282 , 0.0002399 , 0.00025253, 0.00026631, 0.00027779,
		0.0002874 , 0.00027978, 0.00026453, 0.00021912, 0.00019657,
		0.0001781 , 0.0001663 , 0.00015972, 0.00016411, 0.00017705,
		0.00020127, 0.00030065, 0.00037759, 0.0004776 , 0.0005655 ,
		0.00062271, 0.00056237, 0.00047541, 0.00038688, 0.00026098};

		
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
		auto weightsCopy =weights;

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
				
				origins.erase(origins.begin()+std::distance(origins.begin(), itr));
				weights.erase(weights.begin()+ std::distance(origins.begin(), itr));
				
			}
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			if ( rndReplic < (Ndf- int(double(Nfork)/2 + 0.5))*originRate*weightsCopy[indexes[i]] and origin->status==0)
			{

				Replicate(origin);

				
				
				std::vector<int>::iterator itr = std::find(origins.begin(), origins.end(), originsCopy[indexes[i]]);

				origins.erase(origins.begin()+std::distance(origins.begin(), itr));
				weights.erase(weights.begin()+ std::distance(origins.begin(), itr));
				
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
			if ( rndReplic < replicRate and fork->isFork())
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
