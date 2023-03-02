##
##  PolyMSD.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys
import psutil
import math
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import time
import json
from numpy import nanmean
import numba
from numba import jit, int32
from cooler.create import ArrayLoader
import cooler
import warnings
warnings.filterwarnings('ignore')
import pandas as pd

outputDir = sys.argv[1]






exp_car=np.loadtxt("/Users/mariachiara/Desktop/ENS_internship/LatticePoly/LatticePoly/data/CAR_ch4.in")

for folder in os.listdir(outputDir):
	print(folder)
	if(folder!=".DS_Store"):
		for subfolder in os.listdir(outputDir+"/"+folder):
			if(subfolder!=".DS_Store"):
				print(subfolder)
				data=np.loadtxt(outputDir+"/"+folder+"/"+subfolder+ "/cohesion_pattern_cis1.res")
				CHIP=np.zeros(1531)
				for k in data:
					CHIP[int(k)]+=0.5
				CHIP_ren=np.zeros(len(exp_car))
				for i in range(len(exp_car)):
					if(exp_car[i]!=0):
						CHIP_ren[i]=CHIP[i]+CHIP[i+1]+CHIP[i-1]
				np.savetxt(outputDir+"/"+folder+"/"+subfolder+ "/N"+folder+"_d"+subfolder+".res", CHIP_ren)
	
		
	


