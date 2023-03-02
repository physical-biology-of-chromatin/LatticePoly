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
scriptname = sys.argv[2]





@numba.jit()
def computeCHIP(outputDir):
	CHIP=np.zeros(1531)
	for folder in os.listdir(outputDir):
		print(folder)
		if(folder.startswith('.')==False and folder.endswith('.scool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False):
			for file_name in os.listdir(outputDir+'/'+folder):
				data=np.loadtxt(outputDir+"/"+folder+ "/cohesion_pattern_cis1.res")
				for k in data:
					CHIP[int(k)]+=0.5


	return CHIP

CHIP=computeCHIP(outputDir)
exp_car=np.loadtxt("/Users/mariachiara/Desktop/ENS_internship/LatticePoly/LatticePoly/data/CAR_ch4.in")
CHIP_ren=np.zeros(len(exp_car))
for i in range(len(exp_car)):
	if(exp_car[i]!=0):
		CHIP_ren[i]=CHIP[i]+CHIP[i+1]+CHIP[i-1]
plt.figure(figsize=(10,2))
plt.plot(exp_car/sum(exp_car),"b")
#plt.plot(CHIP_ren/sum(CHIP_ren))

np.savetxt(outputDir+"/"+scriptname+".res", CHIP_ren)

plt.savefig(outputDir+"/"+scriptname+".png")
plt.figure(figsize=(10,2))
plt.plot(CHIP_ren/sum(CHIP_ren),"y")
#plt.plot(exp_car/sum(exp_car))
plt.savefig(outputDir+"/"+scriptname+"inv.png")


