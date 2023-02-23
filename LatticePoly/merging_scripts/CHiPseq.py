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
print(CHIP)
np.savetxt(outputDir+"/"+scriptname+".res", CHIP)

