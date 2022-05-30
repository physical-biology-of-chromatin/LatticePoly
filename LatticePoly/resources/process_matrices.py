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

from utils import msdFFT
from vtkReader import vtkReader
import networkx as nx
import time
import json
from numpy import nanmean
import numba
from numba import jit, int32


outputDir = sys.argv[1]
timeframe = sys.argv[2]


cis=[]
trans=[]
allint=[]

@numba.jit()
def merge_matrices(outputDir,timeframe):
	for folder in os.listdir(outputDir):
		if(folder.endswith('.gz')==False):
			for file_name in os.listdir(outputDir+'/'+folder):
				print(file_name)
				if file_name.startswith(str(timeframe+"cont":)):
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					cis.append(np.loadtxt(file_path))
					break;


	finalcis=nanmean(cis,axis=0)
	return finalcis
np.savetxt(outputDir+"/"+str(timeframe)+"finalcis.res", merge_matrices(outputDir,timeframe))
