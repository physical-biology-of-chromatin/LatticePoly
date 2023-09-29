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
from cooler.create import ArrayLoader
import cooler
import warnings
warnings.filterwarnings('ignore')
import pandas as pd

outputDir = sys.argv[1]




@numba.jit()
def merge_matrices(outputDir):
	matrices=[]
	for folder in os.listdir(outputDir):
		if(folder.endswith('.gz')==False and folder.endswith('.res')==False):
			#print(folder)
			for file_name in os.listdir(outputDir+'/'+folder):
				check1=0
				check2=0

				if file_name.endswith('duringRepl.res'):
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					before=np.loadtxt(file_path)
					before=np.array(before).T[1]
					check1=1
					break;
			for file_name in os.listdir(outputDir+'/'+folder):
				if file_name.endswith('afterRepl.res'):
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					after=np.loadtxt(file_path)
					check2=1
					break;
			if(check1==check2 and check1==1 and len(before)+len(after)==1999):
				final=np.concatenate((before,after))
				matrices.append(final)




	rawdata=np.nanmean(matrices,axis=0)
	error=np.nanstd(matrices,axis=0)/len(matrices)**0.5

	return [rawdata,error]

merge=merge_matrices(outputDir)
np.savetxt(outputDir+"/merged_Rcm_mixing.res", merge[0])
np.savetxt(outputDir+"/merged_Rcm_mixing_err.res", merge[1])


