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
def merge_matrices(outputDir):
	matrices=[]
	for folder in os.listdir(outputDir)[:10]:
		if(folder.endswith('.scool')==False  and folder.endswith('.gz')==False and folder.endswith('.res')==False):
			print(folder)
			for file_name in os.listdir(outputDir+'/'+folder):
				check1=0
				file_path = os.path.join(outputDir+'/'+folder, file_name)
				if file_name.endswith(scriptname+".res"):
					matrix=np.loadtxt(file_path)
					if(len(matrix)>1):
						n=np.array(matrix).T[0]
						before=np.array(matrix).T[1]
						check1=1
				if(check1==1):
					stop=0
					while(n[stop]!=1000.0):
						stop+=1
					after=list(before[stop:])
					while(len(after)!=len(before)):
						after.append(np.nan)
					matrices.append(after)
					break


	
	print(np.shape(matrices))
	rawdata=np.nanmean(matrices,axis=0)
	error=np.nanstd(matrices,axis=0)/len(matrices)**0.5


	return [rawdata,error]

merge=merge_matrices(outputDir)
np.savetxt(outputDir+"/"+scriptname+"after_repl.res", merge[0])
np.savetxt(outputDir+"/err_"+scriptname+"after_repl.res", merge[1])

