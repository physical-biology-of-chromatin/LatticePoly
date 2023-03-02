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
scriptname = sys.argv[2]





@numba.jit()
def merge_matrices(outputDir):
	matrices=[]
	for folder in os.listdir(outputDir):
		if(folder.endswith('.gz')==False and folder.endswith('.res')==False):
			for file_name in os.listdir(outputDir+'/'+folder):
				check1=0
				if file_name.endswith(scriptname+".res"):
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					before=np.loadtxt(file_path)
					before=np.array(before).T[0]
					check1=1
				if(check1==1):
					for i in range(len(before)):
						if before[i]==872:
							matrices.append(i)
							break
					break


		
		
	rawdata=np.nanmean(matrices,axis=0)


	return rawdata



for folder in os.listdir(outputDir):
	for subfolder in  os.listdir(outputDir):
		np.loadtxt(outputDir+'/'+folder+"/"+subfolder)
			
		


