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
from iced import normalization

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
final_name = sys.argv[2]




@numba.jit()
def merge_matrices(outputDir,name):
	matrices=[]
	for folder in os.listdir(outputDir):
		if(folder.endswith('.err')==False and folder.endswith('.out')==False and folder.endswith('.vtp')==False and folder.endswith('.scool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False):
			#print(folder)
			for file_name in os.listdir(outputDir+'/'+folder):
				check=0
				if file_name==name:
					#print("file name = "+file_name)
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					matrices.append(np.loadtxt(file_path))
					check=1
					break;
			if(check==0):
				print("missing "+ name+" from "+folder)





	rawdata=np.nansum(matrices,axis=0)

	return [rawdata,len(matrices)]

@numba.jit()
def merge_times(outputDir,name):
	matrices=[]
	for folder in os.listdir(outputDir):
		if(folder.endswith('.err')==False and folder.endswith('.out')==False and folder.endswith('.vtp')==False and folder.endswith('.scool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False):
			#print(folder)
			for file_name in os.listdir(outputDir+'/'+folder):
				if file_name==name:
					#print("file name = "+file_name)
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					matrices.append(np.loadtxt(file_path))
					break;


	rawdata=np.nansum(matrices)/len(matrices)
	return rawdata

@numba.jit()
def merge_copyweight(outputDir,name):
	matrices=[]
	for folder in os.listdir(outputDir):
		print(folder)
		if(folder.endswith('.err')==False and folder.endswith('.out')==False and folder.endswith('.vtp')==False and folder.endswith('.scool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False):
			for file_name in os.listdir(outputDir+'/'+folder):
				if file_name==name:
					#print("file name = "+file_name)
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					matrices.append(np.loadtxt(file_path))
					break;


	rawdata=np.nanmean(matrices,axis=0)
	return rawdata


#Define all matrices names


merge=merge_matrices(outputDir,final_name)
rawdata=merge[0]
traj=merge[1]
np.savetxt(outputDir+"/"+final_name, rawdata)
avtime=merge_times(outputDir,"cycles_"+final_name)
avcopyweight=merge_copyweight(outputDir,"copy_weights_"+final_name)
np.savetxt(outputDir+"/"+"copy_weights_"+final_name, avcopyweight)
np.savetxt(outputDir+"/"+"cycles_"+final_name, avtime)
