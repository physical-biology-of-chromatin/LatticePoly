##1
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
chrom = sys.argv[2]
final_name = sys.argv[3]




@numba.jit()
def merge_matrices(outputDir,name):
	matrices=[]
	for folder in os.listdir(outputDir):
		if(folder.endswith('.gz')==False and folder.endswith('.scool')==False and folder.endswith('.cool')==False and folder.endswith('.scool')==False and folder.endswith('.res')==False):
			file_path = os.path.join(outputDir+'/'+folder, name)
			matrices.append(file_path)

	merged_clr=cooler.merge_coolers(outputDir+"/"+name,matrices,10000)
	return [len(matrices)]

@numba.jit()
def merge_times(outputDir,name):
	matrices=[]
	for folder in os.listdir(outputDir):
		if(folder.endswith('.gz')==False and folder.endswith('.res')==False):
			#print(folder)
			if(folder.endswith('.gz')==False and folder.endswith('.scool')==False and folder.endswith('.cool')==False and folder.endswith('.scool')==False and folder.endswith('.res')==False):
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
		if(folder.endswith('.gz')==False and folder.endswith('.scool')==False and folder.endswith('.cool')==False and folder.endswith('.scool')==False and folder.endswith('.res')==False):
			#print(folder)
			for file_name in os.listdir(outputDir+'/'+folder):
				if file_name==name:
					#print("file name = "+file_name)
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					matrices.append(np.loadtxt(file_path))
					break;


	rawdata=np.nanmean(matrices,axis=0)
	return rawdata


#Define all matrices names
matric_names=[]
for folder in os.listdir(outputDir):
	if(folder.endswith('.gz')==False and folder.endswith('.scool')==False and folder.endswith('.cool')==False and folder.endswith('.scool')==False and folder.endswith('.res')==False):
		for file_name in os.listdir(outputDir+'/'+folder):
			if file_name.endswith('hic.cool')==True and file_name.startswith('cycles')==False and file_name.startswith('copy')==False:
				matric_names.append(file_name)
		break

bins_dict={}
bins_dict2={}
pixels_dict={}

for e in range(len(matric_names)):
	merge=merge_matrices(outputDir,matric_names[e])
	traj=merge[0]
	clr=cooler.Cooler(outputDir+"/"+matric_names[e])
	#print(np.array(clr.matrix(balance=False)[:]))
	avtime=merge_times(outputDir,"cycles_"+matric_names[e][-5:]+".res")
	#avcopyweight=merge_copyweight(outputDir,"copy_weights_"+matric_names[e][-4:]+".res")
	mymatrix = (np.array(clr.matrix(balance=False)[:]))
	#NB matrix must have raw counts: here I multiply by # trajectories and # timestep
	binsize = 1000
	#open a cooler file of experimental to recover information regarding chromosome sizes
	#clr = cooler.Cooler('./GSM4585143_23C-15min.mcool::/resolutions/200')
	#create a series with the chromosome of interest
	ser={str(chrom):1531000}
	chromsizes=pd.Series(ser)
	chromsizes=chromsizes.astype('int64')

	

	bins = cooler.binnify(chromsizes, binsize)
	
	#add  weights
	bins["raw"]=1
	bins["weight"]=1/(avtime*traj)**0.5
	#add copy weights
	pixels = ArrayLoader(bins, mymatrix, chunksize=10000000)
	bins_dict[matric_names[e][:-9]]=bins
	pixels_dict[matric_names[e][:-9]]=pixels
	#create cooler file
#print(bins_dict)
cooler.create_scool(outputDir+"/"+final_name+".scool",bins_dict,pixels_dict)


	
