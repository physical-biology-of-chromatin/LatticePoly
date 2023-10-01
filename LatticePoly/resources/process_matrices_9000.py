







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
	first_matrix=0
	for folder in os.listdir(outputDir)[:]:
		if(folder.startswith(".")==False and folder.endswith('.scool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False):
			print(folder)
			file_path = os.path.join(outputDir+'/'+folder, name)
			if(first_matrix==0):
				matrix=pd.read_csv(file_path,sep=" " , header=None).to_numpy(dtype=np.float)
				first_matrix+=1
			else:
				newmatr=pd.read_csv(file_path,sep=" " , header=None).to_numpy(dtype=np.float)
				if(len(newmatr)==9000):
					first_matrix+=1
					matrix=matrix+newmatr
			print("finished_importing")




	return [matrix,first_matrix]

@numba.jit()
def merge_times(outputDir,name):
	matrices=[]
	for folder in os.listdir(outputDir)[:]:
		if(folder.endswith('.scool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False):
					#print("file name = "+file_name)
			file_path = os.path.join(outputDir+'/'+folder, name)
			matrices.append(np.loadtxt(file_path))


	rawdata=np.nansum(matrices)/len(matrices)
	return rawdata

@numba.jit()
def merge_copyweight(outputDir,name):
	matrices=[]
	for folder in os.listdir(outputDir):
		print(folder)
		if(folder.endswith('.scool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False):
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
	if(folder.endswith('.gz')==False and folder.endswith(".res")==False):
		for file_name in os.listdir(outputDir+'/'+folder):
			if file_name.endswith('hic.res')==True and file_name.startswith('cycles')==False and file_name.startswith('copy')==False:
				print("Found matrix name")
				matric_names.append(file_name)
		break

bins_dict={}
bins_dict2={}
pixels_dict={}

for e in range(len(matric_names)):
	traj=96
	if( matric_names[e]!="r_4.0_G1_hic.res" and matric_names[e]!="r_5.0_G1_hic.res" and matric_names[e]!="r_3.0_G1_hic.res"):
		print(matric_names[e])
		merge=merge_matrices(outputDir,matric_names[e])
		rawdata=merge[0]
		traj=merge[1]
		np.savetxt(outputDir+"/"+matric_names[e], rawdata)
	avtime=1201*1000
	#avcopyweight=merge_copyweight(outputDir,"copy_weights_"+matric_names[e])
	mymatrix = np.loadtxt(outputDir+"/"+matric_names[e])
	#NB matrix must have raw counts: here I multiply by # trajectories and # timestep
	binsize = 2500
	#open a cooler file of experimental to recover information regarding chromosome sizes
	#clr = cooler.Cooler('./LatticePoly/LatticePoly/data/GSM4585143_23C-15min.mcool::/resolutions/3200')
	#clr = cooler.Cooler('./GSM4585143_23C-15min.mcool::/resolutions/200')
	ser={"chr1":str(int(binsize*9_000))}
	chromsizes=pd.Series(ser)
	chromsizes=chromsizes.astype('int64')

	#Check that the experimental and simulated chromsizes match:
	#here for exempleI needed to cut last bin (ideally it should not happen)
	if(round(chromsizes[0]/binsize)>len(mymatrix)):
		while(round(chromsizes[0]/binsize)>len(mymatrix)):
			chromsizes[0]=chromsizes-binsize
	else:
		while(round(chromsizes[0]/binsize)<len(mymatrix)):
			chromsizes[0]=chromsizes-binsize

	bins = cooler.binnify(chromsizes, binsize)
	#find ICE bins
	norm=normalization.ICE_normalization(mymatrix)
	weight_ice=[]
	weight_ice.append((norm[0][0]/mymatrix[0][0])**0.5)
	for i in range(1,len(mymatrix[0])):
		if(mymatrix[0][i]!=0):
			weight_ice.append((norm[0][i]/((weight_ice[0]*mymatrix[0][i]))))
		else:
			k=1
			while(k<i):
				if(mymatrix[k][i]==0):
					k+=1
				else:
					break
			weight_ice.append((norm[k][i]/((weight_ice[k]*mymatrix[k][i]))))
			
	#add  weights
	bins["raw"]=1
	#bins["copyweight"]=1/(avcopyweight*(avtime*traj)**0.5)
	bins["weight"]=1/(avtime*traj)**0.5
	bins["ICE"]=np.array(weight_ice)*1/(avtime*traj)**0.5
	#add copy weights
	pixels = ArrayLoader(bins, mymatrix, chunksize=10000000)
	bins_dict[matric_names[e][:-8]]=bins
	pixels_dict[matric_names[e][:-8]]=pixels
	#create cooler file
#print(bins_dict)
cooler.create_scool(outputDir+"/"+final_name+".scool",bins_dict,pixels_dict)
