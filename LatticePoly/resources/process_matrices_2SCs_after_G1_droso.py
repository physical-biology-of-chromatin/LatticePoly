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
	time=0
	for folder in os.listdir(outputDir)[:]:
		if(folder.endswith('.scool')==False and folder.endswith('.cool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False):
			print(folder)
			for file_name in os.listdir(outputDir+'/'+folder):
				check=0
				if file_name==name:
					#1print("file name = "+file_name)
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					clr=cooler.Cooler(file_path)
					matr= (np.array(clr.matrix(balance=False)[:]))
					if(len(matr)>0):
						matrices.append(matr)
						time+=(matr[0][1])
						check=1
					break;
			if(check==0):
				print("missing "+ name+" from "+folder)





	rawdata=np.nansum(matrices,axis=0)

	return [rawdata,len(matrices),(time)/len(matrices)]

@numba.jit()
def merge_times(outputDir,name):
	matrices=[]
	for folder in os.listdir(outputDir):
		if(folder.endswith('.scool')==False and folder.endswith('.gz')==False and  folder.endswith('.cool')==False and folder.endswith('.res')==False):
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
for folder in os.listdir(outputDir)[10:]:
	print(folder)
	if(folder.endswith("scool")==False and folder.endswith(".res")==False and folder.endswith('.gz')==False):
		for file_name in os.listdir(outputDir+'/'+folder):
			print(file_name)
			if file_name.startswith("r") and file_name.endswith('hic.cool')==True and file_name.startswith('cycles')==False and file_name.startswith('copy')==False:
				print("Found matrix name")
				matric_names.append(file_name)
		print(matric_names)
		break

bins_dict={}
bins_dict2={}
pixels_dict={}

for e in range(len(matric_names)):
	print(matric_names[e])
	merge=merge_matrices(outputDir,matric_names[e])
	rawdata=merge[0]
	traj=merge[1]
	print(traj)
	#np.savetxt(outputDir+"/"+matric_names[e][:-5]+".res", rawdata)
	#print(outputDir+"/"+matric_names[e][:-5]+".res")
	avtime=merge[2]
	#avcopyweight=merge_copyweight(outputDir,"copy_weights_"+matric_names[e])
	#mymatrix = np.loadtxt(outputDir+"/"+matric_names[e][:-5]+".res")
	mymatrix = rawdata
	#NB matrix must have raw counts: here I multiply by # trajectories and # timestep
	binsize = 2000
	#open a cooler file of experimental to recover information regarding chromosome sizes
	#clr = cooler.Cooler('./LatticePoly/LatticePoly/data/GSM4585143_23C-15min.mcool::/resolutions/3200')
	#clr = cooler.Cooler('./GSM4585143_23C-15min.mcool::/resolutions/200')
	#create a series with the chromosome of interest
	ser={"2L":(len(mymatrix))*binsize}
	chromsizes=pd.Series(ser)
	chromsizes=chromsizes.astype('int64')

	#Check that the experimental and simulated chromsizes match:
	#here for exempleI needed to cut last bin (ideally it should not happen)
	bins = cooler.binnify(chromsizes, binsize)
	#find ICE bins
			
	#add  weights
	bins["raw"]=1
	#bins["copyweight"]=1/(avcopyweight*(avtime*traj)**0.5)
	bins["weight"]=1/(avtime*traj)**0.5
	#bins["ICE"]=weight_ice*1/(avtime*traj)**0.5
	#add copy weights
	pixels = ArrayLoader(bins, mymatrix, chunksize=10000000)
	bins_dict[matric_names[e][:-9]]=bins
	pixels_dict[matric_names[e][:-9]]=pixels
	#create cooler file
#print(bins_dict)
cooler.create_scool(outputDir+"/"+final_name+".scool",bins_dict,pixels_dict)

