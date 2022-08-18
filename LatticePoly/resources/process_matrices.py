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

outputDir = sys.argv[1]
chrom = sys.argv[2]




@numba.jit()
def merge_matrices(outputDir,name):
	matrices=[]
	for folder in os.listdir(outputDir):
		if(folder.endswith('.gz')==False):
			print(folder)
			for file_name in os.listdir(outputDir+'/'+folder):
				if file_name==name:
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					matrices.append(np.loadtxt(file_path))
					break;


	rawdata=nansum(cis,axis=0)
	return rawdata

#Define all matrices names
matric_names=[]
for folder in os.listdir(outputDir):
	if(folder.endswith('.gz')==False):
		print(folder)
		for file_name in os.listdir(outputDir+'/'+folder):
			if file_name.endswith('.res'):
				if file_name  in matric_names == False:
					matric_names.append(file_name)
		break

matric_names
bins_dict={}
pixels_dict={}
for name in matric_names:
	np.savetxt(outputDir+"/"+name+"_hic.res", merge_matrices(outputDir,name))

	mymatrix = np.loadtxt(outputDir+"/"+name+"_hic.res")
	#NB matrix must have raw counts: here I multiply by # trajectories and # timestep
	binsize = 1250
	#open a cooler file of experimental to recover information regarding chromosome sizes
	clr = cooler.Cooler('/LatticePoly/LatticePoly/data/GSM4585143_23C-15min.mcool::/resolutions/200')

	#create a series with the chromosome of interest
	ser={str(chrom):clr.chromsizes.loc[str(chrom)]}
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
	#add uniform weights
	bins["weight"]=1/(mymatrix[0][1])
	pixels = ArrayLoader(bins, mymatrix, chunksize=10000000)
	bins_dict[name[:-7]]=bins
	pixels_dict[name[:-7]]=pixels
	#create cooler file

cooler.create_scool(outputDir+"/hic_library.scool",bins_dict,pixel_dict)

