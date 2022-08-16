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
region = sys.argv[2]




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

for name in matric_names:
	np.savetxt(outputDir+"finalcis.res", merge_matrices(outputDir,timeframe))

	mymatrix = np.loadtxt(outputDir+"/finalcis.res")
	#NB matrix must have raw counts: here I multiply by # trajectories and # timestep
	binsize = 1250
	#open a cooler file of experimental to recover information regarding chromosome sizes
	clr = cooler.Cooler('/Volumes/KESU/ENS/PhD_data/GSM4585143_23C-15min.mcool::/resolutions/200')

	#create a series with the chromosome of interest
	ser={"chrVII":clr.chromsizes.loc[region]}
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
	#create cooler file
	cooler.io.create('/Volumes/KESU/ENS/PhD_data/May22/saner/Jsister10_Jf0/test1250bp.cool', bins, pixels)

	#open the generated cool file
	clr = cooler.Cooler(outputDir+"/"+name[:-3]+".mcool")
	