# -*- coding: utf-8 -*-
##
# MonomerDist.py
# LatticePoly
##
# Based on mtortora's  code.
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import sys
import pandas as pd
import numpy as np
import cooler
from scipy.spatial import cKDTree
from cooler.create import ArrayLoader
from vtkReader import vtkReader

from scipy.spatial.distance import pdist, squareform
class MonomerDmap():
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame,readLiq=False, readPoly=True)
		self.contactFile = os.path.join(self.reader.outputDir,"r_"+str(r)+ "_trans_hic.cool")
		self.timeFile = os.path.join(self.reader.outputDir,"cycles_r_"+str(r)+ "_trans_hic.res")

		if os.path.exists(self.contactFile):
			print("Files %s' already exist - aborting" % (self.contactFile))
			sys.exit()
		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1
		timepoint=0 #find frame where replication reach the desired percentage
		
		fullyreplicated=False
		for i in range(self.reader.N):
			if(self.reader.nTad<2*self.Nchain):
				next(self.reader)
				timepoint+=1
			if(self.reader.nTad==2*self.Nchain):
				fullyreplicated=True
				break
	
	
		
		if(fullyreplicated==False):
			print("Chromosome not fully replicated ")
			sys.exit()
		Nframe=self.reader.N-timepoint
		#compute the hic for the minutes
		self.Compute(Nframe)
		np.savetxt(self.timeFile, [Nframe] )

		self.Print()
			






		
	#NB here is not the finalFrame but the number of iterations
	def Compute(self,finalFrame):
		#self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
		self.contactProb = np.zeros((self.Nchain, self.Nchain), dtype=np.float32)
		self.copy_weight = np.zeros(self.Nchain, dtype=np.float32)


		for i in range(0, finalFrame):
			self.ProcessFrame(i)

			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" %
					  (i+1, finalFrame))
#self.contactProb=np.rint(self.contactProb/finalFrame)

	def ProcessFrame(self, i):
		data = next(self.reader)


		tree1	= cKDTree(data.polyPos[:], boxsize = None)
		pairs = tree1.query_pairs(r = r*0.71) # NN distance FCC lattice 1/np.sqrt(2) = 0.71
		for (i,j) in pairs:
			if((i>=self.Nchain and j<self.Nchain) or (i<self.Nchain and j>=self.Nchain)):
				k=i
				z=j
				if(i>=self.Nchain):
					k=data.SisterID[i]
				if(j>=self.Nchain):
					z=data.SisterID[j]
				if(k!=z):
					self.contactProb[k,z] = self.contactProb[k,z] + 1
					self.contactProb[z,k] = self.contactProb[z,k] + 1
				else:
					self.contactProb[k,z] = self.contactProb[k,z] + 1
					

		for i in range(self.Nchain):
			self.contactProb[i,i] = self.contactProb[i,i] + 1


	def Print(self):
		#np.savetxt(self.contactFile, self.contactProb )
		ser={str(chrom):1531000}
		chromsizes=pd.Series(ser)
		chromsizes=chromsizes.astype('int64')	
		bins = cooler.binnify(chromsizes, 1000)
		pixels = ArrayLoader(bins, self.contactProb, chunksize=10000000)
		cooler.create_cooler(self.contactFile,bins,pixels)
		
		
		

		print("\033[1;32mPrinted avg.contact probability to r'%s'\033[0m" %self.contactFile)


if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("\033[1;31mUsage is %s outputDir initFrame r\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	r=float(sys.argv[3])
	chrom = sys.argv[4]



	monom = MonomerDmap(outputDir, initFrame=initFrame)
