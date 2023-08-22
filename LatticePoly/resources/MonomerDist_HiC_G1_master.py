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

import numpy as np

from scipy.spatial import cKDTree

from vtkReader import vtkReader

from scipy.spatial.distance import pdist, squareform


class MonomerDmap():
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame,readLiq=False, readPoly=True)
		self.contactFile = os.path.join(self.reader.outputDir,"r_"+str(r)+ "_G1_hic.res")
		self.timeFile = os.path.join(self.reader.outputDir,"cycles_r_"+str(r)+ "_G1_hic.res")

		if os.path.exists(self.contactFile):
			print("Files %s' already exist - aborting" % (self.contactFile))
			sys.exit()
		self.Nchain=len(self.reader.polyPos)
		self.timepoint=initFrame #find frame where replication reach the desired percentage
		if(self.reader.nTad>self.Nchain):
			print("Chromosome already replicating ")
		self.timepoint=initFrame #find frame where replication starts
		for i in range(self.reader.N):
			if(self.reader.nTad==self.Nchain):
				next(self.reader)
				self.timepoint+=1


		self.reader = vtkReader(outputDir, initFrame,readLiq=False, readPoly=True)
		#restarted vtk reader from middle frame of desired percentage
		#compute the hic for the minutes
		self.Compute(self.timepoint-initFrame)
		np.savetxt(self.timeFile, [self.timepoint-initFrame] )

		self.Print()
			






		
	#NB here is not the finalFrame but the number of iterations
	def Compute(self,finalFrame):
		#self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
		self.contactProb = np.zeros((self.Nchain, self.Nchain), dtype=np.float32)

		

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
			if((i>=self.Nchain or j>=self.Nchain) or (i<self.Nchain or j<self.Nchain)):#always true
				k=i
				z=j
				if(k!=z):
					self.contactProb[k,z] = self.contactProb[k,z] + 1
					self.contactProb[z,k] = self.contactProb[z,k] + 1
				else:
					self.contactProb[k,z] = self.contactProb[k,z] + 1
		for i in range(self.Nchain):
			self.contactProb[i,i] = self.contactProb[i,i] + 1

		for i in range(self.Nchain):
			self.contactProb[i,i] = self.contactProb[i,i] + 1

					



	def Print(self):
		np.savetxt(self.contactFile, self.contactProb )

		print("\033[1;32mPrinted avg.contact probability to '%s'\033[0m" %self.contactFile)


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("\033[1;31mUsage is %s outputDir initFrame r\033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	r=float(sys.argv[3])


monom = MonomerDmap(outputDir, initFrame=initFrame)

