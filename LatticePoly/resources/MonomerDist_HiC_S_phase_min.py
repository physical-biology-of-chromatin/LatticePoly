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
		self.contactFile = os.path.join(self.reader.outputDir, "r_"+str(r)+"_"+str(minutes)+"min_"+str(round(initFrame))+"endmin_hic.res")
		self.timeFile = os.path.join(self.reader.outputDir, "cycles_r_"+str(r)+"_"+str(minutes)+"min_"+str(round(initFrame))+"endmin_hic.res")
		self.copyFile = os.path.join(self.reader.outputDir,"copy_weights_r_"+str(r)+"_"+str(minutes)+"min_"+str(round(initFrame))+"endmin_hic.res")

		self.finalFrame=initFrame
		self.frame_minute=round(200_000/Niter)#200_000 cycles in a minute
		if os.path.exists(self.contactFile):
			print("Files %s' already exist - aborting" % (self.contactFile))
			sys.exit()
		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1
		
		self.reader = vtkReader(outputDir, initFrame-(minutes*self.frame_minute)/2 ,readLiq=False, readPoly=True)
		#compute the hic for the minutes
		self.Compute(round(minutes*self.frame_minute))
		np.savetxt(self.timeFile, [minutes*self.frame_minute] )

		self.Print()
			






		

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
			

		for i in range(self.Nchain):
			self.contactProb[i,i] = self.contactProb[i,i] + 1
		copy_weight = np.ones(self.Nchain, dtype=np.float32)
		for tad in range(self.Nchain):
			if(data.status[tad]!=0):
				copy_weight[tad]+=1
		self.copy_weight+=copy_weight


	def Print(self):
		np.savetxt(self.contactFile, self.contactProb )
		np.savetxt(self.copyFile, self.copy_weight/(minutes*frame_minute) )
		self.copy_weight = np.zeros(self.Nchain, dtype=np.float32)

		print("\033[1;32mPrinted avg.contact probability to '%s'\033[0m" %self.contactFile)


if __name__ == "__main__":
	if len(sys.argv) != 6:
		print("\033[1;31mUsage is %s outputDir initFrame minutes Niter r \033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	minutes = float(sys.argv[3])
	Niter = int(sys.argv[4])
	r = float(sys.argv[5])



	monom = MonomerDmap(outputDir, initFrame=initFrame)
