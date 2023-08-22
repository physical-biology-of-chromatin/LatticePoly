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
		self.contactFile = os.path.join(self.reader.outputDir, "r_"+str(r)+"_"+str(round(minutes))+"min_"+str(round(percentage))+"perc_chr1_hic.res")
		self.timeFile = os.path.join(self.reader.outputDir, "cycles_r_"+str(r)+"_"+str(round(minutes))+"min_"+str(round(percentage))+"perc_chr1_hic.res")
			self.copyFile = os.path.join(self.reader.outputDir,"copy_weights_r_"+str(r)+ "_"+str(initFrame)+"_"+str(FinalFrame)+"perc_chr1_hic.res")

		self.finalFrame=initFrame
		frame_minute=round(100_000/Niter)#100_000 cycles in a minute
		if os.path.exists(self.contactFile):
			print("Files %s' already exist - aborting" % (self.contactFile))
			sys.exit()
		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1
		timepoint=initFrame #find frame where replication reach the desired percentage
		while(round(100*(self.reader.nTad-self.Nchain)/self.Nchain)<percentage):
			next(self.reader)
			timepoint+=1
		end_point=0 #find where replication percentage changes
		while(round(100*(self.reader.nTad-self.Nchain)/self.Nchain)==percentage):
			self.reader=next(self.reader)
			end_point+=1
			if(timepoint==initFrame):
				print("percentage of replication not recovered ")
				sys.exit()
			if(round(timepoint+round(end_point/2)-(minutes*frame_minute)/2)<0):
				print("minute interval too wide for time point")
				sys.exit()
		print(timepoint)	
		print(end_point)
		print(frame_minute)
		#restarted vtk reader from middle frame of desired percentage
		self.reader = vtkReader(outputDir, round(timepoint+round(end_point/2)-(minutes*frame_minute)/2) ,readLiq=False, readPoly=True)
		#compute the hic for the minutes
		self.Compute(round(minutes*frame_minute))
		np.savetxt(self.timeFile, [minutes*frame_minute] )

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


		tree1	= cKDTree(data.polyPos[:self.Nchain], boxsize = None)
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
			




	def Print(self):
		np.savetxt(self.contactFile, self.contactProb )
		np.savetxt(self.copyFile, self.copy_weight )

		print("\033[1;32mPrinted avg.contact probability to '%s'\033[0m" %self.contactFile)


if __name__ == "__main__":
	if len(sys.argv) != 7:
		print("\033[1;31mUsage is %s outputDir initFrame percentage minutes Niter r \033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	percentage = int(sys.argv[3])
	minutes = float(sys.argv[4])
	Niter = int(sys.argv[5])
	r = float(sys.argv[6])



	monom = MonomerDmap(outputDir, initFrame=initFrame)
