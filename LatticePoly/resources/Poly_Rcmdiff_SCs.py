##  LatticePoly
##
##  Created by ddasaro on the model of mtortora script.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##


import os
import sys

import numpy as np

from vtkReader import vtkReader
import time

class Poly_diff_Rcms():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
					
		self.diffRcmFile = os.path.join(self.reader.outputDir, str(time.time())+"diff_rcm.res")

		self.Nchain=0



	def Compute(self):
		self.diff=[]
		self.n_mon=[]
		for i in range(self.reader.N):
			self.ProcessFrame(i)
			
			#if (i+1) % 10 == 0:
			#	print("Processed %d out of %d configurations" % (i+1, self.reader.N))
			
				
	def ProcessFrame(self, i):
		data = next(self.reader)
		
		if(i==0):
			for t in range(self.reader.nTad):
				if(data.status[t]==-1 or data.status[t]==0):
					self.Nchain+=1
		if(data.nTad<2*self.Nchain or 0==0):
			rcm1=[]
			rcm2=[]

			for n in range(data.nTad):
				if(data.status[n]==0):
					rcm1.append(data.polyPos[n])
					rcm2.append(data.polyPos[n])
				if(data.status[n]==-1):
					rcm1.append(data.polyPos[n])
				if(data.status[n]==1):
					rcm2.append(data.polyPos[n])
	

			rcm1=np.array(rcm1)
			rcm2=np.array(rcm2)
			diff_rcm=np.mean(rcm1,axis=0) - np.mean(rcm2,axis=0)
			self.diff.append(np.dot(diff_rcm,diff_rcm.T))
			self.n_mon.append(data.nTad-self.Nchain)


	def Print(self):
		cumul=[self.n_mon,np.array(self.diff)**0.5]
		cumul=np.array(cumul)
		np.savetxt(self.diffRcmFile,cumul.T)
		

		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	rcm_diff = Poly_diff_Rcms(outputDir, initFrame=initFrame)

	rcm_diff.Compute()
	rcm_diff.Print()
