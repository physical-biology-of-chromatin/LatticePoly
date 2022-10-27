##
##  PolyGyration.py
##  LatticePoly
##
##  Created by mtortora on 02/08/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
import os
import sys

import numpy as np

from scipy.spatial import cKDTree

from vtkReader import vtkReader

from scipy.spatial.distance import pdist, squareform
import time

class corr_distances():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
					
		self.corr_distancesFile = os.path.join(self.reader.outputDir, str(time.time())+"correlation_distances.res")

		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1


	def Compute(self):
		self.pdists=[]
		for i in range(self.reader.N):
			self.ProcessFrame(i)
			
			#if (i+1) % 10 == 0:
			#	print("Processed %d out of %d configurations" % (i+1, self.reader.N))
			
				
	def ProcessFrame(self, i):
		data = next(self.reader)
		
		if(data.nTad==self.Nchain and len(self.pdists)==0):
			self.pdists.append(pdist(data.polyPos[:self.Nchain]))
			print(i)
			
		if(data.nTad==2*self.Nchain and len(self.pdists)!=2):
			self.pdists.append(pdist(data.polyPos[:self.Nchain]))
			print(i)


	def Print(self):
		if(len(self.pdists)!=2):
			print("Not find two configuration")
			sys.exit()
		corr=[]
		corr.append(np.corrcoef(self.pdists[0],self.pdists[1])[0][1])
		np.savetxt(self.corr_distancesFile, corr)
		

		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame \033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	corr = corr_distances(outputDir, initFrame=initFrame)

	corr.Compute()
	corr.Print()
