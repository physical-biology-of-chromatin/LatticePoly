##  LatticePoly
##
##  Created by ddasaro on the model of mtortora script.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

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
					
		self.corr_distancesFile = os.path.join(self.reader.outputDir, str(time.time())+"correlation_distances_after.res")

		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1


	def Compute(self):
		self.corr_t=[]
		self.n_mon=[]

		for i in range(self.reader.N):
			self.ProcessFrame(i)
			
			#if (i+1) % 10 == 0:
			#	print("Processed %d out of %d configurations" % (i+1, self.reader.N))
			
				
	def ProcessFrame(self, i):
		data = next(self.reader)
		
		if(data.nTad==2*self.Nchain and len(self.corr_t)==0):
			self.pdist_init=pdist(data.polyPos[:self.Nchain])
			chain1 = data.polyPos[:self.Nchain]
			chain2 =[data.polyPos[data.SisterID[k]] for k in range(self.Nchain)]

			self.corr_t.append((np.corrcoef((self.pdist_init,pdist(chain1)))[0][1]+np.corrcoef(self.pdist_init,pdist(chain2))[0][1])/2)
			self.n_mon.append(data.nTad-self.Nchain)
			
		if(len(self.corr_t)>0):
			chain1 = data.polyPos[:self.Nchain]
			chain2 =[data.polyPos[data.SisterID[k]] for k in range(self.Nchain)]

			self.corr_t.append((np.corrcoef((self.pdist_init,pdist(chain1)))[0][1]+np.corrcoef(self.pdist_init,pdist(chain2))[0][1])/2)
			self.n_mon.append(data.nTad-self.Nchain)



	def Print(self):
		cumul=[self.n_mon,np.array(self.corr_t)]
		cumul=np.array(cumul)
		np.savetxt(self.corr_distancesFile, cumul.T)
		

		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame \033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	corr = corr_distances(outputDir, initFrame=initFrame)

	corr.Compute()
	corr.Print()
