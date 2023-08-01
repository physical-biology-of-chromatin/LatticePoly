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

class replicons():

	def __init__(self, outputDir, initFrame, threshold=0.5, nMax=10, cutoff=1/2**0.5 + 1e-3):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		self.posHist = []
		self.ForkPos =[]
		self.SisterID=[]
		self.Status=[]
		self.nMax = nMax
		self.cutoff = cutoff
		self.threshold = threshold
		

		self.ClusterFile = os.path.join(self.reader.outputDir,str(time.time())+ "Cluster.json")
		self.ForksnumbFile = os.path.join(self.reader.outputDir,str(time.time())+ "Forksnumb.res")
		self.OriginsnumbFile = os.path.join(self.reader.outputDir,str(time.time())+ "Origins.res")
		self.TimingFile = os.path.join(self.reader.outputDir,str(time.time())+ "timing.res")
		
		self.origins=[525,575,628,677]
		self.repltime=[]


	def ReadHist(self):
		self.Nchain=0

		for i in range(self.reader.N):
			data = next(self.reader)
			self.posHist.append(data.polyPos)
			self.SisterID.append(data.SisterID)
			self.Status.append(data.status)
			if (i==0):
				for t in range(self.reader.nTad):
					if(self.reader.status[t]==-1 or self.reader.status[t]==0):
						self.Nchain+=1
			self.dims=data.boxDim

	
	
	
	def computenumber(self):
		for origin in self.origins: 
			for step in range(self.reader.N):
				if(self.SisterID[step][origin]!=-1):
					self.repltime.append(step)
					break
		for i in range(1,len(self.repltime)):
			if(abs(self.repltime[0]-self.repltime[i])<100):
				self.arraysDistance(i)
				



	def arraysDistance(self,i):
		arrayDistance=[]
		midreplicationtime=round((self.repltime[0]+self.repltime[i])/2)
		for step in range(midreplicationtime-100,midreplicationtime+66):
			diff=self.posHist[step][self.origins[0]]-self.posHist[step][self.origins[i]]
			arrayDistance.append(np.sqrt(np.dot(diff.T,diff)))
		np.savetxt(os.path.join(self.reader.outputDir,str(time.time())+ "distance_replicons_0_"+str(i)+".res"),arrayDistance)
										   


	

	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	replicons = replicons(outputDir, initFrame=initFrame)



	if len(sys.argv) == 3:
		replicons.ReadHist()
		replicons.computenumber()

