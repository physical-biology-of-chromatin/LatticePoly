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

class CohesinNumb():

	def __init__(self, outputDir, initFrame, threshold=0.5, nMax=10, cutoff=1/2**0.5 + 1e-3):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		self.posHist = []
		self.ForkPos =[]
		self.SisterID=[]
		self.Status=[]
		self.Cohesin=[]
		self.nMax = nMax
		self.cutoff = cutoff
		self.threshold = threshold
		

		self.cohesinFile = os.path.join(self.reader.outputDir,str(time.time())+ "cohesin.res")
		
		self.origins=[0,7,12,17,36,40,68,76,80,99,110,126,170,185,188,203,205,253,263,281,326,348,355,370,381,387,
					  404,444,454,503,511,515,562,576,598,602,644,676,687,703,719,723,731,737,754,780,813,817,
					  826,846,888,904,927,932,963,992,1021,1042,1046,1082,1103,1108,1123,1158,1163,1169,1189,1199,1202,1204,1219]

	def ReadHist(self):
		self.Nchain=0

		for i in range(self.reader.N):
			data = next(self.reader)
			self.posHist.append(data.polyPos)
			self.ForkPos.append(data.fork)
			self.SisterID.append(data.SisterID)
			self.Status.append(data.status)
			self.Cohesin.append(data.Cohesin)

			if (i==0):
				for t in range(self.reader.nTad):
					if(self.reader.status[t]==-1 or self.reader.status[t]==0):
						self.Nchain+=1
			self.dims=data.boxDim

	def computenumber(self):
		self.cohesinNumb=np.zeros(self.Nchain)

		for step in range(self.reader.N):
			for i in range(len(self.Cohesin[step])):
				if(self.Cohesin[step][i]==1 ):
					#print(i)
					#print(self.SisterID[i])
					if(i<self.Nchain):
						self.cohesinNumb[i]+=1
					else:
						j=self.SisterID[step][i]
						self.cohesinNumb[j]+=1

						


	
	
	def Print(self):
		np.savetxt(self.cohesinFile, self.cohesinNumb)

		

		
		print("\033[1;32mPrinted Forksmove  to '%s'\033[0m" % self.cohesinFile)


	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	cohesinNumb = CohesinNumb(outputDir, initFrame=initFrame)



	if len(sys.argv) == 3:
		cohesinNumb.ReadHist()
		cohesinNumb.computenumber()
		cohesinNumb.Print()

