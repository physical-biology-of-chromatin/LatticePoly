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

class Forksnumber():

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
		self.FiringFile = os.path.join(self.reader.outputDir,str(time.time())+ "firing.res")

		
		self.origins=[0,7,12,17,36,40,68,76,80,99,110,126,135,170,176,185,188,203,205,253,263,281,286,307,326,348,355,370,381,387,404,444,454,488,503,511,515,551,562,576,591,598,602,611,622,632,644,656,666,676,687,703,719,723,731,737,754,763,780,792,813,817,826,846,888,904,908,927,932,963,992,1020,1021,1042,1044,1046,1082,1084,1103,1108,1123,1143,1158,1163,1169,1175,1176,1189,1199,1202,1204,1219,1225]

	def ReadHist(self):
		self.Nchain=0

		for i in range(self.reader.N):
			data = next(self.reader)
			self.posHist.append(data.polyPos)
			self.ForkPos.append(data.Forks)
			self.SisterID.append(data.SisterID)
			self.Status.append(data.Status)
			if (i==0):
				for t in range(self.reader.nTad):
					if(self.reader.Status[t]==-1 or self.reader.Status[t]==0):
						self.Nchain+=1
			self.dims=data.boxDim

	def computenumber(self):
		self.Forksnumber=np.zeros(self.reader.N)
		self.Originsnumber=np.zeros(self.reader.N)
		self.Firing=np.zeros(self.reader.N)


		for step in range(self.reader.N):
			unrepl=0
			for i in range(len(self.ForkPos[step])):
				if(self.Status[step][i]==0):
					unrepl+=1
				if(self.ForkPos[step][i]==-1 or self.ForkPos[step][i]==1 ):
					self.Forksnumber[step]+=1
				if (self.Status[step][i]==0 and i in self.origins):
					self.Originsnumber[step]+=1
			if(unrepl!=0):
				self.Firing[step]=((10-self.Forksnumber[step]/2)*self.Originsnumber[step])/unrepl
		print(len(self.Forksnumber))

	def computeclusters(self):
		clustersdict={}

		for step in range(self.reader.N):
			forksID=[]
			for i in range(len(self.ForkPos[step])):
				if(self.ForkPos[step][i]==-1 or self.ForkPos[step][i]==1):
					forksID.append(i)
			A=np.zeros([len(forksID),len(forksID)])
			for i in range(len(forksID)):
				for j in range(i+1,len(forksID)):
					pDist = 0.
					for k in range(3):
						delta = self.posHist[step][forksID[j]][k] - self.posHist[step][forksID[i]][k]
						while abs(delta) > self.dims[k] / 2.:
							shift = math.copysign(self.dims[k], delta)

							self.posHist[step][forksID[j]][k] -= shift
							delta -=  shift
					
						pDist += delta**2
					if pDist < 3*self.cutoff**2:
						A[i][j]=1
						A[j][i]=1
			G = nx.from_numpy_matrix(A)
			if(np.array([len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]).size>0):
				clustersdict[step]=[len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
			else:
				clustersdict[step]=[]
		json.dump( clustersdict, open( self.ClusterFile, 'w' ) )
						

				
	def computetiming(self):
		timing=np.zeros(self.Nchain)
		for i in range(len(timing)):
			timing[i]=-1
		for step in range(self.reader.N):
			for i in range(self.Nchain):
				if(self.SisterID[step][i]!=-1):
					if(timing[i]<0):
						timing[i]=step
		np.savetxt(self.TimingFile, timing)

					


	
	
	def Print(self):
		np.savetxt(self.ForksnumbFile, self.Forksnumber)
		np.savetxt(self.OriginsnumbFile, self.Originsnumber)
		np.savetxt(self.FiringFile, self.Firing)


		

		
		print("\033[1;32mPrinted Forksmove  to '%s'\033[0m" % self.ForksnumbFile)


	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	forksnumb = Forksnumber(outputDir, initFrame=initFrame)



	if len(sys.argv) == 3:
		forksnumb.ReadHist()
		forksnumb.computenumber()
		forksnumb.computeclusters()
		forksnumb.computetiming()
		forksnumb.Print()

