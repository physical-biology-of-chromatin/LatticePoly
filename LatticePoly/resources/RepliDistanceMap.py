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
from scipy.spatial.distance import pdist, squareform,cdist
from numpy import nanmean



class RepliDistanceMap():

	def __init__(self, outputDir, initFrame, finalFrame, threshold=0.5, nMax=10, cutoff=1/2**0.5 + 1e-3):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		self.posHist = []
		self.SisterID=[]
		self.Status=[]
		self.nMax = nMax
		self.cutoff = cutoff
		self.threshold = threshold
		self.posHistsister=[]

		self.CisMap=[]
		self.TransMap=[]
		self.AllInter=[]
		

		self.CisMapFile = os.path.join(self.reader.outputDir,str(time.time())+"_"+str(initFrame)+"_"+str(finalFrame)+ "CisMap.res")
		self.TransMapFile = os.path.join(self.reader.outputDir,str(time.time())+"_"+str(initFrame)+"_"+str(finalFrame)+ "TransMap.res")
		self.AllFile = os.path.join(self.reader.outputDir,str(time.time())+"_"+str(initFrame)+"_"+str(finalFrame)+ "All.res")
		
		
	def ReadHist(self):
		self.Nchain=0

		for i in range(finalFrame-initFrame):
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
	
	def computesister(self):
		for i in range(len(self.posHist)):
			sispos=[]
			for tad in range(self.Nchain):
				if(self.status[i][tad]!=0):
					sispos.append(self.posHist[i][self.SisterID[i][tad]])
				
				else:
					sispos.append([np.nan,np.nan,np.nan])


			self.posHistsister.append(sispos)

						
	
	
	def Maps(self):
		cismap=[]
		transmap=[]
		all=[]
		tree1 =cKDTree(
		for i in range(len(self.posHist)):
			distscis1 = pdist(self.posHist[i][0:self.Nchain])
			distanceMapcis1 = squareform(distscis)	
			
			
			distscis1 = pdist(self.posHistsister1[i][0:self.Nchain])
			distanceMapcis1 = squareform(distscis1)	
			distscis2 = pdist(self.posHistsister2[i][0:self.Nchain])
			distanceMapcis2 = squareform(distscis2)		
			diststrans = cdist(self.posHistsister1[i][0:self.Nchain],self.posHistsister2[i][0:self.Nchain])
			#distsall=pdist(self.nonRepli[i][0:self.Nchain])

			
			"""tree1   = cKDTree(data.polyPos[domainstart:domainend], boxsize = None)
				pairs = tree1.query_pairs(r = 3.53) # NN distance FCC lattice 1/np.sqrt(2) = 0.71
				for (i,j) in pairs:
				self.contactProb[i,j] = self.contactProb[i,j] + 1"""
			cismap.append(distanceMapcis)
			#cismap.append(distanceMapcis2)		
			#transmap.append(nanmean([diststrans,diststrans.transpose()],axis=0))
			#all.append(distanceMapcis1)
			#all.append(distanceMapcis2)
			#all.append(transmap[-1])
			#all.append(nanmean([distsall,distsall.transpose()],axis=0))
			
		self.CisMap=nanmean(cismap,axis=0)
		#self.transMap=nanmean(transmap,axis=0)
		#self.AllMap=nanmean(all,axis=0)
		
		np.savetxt(self.CisMapFile, self.CisMap)
		#np.savetxt(self.TransMapFile, self.transMap)
		#np.savetxt(self.AllFile, self.AllMap)

	
	
if __name__ == "__main__":
	if len(sys.argv) not in [ 4]:
		print("\033[1;31mUsage is %s outputDir initFrame finalFrame \033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	finalFrame = int(sys.argv[3])
	
	RepliDistanceMap = RepliDistanceMap(outputDir, initFrame=initFrame,finalFrame=finalFrame )



	if len(sys.argv) == 4:
		RepliDistanceMap.ReadHist()
		RepliDistanceMap.computesister()
		RepliDistanceMap.Maps()

