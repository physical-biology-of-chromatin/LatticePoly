##
##  DistanceMap.py
##  LatticePoly
##
##  Created by mtortora on 12/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys
import numba
import psutil
import itertools

import math as m
import numpy as np
import matplotlib.pyplot as plt

from vtkReader import vtkReader

from matplotlib.colors import LogNorm
from scipy.spatial.distance import squareform


class DistanceMap():

	def __init__(self, outputDir, initFrame, nStride, printAllFrames=True, cutoff=1/2**0.5 + 1e-3):
		self.reader = vtkReader(outputDir, initFrame, readPoly=True)

		self.nStride = int(nStride)
		self.cutoffs = np.arange(1, self.nStride+1, dtype=np.float32) * cutoff
		
		self.nBins = self.reader.nTad // self.nStride
		self.nPrune = self.reader.nTad % self.nStride
		
		self.printAllFrames = printAllFrames
		self.typeFile = os.path.join(self.reader.outputDir, "polyType.res")
			
		fontdict = {'family':'serif', 'size':'12'}

		plt.rc('font', **fontdict)
		plt.rcParams.update({'mathtext.fontset':'cm', 'mathtext.rm':'serif'})
		
		if printAllFrames:
			mapDir = os.path.join(self.reader.outputDir, "distanceMaps")
			contactDir = os.path.join(self.reader.outputDir, "contactProbs_%03d")
						
			self.mapFile = os.path.join(mapDir, "map%05d.png")
			self.contactFile = os.path.join(contactDir, "prob%05d.res")

			contactDirs = [contactDir % (i+1) for i in range(self.nStride)]

			for dir in itertools.chain([mapDir], contactDirs):
				if not os.path.exists(dir):
					os.makedirs(dir)
				
		else:
			self.mapFile = os.path.join(self.reader.outputDir, "distanceMap.pdf")
			self.contactFile = os.path.join(self.reader.outputDir, "contactProb_%03d.res")
			
			if os.path.exists(self.mapFile):
				print("File '%s' already exists - aborting" % self.mapFile)
				sys.exit()
		

	def Compute(self):
		self.cumulSqDist = 0.
		self.cumulContHist = 0

		vMem = psutil.virtual_memory()
		nStride_min = int(np.ceil(self.reader.nTad * 2 / vMem.available**0.5))

		if nStride >= nStride_min:
			polyType = self.reader.polyType[:self.reader.nTad-self.nPrune]
			polyType = polyType.reshape((self.nBins, self.nStride)).mean(axis=1)

			self.polyType = np.round(polyType).astype(np.int32)
		
			for i in range(self.reader.N):
				self.ProcessFrame(i)
			
				if (i+1) % 10 == 0:
					print("Processed %d out of %d configurations" % (i+1, self.reader.N))

		else:
			print("Memory overflow likely - increase chosen nStride (minimal value: %d)" % nStride_min)
			sys.exit()
			
			
	def ProcessFrame(self, i):
		data = next(self.reader)
		
		self.sqDist = np.zeros(self.nBins*(self.nBins-1)//2, dtype=np.float32)
		self.contHist = np.zeros((self.nBins-1, 2, self.nStride), dtype=np.float32)
		
		polyPos = data.polyPos[:self.reader.nTad-self.nPrune]
		polyPos = polyPos.reshape((self.nBins, self.nStride, 3)).mean(axis=1)

		self._sqDistPBC(data.boxDim, polyPos, self.polyType, self.cutoffs, self.contHist, self.sqDist)
		
		self.cumulSqDist += self.sqDist
		self.cumulContHist += self.contHist

		if self.printAllFrames:
			self.Print()
	

	def Print(self):
		tadDist = self.sqDist**0.5 if self.printAllFrames else self.cumulSqDist**0.5 / (self.reader.frame-self.reader.initFrame)
		tadContact = self.contHist if self.printAllFrames else self.cumulContHist

		tadMap = squareform(tadDist)
		tadContact /= tadContact.sum(axis=(0,1), keepdims=True)

		if plt.get_fignums():
			plt.clf()
		else:
			fig = plt.figure()

		dMap = plt.imshow(tadMap, extent=(0, self.reader.nTad, 0, self.reader.nTad), origin='lower', norm=LogNorm())
		
		plt.colorbar(dMap)
		
		for d in self.reader.domains:
			if len(d) > 0:
				x = [d[0], d[-1], self.reader.nTad]
				
				y1 = [d[0], d[-1], d[-1]]
				y2 = [d[0], d[0], d[0]]
	
				plt.fill_between(x=x, y1=y1, y2=y2,  color='red', alpha=0.5, lw=0)
				plt.fill_between(x=x[:2], y1=y1[:2], color='red', alpha=0.5, lw=0)

		plt.xlim([0, self.reader.nTad])
		plt.ylim([0, self.reader.nTad])

		if self.printAllFrames:
			mapFile_ = self.mapFile % (self.reader.frame-1)
			plt.savefig(mapFile_, format="png", dpi=300)

			for i in range(self.nStride):
				contactFile_ = self.contactFile % (i+1, self.reader.frame-1)
				np.savetxt(contactFile_, tadContact[:, :, i])
			
			print("\033[1;32mPrinted distance map to '%s'\033[0m" % mapFile_)

		else:
			np.savetxt(self.typeFile, self.polyType)
			plt.savefig(self.mapFile, format="pdf", transparent=True)
			
			print("\033[1;32mPrinted TAD type(s) to '%s'\033[0m" % self.typeFile)
			print("\033[1;32mPrinted distance map to '%s'\033[0m" % self.mapFile)

			for i in range(self.nStride):
				contactFile_ = self.contactFile % (i+1)
				np.savetxt(contactFile_, tadContact[:, :, i])

				
	@staticmethod
	@numba.jit("void(f4[:], f4[:,:], i4[:], f4[:], f4[:,:,:], f4[:])", nopython=True)
	def _sqDistPBC(dims, pts, types, cutoffs, contHist, sqDist):
		cnt = 0
		
		nPoints = pts.shape[0]
		nCuts = cutoffs.shape[0]
				
		for i in range(nPoints-1):
			for j in range(i+1, nPoints):
				pDist = 0.
				
				for k in range(3):
					delta = pts[j, k] - pts[i, k]
					
					while abs(delta) > dims[k] / 2.:
						delta -= m.copysign(dims[k], delta)
					
					pDist += delta**2
				
				for k in range(nCuts):
					if pDist < cutoffs[k]**2:
						for l in range(k, nCuts):
							if (types[i] == 1) & (types[j] == 1):
								contHist[j-i-1, 0, l] += 1
							else:
								contHist[j-i-1, 1, l] += 1
							
						break
					
				sqDist[cnt] = pDist
				cnt += 1
				

if __name__ == "__main__":
	if len(sys.argv) not in [4,5]:
		print("\033[1;31mUsage is %s outputDir initFrame nStride [printAllFrames]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	nStride = int(sys.argv[3])
	printAllFrames = False if len(sys.argv) == 4 else bool(sys.argv[4])
	
	distMap = DistanceMap(outputDir, initFrame=initFrame, nStride=nStride, printAllFrames=printAllFrames)
	
	distMap.Compute()
	distMap.Print()
