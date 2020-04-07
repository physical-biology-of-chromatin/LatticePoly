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

import math as m
import numpy as np
import matplotlib.pyplot as plt

from vtkReader import vtkReader

from matplotlib.colors import LogNorm
from scipy.spatial.distance import squareform


class DistanceMap():

	def __init__(self, outputDir, initFrame, printAllFrames=True):
		self.reader = vtkReader(outputDir, initFrame, readPoly=True)

		self.printAllFrames = printAllFrames
		
		fontdict = {'family':'serif', 'size':'12'}

		plt.rc('font', **fontdict)
		plt.rcParams.update({'mathtext.fontset':'cm',
							 'mathtext.rm':'serif'})
							 
		if printAllFrames:
			mapDir = os.path.join(self.reader.outputDir, "distanceMaps")
			self.mapFile = os.path.join(mapDir, "map%05d.png")

			if not os.path.exists(mapDir):
				os.makedirs(mapDir)
						
		else:
			self.mapFile = os.path.join(self.reader.outputDir, "distanceMap.pdf")
					

	def Compute(self):
		self.cumulDist = 0.

		for i in range(self.reader.N):
			self.ProcessFrame()
			
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))

			
	def ProcessFrame(self):
		data = next(self.reader)

		self.dist = np.sqrt(self._sqDistPBC(data.boxDim, data.polyPos))
		self.cumulDist += self.dist
					
		if self.printAllFrames:
			self.Print()
	

	def Print(self):
		if self.reader.frame == self.reader.initFrame:
			print("Did not process any files - nothing to print")
						
		else:
			tadDist = self.dist if self.printAllFrames else self.cumulDist / (self.reader.frame-self.reader.initFrame)
			tadMap = squareform(tadDist)
			
			tadDomains = np.nonzero(self.reader.polyType)[0]
			tadDomains = np.split(tadDomains, np.where(np.diff(tadDomains) != 1)[0]+1)
		
			if plt.get_fignums():
				plt.clf()
				
			else:
				fig = plt.figure()

			dMap = plt.imshow(tadMap, norm=LogNorm())
			plt.colorbar(dMap)
			
			for domain in tadDomains:
				x = [domain[0], domain[-1], self.reader.nTad]
				
				y1 = [domain[0], domain[-1], domain[-1]]
				y2 = [domain[0], domain[0], domain[0]]
		
				plt.fill_between(x=x, y1=y1, y2=y2,  color='red', alpha=0.25)
				plt.fill_between(x=x[:2], y1=y1[:2], color='red', alpha=0.25)

			plt.xlim([0, self.reader.nTad])
			plt.ylim([0, self.reader.nTad])

			if self.printAllFrames:
				mapFile_ = self.mapFile % (self.reader.frame-1)
				
				plt.savefig(mapFile_, format="png", dpi=300)
				print("\033[1;32mPrinted figure to '%s'\033[0m" % mapFile_)

			else:
				plt.savefig(self.mapFile, format="pdf", transparent=True)
				print("\033[1;32mPrinted figure to '%s'\033[0m" % self.mapFile)

				plt.show()
				
				
	@staticmethod
	@numba.jit("f4[:](f4[:], f4[:,:])", nopython=True)
	def _sqDistPBC(dims, pts):
		cnt = 0
		n = pts.shape[0]
		
		sqDist = np.zeros(n*(n-1)//2, dtype=np.float32)
		
		for i in range(n-1):
			for j in range(i+1, n):
				for k in range(3):
					delta = pts[j, k] - pts[i, k]
					
					while abs(delta) > dims[k] / 2.:
						delta -= m.copysign(dims[k], delta)

					sqDist[cnt] += delta**2
					
				cnt += 1
				
		return sqDist


if __name__ == "__main__":
	if len(sys.argv) not in [3,4]:
		print("\033[1;31mUsage is %s outputDir initFrame [printAllFrames]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	printAllFrames = False if len(sys.argv) == 3 else bool(sys.argv[3])
	distMap = DistanceMap(outputDir, initFrame=initFrame, printAllFrames=printAllFrames)
	
	distMap.Compute()
	distMap.Print()
