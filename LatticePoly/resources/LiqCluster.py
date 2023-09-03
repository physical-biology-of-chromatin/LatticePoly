##
##  LiqCluster.py
##  LatticePoly
##
##  Created by mtortora on 10/11/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

import os
import sys
import numba

import math as m
import numpy as np
import networkx as nx

from vtkReader import vtkReader


class LiqCluster():

	def __init__(self, outputDir, initFrame, threshold=0.5, nMax=10, cutoff=1/2**0.5 + 1e-3):
		self.reader = vtkReader(outputDir, initFrame, readLiq=True, readPoly=False)
		
		
		self.nMax = nMax
		self.cutoff = cutoff
		self.threshold = threshold
		
		self.numFile = os.path.join(self.reader.outputDir, "liqDropNum.res")
		self.radFile = os.path.join(self.reader.outputDir, "liqDropRad.res")
		self.meanFile = os.path.join(self.reader.outputDir, "liqRadii.res")
		self.fractionFile = os.path.join(self.reader.outputDir, "liqFraction.res")

		if os.path.exists(self.numFile) & os.path.exists(self.radFile) & os.path.exists(self.meanFile) & os.path.exists(self.fractionFile):
			print("Files '%s', '%s' and '%s' already exist - aborting" % (self.numFile, self.radFile, self.fractionFile))
			sys.exit()

	def Compute(self):
		self.dropNum = np.zeros(self.reader.N, dtype=np.int32)
		self.dropMean = np.zeros(self.reader.N, dtype=np.float32)
		self.dropRad = np.zeros((self.reader.N, self.nMax), dtype=np.float32)
		self.liqFraction = np.zeros(self.reader.N,dtype=np.float32)

		for i in range(self.reader.N):
			self.ProcessFrame(i)
									
			if (i+1) % 1000 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))

			
	def ProcessFrame(self, i):
		data = next(self.reader)
		nLiq = self.reader.nLiq

		dropPos = data.liqPos[data.liqDens > self.threshold]
		connect = self._connectPBC(data.boxDim, dropPos, self.cutoff)
		
		graph = nx.from_numpy_array(connect)
		clusters = nx.connected_components(graph)
		clusters = sorted(clusters, key=len, reverse=True)
		
		sizes = [len(cluster) for cluster in clusters]
		volumes = np.asarray(sizes) / 4.
		
		
		radii = (3 * volumes / (4*np.pi))**(1/3.)
		n = len(radii)
		m = min(n, self.nMax)
		
		self.liqFraction[i] = np.sum(np.asarray(sizes))/nLiq
		self.dropNum[i] = n
		self.dropMean[i] = radii.mean() if n > 0 else 0.
		
		self.dropRad[i, :m] = radii[:m]

	
	def Print(self):
		np.savetxt(self.radFile, self.dropRad)
		np.savetxt(self.numFile, self.dropNum)
		np.savetxt(self.meanFile, self.dropMean)
		np.savetxt(self.fractionFile, self.liqFraction)
		print("\033[1;32mPrinted droplet numbers to '%s'\033[0m" % self.numFile)
		print("\033[1;32mPrinted individual/mean droplet radii to '%s' and '%s'\033[0m" % (self.radFile, self.meanFile))
		print("\033[1;32mPrinted liquid fraction to '%s'\033[0m" % self.fractionFile)

	@staticmethod
	@numba.jit("i4[:,:](f4[:], f4[:,:], f4)", nopython=True)
	def _connectPBC(dims, pts, cutoff):
		nPoints = pts.shape[0]
		connect = np.zeros((nPoints, nPoints), dtype=np.int32)
		
		for i in range(nPoints):
			for j in range(i+1, nPoints):
				pDist = 0.
						
				for k in range(3):
					delta = pts[j, k] - pts[i, k]
					
					while abs(delta) > dims[k] / 2.:
						shift = m.copysign(dims[k], delta)
					
						pts[j, k] -= shift
						delta -= shift
						
					pDist += delta**2
				
				if pDist < cutoff**2:
					connect[i, j] = 1
					connect[j, i] = 1

		return connect
					
					
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	cluster = LiqCluster(outputDir, initFrame=initFrame)

	cluster.Compute()
	cluster.Print()
