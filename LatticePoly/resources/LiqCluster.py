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
from scipy.spatial import cKDTree


class LiqCluster():

	def __init__(self, outputDir, initFrame, nMax=100, threshold=0.5, cutoff=1/2**0.5, tol=1e-3):
		self.reader = vtkReader(outputDir, initFrame, readLiq=True, readPoly=False)
		
		self.nMax = nMax
		
		self.cutoff = cutoff + tol
		self.threshold = threshold - tol
		
		self.numFile = os.path.join(self.reader.outputDir, "liqNum.res")
		
		self.radFile = os.path.join(self.reader.outputDir, "liqRadii.res")
		self.anisoFile = os.path.join(self.reader.outputDir, "liqAniso.res")

		self.meanRadFile = os.path.join(self.reader.outputDir, "liqMeanRadii.res")
		self.meanAnisoFile = os.path.join(self.reader.outputDir, "liqMeanAniso.res")
		
		self.fileList = [self.numFile, self.radFile, self.anisoFile, self.meanRadFile, self.meanAnisoFile]

		if all(map(os.path.exists, self.fileList)):
			print("Files %s already exist - aborting" % ", ".join(f"'{filePath}'" for filePath in self.fileList))
			sys.exit()


	def Compute(self):
		self.dropNum = np.zeros(self.reader.N, dtype=np.int32)
		
		self.radHist = np.zeros((self.reader.N, self.nMax), dtype=np.float32)
		self.anisoHist = np.zeros((self.reader.N, self.nMax), dtype=np.float32)

		self.meanRad = np.zeros(self.reader.N, dtype=np.float32)
		self.meanAniso = np.zeros(self.reader.N, dtype=np.float32)
	
		for i in range(self.reader.N):
			self.ProcessFrame(i)
									
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))

			
	def ProcessFrame(self, i):
		data = next(self.reader)

		dropPos = data.liqPos[data.liqDens >= self.threshold]
		liqTree = cKDTree(dropPos, boxsize=data.boxDim)
		
		pairs = liqTree.query_pairs(r=self.cutoff)
		
		graph = nx.Graph(pairs)
		clusters = nx.connected_components(graph)
		
		clusters = sorted(clusters, key=len, reverse=True)
		cluster_ids = [np.asarray(list(cluster), dtype=np.int32) for cluster in clusters]
				
		sizes = [len(cluster) for cluster in clusters]
		
		volumes = np.asarray(sizes, dtype=np.float32) / 4.
		radii = (3 * volumes / (4*np.pi))**(1/3.)

		aniso = np.zeros_like(volumes)
		
		for j, ids in enumerate(cluster_ids):
			if sizes[j] > 2:
				clusterPBC = dropPos[ids]
				
				self._fixClusterPBC(data.boxDim, clusterPBC)
				
				clusterPBC -= clusterPBC.mean(axis=0, keepdims=True)
				eig = np.linalg.svd(clusterPBC, compute_uv=False)
				
				aniso[j] = 3/2.*(eig**4).sum(axis=-1)/(eig**2).sum(axis=-1)**2 - 1/2.
				
		n = len(radii)
		m = min(n, self.nMax)
		
		self.dropNum[i] = n

		self.radHist[i, :m] = radii[:m]
		self.anisoHist[i, :m] = aniso[:m]
		
		self.meanRad[i] = radii.mean() if n > 0 else 0.
		self.meanAniso[i] = np.average(aniso, weights=volumes) if n > 0 else 0.

	
	def Print(self):
		np.savetxt(self.numFile, self.dropNum)

		np.savetxt(self.radFile, self.radHist)
		np.savetxt(self.anisoFile, self.anisoHist)

		np.savetxt(self.meanRadFile, self.meanRad)
		np.savetxt(self.meanAnisoFile, self.meanAniso)

		print("\033[1;32mPrinted clustering analysis to %s\033[0m" % ", ".join(f"'{filePath}'" for filePath in self.fileList))
		
	
	@staticmethod
	@numba.jit("void(f4[:], f4[:,:])", nopython=True)
	def _fixClusterPBC(dims, pts):
		nPoints = pts.shape[0]
		
		for i in range(1, nPoints):
			for j in range(3):
				delta = pts[i, j] - pts[0, j]
					
				while abs(delta) > dims[j] / 2.:
					shift = m.copysign(dims[j], delta)
					
					pts[i, j] -= shift
					delta -= shift
											
					
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	cluster = LiqCluster(outputDir, initFrame=initFrame)

	cluster.Compute()
	cluster.Print()
