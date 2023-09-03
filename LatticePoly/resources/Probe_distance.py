##
##  Probe_distance.py
##  LatticePoly
##
##  Created by mtortora on 10/11/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

import os
import sys
import numba
import pandas as pd

import math as m
import numpy as np
import networkx as nx

from vtkReader import vtkReader


class Probe_distance():

	def __init__(self, outputDir, initFrame,gene, threshold=0.5, nMax=10, cutoff=1/2**0.5 + 1e-3):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		
		self.Ntot = self.reader.nTad
		self.gene = gene

		if self.gene == "BX":
			gene_start, gene_end = 16300000, 17700000
			probe_localization = pd.read_excel('/home/ppuel/Painter_Clone5/LatticePoly/LatticePoly/data/BX_probe_localization.ods')
		elif self.gene == 'ANT':
			gene_start, gene_end = 4500000, 8000000
			probe_localization = pd.read_excel('/home/ppuel/Painter_Clone5/LatticePoly/LatticePoly/data/ANT_probe_localization.ods')

		bp_per_tad = 800

		shift = int((self.Ntot - (gene_end - gene_start)/bp_per_tad - 1 ))//2 

		self.probe_size = 1e-26

		self.probe_N = len(probe_localization.keys())

		for key in probe_localization.keys():
			probe_start, probe_end = probe_localization[key]
			self.probe_size = max(self.probe_size,(probe_end - probe_start)//bp_per_tad )


		self.probe_tad = np.zeros((self.probe_N,self.probe_size),dtype=np.int64)

		for i in range(self.probe_N):
			probe_start, probe_end = probe_localization[probe_localization.keys()[i]]
			self.probe_tad[i] = np.array([int(probe_start-gene_start)//bp_per_tad+shift+j for j in range(self.probe_size)], dtype = np.int64)


		self.nMax = nMax
		self.cutoff = cutoff
		self.threshold = threshold
		
		self.distFile = os.path.join(self.reader.outputDir, "Probe_distance.res")

		if os.path.exists(self.distFile):
			print("Files '%s' already exist - aborting" % (self.distFile))
			sys.exit()

	def Compute(self):
		self.distProbe = np.zeros((self.reader.N, self.probe_N)) 		# /!\ This part of the script only work if there is 3 probes

		for i in range(self.reader.N):
			self.ProcessFrame(i)
									
			if (i+1) % 1000 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))

			
	def ProcessFrame(self, i):
		data = next(self.reader)
		
		probe_pos = np.zeros((self.probe_N, self.probe_size, 3))

		for j in range(self.probe_N):
			probe_pos[j] = data.polyPos[self.probe_tad[j],:]
		
		probe_center = np.mean(probe_pos,1)

		# /!\ This part of the script only work if there is 3 probes

		self.distProbe[i,0] = np.sqrt((probe_center[0,0] - probe_center[1,0])**2+(probe_center[0,1] - probe_center[1,1])**2+(probe_center[0,0] - probe_center[1,2])**2)
		self.distProbe[i,1] = np.sqrt((probe_center[2,0] - probe_center[1,0])**2+(probe_center[2,1] - probe_center[1,1])**2+(probe_center[2,0] - probe_center[1,2])**2)
		self.distProbe[i,2] = np.sqrt((probe_center[0,0] - probe_center[2,0])**2+(probe_center[0,1] - probe_center[2,1])**2+(probe_center[0,0] - probe_center[2,2])**2)


	
	def Print(self):
		np.savetxt(self.distFile, self.distProbe)
		print("\033[1;32mPrinted distance Probe to '%s'\033[0m" % self.distFile)

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
	if len(sys.argv) != 4:
		print("\033[1;31mUsage is %s outputDir initFrame gene\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	gene = sys.argv[3]

	cluster = Probe_distance(outputDir, initFrame=initFrame, gene=gene)

	cluster.Compute()
	cluster.Print()
