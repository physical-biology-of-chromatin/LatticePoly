##
##  LiqPolyContact.py
##  LatticePoly
##
##  Created by mtortora on 17/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from vtkReader import vtkReader
from scipy.spatial import cKDTree


class LiqPolyContact():
	
	def __init__(self, outputDir, initFrame, cutoff=1/2**0.5 + 1e-3):
		self.reader = vtkReader(outputDir, initFrame, readLiq=True, readPoly=True, backInBox=True)
			
		self.cutoff = cutoff
		
		self.liqFile = os.path.join(self.reader.outputDir, "liqContact.res")
		self.polyFile = os.path.join(self.reader.outputDir, "polyContact.res")
		
		if os.path.exists(self.liqFile) & os.path.exists(self.polyFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.liqFile, self.polyFile))
			sys.exit()


	def Compute(self):
		self.liqCont = np.zeros(self.reader.N, dtype=np.float32)
		self.polyCont = np.zeros(self.reader.N, dtype=np.float32)

		for i in range(self.reader.N):
			self.ProcessFrame(i)
			
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))
				
				
	def ProcessFrame(self, i):
		data = next(self.reader)
		
		if data.nHet > 0:
			liqTree = cKDTree(data.liqPos, boxsize=data.boxDim)
			polyTree = cKDTree(data.polyPos[data.polyType == 1], boxsize=data.boxDim)
			
			liqPolyIds = liqTree.query_ball_tree(polyTree, self.cutoff)
			liqPolyIds = list(filter(None, liqPolyIds))
			
			polyIds = [id for ids in liqPolyIds for id in ids]
			polyIds = set(polyIds)
			
			numLiqCont = len(liqPolyIds)
			numPolyCont = len(polyIds)
		
			self.liqCont[i] = numLiqCont / float(data.nLiq)
			self.polyCont[i] = numPolyCont / float(data.nHet)
			
		else:
			self.liqCont[i] = 0.
			self.polyCont[i] = 0.


	def Print(self):
		np.savetxt(self.liqFile, self.liqCont)
		np.savetxt(self.polyFile, self.polyCont)
		
		print("\033[1;32mPrinted mean liquid contact fractions to '%s'\033[0m" % self.liqFile)
		print("\033[1;32mPrinted mean polymer contact fractions to '%s'\033[0m" % self.polyFile)

		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	contact = LiqPolyContact(outputDir, initFrame=initFrame)

	contact.Compute()
	contact.Print()
