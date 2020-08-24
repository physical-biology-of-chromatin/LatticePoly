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
from scipy.spatial.distance import cdist


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
		data = next(self.reader) if i > 0 else self.reader
		liqPolyDist = cdist(data.liqPos, data.polyPos[data.polyType == 1])
		
		liqDist = np.min(liqPolyDist, axis=1)
		polyDist = np.min(liqPolyDist, axis=0)
		
		numLiqCont = np.count_nonzero(liqDist < self.cutoff)
		numPolyCont = np.count_nonzero(polyDist < self.cutoff)
		
		self.liqCont[i] = numLiqCont / float(self.reader.nLiq)
		self.polyCont[i] = numPolyCont / float(self.reader.nHet)


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
