##
##  LiqDensity.py
##  LatticePoly
##
##  Created by mtortora on 29/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import sys

import numpy as np

from utils import getInputParam
from vtkReader import vtkReader


class LiqDensity(vtkReader):

	def __init__(self, outputDir, initFrame):
		vtkReader.__init__(self, outputDir)

		self.InitReader(initFrame, readLiq=True)

		self.inputFile = "%s/.input.cfg" % self.outputDir
		
		self.densFile = "%s/liqDensity.res" % self.outputDir
		self.fracFile = "%s/liqFraction.res" % self.outputDir


	def Compute(self):
		self.densHist = np.zeros(self.N, dtype=np.float32)

		for i in range(self.N):
			self.ProcessFrame()
			
			densTot = self.liqDens.sum()
			latSize = 4*self.boxDim.prod()
			
			self.densHist[i] = densTot / latSize
			
			if (self.frame-self.initFrame) % 10 == 0:
				print("Processed %d out of %d configurations" % (self.frame-self.initFrame, self.N))

			
	def ProcessFrame(self):
		self.ReadLiqFrame()
	
		self.frame += 1

	
	def Print(self):
		Ldens = float(getInputParam('Ldens', self.inputFile))

		np.savetxt(self.densFile, self.densHist)
		np.savetxt(self.fracFile, self.densHist / Ldens)

		print("\033[1;32mPrinted liquid density to '%s'\033[0m" % self.densFile)
		print("\033[1;32mPrinted liquid dense fraction to '%s'\033[0m" % self.fracFile)


if len(sys.argv) != 3:
	print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])

	sys.exit()

outputDir = sys.argv[1]
initFrame = int(sys.argv[2])

density = LiqDensity(outputDir, initFrame=initFrame)

density.Compute()
density.Print()
