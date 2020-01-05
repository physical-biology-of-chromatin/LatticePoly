##
##  LiqDensity.py
##  LatticePoly
##
##  Created by mtortora on 29/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from vtkReader import vtkReader


class LiqDensity(vtkReader):

	def __init__(self, outputDir, initFrame):
		vtkReader.__init__(self, outputDir)

		self.InitReader(initFrame, readLiq=True)
		
		self.meanFile = os.path.join(self.outputDir, "liqMean.res")
		self.stdFile = os.path.join(self.outputDir, "liqSTD.res")


	def Compute(self):
		self.meanHist = np.zeros(self.N, dtype=np.float32)
		self.stdHist = np.zeros(self.N, dtype=np.float32)

		for _ in range(self.N):
			self.ProcessFrame()
									
			if (self.frame-self.initFrame) % 10 == 0:
				print("Processed %d out of %d configurations" % (self.frame-self.initFrame, self.N))

			
	def ProcessFrame(self):
		self.ReadLiqFrame()
				
		idx = self.frame-self.initFrame
		
		meanDens = self.liqDens.sum()
		stdDens = np.square(self.liqDens-self.liqDens.mean()).sum()
		
		self.meanHist[idx] = meanDens
		self.stdHist[idx] = stdDens

		self.frame += 1

	
	def Print(self):
		np.savetxt(self.meanFile, self.meanHist / self.nLiq)
		np.savetxt(self.stdFile, np.sqrt(self.stdHist / self.nLiq))

		print("\033[1;32mPrinted liquid mean densities to '%s'\033[0m" % self.meanFile)
		print("\033[1;32mPrinted liquid density STDs fraction to '%s'\033[0m" % self.stdFile)


if len(sys.argv) != 3:
	print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])

	sys.exit()

outputDir = sys.argv[1]
initFrame = int(sys.argv[2])

density = LiqDensity(outputDir, initFrame=initFrame)

density.Compute()
density.Print()
