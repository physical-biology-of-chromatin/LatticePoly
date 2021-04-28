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


class LiqDensity():

	def __init__(self, outputDir, initFrame, threshold=0.5):
		self.reader = vtkReader(outputDir, initFrame, readLiq=True, readPoly=False)
		
		self.threshold = threshold

		self.meanFile = os.path.join(self.reader.outputDir, "liqMean.res")
		self.stdFile = os.path.join(self.reader.outputDir, "liqSTD.res")
		
		if os.path.exists(self.meanFile) & os.path.exists(self.stdFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.meanFile, self.stdFile))
			sys.exit()


	def Compute(self):
		self.meanHist = np.zeros(self.reader.N, dtype=np.float32)
		self.stdHist = np.zeros(self.reader.N, dtype=np.float32)

		for i in range(self.reader.N):
			self.ProcessFrame(i)
									
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))

			
	def ProcessFrame(self, i):
		data = next(self.reader)

		meanDens = data.liqDens.sum()
		stdDens = np.square(data.liqDens - data.liqDens.mean()).sum()
		
		self.meanHist[i] = np.count_nonzero(data.liqDens > self.threshold)
		self.stdHist[i] = stdDens

	
	def Print(self):
		np.savetxt(self.meanFile, self.meanHist / self.reader.nLiq)
		np.savetxt(self.stdFile, np.sqrt(self.stdHist / self.reader.nLiq))

		print("\033[1;32mPrinted liquid mean densities to '%s'\033[0m" % self.meanFile)
		print("\033[1;32mPrinted liquid density STDs fraction to '%s'\033[0m" % self.stdFile)


if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	density = LiqDensity(outputDir, initFrame=initFrame)

	density.Compute()
	density.Print()
