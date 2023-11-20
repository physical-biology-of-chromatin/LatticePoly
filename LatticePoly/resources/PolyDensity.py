##
##  PolyGyration.py
##  LatticePoly
##
##  Created by amithzafal on 02/08/2020.
##  Copyright © 2023 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from vtkReader import vtkReader


class PolyGyration():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
					
		self.anisoFile = os.path.join(self.reader.outputDir, "homDensity.res")
		self.gyrationFile = os.path.join(self.reader.outputDir, "polyDensity.res")

		if os.path.exists(self.anisoFile) & os.path.exists(self.gyrationFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.anisoFile, self.gyrationFile))
			sys.exit()


	def Compute(self):
		self.polyAniso = np.zeros((self.reader.N, self.reader.nTad), dtype=np.float32)
		self.polyGyration = np.zeros((self.reader.N, self.reader.nTad), dtype=np.float32)

		for i in range(self.reader.N):
			self.ProcessFrame(i)
			
			if (i+1) % 1000 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))
				
				
	def ProcessFrame(self, i):
		data = next(self.reader)					
		self.polyAniso[i]     = data.homDensity
		self.polyGyration[i]  = data.density


	def Print(self):
		np.savetxt(self.anisoFile, self.polyAniso)
		np.savetxt(self.gyrationFile, self.polyGyration)
		
		print("\033[1;32mPrinted local density to '%s'\033[0m" % self.gyrationFile)
		print("\033[1;32mPrinted local chain specific density factors to '%s'\033[0m" % self.anisoFile)

		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	gyr = PolyGyration(outputDir, initFrame=initFrame)

	gyr.Compute()
	gyr.Print()
