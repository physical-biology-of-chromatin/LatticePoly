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
	
	def __init__(self, outputDir, chrom1, chrom2, initFrame):
		self.reader1 = vtkReader(outputDir, chrom1, initFrame, readLiq=False, readPoly=True)
		self.reader2 = vtkReader(outputDir, chrom2, initFrame, readLiq=False, readPoly=True)
					
		self.anisoFile = os.path.join(self.reader1.outputDir, "homDensity.res")
		self.gyrationFile = os.path.join(self.reader1.outputDir, "polyDensity.res")

		if os.path.exists(self.anisoFile) & os.path.exists(self.gyrationFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.anisoFile, self.gyrationFile))
			sys.exit()


	def Compute(self):
		self.polyAniso = np.zeros((self.reader1.N, 2*self.reader1.nTad), dtype=np.float32)
		self.polyGyration = np.zeros((self.reader1.N, 2*self.reader1.nTad), dtype=np.float32)

		for i in range(self.reader1.N):
			self.ProcessFrame(i)
			
			if (i+1) % 1000 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))
				
				
	def ProcessFrame(self, i):
		data1 = next(self.reader1)	
		data2 = next(self.reader2)				
		self.polyAniso[i]     = np.append(data1.homDensity,data2.homDensity)
		self.polyGyration[i]  = np.append(data1.density,data2.density)


	def Print(self):
		np.savetxt(self.anisoFile, self.polyAniso)
		np.savetxt(self.gyrationFile, self.polyGyration)
		
		print("\033[1;32mPrinted local density to '%s'\033[0m" % self.gyrationFile)
		print("\033[1;32mPrinted local chain specific density factors to '%s'\033[0m" % self.anisoFile)

		
if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("\033[1;31mUsage is %s outputDir chrom1 chrom2 initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	chrom1 = sys.argv[2]
	chrom2 = sys.argv[3]
	initFrame = int(sys.argv[4])

	gyr = PolyGyration(outputDir, chrom1, chrom2, initFrame=initFrame)

	gyr.Compute()
	gyr.Print()
