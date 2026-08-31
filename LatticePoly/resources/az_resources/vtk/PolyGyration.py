##
##  PolyGyration.py
##  LatticePoly
##
##  Created by mtortora on 02/08/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from vtkReader import vtkReader


class PolyGyration():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
					
		self.anisoFile = os.path.join(self.reader.outputDir, "polyAniso.res")
		self.gyrationFile = os.path.join(self.reader.outputDir, "polyGyration.res")

		if os.path.exists(self.anisoFile) & os.path.exists(self.gyrationFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.anisoFile, self.gyrationFile))
			sys.exit()


	def Compute(self):
		self.polyAniso = np.zeros(self.reader.N, dtype=np.float32)
		self.polyGyration = np.zeros(self.reader.N, dtype=np.float32)

		for i in range(self.reader.N):
			self.ProcessFrame(i)
			
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))
				
				
	def ProcessFrame(self, i):
		data = next(self.reader)
		
		norm = 0
		
		for id, d in enumerate(self.reader.domains):
			if d.size > 2:
				pos = data.polyPos[d]
				pos -= pos.mean(axis=0, keepdims=True)
			
				diag = np.linalg.svd(pos, compute_uv=False) * np.sqrt(12)/d.size
				r2_gyr = np.square(diag).sum(axis=-1)
			
				r_gyr = np.sqrt(r2_gyr)
				aniso = 3/2.*(diag**4).sum(axis=-1)/r2_gyr**2 - 1/2.
			
				norm += d.size
				
				self.polyAniso[i] += aniso * d.size
				self.polyGyration[i] += r_gyr * d.size
									
		self.polyAniso[i] /= norm if norm > 0 else 1
		self.polyGyration[i] /= norm if norm > 0 else 1


	def Print(self):
		np.savetxt(self.anisoFile, self.polyAniso)
		np.savetxt(self.gyrationFile, self.polyGyration)
		
		print("\033[1;32mPrinted domain gyration radii to '%s'\033[0m" % self.gyrationFile)
		print("\033[1;32mPrinted domain anisotropy factors to '%s'\033[0m" % self.anisoFile)

		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	gyr = PolyGyration(outputDir, initFrame=initFrame)

	gyr.Compute()
	gyr.Print()
