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
import time

class PolyGyration():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
			
		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1		
		self.anisoFile = os.path.join(self.reader.outputDir, str(time.time())+"polyAniso_chr1.res")
		self.gyrationFile = os.path.join(self.reader.outputDir, str(time.time())+"polyGyration_chr1.res")

		if os.path.exists(self.anisoFile) & os.path.exists(self.gyrationFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.anisoFile, self.gyrationFile))
			sys.exit()


	def Compute(self):
		self.polyAniso = np.zeros((self.reader.N,2))
		self.polyGyration = np.zeros((self.reader.N,2))

		for i in range(self.reader.N):
			self.ProcessFrame(i)
			
			#if (i+1) % 10 == 0:
			#	print("Processed %d out of %d configurations" % (i+1, self.reader.N))
				
				
	def ProcessFrame(self, i):
		data = next(self.reader)

		
		
		pos = data.polyPos[:self.Nchain]
		pos -= pos.mean(axis=0, keepdims=True)
		pos=pos/self.Nchain**0.5
		diag = np.linalg.svd(pos, compute_uv=False)
		r2_gyr = np.square(diag).sum(axis=-1)
			
		r_gyr = np.sqrt(r2_gyr)
		aniso = 3/2.*(diag**4).sum(axis=-1)/r2_gyr**2 - 1/2.
			
		self.polyAniso[i][0] += data.nTad- self.Nchain
		self.polyGyration[i][0] += data.nTad- self.Nchain 
		self.polyAniso[i][1] += aniso 
		self.polyGyration[i][1] += r_gyr 

									


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
