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

class Gyration_moments():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
			
		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1		
		self.lambdaxFile = os.path.join(self.reader.outputDir, "lambdax.res")
		self.lambdayFile = os.path.join(self.reader.outputDir, "lambday.res")
		self.lambdazFile = os.path.join(self.reader.outputDir, "lambdaz.res")


		if os.path.exists(self.lambdaxFile) & os.path.exists(self.lambdayFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.anisoFile, self.gyrationFile))
			sys.exit()


	def Compute(self):
		self.lambdax = np.zeros((self.reader.N,2))
		self.lambday = np.zeros((self.reader.N,2))
		self.lambdaz = np.zeros((self.reader.N,2))

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
			
			
		self.lambdax[i][0] += data.nTad- self.Nchain
		self.lambday[i][0] += data.nTad- self.Nchain
		self.lambdaz[i][0] += data.nTad- self.Nchain
		
		self.lambdax[i][1] += diag[0]
		self.lambday[i][1] += diag[1]
		self.lambdaz[i][1] += diag[2]

									


	def Print(self):
		np.savetxt(self.lambdaxFile, self.lambdax)
		np.savetxt(self.lambdayFile, self.lambday)
		np.savetxt(self.lambdazFile, self.lambdaz)
		


		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	lamdamoments = Gyration_moments(outputDir, initFrame=initFrame)

	lamdamoments.Compute()
	lamdamoments.Print()
