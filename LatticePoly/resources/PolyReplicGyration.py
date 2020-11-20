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
					
		self.anisoFileTempl = os.path.join(self.reader.outputDir, "polyAnisoTempl.res")
		self.gyrationFileTempl = os.path.join(self.reader.outputDir, "polyGyrationTempl.res")
		self.anisoFileRepl = os.path.join(self.reader.outputDir, "polyAnisoRepl.res")
		self.gyrationFileRepl = os.path.join(self.reader.outputDir, "polyGyrationRepl.res")

		if os.path.exists(self.anisoFileRepl) & os.path.exists(self.gyrationFileRepl):
			print("Files '%s' and '%s' already exist - aborting" % (self.anisoFileTemp, self.gyrationFileTempl,self.anisoFileRepl,self.gyrationFileRepl))
			sys.exit()


	def Compute(self):
		self.InitTadNumber=self.reader.nTad
		self.MaxTadNumber=self.reader.InitReader(self.reader.N)
		self.polyAnisoTempl = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
		self.polyGyrationTempl = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
		self.polyAnisoRepl = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
		self.polyGyrationRepl = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)

		for i in range(self.reader.N):
			self.ProcessFrame(i)
			
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))
				
				
	def ProcessFrame(self, i):
		data = next(self.reader)
		if(i<self.InitTadNumber):
			for id, d in enumerate(self.reader.domains):
				if(len(self.reader.domains)==1):
					d=np.arange(0,self.InitTadNumber)

				pos = data.polyPos[d]
				pos -= pos.mean(axis=0, keepdims=True)
				diag = np.linalg.svd(pos, compute_uv=False) / d.size**0.5
				r2_gyr = np.square(diag).sum(axis=-1)
				
				r_gyr = np.sqrt(r2_gyr)
				aniso = 3/2.*(diag**4).sum(axis=-1)/r2_gyr**2 - 1/2.
				
				self.polyAnisoTempl[i, id] = aniso
				self.polyGyrationTempl[i, id] = r_gyr
		else:
			for id, d in enumerate(self.reader.domains):
				if(len(self.reader.domains)==1):
					d=np.arange(self.InitTadNumber,self.MaxTadNumber)

				pos = data.polyPos[d]
				pos -= pos.mean(axis=0, keepdims=True)
				diag = np.linalg.svd(pos, compute_uv=False) / d.size**0.5
				r2_gyr = np.square(diag).sum(axis=-1)
				
				r_gyr = np.sqrt(r2_gyr)
				aniso = 3/2.*(diag**4).sum(axis=-1)/r2_gyr**2 - 1/2.
				
				self.polyAnisoRepl[i, id] = aniso
				self.polyGyrationTempl[i, id] = r_gyr
			
					

	def Print(self):
		np.savetxt(self.anisoFileTempl, self.polyAnisoTempl)
		np.savetxt(self.gyrationFileTempl, self.polyGyrationTempl)
		np.savetxt(self.anisoFileRepl, self.polyAnisoRepl)
		np.savetxt(self.gyrationFileRepl, self.polyGyrationRepl)
		
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
