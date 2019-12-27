##
##  PolyMSD.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import sys

import numpy as np

from msdTools import msdFFT
from vtkReader import vtkReader


class PolyMSD(vtkReader):

	def __init__(self, outputDir, initFrame):
		vtkReader.__init__(self, outputDir)
		
		self.InitReader(initFrame, readPoly=True)
		
		self.homFile = "%s/msdPolyHom.res" % self.outputDir
		self.hetFile = "%s/msdPolyHet.res" % self.outputDir
			
		self.cumulDistHom = np.zeros(self.N, dtype=np.float32)
		self.cumulDistHet = np.zeros(self.N, dtype=np.float32)
		
		self.tadPosHist = np.zeros((self.N, 3), dtype=np.float32)


	def Compute(self):
		for idxTad in range(self.nLoc):
			self.ProcessTad(idxTad)
				
			if self.polyType[idxTad] == 0:
				self.cumulDistHom += msdFFT(self.tadPosHist)
					
			else:
				self.cumulDistHet += msdFFT(self.tadPosHist)
				
			if (idxTad+1) % 10 == 0:
				print("Processed %d out of %d TADs" % (idxTad+1, self.nLoc))


	def ProcessTad(self, idxTad):
		self.frame = self.initFrame

		for i in range(self.N):
			self.ReadPolyFrame()

			self.tadPosHist[i] = self.polyPos[idxTad]			
			self.frame += 1
			

	def Print(self):
		msdHom = self.cumulDistHom / max(1, self.nHom)
		msdHet = self.cumulDistHet / max(1, self.nHet)

		np.savetxt(self.homFile, msdHom)
		np.savetxt(self.hetFile, msdHet)
		
		print("\033[1;32mPrinted polymer MSDs to '%s' and '%s'\033[0m" % (self.homFile, self.hetFile))


	def PrintTad(self, idxTad):
		msd = msdFFT(self.tadPosHist)
		msdFile = self.outputDir + "/msdTAD%05d.res" % idxTad
		
		np.savetxt(msdFile, msd)
		
		print("\033[1;32mPrinted MSD to '%s'\033[0m" % msdFile)
	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])

		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	msd = PolyMSD(outputDir, initFrame=initFrame)

	if len(sys.argv) == 3:
		msd.Compute()
		msd.Print()
		
	elif len(sys.argv) == 4:
		idxTad = int(sys.argv[3])
		
		msd.ProcessTad(idxTad)
		msd.PrintTad(idxTad)
