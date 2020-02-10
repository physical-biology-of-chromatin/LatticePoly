##
##  PolyMSD.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys
import psutil

import numpy as np

from utils import msdFFT
from vtkReader import vtkReader


class PolyMSD(vtkReader):

	def __init__(self, outputDir, initFrame):
		vtkReader.__init__(self, outputDir, initFrame, readPoly=True)
		
		self.msdHetFile = os.path.join(self.outputDir, "polyHetMSD.res")
		self.msdHomFile = os.path.join(self.outputDir, "polyHomMSD.res")


	def Compute(self):
		vMem = psutil.virtual_memory()
		sizeTot = self.N * self.polyPos.nbytes
		
		if sizeTot < vMem.available:
			self.cumulDistHet = 0
			self.cumulDistHom = 0
		
			posHist = self.ReadHist()
			
			for idxTad in range(self.nLoc):
				if self.polyType[idxTad] == 1:
					self.cumulDistHet += msdFFT(posHist[:, idxTad])
								
				else:
					self.cumulDistHom += msdFFT(posHist[:, idxTad])
							
				if (idxTad+1) % 10 == 0:
					print("Processed %d out of %d TADs" % (idxTad+1, self.nLoc))
					
		else:
			print("Memory overflow - reduce chosen number of frames")
			
			sys.exit()


	def ComputeTad(self, idxTad):
		tadPosHist = np.zeros((self.N, 3), dtype=np.float32)

		for i in range(self.N):
			self.ReadPolyFrame(readAttr=False)

			tadPosHist[i] = self.polyPos[idxTad]
			self.frame += 1
			
		self.distTad = msdFFT(tadPosHist)
				

	def ReadHist(self):
		posHist = np.zeros((self.N, self.nLoc, 3), dtype=np.float32)
		
		for i in range(self.N):
			self.ReadPolyFrame()

			posHist[i] = self.polyPos
			self.frame += 1
			
		return posHist
	
	
	def Print(self):
		if self.nHet > 0:
			msdHet = self.cumulDistHet /  self.nHet
			
			np.savetxt(self.msdHetFile, msdHet)
			print("\033[1;32mPrinted heterochromatic MSDs to '%s'\033[0m" % self.msdHetFile)
			
		if self.nHom > 0:
			msdHom = self.cumulDistHom / self.nHom
		
			np.savetxt(self.msdHomFile, msdHom)
			print("\033[1;32mPrinted euchromatic MSDs to '%s'\033[0m" % self.msdHomFile)

	
	def PrintTad(self, idxTad):
		msdFile = self.outputDir + "/msdTad%05d.res" % idxTad
		np.savetxt(msdFile, self.distTad)
				
		print("\033[1;32mPrinted TAD MSD to '%s'\033[0m" % msdFile)
	
	
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
	
		msd.ComputeTad(idxTad)
		msd.PrintTad(idxTad)
