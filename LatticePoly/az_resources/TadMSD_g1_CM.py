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


class PolyMSD():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		
		self.msdHetFile = os.path.join(self.reader.outputDir, "g1_160kb_DD.res")
		self.msdHomFile = os.path.join(self.reader.outputDir, "g3_160kb_DD.res")

		if os.path.exists(self.msdHetFile) & os.path.exists(self.msdHomFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.msdHetFile, self.msdHomFile))
			sys.exit()


	def Compute(self, idxTad, idxTadend):
		posHist, posCM = self.ReadHist(idxTad, idxTadend)
		self.cumulDistHet = 0	
		self.msdCM = 0
		for i in range(idxTad, idxTadend):
			self.cumulDistHet += msdFFT(posHist[:, i-idxTad])
					
			if (i+1) % 100 == 0:
				print("Processed %d out of %d TADs" % (i+1-idxTad, idxTadend-idxTad))
		self.msdCM = msdFFT(posCM)			

	def ComputeTad(self, idxTad):
		tadPosHist = np.zeros((self.reader.N, 3), dtype=np.float32)

		for i in range(self.reader.N):
			data = next(self.reader)
			tadPosHist[i] = data.polyPos[idxTad]


		self.distTad = msdFFT(tadPosHist)
				

	def ReadHist(self, idxTad, idxTadend):
		posHist = np.zeros((self.reader.N,idxTadend-idxTad, 3), dtype=np.float32)
		posCM   = np.zeros((self.reader.N, 3), dtype=np.float32)
		for i in range(self.reader.N):
			data = next(self.reader)
			posHist[i] = data.polyPos[idxTad:idxTadend]
			
			for j in range(idxTad,idxTadend):
				posCM[i] += data.polyPos[j] 

		return posHist, posCM/(idxTadend - idxTad)
	
	
	def Print(self, idxTad, idxTadend):

			msdHet = self.cumulDistHet /  (idxTadend -idxTad)
			np.savetxt(self.msdHetFile, msdHet)
			np.savetxt(self.msdHomFile, self.msdCM)

			print("\033[1;32mPrinted domain monomer MSDs to '%s'\033[0m" % self.msdHetFile)
			print("\033[1;32mPrinted domain CM MSDs to '%s'\033[0m" % self.msdHomFile)		
			

	
	def PrintTad(self, idxTad):
		msdFile = self.reader.outputDir + "/msdTad%05d.res" % idxTad
		np.savetxt(msdFile, self.distTad)
		
		print("\033[1;32mPrinted TAD MSD to '%s'\033[0m" % msdFile)
	
	
if __name__ == "__main__":
	if len(sys.argv) not in [5]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	msd = PolyMSD(outputDir, initFrame=initFrame)
	
	idxTad = int(sys.argv[3])
	idxTadend = int(sys.argv[4])

	msd.Compute(idxTad,idxTadend)
	msd.Print(idxTad,idxTadend)
