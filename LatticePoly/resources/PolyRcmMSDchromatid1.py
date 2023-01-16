##
##  LatticePoly
##
##  Created by ddasaro on the model of mtortora script.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##


import os
import sys
import psutil

import numpy as np

from utils import msdFFT
from vtkReader import vtkReader
import time


class PolyRcmMSDchromatid1():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		self.Nchain=0
		for t in range(self.reader.nTad):
				if(self.reader.status[t]==-1 or self.reader.status[t]==0):
					self.Nchain+=1
		self.msdHomFile = os.path.join(self.reader.outputDir,str(time.time())+ "polyRcmMSDchromatid1.res")

		if os.path.exists(self.msdHomFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.msdHomFile))
			sys.exit()


	def Compute(self):
		vMem = psutil.virtual_memory()
		sizeTot = self.reader.N * self.reader.polyPos.nbytes
		
		if sizeTot < vMem.available:
			self.cumulDistHet = 0
			self.cumulDistHom = 0

			posHist = self.ReadHist()

			Rg=np.zeros((len(posHist),3), dtype=np.float32)
			for i in range(len(posHist)):
				Rg[i]=sum(np.array(posHist[i]))/len(posHist[0])
		
			self.cumulDistHom= msdFFT(Rg)
							

		else:
			print("Memory overflow - reduce chosen number of frames")
			sys.exit()


	def ComputeTad(self, idxTad):
		tadPosHist = np.zeros((self.reader.N, 3), dtype=np.float32)

		for i in range(self.reader.N):
			data = next(self.reader)
			tadPosHist[i] = data.polyPos[idxTad]
			
		self.distTad = msdFFT(tadPosHist)
				

	def ReadHist(self):
		posHist = np.zeros((self.reader.N, self.Nchain, 3), dtype=np.float32)
		
		for i in range(self.reader.N):
			data = next(self.reader)
			posHist[i] = data.polyPos[:self.Nchain]
			
		return posHist
	
	
	def Print(self):

		msdHom = self.cumulDistHom
		np.savetxt(self.msdHomFile, msdHom)
		
		print("\033[1;32mPrinted msdrcm1 to '%s'\033[0m" % self.msdHomFile)

	
	def PrintTad(self, idxTad):
		msdFile = self.reader.outputDir + "/msdTad%05d.res" % idxTad
		np.savetxt(msdFile, self.distTad)
		
		print("\033[1;32mPrinted TAD msdrcm1 to '%s'\033[0m" % msdFile)
	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	msdrcm1 = PolyRcmMSDchromatid1(outputDir, initFrame=initFrame)

	if len(sys.argv) == 3:
		msdrcm1.Compute()
		msdrcm1.Print()
		
	elif len(sys.argv) == 4:
		idxTad = int(sys.argv[3])
	
		msdrcm1.ComputeTad(idxTad)
		msdrcm1.PrintTad(idxTad)
