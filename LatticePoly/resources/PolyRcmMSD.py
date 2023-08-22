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

class polyRcmMSD():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		
		self.msdHomFile = os.path.join(self.reader.outputDir, str(time.time())+"polyRcmMSD.res")

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
		posHist = np.zeros((self.reader.N, self.reader.nTad, 3), dtype=np.float32)
		
		for i in range(self.reader.N):
			data = next(self.reader)
			posHist[i] = data.polyPos
			
		return posHist
	
	
	def Print(self):

		msdHom = self.cumulDistHom
		np.savetxt(self.msdHomFile, msdHom)
		
		print("\033[1;32mPrinted msdrg to '%s'\033[0m" % self.msdHomFile)

	
	def PrintTad(self, idxTad):
		msdFile = self.reader.outputDir + "/msdTad%05d.res" % idxTad
		np.savetxt(msdFile, self.distTad)
		
		print("\033[1;32mPrinted TAD msdrg to '%s'\033[0m" % msdFile)
	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	msdrcm = polyRcmMSD(outputDir, initFrame=initFrame)

	if len(sys.argv) == 3:
		msdrcm.Compute()
		msdrcm.Print()
		
	elif len(sys.argv) == 4:
		idxTad = int(sys.argv[3])
	
		msdrcm.ComputeTad(idxTad)
		msdrcm.PrintTad(idxTad)
