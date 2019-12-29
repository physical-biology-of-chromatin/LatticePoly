##
##  LiqMSD.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import sys
import psutil

import numpy as np

from utils import msdFFT
from vtkReader import vtkReader


class LiqMSD(vtkReader):

	def __init__(self, outputDir, initFrame):
		vtkReader.__init__(self, outputDir)

		self.InitReader(initFrame, readLiq=True)
		
		self.msdFile = "%s/liqMSD.res" % self.outputDir


	def Compute(self):
		vMem = psutil.virtual_memory()
		sizeTot = self.N * self.liqPos.nbytes
		
		if sizeTot < vMem.available:
			self.cumulDist = 0
			self.liqPosInit = self.liqPos
			
			posHist = self.GetHist()

			for idxSpin in range(self.nLiq):
				self.cumulDist += msdFFT(posHist[:, idxSpin])
										
				if (idxSpin+1) % 10 == 0:
					print("Processed %d out of %d spins" % (idxSpin+1, self.nLiq))
											
		else:
			print("Memory overflow - reduce chosen number of frames")
			
			sys.exit()


	def GetHist(self):
		posHist = np.zeros((self.N, self.nLiq, 3), dtype=np.float32)
		
		for i in range(self.N):
			self.ReadLiqFrame(readAttr=False)
						
			posHist[i] = self.liqPosInit + self.liqDisp
			self.frame += 1
			
		return posHist
		

	def Print(self):
		msdLiq = self.cumulDist / self.nLiq

		np.savetxt(self.msdFile, msdLiq)
		
		print("\033[1;32mPrinted liquid MSDs to '%s'\033[0m" % self.msdFile)
	
	
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])

		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	msd = LiqMSD(outputDir, initFrame=initFrame)

	msd.Compute()
	msd.Print()
