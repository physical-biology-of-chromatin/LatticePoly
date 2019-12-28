##
##  LiqMSD.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import sys

import numpy as np

from msdTools import msdFFT
from vtkReader import vtkReader


class LiqMSD(vtkReader):

	def __init__(self, outputDir, initFrame):
		vtkReader.__init__(self, outputDir)

		self.InitReader(initFrame, readLiq=True)
		
		self.msdFile = "%s/msdLiq.res" % self.outputDir
		
		self.liqPosInit = self.liqPos
		self.cumulDist = np.zeros(self.N, dtype=np.float32)


	def Compute(self, sizeMax=1e9):
		sizeTot = self.N*self.nLiq
		
		if sizeTot < sizeMax:
			posHist = self.GetHist()

			for idxSpin in range(self.nLiq):
				self.cumulDist += msdFFT(posHist[:, idxSpin])
										
				if (idxSpin+1) % 10 == 0:
					print("Processed %d out of %d spins" % (idxSpin+1, self.nLiq))
											
		else:
			print("Likely memory overflow - aborting")
			
			sys.exit()


	def GetHist(self):
		print("\033[1;36mBuilding position histogram...\033[0m")

		posHist = np.zeros((self.N, self.nLiq, 3), dtype=np.float32)
		
		for i in range(self.N):
			self.ReadLiqFrame()
						
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
