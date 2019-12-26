import sys

import numpy as np

from msdTools import msdFFT
from vtkReader import vtkReader


class LiqMSD(vtkReader):

	def __init__(self, outputDir, frameInit):
		vtkReader.__init__(self, outputDir, frameInit)
	
		if frameInit > self.maxFrameLiq:
			print("Trajectory too short for chosen initial frame (max. frame: %d)" % self.maxFrameLiq)
		
			sys.exit()
		
		self.ReadBox(readLiq=True)
		
		self.msdFile = "%s/msdLiq.res" % self.outputDir
		
		self.frameInit = frameInit
		self._N = self.maxFrameLiq - frameInit + 1
		
		self.liqPosInit = self.liqPos
		
		self.cumulDist = np.zeros(self._N, dtype=np.float32)
		self.spinPosHist = np.zeros((self._N, 3), dtype=np.float32)


	def Compute(self):
		for idxSpin in range(self.nLiq):
			self.ProcessSpin(idxSpin)

			self.cumulDist += msdFFT(self.spinPosHist)
									
			if (idxSpin+1) % 10 == 0:
				print("Processed %d spins" % (idxSpin+1))


	def ProcessSpin(self, idxSpin):
		self.frame = self.frameInit
		
		for i in range(self._N):
			self.ReadLiqFrame()
				
			spinPos = self.liqPosInit[idxSpin] + self.liqDisp[idxSpin]
			self.spinPosHist[i] = spinPos
			
			self.frame += 1
				

	def Print(self):
		msdLiq = self.cumulDist / self.nLiq

		np.savetxt(self.msdFile, msdLiq)
		
		print("\033[1;32mPrinted liquid MSDs to '%s'\033[0m" % self.msdFile)
	
	
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir Neq\033[0m" % sys.argv[0])

		sys.exit()

	outputDir = sys.argv[1]
	Neq = int(sys.argv[2])

	msd = LiqMSD(outputDir, frameInit=Neq)

	msd.Compute()
	msd.Print()
