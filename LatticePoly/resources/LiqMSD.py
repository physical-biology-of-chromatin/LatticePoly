import sys

import numpy as np
import fileseq as fs

from msdTools import msdFFT
from vtkReader import vtkReader


class LiqMSD(vtkReader):

	def __init__(self, outputDir, frameInit):
		vtkReader.__init__(self, outputDir, frameInit)
	
		self.ReadBox(readLiq=True)
		
		self.msdFile = "%s/msdLiq.res" % self.outputDir

		frameFinal = fs.findSequenceOnDisk(self.outputDir + '/liq@.vtp').end()
		self.frameInit = frameInit
		
		self.N = frameFinal - frameInit + 1
		
		self.liqPosInit = self.liqPos
		
		self.cumulDist = np.zeros(self.N, dtype=np.float32)
		self.spinPosHist = np.zeros((self.N, 3), dtype=np.float32)


	def Compute(self):
		for idxSpin in range(self.nLiq):
			self.ProcessSpin(idxSpin)

			self.cumulDist += msdFFT(self.spinPosHist)
									
			if (idxSpin+1) % 10 == 0:
				print("Processed %d spins" % (idxSpin+1))


	def ProcessSpin(self, idxSpin):
		self.frame = self.frameInit
		
		while True:
			try:
				self.ReadLiqFrame()
				
				i = self.frame - self.frameInit
				spinPos = self.liqPosInit[idxSpin] + self.liqDisp[idxSpin]
									
				self.spinPosHist[i] = spinPos
				self.frame += 1
				
			except IOError:
				break
				

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
