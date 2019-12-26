import sys

import numpy as np

from msdTools import msdFFT
from vtkReader import vtkReader


class PolyMSD(vtkReader):

	def __init__(self, outputDir, frameInit):
		vtkReader.__init__(self, outputDir, frameInit)
	
		if frameInit > self.maxFramePoly:
			print("Trajectory too short for chosen initial frame (max. frame: %d)" % self.maxFramePoly)
		
			sys.exit()
		
		self.ReadBox(readPoly=True)
		
		self.homFile = "%s/msdPolyHom.res" % self.outputDir
		self.hetFile = "%s/msdPolyHet.res" % self.outputDir
			
		self.frameInit = frameInit
		self._N = self.maxFramePoly - frameInit + 1
		
		self.cumulDistHom = np.zeros(self._N, dtype=np.float32)
		self.cumulDistHet = np.zeros(self._N, dtype=np.float32)
		
		self.tadPosHist = np.zeros((self._N, 3), dtype=np.float32)


	def Compute(self):
		for idxTad in range(self.nLoc):
			self.ProcessTad(idxTad)
				
			if self.polyType[idxTad] == 0:
				self.cumulDistHom += msdFFT(self.tadPosHist)
					
			else:
				self.cumulDistHet += msdFFT(self.tadPosHist)
				
			if (idxTad+1) % 10 == 0:
				print("Processed %d TADs" % (idxTad+1))


	def ProcessTad(self, idxTad):
		self.frame = self.frameInit

		for i in range(self._N):
			self.ReadPolyFrame()

			tadPos = self.polyPos[idxTad]
			self.tadPosHist[i] = tadPos
			

	def Print(self):
		msdHom = self.cumulDistHom / max(1, self.nHom)
		msdHet = self.cumulDistHet / max(1, self.nHet)

		np.savetxt(self.homFile, msdHom)
		np.savetxt(self.hetFile, msdHet)
		
		print("\033[1;32mPrinted polymer MSDs to '%s' and '%s'\033[0m" % (self.homFile, self.hetFile))


	def PrintTad(self, idxTad):
		msd = msdFFT(self.tadPosHist)
		msdFile = self.outputDir + "/msd%04d.res" % idxTad
		
		np.savetxt(msdFile, msd)
		
		print("\033[1;32mPrinted MSD to '%s'\033[0m" % msdFile)
	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir Neq [idxTad]\033[0m" % sys.argv[0])

		sys.exit()

	outputDir = sys.argv[1]
	Neq = int(sys.argv[2])

	msd = PolyMSD(outputDir, frameInit=Neq)

	if len(sys.argv) == 3:
		msd.Compute()
		msd.Print()
		
	elif len(sys.argv) == 4:
		idxTad = int(sys.argv[3])
		
		msd.ProcessTad(idxTad)
		msd.PrintTad(idxTad)
