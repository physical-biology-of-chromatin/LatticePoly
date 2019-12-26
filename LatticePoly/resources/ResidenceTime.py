import sys

import numpy as np
import matplotlib.pyplot as plt

from vtkReader import vtkReader
from scipy.spatial.distance import cdist


class ResidenceTime(vtkReader):
	
	def __init__(self, outputDir, frameInit, cutoff=0.1):
		vtkReader.__init__(self, outputDir, frameInit)
		
		frameMax = min(self.maxFrameLiq, self.maxFramePoly)
		
		if frameInit > frameMax:
			print("Trajectory too short for chosen initial frame (max. frame: %d)" % frameMax)
		
			sys.exit()
			
		self.ReadBox(readLiq=True, readPoly=True)

		self.cutoff = cutoff
		self.histFile = self.outputDir + "/contactHist.pdf"

		self.frameInit = frameInit
		self._N = frameMax - frameInit + 1
				
		self.contNum = np.zeros(self.nLiq*self.nLoc, dtype=np.float32)
		self.contHist = np.zeros(self.nLiq*self.nLoc, dtype=np.float32)
		

	def Compute(self):
		self.idContOld = -1

		for _ in range(self._N):
			self.ProcessFrame()
			
			if (self.frame-self.frameInit) % 10 == 0:
				print("Processed %d configurations" % (self.frame-self.frameInit))
				
				
	def ProcessFrame(self):
		self.ReadLiqFrame()
		self.ReadPolyFrame(backInBox=True)

		liqPolyDist = cdist(self.liqPos, self.polyPos)
		
		idSpin, idTad = np.nonzero(liqPolyDist < self.cutoff)
		idCont = self.nLoc*idSpin + idTad
		
		newCont = np.setdiff1d(idCont, self.idContOld, assume_unique=True)
		commonCont = np.intersect1d(idCont, self.idContOld, assume_unique=True)
		
		self.idContOld = idCont.copy()

		self.contNum[newCont] += 1.
		self.contHist[commonCont] += 1.
				
		self.frame += 1


	def Print(self, span=20):
		contNum = self.contNum.reshape((self.nLiq, self.nLoc))
		contHist = self.contHist.reshape((self.nLiq, self.nLoc))

		contHetNum = contNum[:, self.polyType.astype(bool)]
		contHomNum = contNum[:,~self.polyType.astype(bool)]

		contHetHist = contHist[:, self.polyType.astype(bool)]
		contHomHist = contHist[:,~self.polyType.astype(bool)]

		normHet = np.where(contHetNum > 0., contHetNum, 1.)
		normHom = np.where(contHomNum > 0., contHomNum, 1.)

		contHetAveTime = contHetHist / normHet
		contHomAveTime = contHomHist / normHom
		
		contHetAveTime = contHetAveTime[np.nonzero(contHetAveTime)].flatten()
		contHomAveTime = contHomAveTime[np.nonzero(contHomAveTime)].flatten()

		meanHet = contHetAveTime.mean()
		meanHom = contHomAveTime.mean()

		print("Mean homochromatic contact time: %.3f MC frames" % meanHom)
		print("Mean heterochromatic contact time: %.3f MC frames" % meanHet)
		
		fig = plt.figure()

		ctrHet = max(int(meanHet), 1)
		ctrHom = max(int(meanHom), 1)

		plt.hist(contHetAveTime, bins=np.linspace(0.5, span*ctrHet+0.5, num=span*ctrHet+1))
		plt.hist(contHomAveTime, bins=np.linspace(0.5, span*ctrHom+0.5, num=span*ctrHom+1), alpha=0.5)

		plt.savefig(self.histFile, format="pdf", transparent=True)
		print("\033[1;32mPrinted figure to '%s'\033[0m" % self.histFile)

		plt.show()
				
		
if __name__ == "__main__":
	if len(sys.argv) not in [3,4]:
		print("\033[1;31mUsage is %s outputDir Neq [cutoff]\033[0m" % sys.argv[0])

		sys.exit()

	outputDir = sys.argv[1]
	Neq = int(sys.argv[2])

	if len(sys.argv) == 3:
		res = ResidenceTime(outputDir, frameInit=Neq)
		
	elif len(sys.argv) == 4:
		cutoff = float(sys.argv[3])
		res = ResidenceTime(outputDir, frameInit=Neq, cutoff=cutoff)

	res.Compute()
	res.Print()
