##
##  ResidenceTime.py
##  LatticePoly
##
##  Created by mtortora on 17/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

from vtkReader import vtkReader
from scipy.spatial.distance import cdist


class ResidenceTime(vtkReader):
	
	def __init__(self, outputDir, initFrame, cutoff=0.1):
		vtkReader.__init__(self, outputDir)
			
		self.InitReader(initFrame, readLiq=True, readPoly=True)

		self.cutoff = cutoff
		self.residFile = os.path.join(self.outputDir, "residenceTime.pdf")
		
		fontdict = {'family':'serif', 'size':'12'}

		plt.rc('font', **fontdict)
		plt.rcParams.update({'mathtext.fontset':'cm',
							 'mathtext.rm':'serif'})
		

	def Compute(self):
		self.idContOld = -1

		self.contNum = np.zeros(self.nLiq*self.nLoc, dtype=np.float32)
		self.contHist = np.zeros(self.nLiq*self.nLoc, dtype=np.float32)
		
		for _ in range(self.N):
			self.ProcessFrame()
			
			if (self.frame-self.initFrame) % 10 == 0:
				print("Processed %d out of %d configurations" % (self.frame-self.initFrame, self.N))
				
				
	def ProcessFrame(self):
		self.ReadLiqFrame(readAttr=False)
		self.ReadPolyFrame(readAttr=False, backInBox=True)

		liqPolyDist = cdist(self.liqPos, self.polyPos)
		
		idSpin, idTad = np.nonzero(liqPolyDist < self.cutoff)
		idCont = self.nLoc*idSpin + idTad
		
		newCont = np.setdiff1d(idCont, self.idContOld, assume_unique=True)
		commonCont = np.intersect1d(idCont, self.idContOld, assume_unique=True)
		
		self.idContOld = idCont.copy()

		self.contNum[newCont] += 1.
		self.contHist[commonCont] += 1.
				
		self.frame += 1


	def Print(self):
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
		
		contHetAveTime = contHetAveTime[contHetAveTime > 0.]
		contHomAveTime = contHomAveTime[contHomAveTime > 0.]

		meanHet = contHetAveTime.mean()
		meanHom = contHomAveTime.mean()

		print("Mean homochromatic contact time: %.3f MC frames" % meanHom)
		print("Mean heterochromatic contact time: %.3f MC frames" % meanHet)
							 
		fig = plt.figure()

		plt.hist(contHetAveTime, bins=np.linspace(0.5, self.N+0.5, num=self.N+1), label=r'${\rm Heterochromatic}$')
		plt.hist(contHomAveTime, bins=np.linspace(0.5, self.N+0.5, num=self.N+1), alpha=0.5, label=r'${\rm Euchromatic}$')
		
		plt.xlabel(r'$\tau_{\rm res}$', size=16)
		plt.legend(loc='upper right', fontsize=16)
		
		plt.savefig(self.residFile, format="pdf", transparent=True)
		print("\033[1;32mPrinted figure to '%s'\033[0m" % self.residFile)

		plt.show()
				
		
if __name__ == "__main__":
	if len(sys.argv) not in [3,4]:
		print("\033[1;31mUsage is %s outputDir initFrame [cutoff]\033[0m" % sys.argv[0])

		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	if len(sys.argv) == 3:
		res = ResidenceTime(outputDir, initFrame=initFrame)
		
	elif len(sys.argv) == 4:
		cutoff = float(sys.argv[3])
		res = ResidenceTime(outputDir, initFrame=initFrame, cutoff=cutoff)

	res.Compute()
	res.Print()
