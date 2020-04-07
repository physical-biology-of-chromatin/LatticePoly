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


class ResidenceTime():
	
	def __init__(self, outputDir, initFrame, cutoff=0.1):
		self.reader = vtkReader(outputDir, initFrame, readLiq=True, readPoly=True, backInBox=True)
			
		self.cutoff = cutoff
		self.residFile = os.path.join(self.reader.outputDir, "residenceTime.pdf")
		
		fontdict = {'family':'serif', 'size':'12'}

		plt.rc('font', **fontdict)
		plt.rcParams.update({'mathtext.fontset':'cm',
							 'mathtext.rm':'serif'})
		

	def Compute(self):
		self.idContOld = -1

		self.contNum = np.zeros(self.reader.nLiq*self.reader.nTad, dtype=np.float32)
		self.contHist = np.zeros(self.reader.nLiq*self.reader.nTad, dtype=np.float32)
		
		for i in range(self.reader.N):
			self.ProcessFrame()
			
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))
				
				
	def ProcessFrame(self):
		data = next(self.reader)
		liqPolyDist = cdist(data.liqPos, data.polyPos)
		
		idSpin, idTad = np.nonzero(liqPolyDist < self.cutoff)
		idCont = data.nTad*idSpin + idTad
		
		newCont = np.setdiff1d(idCont, self.idContOld, assume_unique=True)
		commonCont = np.intersect1d(idCont, self.idContOld, assume_unique=True)
		
		self.idContOld = idCont.copy()

		self.contNum[newCont] += 1.
		self.contHist[commonCont] += 1.
		

	def Print(self):
		contNum = self.contNum.reshape((self.reader.nLiq, self.reader.nTad))
		contHist = self.contHist.reshape((self.reader.nLiq, self.reader.nTad))

		contHetNum = contNum[:, self.reader.polyType.astype(bool)]
		contHomNum = contNum[:,~self.reader.polyType.astype(bool)]

		contHetHist = contHist[:, self.reader.polyType.astype(bool)]
		contHomHist = contHist[:,~self.reader.polyType.astype(bool)]

		normHet = np.where(contHetNum > 0., contHetNum, 1.)
		normHom = np.where(contHomNum > 0., contHomNum, 1.)

		contHetAveTime = contHetHist / normHet
		contHomAveTime = contHomHist / normHom
		
		contHetAveTime = contHetAveTime[contHetAveTime > 0.]
		contHomAveTime = contHomAveTime[contHomAveTime > 0.]

		meanHet = contHetAveTime.mean()
		meanEuc = contHomAveTime.mean()

		print("Mean homochromatic contact time: %.3f MC frames" % meanEuc)
		print("Mean heterochromatic contact time: %.3f MC frames" % meanHet)
							 
		fig = plt.figure()

		plt.hist(contHetAveTime, bins=np.linspace(0.5, self.reader.N+0.5, num=self.reader.N+1), label=r'${\rm Heterochromatic}$')
		plt.hist(contHomAveTime, bins=np.linspace(0.5, self.reader.N+0.5, num=self.reader.N+1), alpha=0.5, label=r'${\rm Euchromatic}$')
		
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
