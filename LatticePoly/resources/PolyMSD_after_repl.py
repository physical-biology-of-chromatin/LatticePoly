##
##  PolyMSD.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys
import psutil

import numpy as np

from utils import msdFFT
from vtkReader import vtkReader
import time

class PolyMSD():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		
		self.msdHetFile = os.path.join(self.reader.outputDir, "polyHetMSD.res")
		self.msdHomFile = os.path.join(self.reader.outputDir, str(time.time())+"polyHomMSD.res")

		if os.path.exists(self.msdHetFile) & os.path.exists(self.msdHomFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.msdHetFile, self.msdHomFile))
			sys.exit()
		
		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1
		self.timepoint=0 #find frame where replication reach the desired percentage
		
		fullyreplicated=False
		for i in range(self.reader.N):
			if(self.reader.nTad<2*self.Nchain):
				next(self.reader)
				self.timepoint+=1
			if(self.reader.nTad==2*self.Nchain):
				fullyreplicated=True
				break
	
	
		
		if(fullyreplicated==False):
			print("Chromosome not fully replicated ")
			sys.exit()
			
			

	def Compute(self):
			self.cumulDistHet = 0
			self.cumulDistHom = 0
			posHist = self.ReadHist()
			for idxTad in range(self.reader.nTad):
				if self.reader.polyType[idxTad] == 1:
					self.cumulDistHet += msdFFT(posHist[:, idxTad])
				else:
					self.cumulDistHom += msdFFT(posHist[:, idxTad])
							
				if (idxTad+1) % 1000 == 0:
					print("Processed %d out of %d TADs" % (idxTad+1, self.reader.nTad))
					


	
	def ReadHist(self):
		posHist = np.zeros((self.reader.N-self.timepoint, self.reader.nTad, 3), dtype=np.float32)
		
		for i in range(self.reader.N):
			data = next(self.reader)
			posHist[i] = data.polyPos
			
		return posHist
	
	
	def Print(self):
		if self.reader.nHet > 0:
			msdHet = self.cumulDistHet /  self.reader.nHet
			np.savetxt(self.msdHetFile, msdHet)
			
			print("\033[1;32mPrinted heterochromatic MSDs to '%s'\033[0m" % self.msdHetFile)
			
		if self.reader.nEuc > 0:
			msdHom = self.cumulDistHom / self.reader.nEuc
			np.savetxt(self.msdHomFile, msdHom)
			
			print("\033[1;32mPrinted euchromatic MSDs to '%s'\033[0m" % self.msdHomFile)

	

	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	msd = PolyMSD(outputDir, initFrame=initFrame)



	if len(sys.argv) == 3:
		msd.Compute()
		msd.Print()
		
