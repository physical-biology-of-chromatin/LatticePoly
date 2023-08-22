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
from itertools import zip_longest




class ReplicPolyMSD():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		
		self.msdHetFile = os.path.join(self.reader.outputDir,time.time())+ "polyReplicHetMSD.res")
		self.msdHomFile = os.path.join(self.reader.outputDir, time.time())+"polyReplicHomMSD.res")
		
		if os.path.exists(self.msdHetFile) & os.path.exists(self.msdHomFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.msdHetFile, self.msdHomFile))
			sys.exit()


	def Compute(self):
		vMem = psutil.virtual_memory()
		sizeTot = 2*self.reader.N * self.reader.polyPos.nbytes
		
		if sizeTot < vMem.available:
			self.cumulDistHet = 0
			self.cumulDistHom = 0
			self.distTad=[]
		
		#self.reader.InitReader(self.reader.N-1)
		#	data = next(self.reader)
		#	FinalNumberMonomers=data.nTad
			


#	for idxTad in range(FinalNumberMonomers):
#				self.reader.InitReader(0)
#				self.ComputeTad(idxTad)

			posHist=self.ReadHist()
			StartingMonomer=0
			MaxNumberofMonomers=len(posHist[len(posHist)-1])
			print("\033[1;32mMax number of Monomer '%s'\033[0m" % MaxNumberofMonomers)

			for timestep in range(0,self.reader.N):
				if StartingMonomer!=MaxNumberofMonomers:
					monomer_at_timestep=len(posHist[timestep])
					for idxTad in range(StartingMonomer,monomer_at_timestep):
						tadPosHist = np.zeros((self.reader.N-timestep, 3), dtype=np.float32)
						for i in range(0,self.reader.N-timestep):
							tadPosHist[i] = posHist[timestep+i][idxTad]
						self.distTad.append(msdFFT(tadPosHist))
					StartingMonomer=len(posHist[timestep])
			TransposeDistTad = [list(filter(None,i)) for i in zip_longest(*self.distTad)]
			self.msdHom=np.zeros(self.reader.N, dtype=np.float32)
			for i in range(len(TransposeDistTad)):
				self.msdHom[i]=np.mean(TransposeDistTad[i])

	

	

							


	def ComputeTad(self, idxTad):
		tadPosHist = np.zeros((self.reader.N, 3), dtype=np.float32)
		
		for i in range(self.reader.N):
			data = next(self.reader)
			tadPosHist[i] = data.polyPos[idxTad]
			
			self.distTad = msdFFT(tadPosHist)

	
	def ReadHist(self):
		posHist = []
		
		for i in range(self.reader.N):
			data = next(self.reader)
			posHist.append(data.polyPos)
		
		return posHist
	
	
	def Print(self):
		if self.reader.nHet > 0:
			msdHet = self.cumulDistHet /  self.reader.nHet
			np.savetxt(self.msdHetFile, msdHet)
			
			print("\033[1;32mPrinted heterochromatic MSDs to '%s'\033[0m" % self.msdHetFile)
		
		if self.reader.nEuc > 0:
			np.savetxt(self.msdHomFile, self.msdHom)
			
			print("\033[1;32mPrinted euchromatic MSDs to '%s'\033[0m" % self.msdHomFile)


	def PrintTad(self, idxTad):
		msdFile = self.reader.outputDir + "/msdTad%05d.res" % idxTad
		np.savetxt(msdFile, self.distTad)
		
		print("\033[1;32mPrinted TAD MSD to '%s'\033[0m" % msdFile)


if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	msd = ReplicPolyMSD(outputDir, initFrame=initFrame)

if len(sys.argv) == 3:
	msd.Compute()
	msd.Print()
	
elif len(sys.argv) == 4:
	idxTad = int(sys.argv[3])
	
	msd.ComputeTad(idxTad)

