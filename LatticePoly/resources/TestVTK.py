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

# outputDir = "/home/adminppuel/Documents/cpp/data/Experiment11_Jll0.4_ldens0.001_2PRE_RELAX/Jlp/0.0000/Jlpp/0.0000/"
# initFrame = -1

# reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)

# posHist = np.zeros((reader.N, reader.nTad, 3), dtype=np.float32)
# for i in range(reader.N):
# 	data = next(reader)
# 	posHist[i] = data.polyPos

# polyType_idx_list = []
# for idxTad in range(reader.nTad):
# 	if reader.polyPainter[idxTad] == 1:
# 		polyType_idx_list.append(idxTad)



# def distance_box(posHist,frame,idx1,idx2):
# 	return(np.sqrt((posHist[frame,idx1,0]-posHist[frame,idx2,0])**2+(posHist[frame,idx1,1]-posHist[frame,idx2,1])**2+(posHist[frame,idx1,2]-posHist[frame,idx2,2])**2))

# print(distance_box(posHist,0,polyType_idx_list[0],polyType_idx_list[-1]))

class PolyMSD():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		
		self.msdHetFile = os.path.join(self.reader.outputDir, "polyHetMSD.res")
		self.msdHomFile = os.path.join(self.reader.outputDir, "polyHomMSD.res")
		self.msdPREFile = os.path.join(self.reader.outputDir, "polyPREMSD.res")

		if os.path.exists(self.msdHetFile) & os.path.exists(self.msdHomFile) & os.path.exists(self.msdPREFile):
			print("Files '%s', '%s' and '%s' already exist - aborting" % (self.msdHetFile, self.msdHomFile, self.msdPREFile))
			sys.exit()


	def Compute(self):
		vMem = psutil.virtual_memory()
		sizeTot = self.reader.N * self.reader.polyPos.nbytes
		
		if sizeTot < vMem.available:
			self.cumulDistHet = 0
			self.cumulDistHom = 0
			self.cumulDistPRE = 0

			posHist = self.ReadHist()
			
			for idxTad in range(self.reader.nTad):
				if self.reader.polyPainter[idxTad] == 1:
					self.cumulDistPRE += msdFFT(posHist[:, idxTad])
				# elif self.reader.polyType[idxTad] == 1:
				# 	self.cumulDistHet += msdFFT(posHist[:, idxTad])
				# else:
				# 	self.cumulDistHom += msdFFT(posHist[:, idxTad])
							
				if (idxTad+1) % 1000 == 0:
					print("Processed %d out of %d TADs" % (idxTad+1, self.reader.nTad))
					
		else:
			print("Memory overflow likely - reduce chosen number of frames")
			sys.exit()


	def ComputeTad(self, idxTad):
		tadPosHist = np.zeros((self.reader.N, 3), dtype=np.float32)

		for i in range(self.reader.N):
			data = next(self.reader)
			tadPosHist[i] = data.polyPos[idxTad]
			
		self.distTad = msdFFT(tadPosHist)
				

	def ReadHist(self):
		posHist = np.zeros((self.reader.N, self.reader.nTad, 3), dtype=np.float32)
		
		for i in range(self.reader.N):
			data = next(self.reader)
			posHist[i] = data.polyPos
			
		return posHist
	
	
	def Print(self):
		if np.count_nonzero(self.reader.polyPainter == 1) > 0:
			msdPRE = self.cumulDistPRE / np.count_nonzero(self.reader.polyPainter == 1) 
			np.savetxt(self.msdPREFile, msdPRE)
			
			print("\033[1;32mPrinted PRE MSDs to '%s'\033[0m" % self.msdPREFile)

		# if self.reader.nHet > 0:
		# 	msdHet = self.cumulDistHet /  self.reader.nHet
		# 	np.savetxt(self.msdHetFile, msdHet)
			
		# 	print("\033[1;32mPrinted heterochromatic MSDs to '%s'\033[0m" % self.msdHetFile)
			
		# if self.reader.nEuc > 0:
		# 	msdHom = self.cumulDistHom / self.reader.nEuc
		# 	np.savetxt(self.msdHomFile, msdHom)
			
		# 	print("\033[1;32mPrinted euchromatic MSDs to '%s'\033[0m" % self.msdHomFile)

	
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

	msd = PolyMSD(outputDir, initFrame=initFrame)

	if len(sys.argv) == 3:
		msd.Compute()
		msd.Print()
		
	elif len(sys.argv) == 4:
		idxTad = int(sys.argv[3])
	
		msd.ComputeTad(idxTad)
		msd.PrintTad(idxTad)
