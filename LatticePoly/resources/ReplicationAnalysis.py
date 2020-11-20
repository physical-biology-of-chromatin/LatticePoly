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


class ReplicationAnalysis():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		self.posHist = []
		self.ForkPos =[]
		self.SisterID=[]
		self.Status=[]
		self.Rg1=[]
		self.Rg2=[]
		
		self.Nchain=0
		self.Origin=50
		
		self.ForkDistanceFile = os.path.join(self.reader.outputDir, "polyForkDistance.res")
		self.OriginDistanceFile = os.path.join(self.reader.outputDir, "polyOriginDistance.res")
		self.ForkMSDFile = os.path.join(self.reader.outputDir, "polyForkMSD.res")
		self.OriginMSDFile = os.path.join(self.reader.outputDir, "polyOriginMSD.res")
		self.OriginMatrixFile = os.path.join(self.reader.outputDir, "Matrix")
		self.Gyr1File = os.path.join(self.reader.outputDir, "Gyr1.res")
		self.Gyr2File = os.path.join(self.reader.outputDir, "Gyr2.res")





#if os.path.exists(self.ForkDistanceFile) & os.path.exists(self.OriginDistanceFile):
#			print("Files '%s' and '%s' already exist - aborting" % (self.ForkDistanceFile, self.OriginDistanceFile))
#			sys.exit()






	def ComputeForkDistance(self):
		self.ForkDistance=[]
		self.ForkVector=[]
		for step in range(self.reader.N):
			active=False
			for i in range(len(self.ForkPos[0])):
				if(self.ForkPos[step][i]==-1):
					r1=self.posHist[step][i]
					active=True
				if(self.ForkPos[step][i]==+1 and active==True ):
					diff=self.posHist[step][i]-r1
					self.ForkDistance.append(np.sqrt(np.dot(diff.T,diff)))
					self.ForkVector.append(diff)

	def ComputeOriginDistance(self):
		self.OriginDistance=[]
		self.OriginVector=[]

		for step in range(self.reader.N):
			if(self.SisterID[step][self.Origin]==-1):
				self.OriginVector.append(np.zeros(3))
			else:
				diff=self.posHist[step][self.Origin]-self.posHist[step][self.SisterID[step][self.Origin]]
				self.OriginDistance.append(np.sqrt(np.dot(diff.T,diff)))
				self.OriginVector.append(diff)


	def ComputeOriginMSD(self):
		self.OriginMSD = msdFFT(np.array(self.OriginVector))
	
	
	def ComputeForkMSD(self):
		self.ForkMSD = msdFFT(np.array(self.ForkVector))
	
	def DistanceMatrix(self,step):
		DistanceMatrix=np.zeros((self.Nchain,self.Nchain))
		for i in range(self.Nchain):
			sister=self.SisterID[step][i]
			if(sister!=-1):
				for j in range(self.Nchain):
					diff=self.posHist[step][j]-self.posHist[step][sister]
					DistanceMatrix[i][j]=np.sqrt(np.dot(diff.T,diff))
					DistanceMatrix[j][i]=np.sqrt(np.dot(diff.T,diff))
		return DistanceMatrix
	
	def PrintMatrices(self):
		for step in range(self.reader.N):
			np.savetxt(self.OriginMatrixFile+str(step)+'.txt',self.DistanceMatrix(step))


	def CromatidGyration(self):
		for step in range(self.reader.N):
			self.Gyration1pos=[]
			self.Gyration2pos=[]
			for tad in range(len(self.posHist[step])):
				if(self.Status[step][tad]==1):
					self.Gyration1pos.append(self.posHist[step][tad])
				if(self.Status[step][tad]==-1):
					self.Gyration2pos.append(self.posHist[step][tad])

			if(len(self.Gyration1pos)>0):
				self.Gyration1pos=np.array(self.Gyration1pos)
				self.Gyration1pos-= self.Gyration1pos.mean(axis=0, keepdims=True)
				diag1 = np.linalg.svd(self.Gyration1pos, compute_uv=False) / len(self.Gyration1pos)**0.5
				r2_gyr1 = np.square(diag1).sum(axis=-1)
				r_gyr1 = np.sqrt(r2_gyr1)
				self.Rg1.append(r_gyr1)
			if(len(self.Gyration2pos)>0):
				self.Gyration2pos=np.array(self.Gyration2pos)
				self.Gyration2pos-= self.Gyration2pos.mean(axis=0, keepdims=True)
				diag2 = np.linalg.svd(self.Gyration2pos, compute_uv=False) / len(self.Gyration2pos)**0.5
				r2_gyr2 = np.square(diag2).sum(axis=-1)
				r_gyr2 = np.sqrt(r2_gyr2)
				self.Rg2.append(r_gyr2)


	
	


	def ReadHist(self):
		for i in range(self.reader.N):
			data = next(self.reader)
			self.posHist.append(data.polyPos)
			self.ForkPos.append(data.Forks)
			self.SisterID.append(data.SisterID)
			self.Status.append(data.Status)
			if (i==0):
				for t in range(self.reader.nTad):
					if(self.reader.Status[t]==-1 or self.reader.Status[t]==0):
						self.Nchain+=1
			

	def PrintGyration(self):
		np.savetxt(self.Gyr1File,self.Rg1)
		np.savetxt(self.Gyr2File,self.Rg2)



	
	
	def Print(self):
		
		np.savetxt(self.ForkDistanceFile, self.ForkDistance)
		print("\033[1;32mPrinted ForkDistance to '%s'\033[0m" % self.ForkDistanceFile)

		np.savetxt(self.OriginDistanceFile, self.OriginDistance)
		print("\033[1;32mPrinted OriginDistance to '%s'\033[0m" % self.OriginDistanceFile)

		np.savetxt(self.ForkMSDFile, self.ForkMSD)
		print("\033[1;32mPrinted ForkMSD to '%s'\033[0m" % self.ForkMSDFile)
	
		np.savetxt(self.OriginMSDFile, self.OriginMSD)
		print("\033[1;32mPrinted OriginMSD to '%s'\033[0m" % self.OriginMSDFile)

	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	ReplicationAnalysis = ReplicationAnalysis(outputDir, initFrame=initFrame)



	if len(sys.argv) == 3:
		ReplicationAnalysis.ReadHist()
		ReplicationAnalysis.ComputeOriginDistance()
		ReplicationAnalysis.ComputeForkDistance()
		ReplicationAnalysis.ComputeOriginMSD()
		ReplicationAnalysis.ComputeForkMSD()
		ReplicationAnalysis.Print()
		#ReplicationAnalysis.CromatidGyration()
		#ReplicationAnalysis.PrintGyration()


		
	elif len(sys.argv) == 4:
		idxTad = int(sys.argv[3])
	
		ReplicationAnalysis.ComputeTad(idxTad)
		ReplicationAnalysis.PrintTad(idxTad)
