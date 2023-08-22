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
import math
import time

class ReplicationAnalysis():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		self.posHist = []
		self.ForkPos =[]
		self.SisterID=[]
		self.Status=[]
		self.Rg1=[]
		self.Rg2=[]
		self.sisterdist=[]

		self.RcmDistance=[]
		
		self.Nchain=0
		self.Origin=500
		
		self.ForkDistanceFile = os.path.join(self.reader.outputDir, str(time.time())+"polyForkDistance.res")
		self.OriginDistanceFile = os.path.join(self.reader.outputDir,str(time.time())+ "polyOriginDistance.res")
		self.ForkMSDFile = os.path.join(self.reader.outputDir, "polyForkMSD.res")
		self.ForkMSDtoFile = os.path.join(self.reader.outputDir, "polyForkMSDto.res")
		self.OriginMSDFile = os.path.join(self.reader.outputDir, "polyOriginMSD.res")
		self.OriginMSDtoFile = os.path.join(self.reader.outputDir, "polyOriginMSDto.res")
		self.sisterdistFile = os.path.join(self.reader.outputDir, "sisterdist.res")
		self.sisterdisterrFile = os.path.join(self.reader.outputDir, "sisterdisterr.res")

		
		self.OriginMatrixcisFile = os.path.join(self.reader.outputDir, "Matrixcis")
		self.OriginMatrixtransFile = os.path.join(self.reader.outputDir, "Matrixtrans")
		self.OriginMatrixFile = os.path.join(self.reader.outputDir, "Matrix")
		self.Gyr1File = os.path.join(self.reader.outputDir,self.reader.outputDir,str(time.time())+ "Gyr1.res")
		self.Gyr2File = os.path.join(self.reader.outputDir,self.reader.outputDir,str(time.time())+ "Gyr2.res")
		self.RcmDistanceFile = os.path.join(self.reader.outputDir, "RcmDistance.res")





#if os.path.exists(self.ForkDistanceFile) & os.path.exists(self.OriginDistanceFile):
#			print("Files '%s' and '%s' already exist - aborting" % (self.ForkDistanceFile, self.OriginDistanceFile))
#			sys.exit()






	def ComputeForkDistance(self):
		self.ForkDistance=[]
		self.ForkVector=[]
		self.OriginDistance=[]

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
					diff=self.posHist[step][self.Origin]-self.posHist[step][self.SisterID[step][self.Origin]]
					self.OriginDistance.append(np.sqrt(np.dot(diff.T,diff)))

		np.savetxt(self.ForkDistanceFile, self.ForkDistance)
		print("\033[1;32mPrinted ForkDistance to '%s'\033[0m" % self.ForkDistanceFile)
		np.savetxt(self.OriginDistanceFile, self.OriginDistance)
		print("\033[1;32mPrinted OriginsDistance to '%s'\033[0m" % self.OriginDistanceFile)


	def ComputeOriginDistance(self):
		self.OriginDistance=[]
		self.OriginVector=[]

		for step in range(self.reader.N):
			if(self.SisterID[step][self.Origin]!=-1 and self.SisterID[step][self.Origin]!=self.Origin and len(self.posHist[step][self.Origin])<=self.Nchain*2):
				diff=self.posHist[step][self.Origin]-self.posHist[step][self.SisterID[step][self.Origin]]
				self.OriginDistance.append(np.sqrt(np.dot(diff.T,diff)))
				self.OriginVector.append(diff)
		np.savetxt(self.OriginDistanceFile, self.OriginDistance)
		print("\033[1;32mPrinted OriginsDistance to '%s'\033[0m" % self.OriginDistanceFile)


	def ComputeOriginMSDto(self,i):
		self.OriginMSDto=[]
		for step in range(i,len(self.OriginVector)):
			diff=self.OriginVector[step]-self.OriginVector[i]
			self.OriginMSDto.append(np.sqrt(np.dot(diff.T,diff)))
		return self.OriginMSDto

	def ComputeForkMSDto(self,i):
		self.ForkMSDto=[]
		for step in range(i,len(self.ForkVector)):
			diff=self.ForkVector[step]-self.ForkVector[i]
			self.ForkMSDto.append(np.sqrt(np.dot(diff.T,diff)))
		return self.ForkMSDto

	def ComputepointMSDto(self,i):
		self.PointMSDto=[]
		for step in range(i,self.reader.N):
			diff=self.posHist[step][50]-self.posHist[i][50]
			self.PointMSDto.append(np.sqrt(np.dot(diff.T,diff)))
		return self.PointMSDto

	def ComputedistanceMSDto(self,i):
		self.ForkMSDto=[]
		for step in range(i,self.reader.N):
			diff=(self.posHist[step][20]-self.posHist[step][50])-(self.posHist[i][20]-self.posHist[i][50])
			self.ForkMSDto.append(np.sqrt(np.dot(diff.T,diff)))
		return self.ForkMSDto

	
	def ComputeOriginMSD(self):
		self.OriginMSD = msdFFT(np.array(self.OriginVector))
	
	
	def ComputeForkMSD(self):
		self.ForkMSD = msdFFT(np.array(self.ForkVector))
	
	def DistanceMatrixtrans(self,step):
		DistanceMatrix=np.zeros((self.Nchain,self.Nchain))
		for i in range(self.Nchain):
			sister=self.SisterID[step][i]
			if(sister!=-1):
				for j in range(self.Nchain):
					if(self.Status[step][j]!=0):
						diff=self.posHist[step][j]-self.posHist[step][sister]
						DistanceMatrix[i][j]=np.sqrt(np.dot(diff.T,diff))
						DistanceMatrix[j][i]=np.sqrt(np.dot(diff.T,diff))
		return DistanceMatrix
					
	def DistanceMatrixcis(self,step):
		DistanceMatrix=np.zeros((self.Nchain,self.Nchain))
		for i in range(self.Nchain):
			sister=self.SisterID[step][i]
			if(sister!=-1):
				for j in range(self.Nchain):
					if(self.Status[step][j]!=0):
						diff=self.posHist[step][j]-self.posHist[step][i]
						DistanceMatrix[i][j]=np.sqrt(np.dot(diff.T,diff))
						DistanceMatrix[j][i]=np.sqrt(np.dot(diff.T,diff))
		return DistanceMatrix


	def DistanceMatrix(self,step):
		DistanceMatrix=np.zeros((self.Nchain,self.Nchain))
		for i in range(self.Nchain):
			for j in range(self.Nchain):
				diff=self.posHist[step][j]-self.posHist[step][i]
				DistanceMatrix[i][j]=np.sqrt(np.dot(diff.T,diff))
				DistanceMatrix[j][i]=np.sqrt(np.dot(diff.T,diff))
		return DistanceMatrix


	
	
	def PrintMatrices(self):
		for step in range(self.reader.N):
	#np.savetxt(self.OriginMatrixtransFile+str(step)+'.txt',self.DistanceMatrix(step))
			np.savetxt(self.OriginMatrixtransFile+str(step)+'.txt',self.DistanceMatrixtrans(step))
	#np.savetxt(self.OriginMatrixcisFile+str(step)+'.txt',self.DistanceMatrixcis(step))


	def PrintMoverate(self):
		repl=[]
		unrepl=[]
		fork=[]
		for step in range(1,self.reader.N):
			replcount=0
			replacc=0
			unreplcount=0
			unreplacc=0
			forkcont=0
			forkacc=0
			for tad in range(0,len(self.posHist[0])):
				if(self.ForkPos[step][tad]==0 and self.Status[step][tad]!=0):
					replcount+=1
					if(np.sum(self.posHist[step][tad]-self.posHist[step-1][tad])!=0):
						replacc+=1
				if(self.ForkPos[step][tad]!=0):
					forkcont+=1
					if(np.sum(self.posHist[step][tad]-self.posHist[step-1][tad])!=0):
						forkacc+=1
				if(self.ForkPos[step][tad]==0 and self.Status[step][tad]==0):
					unreplcount+=1
					if(np.sum(self.posHist[step][tad]-self.posHist[step-1][tad])!=0):
						unreplacc+=1
			repl.append(replacc/replcount)
			unrepl.append(unreplacc/unreplcount)
			fork.append(forkacc/forkcont)
		np.savetxt(os.path.join(self.reader.outputDir, "polyForkAcc.res"),fork)
		np.savetxt(os.path.join(self.reader.outputDir, "polyReplAcc.res"),unrepl)
		np.savetxt(os.path.join(self.reader.outputDir, "polyUnreplAcc.res"),repl)






	def CromatidGyration(self):
		for step in range(self.reader.N):
			self.Chr1pos=[]
			self.Chr2pos=[]
			for tad in range(len(self.posHist[step])):
				if(self.Status[step][tad]==1):
					self.Chr1pos.append(self.posHist[step][tad])
				if(self.Status[step][tad]==-1):
					self.Chr2pos.append(self.posHist[step][tad])
			
			if(len(self.Chr1pos)>0):
				#Chr1 and Chr2 Rcm distance
				diff=np.array(self.Chr1pos).mean(axis=0)-np.array(self.Chr2pos).mean(axis=0)
				self.RcmDistance.append(np.sqrt(np.dot(diff.T,diff)))
				#GyrationRadius
				self.Chr1pos=np.array(self.Chr1pos)
				self.Chr1pos-= self.Chr1pos.mean(axis=0, keepdims=True)
				diag1 = np.linalg.svd(self.Chr1pos, compute_uv=False) / len(self.Chr1pos)**0.5
				r2_gyr1 = np.square(diag1).sum(axis=-1)
				r_gyr1 = np.sqrt(r2_gyr1)
				self.Rg1.append(r_gyr1)
			if(len(self.Chr2pos)>0):
				self.Chr2pos=np.array(self.Chr2pos)
				self.Chr2pos-= self.Chr2pos.mean(axis=0, keepdims=True)
				diag2 = np.linalg.svd(self.Chr2pos, compute_uv=False) / len(self.Chr2pos)**0.5
				r2_gyr2 = np.square(diag2).sum(axis=-1)
				r_gyr2 = np.sqrt(r2_gyr2)
				self.Rg2.append(r_gyr2)


	def PrintSisterDist(self):
		#averagedistance=np.zeros(1)
		#error=np.zeros(1)

		#averagedistance[0]=sum(self.sisterdist)/len(self.sisterdist)
		#error[0]=np.array(self.sisterdist).std(axis=0)/len(self.sisterdist)**0.5

		np.savetxt(self.sisterdistFile,self.sisterdist)
		#np.savetxt(self.sisterdisterrFile,error)
	
	def computesisterdist(self):
		for step in range(self.reader.N):
			stepsister=[]
			for i in range(self.Nchain):
				sister=self.SisterID[step][i]
				if(sister!=-1):
					diff=self.posHist[step][i]-self.posHist[step][sister]
					stepsister.append(np.sqrt(np.dot(diff.T,diff)))
			if(len(stepsister)==0):
				self.sisterdist.append(0)
			else:
				self.sisterdist.append(sum(stepsister)/len(stepsister))

	def ReadHist(self):
		for i in range(self.reader.N):
			data = next(self.reader)
			self.posHist.append(data.polyPos)
			self.ForkPos.append(data.fork)
			self.SisterID.append(data.SisterID)
			self.Status.append(data.status)
			if (i==0):
				for t in range(self.reader.nTad):
					if(self.reader.status[t]==-1 or self.reader.status[t]==0):
						self.Nchain+=1
			
	def comulativedistance(self):
		self.comulativedistancearray=[]
		cumldiff=0;
		for step in range(self.reader.N-1):
			diff=self.posHist[step+1][idxTad]-self.posHist[step][idxTad]
			cumldiff+=np.sqrt(np.dot(diff.T,diff))
			self.comulativedistancearray.append(np.sqrt(cumldiff))
		np.savetxt(os.path.join(self.reader.outputDir, "comulativedistance"+str(idxTad)+".res"),self.comulativedistancearray)

				
			

	
						
	def PrintGyration(self):
		np.savetxt(self.Gyr1File,self.Rg1)
		np.savetxt(self.Gyr2File,self.Rg2)


	def PrintMSDto(self):
		to=[1,10,20,40,80,160,310]
		for i in range(len(to)):
			np.savetxt(os.path.join(self.reader.outputDir, "polyForkMSDt"+str(int(to[i]))+".res"),self.ComputeForkMSDto(int(to[i])))
			np.savetxt(os.path.join(self.reader.outputDir, "polyOriginMSDt"+str(int(to[i]))+".res"),self.ComputeOriginMSDto(int(to[i])))

	def PrintMSDtotrial(self):
		to=[1,10,20,40,80,160,310]
		for i in range(len(to)):
			np.savetxt(os.path.join(self.reader.outputDir, "polydistanceMSDtrial"+str(int(to[i]))+".res"),self.ComputedistanceMSDto(int(to[i])))
			np.savetxt(os.path.join(self.reader.outputDir, "polyPointMSDtrial"+str(int(to[i]))+".res"),self.ComputepointMSDto(int(to[i])))
	
	def Print(self):
		
		np.savetxt(self.ForkDistanceFile, self.ForkDistance)
		print("\033[1;32mPrinted ForkDistance to '%s'\033[0m" % self.ForkDistanceFile)

		np.savetxt(self.OriginDistanceFile, self.OriginDistance)
		print("\033[1;32mPrinted OriginDistance to '%s'\033[0m" % self.OriginDistanceFile)

		np.savetxt(self.ForkMSDFile, self.ForkMSD)
		print("\033[1;32mPrinted ForkMSD to '%s'\033[0m" % self.ForkMSDFile)
	
		np.savetxt(self.OriginMSDFile, self.OriginMSD)
		print("\033[1;32mPrinted OriginMSD to '%s'\033[0m" % self.OriginMSDFile)

		np.savetxt(self.RcmDistanceFile, self.RcmDistance)
		print("\033[1;32mPrinted OriginMSD to '%s'\033[0m" % self.RcmDistanceFile)

	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	ReplicationAnalysis = ReplicationAnalysis(outputDir, initFrame=initFrame)



	if len(sys.argv) == 3:
		ReplicationAnalysis.ReadHist()
		#ReplicationAnalysis.PrintMoverate()
		#ReplicationAnalysis.ComputeOriginDistance()
		ReplicationAnalysis.ComputeForkDistance()
		#ReplicationAnalysis.ComputeOriginMSD()
		#ReplicationAnalysis.ComputeForkMSD()
		#ReplicationAnalysis.CromatidGyration()
		#ReplicationAnalysis.PrintGyration()
		#ReplicationAnalysis.PrintMatrices()
		#ReplicationAnalysis.Print()
		#ReplicationAnalysis.PrintMSDto()
		#ReplicationAnalysis.PrintMSDtotrial()
		#ReplicationAnalysis.computesisterdist()
		#ReplicationAnalysis.PrintSisterDist()

		

		



		
	elif len(sys.argv) == 4:
		idxTad = int(sys.argv[3])
	
		ReplicationAnalysis.ReadHist()
		ReplicationAnalysis.comulativedistance()
