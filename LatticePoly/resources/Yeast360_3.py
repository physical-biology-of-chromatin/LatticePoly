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

class Yeast360():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		self.posHist = []
		self.ForkPos =[]
		self.SisterID=[]
		self.Status=[]
		self.Rg1=[]
		self.Rg2=[]
		self.RcmDistance=[]
		
		self.Nchain=0
		self.Origin=260
		
		self.ForkDistanceFile = os.path.join(self.reader.outputDir, str(time.time())+"polyForkDistance.res")
		self.OriginDistanceFile = os.path.join(self.reader.outputDir, str(time.time())+"polyOriginDistance.res")
		self.ForkMSDFile = os.path.join(self.reader.outputDir,str(time.time())+ "polyForkMSD.res")
		self.ForkMSDtoFile = os.path.join(self.reader.outputDir,str(time.time())+ "polyForkMSDto.res")
		self.OriginMSDFile = os.path.join(self.reader.outputDir, str(time.time())+"polyOriginMSD.res")
		self.OriginMSDtoFile = os.path.join(self.reader.outputDir,str(time.time())+ "polyOriginMSDto.res")
		self.OriginMatrixFile = os.path.join(self.reader.outputDir, "Matrix")
		self.OriginMatrixcisFile = os.path.join(self.reader.outputDir, "Matrixcis")
		self.OriginMatrixtransFile = os.path.join(self.reader.outputDir, "Matrixtrans")
		self.OriginMatrixFile = os.path.join(self.reader.outputDir, "Matrix")
		self.Gyr1File = os.path.join(self.reader.outputDir,str(time.time())+ "Gyr1.res")
		self.Gyr2File = os.path.join(self.reader.outputDir, str(time.time())+"Gyr2.res")
		self.RcmDistanceFile = os.path.join(self.reader.outputDir, str(time.time())+"RcmDistance.res")





#if os.path.exists(self.ForkDistanceFile) & os.path.exists(self.OriginDistanceFile):
#			print("Files '%s' and '%s' already exist - aborting" % (self.ForkDistanceFile, self.OriginDistanceFile))
#			sys.exit()


	def signal(self):
		self.signal=[]
		step1=0
		step2=0
		neigh1=0
		neigh2=0
		
		self.repltime=[]
		for step in range(self.reader.N):
			if(self.SisterID[step][round(660/1.25)]!=-1):
				if(step1==0):
					step1=step
			if(self.SisterID[step][round(844/1.25)]!=-1):
				if(step2==0):
					step2=step
			if(self.SisterID[step][680]!=-1):
				if(neigh1==0):
					neigh1=step
			if(self.SisterID[step][263+25]!=-1):
					if(neigh2==0):
						neigh2=step

		#if( neigh1<step1 or neigh2<step2):
		#	sys.exit()

		self.repltime.append(((step1+step2)/2))
		self.repltime.append(abs(step1-step2))


		
		#np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"signal.res"),self.signal)
		np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"repltime_3.res"),self.repltime)



	
	

	def arraysDistance(self):
		self.arrayDistance=[]
		self.sister1=[]
		self.sister2=[]

		for step in range(self.reader.N):
			diff=self.posHist[step][round(660/1.25)]-self.posHist[step][round(844/1.25)]
			self.arrayDistance.append(np.sqrt(np.dot(diff.T,diff)))

					   
		np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"arrayDistance_3.res"),self.arrayDistance)
		#np.savetxt(os.path.join(self.reader.outputDir,str(time.time())+ "sister1.res"),self.sister1)
		#np.savetxt(os.path.join(self.reader.outputDir,str(time.time())+ "sister2.res"),self.sister2)


	def arraysDistance1(self):
		self.arrayDistance=[]
		for step in range(self.reader.N):
			diff=self.posHist[step][284]-self.posHist[step][284+48]
			self.arrayDistance.append(np.sqrt(np.dot(diff.T,diff)))

		np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"arrayDistance1.res"),self.arrayDistance)
	
	def arraysDistancemedium(self):
		self.arrayDistance=[]


		for step in range(self.reader.N):
			pos1=np.zeros(3)
			pos2=np.zeros(3)
			for tad in range(236-4,236+5):
				pos1+=self.posHist[step][tad]
			pos1=pos1/9
			for tad in range(284-4,284+5):
				pos2+=self.posHist[step][tad]
			pos2=pos2/9
			
			diff=pos2-pos1
			self.arrayDistance.append(np.sqrt(np.dot(diff.T,diff)))



		np.savetxt(os.path.join(self.reader.outputDir,str(time.time())+ "arrayDistancemedium.res"),self.arrayDistance)
	

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
			np.savetxt(self.OriginMatrixcisFile+str(step)+'.txt',self.DistanceMatrix(step))
#np.savetxt(self.OriginMatrixtransFile+str(step)+'.txt',self.DistanceMatrixtrans(step))
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

		'''
		np.savetxt(self.ForkMSDFile, self.ForkMSD)
		print("\033[1;32mPrinted ForkMSD to '%s'\033[0m" % self.ForkMSDFile)
	
		np.savetxt(self.OriginMSDFile, self.OriginMSD)
		print("\033[1;32mPrinted OriginMSD to '%s'\033[0m" % self.OriginMSDFile)

		np.savetxt(self.RcmDistanceFile, self.RcmDistance)
		print("\033[1;32mPrinted OriginMSD to '%s'\033[0m" % self.RcmDistanceFile)
		'''
	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	Yeast360 = Yeast360(outputDir, initFrame=initFrame)



	if len(sys.argv) == 3:
		Yeast360.ReadHist()
		#ReplicationAnalysis.PrintMoverate()
		#Yeast360.ComputeOriginDistance()
		#Yeast360.ComputeForkDistance()
		#Yeast360.Print()

		#ReplicationAnalysis.ComputeOriginMSD()
		#ReplicationAnalysis.ComputeForkMSD()
		#ReplicationAnalysis.CromatidGyration()
		#ReplicationAnalysis.PrintGyration()
		#ReplicationAnalysis.PrintMatrices()
		#Yeast360.PrintMatrices()
		#ReplicationAnalysis.PrintMSDto()
		#ReplicationAnalysis.PrintMSDtotrial()
		Yeast360.signal()

		Yeast360.arraysDistance()
		#Yeast360.arraysDistance1()
		#Yeast360.arraysDistancemedium()




		
	elif len(sys.argv) == 4:
		idxTad = int(sys.argv[3])
	
		ReplicationAnalysis.ComputeTad(idxTad)
		ReplicationAnalysis.PrintTad(idxTad)
