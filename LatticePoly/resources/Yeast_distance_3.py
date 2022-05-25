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
			if(self.SisterID[step][527]!=-1 and self.SisterID[step][527]!=527):
				if(step1==0):
					step1=step
			if(self.SisterID[step][675]!=-1 and self.SisterID[step][[675]]!=[675]):
				if(step2==0):
					step2=step
			if(self.SisterID[step][263-25]!=-1):
				if(neigh1==0):
					neigh1=step
			if(self.SisterID[step][263+25]!=-1):
					if(neigh2==0):
						neigh2=step



		self.repltime.append(((step1+step2)/2))
		self.repltime.append(abs(step1-step2))


		
		#np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"signal.res"),self.signal)
		#np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"repltime_3.res"),self.repltime)



	
	

	def arraysDistance(self):
		self.arrayDistance=[]
		self.sister1=[]
		self.sister2=[]
		self.arrayDistance.append(self.repltime[0])
		self.arrayDistance.append(self.repltime[1])
		for step in range(self.reader.N):
			diff=self.posHist[step][527]-self.posHist[step][675]
			self.arrayDistance.append(np.sqrt(np.dot(diff.T,diff)))

					   
		np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"arrayDistance_3.res"),self.arrayDistance)
		#np.savetxt(os.path.join(self.reader.outputDir,str(time.time())+ "sister1.res"),self.sister1)
		#np.savetxt(os.path.join(self.reader.outputDir,str(time.time())+ "sister2.res"),self.sister2)




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
			


	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	Yeast360 = Yeast360(outputDir, initFrame=initFrame)



	if len(sys.argv) == 3:
		Yeast360.ReadHist()

		Yeast360.signal()

		Yeast360.arraysDistance()




		

