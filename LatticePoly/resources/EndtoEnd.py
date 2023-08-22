##  LatticePoly
##
##  Created by ddasaro on the model of mtortora script.
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

class Distance_mon():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1

#if os.path.exists(self.ForkDistanceFile) & os.path.exists(self.OriginDistanceFile):
#			print("Files '%s' and '%s' already exist - aborting" % (self.ForkDistanceFile, self.OriginDistanceFile))
#			sys.exit()




		
		#np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"signal.res"),self.signal)



	
	

	def arraysDistance(self):
		self.arrayDistance=[]
		self.mon=[]
		mon2=self.Nchain-1
		mon1=0
		for step in range(self.reader.N):
			data=next(self.reader)
			diff=data.polyPos[mon1]-data.polyPos[mon2]
			self.arrayDistance.append(np.sqrt(np.dot(diff.T,diff)))
			self.mon.append(data.nTad-self.Nchain)

		np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"_endtoend.res"),np.array([self.mon,self.arrayDistance]).T)
		#np.savetxt(os.path.join(self.reader.outputDir,str(time.time())+ "sister1.res"),self.sister1)
		#np.savetxt(os.path.join(self.reader.outputDir,str(time.time())+ "sister2.res"),self.sister2)


	
if __name__ == "__main__":
	if len(sys.argv) not in [3]:
		print("\033[1;31mUsage is %s outputDir initFrame center distance \033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])




	
	Distance_mon = Distance_mon(outputDir, initFrame=initFrame)
	Distance_mon.arraysDistance()




	
