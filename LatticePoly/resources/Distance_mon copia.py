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


#if os.path.exists(self.ForkDistanceFile) & os.path.exists(self.OriginDistanceFile):
#			print("Files '%s' and '%s' already exist - aborting" % (self.ForkDistanceFile, self.OriginDistanceFile))
#			sys.exit()




		
		#np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"signal.res"),self.signal)



	
	

	def arraysDistance(self):
		self.arrayDistance=[]
		mon2=int(center+distance/2)
		mon1=int(center-distance/2)
		repltime1=-1
		repltime2=-1
		for step in range(self.reader.N):
			data=next(self.reader)
			diff=data.polyPos[mon1]-data.polyPos[mon2]
			if(data.status[mon1]!=0 and repltime1==-1):
				repltime1=step
			if(data.status[mon2]!=0 and repltime2==-1):
				repltime2=step
			self.arrayDistance.append(np.sqrt(np.dot(diff.T,diff)))

		self.arrayDistance.insert(0,(repltime2))
		self.arrayDistance.insert(0,(repltime1))
		np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"_"+str(distance)+"_mon_dist.res"),self.arrayDistance)
		#np.savetxt(os.path.join(self.reader.outputDir,str(time.time())+ "sister1.res"),self.sister1)
		#np.savetxt(os.path.join(self.reader.outputDir,str(time.time())+ "sister2.res"),self.sister2)


	
if __name__ == "__main__":
	if len(sys.argv) not in [5]:
		print("\033[1;31mUsage is %s outputDir initFrame center distance \033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	center = int(sys.argv[3])
	distance = int(sys.argv[4])



	
	Distance_mon = Distance_mon(outputDir, initFrame=initFrame)
	Distance_mon.arraysDistance()




	