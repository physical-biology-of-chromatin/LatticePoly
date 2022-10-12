##
##  PolyGyration.py
##  LatticePoly
##
##  Created by mtortora on 02/08/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from vtkReader import vtkReader
import time

class PolyGyration():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
					
		self.diffRcmFile = os.path.join(self.reader.outputDir, str(time.time())+"diff_rcm.res")

		self.Nchain=0



	def Compute(self):
		self.diff=np.zeros(self.reader.N)
		for i in range(self.reader.N):
			self.ProcessFrame(i)
			
			#if (i+1) % 10 == 0:
			#	print("Processed %d out of %d configurations" % (i+1, self.reader.N))
			
				
	def ProcessFrame(self, i):
		data = next(self.reader)

		if(i==0):
			for t in range(self.reader.nTad):
				if(data.status[t]==-1 or data.status[t]==0):
					self.Nchain+=1
		rcm1=np.zeros((self.Nchain,3))
		rcm2=np.zeros((self.Nchain,3))

		for n in range(self.Nchain):
			if(data.status[n]==0):
				rcm1[n]=data.polyPos[n]
				rcm2[n]=data.polyPos[n]
			else:
				rcm1[n]=data.polyPos[n]
				rcm2[n]=data.polyPos[data.SisterID[n]]


									
		diff_rcm=np.mean(rcm1) - np.mean(rcm2)
		self.diff[i] = np.dot(diff_rcm,diff_rcm.T)


	def Print(self):
		np.savetxt(self.diffRcmFile, self.diff**0.5)
		

		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	gyr = PolyGyration(outputDir, initFrame=initFrame)

	gyr.Compute()
	gyr.Print()
