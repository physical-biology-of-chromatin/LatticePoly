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

from vtkReader_multi import vtkReader
import time

class PolyGyration():
	
	def __init__(self, outputDir, initFrame):
		self.reader1 = vtkReader(outputDir,1, initFrame, readLiq=False, readPoly=True)
		self.reader2 = vtkReader(outputDir,0, initFrame, readLiq=False, readPoly=True)

					
		self.diffRcmFile = os.path.join(self.reader1.outputDir, "diff_rcm.res")

		self.Nchain=0



	def Compute(self):
		self.diff=[]
		for i in range(self.reader1.N):
			self.ProcessFrame(i)
			
			#if (i+1) % 10 == 0:
			#	print("Processed %d out of %d configurations" % (i+1, self.reader.N))
			
				
	def ProcessFrame(self, i):
		data = next(self.reader1)
		data2 = next(self.reader2)

		
		


		rcm1=np.array(data.polyPos[:])
		rcm2=np.array(data2.polyPos[:])
		diff_rcm=np.mean(rcm1,axis=0) - np.mean(rcm2,axis=0)
		self.diff.append(np.dot(diff_rcm,diff_rcm.T))


	def Print(self):
		np.savetxt(self.diffRcmFile, np.array(self.diff)**0.5)
		

		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	gyr = PolyGyration(outputDir, initFrame=initFrame)

	gyr.Compute()
	gyr.Print()
