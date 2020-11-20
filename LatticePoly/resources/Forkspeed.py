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


class ForkSpeed():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		
		self.ForkSpeedFile = os.path.join(self.reader.outputDir, "ForkSpeed.res")

		if os.path.exists(self.ForkSpeedFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.ForkSpeedFile))
			sys.exit()


	def Compute(self):


		
		MonomerPerStep=np.zeros(self.reader.N)
		for step in range(self.reader.N):
			data=next(self.reader)
			MonomerPerStep[step]=data.nTad
		self.Forksmoves=np.zeros(self.reader.N)
		for i in range(len(self.Forksmoves)):
			if (i==0):
				self.Forksmoves[i]=MonomerPerStep[i]
			else:
				self.Forksmoves[i]=MonomerPerStep[i]

	
	
	def Print(self):
		Forksmoves = self.Forksmoves
		np.savetxt(self.ForkSpeedFile, Forksmoves)
		
		print("\033[1;32mPrinted Forksmove  to '%s'\033[0m" % self.ForkSpeedFile)


	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	forkspeed = ForkSpeed(outputDir, initFrame=initFrame)



	if len(sys.argv) == 3:
		forkspeed.Compute()
		forkspeed.Print()

