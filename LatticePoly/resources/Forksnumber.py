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


class Forksnumber():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		self.posHist = []
		self.ForkPos =[]
		self.SisterID=[]
		self.Status=[]

		
		self.ForksnumbFile = os.path.join(self.reader.outputDir, "Forksnumb.res")

	def ReadHist(self):
		self.Nchain=0

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

	def compute(self):
		self.Forksnumber=np.zeros(self.reader.N)
		for step in range(self.reader.N):
			for i in range(len(self.ForkPos[0])):
				if(self.ForkPos[step][i]==-1 or self.ForkPos[step][i]==1 ):
					self.Forksnumber[step]+=1



	
	
	def Print(self):
		np.savetxt(self.ForksnumbFile, self.Forksnumber)
		
		print("\033[1;32mPrinted Forksmove  to '%s'\033[0m" % self.ForksnumbFile)


	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	forksnumb = Forksnumber(outputDir, initFrame=initFrame)



	if len(sys.argv) == 3:
		forksnumb.ReadHist()
		forksnumb.compute()
		forksnumb.Print()

