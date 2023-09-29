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

class stored_lenght():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)


#if os.path.exists(self.ForkDistanceFile) & os.path.exists(self.OriginDistanceFile):
#			print("Files '%s' and '%s' already exist - aborting" % (self.ForkDistanceFile, self.OriginDistanceFile))
#			sys.exit()




		
		#np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"signal.res"),self.signal)



	
	

	def stored_lenght(self):
		mon2=int(center+distance/2)
		mon1=int(center-distance/2)
		self.stored_lenght=[]
		for step in range(self.reader.N):
			data=next(self.reader)
			curr_stored_lenght=0
			for i in range(mon1,mon2):
				if(data.polyPos[i][0]==data.polyPos[i+1][0] and data.polyPos[i][1]==data.polyPos[i+1][1] and data.polyPos[i][2]==data.polyPos[i+1][2]):
					curr_stored_lenght+=1
			self.stored_lenght.append(curr_stored_lenght)

		np.savetxt(os.path.join(self.reader.outputDir, str(time.time())+"_"+str(distance)+"_stored_lenght.res"),self.stored_lenght)
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



	
	stored_mon = stored_lenght(outputDir, initFrame=initFrame)
	stored_mon.stored_lenght()




	

