##
# MonomerDist.py
# LatticePoly
##
##

import os
import sys
import pandas as pd
import numpy as np
from vtkReader_multi import vtkReader
from utils import msdFFT




class Poly_Forks_number():
	def __init__(self, outputDir, initFrame):
		self.readers=[]
		self.N_chain=[]
		for i in range(0,17):
			self.readers.append(vtkReader(outputDir, i,initFrame,readLiq=False, readPoly=True))
			self.N_chain.append(np.sum([(self.readers[i].status==-1)+(self.readers[i].status==0)]))
		self.NForks_n_File = os.path.join(outputDir,"N_forks_full_genome.res")
		self.percentage_File = os.path.join(outputDir,"percentage_full_genome.res")


		self.Compute(self.readers[0].N)
		self.Print()
			






		
	#NB here is not the finalFrame but the number of iterations
	def Compute(self,finalFrame):
		#position at to
		self.percentage = np.full(finalFrame,0)
		self.forks = np.full(finalFrame,0)


		for i in range(0, finalFrame):
			self.ProcessFrame(i)
			if (i+1) % 100 == 0:
				print("Processed %d out of %d configurations" %
					  (i+1, finalFrame))



	def ProcessFrame(self, i):
		for chrom,reader in enumerate(self.readers):
			data = next(reader)
			self.forks[i]+=np.sum(data.fork!=0)
			self.percentage[i]+= len(data.fork)
		#full_status=np.concatenate(full_status)
		#self.Repl_time[i]=full_status



	

	def Print(self):
		np.savetxt(self.NForks_n_File,self.forks)
		np.savetxt(self.percentage_File,self.percentage)

		


if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame \033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])


Poly_Forks_number = Poly_Forks_number(outputDir, initFrame=initFrame)

