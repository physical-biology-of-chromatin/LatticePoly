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




class PolyMSD():
	def __init__(self, outputDir, initFrame,interval):
		self.readers=[]
		self.N_chain=[]
		for i in range(0,17):
			self.readers.append(vtkReader(outputDir, i,initFrame,readLiq=False, readPoly=True))
			self.N_chain.append(np.sum([(self.readers[i].status==-1)+(self.readers[i].status==0)]))
		self.PolyMSD_File = os.path.join(outputDir,"MSD_"+str(initFrame)+ "_"+str(interval)+"full_genome_MSD.res")

		self.Compute(interval)
		self.Print()
			






		
	#NB here is not the finalFrame but the number of iterations
	def Compute(self,finalFrame):
		#position at to
		self.Pos = np.full((finalFrame,np.sum(self.N_chain),3),np.nan)
		#replication_time
		#self.Repl_time = np.full((finalFrame,np.sum(self.N_chain)),np.nan)


		for i in range(0, finalFrame):
			self.ProcessFrame(i)
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" %
					  (i+1, finalFrame))



	def ProcessFrame(self, i):
		full_positions=[]
		full_status=[]

		for chrom,reader in enumerate(self.readers):
			data = next(reader)
			full_positions.append(data.polyPos[:self.N_chain[chrom]])
			#full_status.append(data.status[:self.N_chain[chrom]])

				 
		full_positions=np.concatenate(full_positions)
		self.Pos[i]=full_positions

		#full_status=np.concatenate(full_status)
		#self.Repl_time[i]=full_status



	

	def Print(self):
		MSD=np.zeros((np.sum(self.N_chain),interval))
		MSD[0]=msdFFT(self.Pos[:, 0])
		for i in range(np.sum(self.N_chain)):
			MSD[i]=msdFFT(self.Pos[:, i])
		np.savetxt(self.PolyMSD_File,MSD)
		


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("\033[1;31mUsage is %s outputDir initFrame interval \033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	interval=int(sys.argv[3])


PolyMSD = PolyMSD(outputDir, initFrame=initFrame, interval=interval)

