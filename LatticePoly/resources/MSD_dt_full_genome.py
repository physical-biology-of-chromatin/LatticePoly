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



class MSD_dt():
	def __init__(self, outputDir, initFrame,dt,interval):
		self.readers=[]
		self.N_chain=[]
		for i in range(0,17):
			self.readers.append(vtkReader(outputDir, i,initFrame,readLiq=False, readPoly=True))
			self.N_chain.append(np.sum([(self.readers[i].status==-1)+(self.readers[i].status==0)]))
		self.MSD_dt_File = os.path.join(outputDir,"dt_"+str(dt)+"_"+str(initFrame)+ "_"+str(interval)+"full_genome_MSD.res")
		self.Repl_time_File = os.path.join(outputDir,"Monomer_Repl_time_full_genome.res")


		self.Compute(interval)
		self.Print()
			






		
	#NB here is not the finalFrame but the number of iterations
	def Compute(self,finalFrame):
		#position at to
		self.Pos = np.full((finalFrame,np.sum(self.N_chain),3),np.nan)
		#replication_time
		self.Repl_time = np.full((finalFrame,np.sum(self.N_chain)),np.nan)


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
			full_status.append(data.status[:self.N_chain[chrom]])

				 
		full_positions=np.concatenate(full_positions)
		self.Pos[i]=full_positions
		full_status=np.concatenate(full_status)
		self.Repl_time[i]=full_status



	

	def Print(self):
		diff=self.Pos[:-dt]-self.Pos[dt:]
		squared_displacement=np.nansum(diff**2,axis=2)
		#average over all to and obtain array of the mean squared displacement in dt of each monomer
		#print(np.shape(np.nanmean(squared_displacement,axis=0)))
		np.savetxt(self.MSD_dt_File,squared_displacement)
		np.savetxt(self.Repl_time_File,self.Repl_time)

		


if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("\033[1;31mUsage is %s outputDir initFrame dt interval \033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	dt=int(sys.argv[3])
	interval=int(sys.argv[4])

	if(interval<=dt):
		print("Error! interval<=dt")
		sys.exit()
	
	#init_time=int(sys.argv[4])


MSD_dt = MSD_dt(outputDir, initFrame=initFrame, dt=dt, interval=interval)

