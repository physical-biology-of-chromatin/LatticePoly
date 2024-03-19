# -*- coding: utf-8 -*-
##
# MonomerDist.py
# LatticePoly
##
# Based on mtortora's  code.
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import sys
import pandas as pd
import numpy as np
import cooler
from scipy.spatial import cKDTree
from cooler.create import ArrayLoader
from vtkReader_multi import vtkReader

from scipy.spatial.distance import pdist, squareform


class MonomerDmap():
	def __init__(self, outputDir, initFrame):
		self.readers=[]
		self.N_chain=[]
		for i in range(0,16):
			self.readers.append(vtkReader(outputDir, i,initFrame,readLiq=False, readPoly=True))
			self.N_chain.append(np.sum([(self.readers[i].status==-1)+(self.readers[i].status==0)]))
		self.contactFile = os.path.join(outputDir,"r_"+str(r)+"_"+str(initFrame)+ "_"+str(interval)+"full_genome.cool")


		self.Compute(interval)

		self.Print()
			






		
	#NB here is not the finalFrame but the number of iterations
	def Compute(self,finalFrame):
		#self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
		n_bins=2*np.sum(self.N_chain)
		self.contactProb = np.zeros((n_bins, n_bins), dtype=np.float32)
		
		

		for i in range(0, finalFrame):
			self.ProcessFrame(i)
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" %
					  (i+1, finalFrame))


#self.contactProb=np.rint(self.contactProb/finalFrame)

	def ProcessFrame(self, i):
		full_positions=[]
		for chrom,reader in enumerate(self.readers):
			data = next(reader)
			positions_sister=np.full((2*self.N_chain[chrom],3),np.nan)
			for pos_id,pos in enumerate(data.polyPos):
				if(pos_id<self.N_chain[chrom]):
					positions_sister[pos_id]=pos
				else:
					positions_sister[self.N_chain[chrom]+data.SisterID[pos_id]]=pos
					
			full_positions.append(positions_sister)

				 
		tree1	= cKDTree(np.concatenate(full_positions), boxsize = None)
		pairs = tree1.query_pairs(r = r*0.71) # NN distance FCC lattice 1/np.sqrt(2) = 0.71
		for (i,j) in pairs:
			self.contactProb[i,j] = self.contactProb[i,j] + 1
			self.contactProb[j,i] = self.contactProb[j,i] + 1
		
		
			

					



	def Print(self):
		final_bin_size=int(len(self.contactProb)/2)
		final_map = np.zeros((final_bin_size, final_bin_size), dtype=np.float32)
		for chrom_id1 in range(len(self.N_chain)):
			for chrom_id2 in range(chrom_id1,len(self.N_chain)):
				#define a sub-matrix
				chrom_start_bin1=int(2*np.sum(self.N_chain[:chrom_id1]))
				chrom_end_bin1=int(chrom_start_bin1+2*self.N_chain[chrom_id1])
				chrom_start_bin2=int(2*np.sum(self.N_chain[:chrom_id2]))
				chrom_end_bin2=int(chrom_start_bin2+2*self.N_chain[chrom_id2])
				
				current_matrix_block=self.contactProb[chrom_start_bin1:chrom_end_bin1,chrom_start_bin2:chrom_end_bin2]
				#divide matrix in 4 block and sum them
				print(np.shape(current_matrix_block))

				side1=int(np.shape(current_matrix_block)[0]/2)
				side2=int(np.shape(current_matrix_block)[1]/2)
				block1=current_matrix_block[0:side1,0:side2]
				block2=current_matrix_block[side1:,0:side2]
				block3=current_matrix_block[0:side1,side2:]
				block4=current_matrix_block[side1:,side2:]
				print(np.shape(current_matrix_block))
				print(np.shape(block1))
				print(np.shape(block2))
				print(np.shape(block3))
				print(np.shape(block4))



				
				fchrom_start_bin1=int(chrom_start_bin1/2)
				fchrom_start_bin2=int(chrom_start_bin2/2)
				fchrom_end_bin1=int(chrom_end_bin1/2)
				fchrom_end_bin2=int(chrom_end_bin2/2)

				final_map[fchrom_start_bin1:fchrom_end_bin1,fchrom_start_bin2:fchrom_end_bin2]+=block1
				final_map[fchrom_start_bin1:fchrom_end_bin1,fchrom_start_bin2:fchrom_end_bin2]+=block2
				final_map[fchrom_start_bin1:fchrom_end_bin1,fchrom_start_bin2:fchrom_end_bin2]+=block3
				final_map[fchrom_start_bin1:fchrom_end_bin1,fchrom_start_bin2:fchrom_end_bin2]+=block4



		
		#np.savetxt(outputDir+"full_contactmap.res",self.contactProb)
		
		chrom_names=["chr" + str(i+1) for i in range(len(self.N_chain))]
		chromsizes=pd.Series(1000*np.array(self.N_chain))
		print(chromsizes)
		chromsizes=chromsizes.rename(lambda x: chrom_names[x])
		chromsizes=chromsizes.astype('int64')	
		bins = cooler.binnify(chromsizes, 1000)
		pixels = ArrayLoader(bins, final_map, chunksize=10000000)
		cooler.create_cooler(self.contactFile,bins,pixels)
		
		
		print("\033[1;32mPrinted avg.contact probability to '%s'\033[0m" %self.contactFile)


if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("\033[1;31mUsage is %s outputDir initFrame r interval \033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	r=float(sys.argv[3])
	interval=int(sys.argv[4])
	
	#init_time=int(sys.argv[4])


monom = MonomerDmap(outputDir, initFrame=initFrame)

