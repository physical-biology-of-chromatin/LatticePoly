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
import networkx as nx
from scipy.spatial.distance import pdist, squareform


class cluster_full_genome():
	def __init__(self, outputDir, initFrame):
		self.readers=[]
		self.N_chain=[]
		for i in range(0,17):
			self.readers.append(vtkReader(outputDir, i,initFrame,readLiq=False, readPoly=True))
			self.N_chain.append(np.sum([(self.readers[i].status==-1)+(self.readers[i].status==0)]))
		self.clusterFile = os.path.join(outputDir,str(r)+"_"+str(initFrame)+"_"+str(interval)+"largest_component_full_genome.res")
		#self.clusterFile_cis = os.path.join(outputDir,str(r)+"_"+str(initFrame)+"_"+str(interval)+"total_cluster_full_genome_cis.res")
		#self.ForksFile = os.path.join(outputDir,str(initFrame)+"_"+str(interval)+"forks_full_genome.res")
		#self.Forks_per_clusterFile = os.path.join(outputDir,str(r)+"_"+str(initFrame)+"_"+str(interval)+"_forks_per_cluster_full_genome.res")
		#self.Forks_per_cluster_cisFile = os.path.join(outputDir,str(r)+"_"+str(initFrame)+"_"+str(interval)+"_forks_per_cluster_full_genome_cis.res")




		self.Compute(interval)

		self.Print()
			






		
	def Compute(self,finalFrame):
		self.cluster_number=[]
		#self.forks_per_cluster=np.zeros(500)
		#self.forks_per_cluster_cis=np.zeros(500)
		for i in range(0, finalFrame):
			self.ProcessFrame(i)
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" %
					  (i+1, finalFrame))



	def ProcessFrame(self, i):
		full_positions=[]
		forks_for_chrom=[]
		for chrom,reader in enumerate(self.readers):
			data = next(reader)
			fork_mask=np.array((data.fork==1)+(data.fork==-1))
			if(np.sum(fork_mask)>0):
				full_positions.append(data.polyPos[fork_mask])
			forks_for_chrom.append(np.sum(fork_mask))	 
		tree1	= cKDTree(np.concatenate(full_positions), boxsize = None)
		pairs = tree1.query_pairs(r = r*0.71) # NN distance FCC lattice 1/np.sqrt(2) = 0.71
		A=np.zeros((len(np.concatenate(full_positions)),len(np.concatenate(full_positions))))
		for (i,j) in pairs:
			A[i,j] = A[i,j] + 1
			A[j,i] = A[j,i] + 1

		G = nx.from_numpy_matrix(A)
		self.cluster_number.append(np.max([len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]))

				



		
			

					



	def Print(self):
		#print(self.cluster_number)
		np.savetxt(self.clusterFile, self.cluster_number)


if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("\033[1;31mUsage is %s outputDir initFrame r interval \033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	r=float(sys.argv[3])
	interval=int(sys.argv[4])
	
	#init_time=int(sys.argv[4])


monom = cluster_full_genome(outputDir, initFrame=initFrame)

