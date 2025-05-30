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


def round_to_even(num):
    rounded = round(num)
    # If the rounded number is odd, adjust to the nearest even
    if rounded % 2 != 0:
        if rounded > num:
            return rounded - 1
        else:
            return rounded + 1
    return rounded


class cluster_full_genome():
	def __init__(self, outputDir, initFrame):
		self.readers=[]
		self.N_chain=[]
		for i in range(0,17):
			self.readers.append(vtkReader(outputDir, i,initFrame,readLiq=False, readPoly=True))
			self.N_chain.append(np.sum([(self.readers[i].status==-1)+(self.readers[i].status==0)]))
		self.Forks_per_clusterFile = os.path.join(outputDir,str(r)+"_"+str(initFrame)+"_"+str(interval)+"_forks_"+str(forks_n)+"percentage_missed_links_"+str(percentage)+"_cluster_even_full_genome.res")
		self.Compute(interval)
		self.Print()






	def Compute(self,finalFrame):
		self.forks_per_cluster=np.zeros(500)
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
			if(100*np.random.rand(1)[0]>percentage):
				A[i,j] = A[i,j] + 1
				A[j,i] = A[j,i] + 1

		G = nx.from_numpy_matrix(A)
		tot_forks=np.nansum(np.array([len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]))
		for e in [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]:
			self.forks_per_cluster[round_to_even(forks_n*e/tot_forks)]+=1
		#print(sorted(nx.connected_components(G), key=len, reverse=True))




		
			

					



	def Print(self):
		np.savetxt(self.Forks_per_clusterFile, self.forks_per_cluster)


if __name__ == "__main__":
	if len(sys.argv) != 7:
		print("\033[1;31mUsage is %s outputDir initFrame r interval fork_n percentage(of missed links) \033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	r=float(sys.argv[3])
	interval=int(sys.argv[4])
	forks_n=int(sys.argv[5])
	percentage=int(sys.argv[6])


monom = cluster_full_genome(outputDir, initFrame=initFrame)

