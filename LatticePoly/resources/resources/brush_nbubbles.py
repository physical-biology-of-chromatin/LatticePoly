##  LatticePoly
##
##  Created by ddasaro on the model of mtortora script.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from scipy.spatial import cKDTree

from vtkReader import vtkReader

from scipy.spatial.distance import pdist, squareform ,cdist
import time

class brush_nbubbles():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
					
		self.bubbleFile = os.path.join(self.reader.outputDir, "brush_full_bubbles.res")

		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1


	def Compute(self):
		self.final_data=[]
		self.break_it=False
		self.replicons=0
		for i in range(self.reader.N):
			self.ProcessFrame(i)
			if(self.break_it==True):
				break
			
			
				
	def ProcessFrame(self, i):

		data = next(self.reader)
		current_n_replicons=0
		left_forks=[]
		right_forks=[]
		for i in range(0,len(data.fork)):
			if(data.fork[i]==-1):
				left_forks.append(i)
			if(data.fork[i]==1):
				right_forks.append(i)
		if(len(right_forks)==0):
			return
		if(self.replicons==0 and len(right_forks)!=0):
			self.replicons=len(right_forks)
		if(self.replicons!=len(right_forks)):
			self.break_it=True
			return
		
		
		array_replicons=[]
		array_links=[]
		n_replicon=[]
		n_links=[]
		for i in range(len(left_forks)):
			repl=np.arange(left_forks[i]+1,right_forks[i])
			n_replicon.append(len(repl))
			array_replicons.append([data.polyPos[k] for k in repl])
			repl_hom=[data.SisterID[k] for k in repl]
			array_replicons.append([data.polyPos[k] for k in repl_hom])
			

			if(i!=0):
				link=np.arange(right_forks[i-1],left_forks[i])
				array_links.append([data.polyPos[k] for k in link])
				n_links.append(len(link))

		step=np.full(5,0.0)
		step[0]=np.mean(n_replicon)
		step[1]=np.mean(n_links)
		
		Rg_n=[]
		Rg_m=[]
		d=[]
		for id,link in enumerate(array_links):
			pos = np.array(link)
			pos -= pos.mean(axis=0, keepdims=True)
			pos=pos/len(link)**0.5
			diag = np.linalg.svd(pos, compute_uv=False)
			r2_gyr = np.square(diag).sum(axis=-1)
			Rg_m.append(r2_gyr)
			diff=link[0]-link[-1]
			d.append(np.sqrt(np.dot(diff.T,diff)))
			
		for id,replicon in enumerate(array_replicons):
			pos = np.array(replicon)
			pos -= pos.mean(axis=0, keepdims=True)
			pos=pos/len(replicon)**0.5
			diag = np.linalg.svd(pos, compute_uv=False)
			r2_gyr = np.square(diag).sum(axis=-1)
			Rg_n.append(r2_gyr)
			
		step[2]=np.nanmean(d)
		step[3]=np.nanmean(Rg_n)
		step[4]=np.nanmean(Rg_m)
		self.final_data.append(step)
				
		
		

	def Print(self):
		while(len(self.final_data)!=self.reader.N):
			self.final_data.append(np.full(5,np.nan))

			
		np.savetxt(self.bubbleFile,self.final_data)
		

		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	brush = brush_nbubbles(outputDir, initFrame=initFrame)

	brush.Compute()
	brush.Print()
