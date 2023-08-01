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

class Mixing_nbubbles():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
					
		self.bubbleFile = os.path.join(self.reader.outputDir, "mixing_full_bubbles.res")

		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1


	def Compute(self):
		self.final_data=[]

		self.replicons=0
		for i in range(self.reader.N):
			self.ProcessFrame(i)
			
			#if (i+1) % 10 == 0:
			#	print("Processed %d out of %d configurations" % (i+1, self.reader.N))
			
				
	def ProcessFrame(self, i):
		print("step "+str(i) )

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
		#print(self.replicons)
		if(self.replicons!=len(right_forks)):
			return
		
		
		array_replicons=[]
		array_replicons2=[]
		for i in range(len(left_forks)):
			array_replicons.append(np.arange(left_forks[i]+1,right_forks[i]))
		for i in range(len(array_replicons)):
			current_replicon2=[]
			for k in range(len(array_replicons[i])):
				current_replicon2.append(data.SisterID[array_replicons[i][k]])
			array_replicons2.append(current_replicon2)

		for id,replicon in enumerate(array_replicons):
			step=np.full(5,0.0)
			step[0]=len(replicon)
			intra_bubble=[data.polyPos[k] for k in replicon]
			inter_bubble=[data.polyPos[k] for k in array_replicons2[id]] 
			#print(intra_bubble)
			#print(inter_bubble)

			cis_bubbles_id=[]
			for el in array_replicons:
				if((el[0]!=replicon[0])):
					cis_bubbles_id.extend(el)
			trans_bubbles_id=[]
			for el in array_replicons2:
				if(((el[0]!=array_replicons2[id][0]))):
					trans_bubbles_id.extend(el)

			#print(cis_bubbles_id)
			#print(trans_bubbles_id)
			
			trans_bubbles=[data.polyPos[k] for k in trans_bubbles_id] 
			cis_bubbles=[data.polyPos[k] for k in cis_bubbles_id] 
			#print(cis_bubbles)
			#print(trans_bubbles)


			#IntraLoop
			matr=pdist(intra_bubble)
			step[1]=np.sum(np.sum((matr < r*0.71) * (matr > 0 )))
			#InterLoop
			matr=cdist(intra_bubble,inter_bubble)
			step[2]=np.sum(np.sum((matr < r*0.71) * (matr > 0 )))
			#cisLoop
			matr=cdist(intra_bubble,cis_bubbles)
			step[3]=np.sum(np.sum((matr < r*0.71) * (matr > 0 )))
			#cisLoop
			matr=cdist(intra_bubble,trans_bubbles)
			step[4]=np.sum(np.sum((matr < r*0.71) * (matr > 0 )))

			
			

			#step[1]=step[1]/(step[0]*(step[0]-1))
			#step[2]=step[2]/(2*step[0])**2
			#step[3]=step[3]/(step[0]*(len(cis_bubbles)-step[0]))
			#step[4]=step[4]/(step[0]*(len(cis_bubbles)-step[0]))
			self.final_data.append(step)
				
		
		

	def Print(self):

		np.savetxt(self.bubbleFile,self.final_data)
		

		
if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("\033[1;31mUsage is %s outputDir initFrame r\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	r=int(sys.argv[3])

	mix = Mixing_nbubbles(outputDir, initFrame=initFrame)

	mix.Compute()
	mix.Print()
