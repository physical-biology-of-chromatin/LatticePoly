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

class Writhe():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
					
		self.bubbleFile = os.path.join(self.reader.outputDir, "Writhe.res")

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

		print("")
		print("")
		print("")
		print("")
		print("")
		print("")

		print("step "+str(i) )

		data = next(self.reader)
		bubbles=[]
		branches1=[]
		branches2=[]
		
		for i in range(0,self.Nchain):
			if(data.fork[i]==-1):
				print("start bubble at "+ str(i))
				current_bubble=[]
				current_branch1=[]
				current_branch2=[]

				current_bubble.append(data.polyPos[i])
				k=i+1
				while(data.fork[k]!=1 and k<self.Nchain):
					current_bubble.append(data.polyPos[k])
					current_branch1.append(data.polyPos[k])
					current_branch2.append(data.polyPos[data.SisterID[k]])
					
					k+=1
				if(k<self.Nchain):
					branches1.append(current_branch1)
					branches2.append(current_branch2)

					current_bubble.append(data.polyPos[k])
					print("end bubble at "+ str(k))
					k-=1		
					while(k!=i):
						current_bubble.append(data.polyPos[data.SisterID[k]])
						k-=1
					bubbles.append(current_bubble)
					print("bubble size "+ str(len(current_bubble)))


		print( "I have " +str(len(bubbles))+ "bubble")
			
		
		
		for bubble in bubbles:
			for i in range(1,len(bubble)):
				diff=bubble[i]-bubble[i-1]
				if(np.sqrt(np.dot(diff.T,diff))>0.71):
					print("ERROR check")
		
			
		#create a segment array
		for id,bubble in enumerate(bubbles):
			segments=[]
			bubble_corr=[]
			chain=[]
			bubble_corr.append(bubble[0])
			for i in range(1,len(bubble)):
				diff=bubble[i]-bubble[i-1]
				if(np.sqrt(np.dot(diff.T,diff))!=0):
					bubble_corr.append(bubble[i])
					if(i<=len(bubble)/2):
						chain.append(-1)
					else:
						chain.append(1)
					
			for i in range(len(bubble_corr)-1):
				segments.append([bubble_corr[i],bubble_corr[i+1]])
			segments.append([bubble_corr[-1],bubble_corr[0]])

			branch1=[]
			branch2=[]
			
			
			for i in range(len(branches1[id])-1):
				branch1.append([branches1[id][i],branches1[id][i+1]])
			for i in range(len(branches2[id])-1):
				branch2.append([branches2[id][i],branches2[id][i+1]])
				
			
			
			
			

			
			
			#print(len(segments))
			writhe=0
			for i in range(2,len(segments)):
				for j in range(i):
					if(0==1):
						r13= -np.array(segments[i][0]-segments[j][0])
						r14= -np.array(segments[i][0]-segments[j][1])
						r24= -np.array(segments[i][1]-segments[j][1])
						r23= -np.array(segments[i][1]-segments[j][0])
						r12= np.array(segments[i][1]-segments[i][0])
						r34= np.array(segments[j][1]-segments[j][0])
				
						n1=np.cross(r13,r14)/np.sqrt(np.dot(np.cross(r13,r14).T,np.cross(r13,r14)))
						n2=np.cross(r14,r24)/np.sqrt(np.dot(np.cross(r14,r24).T,np.cross(r14,r24)))
						n3=np.cross(r24,r23)/np.sqrt(np.dot(np.cross(r24,r23).T,np.cross(r24,r23)))
						n4=np.cross(r23,r13)/np.sqrt(np.dot(np.cross(r23,r13).T,np.cross(r23,r13)))
					
		
						omega_star=np.nansum([np.arcsin(np.dot(n1,n2)),np.arcsin(np.dot(n2,n3)),np.arcsin(np.dot(n3,n4)),np.arcsin(np.dot(n4,n1))])
						omega=omega_star*np.sign(np.dot(np.cross(r34,r12),r13))/(4*np.pi)
						writhe+=omega


			writhe=2*writhe
			print("writhe is " +str(writhe))

			Lk=0
			for i in range(1,len(branch1)):
				for j in range(len(branch2)-1):
					if(0==1):
						r13= -np.array(branch1[i][0]-branch2[j][0])
						r14= -np.array(branch1[i][0]-branch2[j][1])
						r24= -np.array(branch1[i][1]-branch2[j][1])
						r23= -np.array(branch1[i][1]-branch2[j][0])
						r12= np.array(branch1[i][1]-branch1[i][0])
						r34= np.array(branch2[j][1]-branch2[j][0])
				
						n1=np.cross(r13,r14)/np.sqrt(np.dot(np.cross(r13,r14).T,np.cross(r13,r14)))
						n2=np.cross(r14,r24)/np.sqrt(np.dot(np.cross(r14,r24).T,np.cross(r14,r24)))
						n3=np.cross(r24,r23)/np.sqrt(np.dot(np.cross(r24,r23).T,np.cross(r24,r23)))
						n4=np.cross(r23,r13)/np.sqrt(np.dot(np.cross(r23,r13).T,np.cross(r23,r13)))
					
		
						omega_star=np.nansum([np.arcsin(np.dot(n1,n2)),np.arcsin(np.dot(n2,n3)),np.arcsin(np.dot(n3,n4)),np.arcsin(np.dot(n4,n1))])

						omega=omega_star*np.sign(np.dot(np.cross(r34,r12),r13))/(4*np.pi)

						Lk+=omega


			#print("Lk is " +str(Lk))
				
			#print(len(branch1))
			#print(len(branch2))

			Lk=0
			for i in range(len(branch1)):
				for j in range(len(branch2)):
					if(0==0):
						dr1=np.array(branch1[i][1]-branch1[i][0])
						dr2=np.array(branch2[j][1]-branch2[j][0])
						r12=-np.array(branch1[i][0]-branch2[j][0])
						omega=(1/(4*np.pi))*np.dot(np.cross(dr1,dr2),r12)/(np.dot(r12.T,r12))**3/2

						#if((np.dot(r12.T,r12))!=0):
						Lk+=omega
			
			
			print("Lk is " +str(Lk))

			
	
			#step[4]=step[4]/(step[0]*(len(cis_bubbles)-step[0]))
						
		
		

	def Print(self):

		np.savetxt(self.bubbleFile,self.final_data)
		

		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame r\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	writhe = Writhe(outputDir, initFrame=initFrame)

	writhe.Compute()
	writhe.Print()
