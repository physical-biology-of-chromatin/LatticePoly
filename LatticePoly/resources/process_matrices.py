import os
import vtkReader
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import math

"""
vmax=0
averagematrices1=[]
for step in range(0,1001):
	print(step)
	matrices_step=np.zeros((100,100))
	for i in range(0,1):
		data = np.loadtxt('/Volumes/My Passport/Yeast360/pairingtrial/ref/Matrixtrans'+str(step)+'.txt')
		matrices_step+=data
	vmax=np.max(matrices_step)
	averagematrices1.append(matrices_step)
#averagematrices1.append(matrices_step/50)

os.makedirs('/Volumes/My Passport/Yeast360/pairingtrial/ref/averagetrans/')

for i in range(len(averagematrices1)):
	fig, ax = plt.subplots()
	im = ax.imshow(np.log(averagematrices1[i]+0.1),vmin=math.log(0.1), vmax=math.log(vmax))
	fig.colorbar(im)
	plt.suptitle('Simulation step = '+str(i),size=15)
	plt.savefig('/Volumes/My Passport/Yeast360/pairingtrial/ref/averagetrans/%04d.png' %i)
"""
vmax=0
averagematrices1=[]
for step in range(0,1001):
	print(step)
	matrices_step=np.zeros((100,100))
	for i in range(0,1):
		data = np.loadtxt('/Volumes/My Passport/Yeast360/pairingtrial/10kT/Matrixtrans'+str(step)+'.txt')
		matrices_step+=data
	if(vmax<np.max(matrices_step)):
		vmax=np.max(matrices_step)
	averagematrices1.append(matrices_step)
#averagematrices1.append(matrices_step/50)

os.makedirs('/Volumes/My Passport/Yeast360/pairingtrial/10kT/averagetrans/')

for i in range(len(averagematrices1)):
	fig, ax = plt.subplots()
	im = ax.imshow(np.log(averagematrices1[i]+0.1),vmin=math.log(0.1), vmax=math.log(vmax))
	fig.colorbar(im)
	plt.suptitle('Simulation step = '+str(i),size=15)
	plt.savefig('/Volumes/My Passport/Yeast360/pairingtrial/10kT/averagetrans/%04d.png' %i)

averagematrices1=[]
for step in range(0,1001):
	print(step)
	matrices_step=np.zeros((100,100))
	for i in range(0,1):
		data = np.loadtxt('/Volumes/My Passport/Yeast360/pairingtrial/10mon 10kT/Matrixtrans'+str(step)+'.txt')
		matrices_step+=data
	if(vmax<np.max(matrices_step)):
		vmax=np.max(matrices_step)
	averagematrices1.append(matrices_step)
#averagematrices1.append(matrices_step/50)

os.makedirs('/Volumes/My Passport/Yeast360/pairingtrial/10mon 10kT/averagetrans/')

for i in range(len(averagematrices1)):
	fig, ax = plt.subplots()
	im = ax.imshow(np.log(averagematrices1[i]+0.1),vmin=math.log(0.1), vmax=math.log(vmax))
	fig.colorbar(im)
	plt.suptitle('Simulation step = '+str(i),size=15)
	plt.savefig('/Volumes/My Passport/Yeast360/pairingtrial/10mon 10kT/averagetrans/%04d.png' %i)


