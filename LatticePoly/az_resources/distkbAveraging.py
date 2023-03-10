##
# distkbAveraging.py
# LatticePoly
##
# 
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import numpy as np




def fetchFiles(dirPath, fileName):
	fileList = []

	for root, dirs, files in os.walk(dirPath):
		if not dirs:
			if fileName in files:
				filePath = os.path.join(root, fileName)
				fileList.append(filePath)
	
	return fileList

fileList_1 = fetchFiles('/home/adminamith/3D/LatticePoly/LatticePoly/data/output/Living/Nucleation/', 'contactProb_z.res')
Nmeas    = 200
Nchain   = 120
dist_avg = np.zeros((Nchain,Nchain))
for i in range(0,len(fileList_1)):
    dist_avg = dist_avg + np.loadtxt(fileList_1[i])

dist_avg = dist_avg/len(fileList_1)   

dist_kb  = np.zeros(Nchain) 
for i in range(0, len(dist_avg)):
    for j in range(i, len(dist_avg)):
        dist_kb[j-i] = dist_kb[j-i] + dist_avg[i,j]

for i in range(0,len(dist_kb)):
    dist_kb[i] = dist_kb[i]/(Nchain-i)

np.savetxt('/home/adminamith/3D/LatticePoly/LatticePoly/data/output/Living/Nucleation/CProb_NucF.res', dist_avg)
np.savetxt('/home/adminamith/3D/LatticePoly/LatticePoly/data/output/Living/Nucleation/Cdecay_NucF.res', dist_kb)
