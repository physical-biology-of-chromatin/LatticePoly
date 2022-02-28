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

fileList_1 = fetchFiles('/home/aabdulla/3D/LatticePoly/LatticePoly/data/201kb_20MB/0.0000/', 'Monomdist.res')
Nmeas    = 1000
Nchain   = 201
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


np.savetxt('/home/aabdulla/3D/LatticePoly/LatticePoly/data/201kb_20MB/0.0000/dist_kb.res', dist_kb)
