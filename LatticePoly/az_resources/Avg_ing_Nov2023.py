##
# radavg.py
# LatticePoly
##
# 
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import sys
import numpy as np



dpath = '/home/aabdulla/3D/Nov2023/LatticePoly/LatticePoly/data/NExtruders/'
Nchain = 10240
Nmeas = 1000

def fetchFiles(dirPath, fileName):
	fileList = []

	for root, dirs, files in os.walk(dirPath):
		if not dirs:
			if fileName in files:
				filePath = os.path.join(root, fileName)
				fileList.append(filePath)
	
	return fileList

fileList_1 = fetchFiles(dpath, 'homDensity.res')
fileList_2 = fetchFiles(dpath, 'polyDensity.res')

fileList_1 = sorted(fileList_1)
fileList_2 = sorted(fileList_2)

np.savetxt('/home/aabdulla/3D/Nov2023/LatticePoly/LatticePoly/data/NExtruders/filelist1.txt', fileList_1, fmt='%s')
np.savetxt('/home/aabdulla/3D/Nov2023/LatticePoly/LatticePoly/data/NExtruders/filelist2.txt', fileList_2, fmt='%s')

# Create a function to read .res file Nchain by Nmeas numpy matrix, average the rows, and return a Nmeas numpy array
def readResFile(filePath):
    data = np.loadtxt(filePath)
    data = np.mean(data, axis=0)
    return data

def SubResFile(filePath1, filePath2):
    data1 = np.loadtxt(filePath1)
    data2 = np.loadtxt(filePath2)
    data = np.mean(data2-data1, axis=0)
    return data

for i in range(0,len(fileList_1)):
   endtoendFile1  = os.path.join(dpath, "avgHomogeneity_%d.res" %i)
   endtoendFile2 = os.path.join(dpath, "avgDensity_%d.res" %i)
   endtoendFile3 = os.path.join(dpath, "avgHeterogeneity_%d.res" %i)
   dt = readResFile(fileList_1[i])
   du = readResFile(fileList_2[i])
   dv = SubResFile(fileList_1[i], fileList_2[i])
   np.savetxt(endtoendFile1, dt)
   np.savetxt(endtoendFile2, du)
   np.savetxt(endtoendFile3, dv)
   print(i)

print("\033[1;32mChain1:Printed avg homogeneity '%s'\033[0m" %endtoendFile1)

