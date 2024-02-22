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



dpath = '/Users/amith/Topo/ExtrusionData/Phase/'
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

def fetchFilesStart(dirPath, fileName):
    fileList = []

    for root, dirs, files in os.walk(dirPath):
        if not dirs:
            for file in files:   
                if file.startswith(fileName):
                    filePath = os.path.join(root, file)
                    fileList.append(filePath)
    
    return fileList


fileList_1 = fetchFilesStart(dpath, 'avgDensity')
fileList_2 = fetchFilesStart(dpath, 'avgHomogeneity')
fileList_3 = fetchFilesStart(dpath, 'avgHeterogeneity')

def extract_number(filename):
    return int(''.join(filter(str.isdigit, filename)))

fileList_1 = sorted(fileList_1, key=extract_number)
fileList_2 = sorted(fileList_2, key=extract_number)
fileList_3 = sorted(fileList_3, key=extract_number)


np.savetxt('/Users/amith/Topo/ExtrusionData/Phase/fileList_1.res', fileList_1, fmt='%s')

matrixDensity       = np.zeros((11, 11), dtype=np.float32)
matrixHeterogeneity = np.zeros((11, 11), dtype=np.float32)
matrixHomogeneity   = np.zeros((11, 11), dtype=np.float32)

k = 0
for i in range(0,11):
    for j in range(0,11):
        data1 = np.loadtxt(fileList_1[k])
        data2 = np.loadtxt(fileList_2[k])
        data3 = np.loadtxt(fileList_3[k])
        matrixDensity[i,j] = data1[Nmeas]
        matrixHomogeneity[i,j] = data2[Nmeas]
        matrixHeterogeneity[i,j] = data3[Nmeas]
        k = k+1

np.savetxt('/Users/amith/Topo/ExtrusionData/Phase/matrixDensity.res', matrixDensity)  
np.savetxt('/Users/amith/Topo/ExtrusionData/Phase/matrixHeterogeneity.res', matrixHeterogeneity)  
np.savetxt('/Users/amith/Topo/ExtrusionData/Phase/matrixHomogeneity.res', matrixHomogeneity)
