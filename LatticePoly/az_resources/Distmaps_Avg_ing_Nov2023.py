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



dpath = '/home/aabdulla/3D/Nov2023/LatticePoly/LatticePoly/data/output/'
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

fileList_1 = fetchFiles(dpath, 'r950_cProbTopo.res')
fileList_2 = fetchFiles(dpath, 'r950_distmatTopo.res')

fileList_1 = sorted(fileList_1)
fileList_2 = sorted(fileList_2)

np.savetxt('/home/aabdulla/3D/Nov2023/LatticePoly/LatticePoly/data/listofFiles.txt', fileList_1, fmt='%s')


def readmat064(fileList):
    matrix = np.zeros((Nchain, Nchain), dtype=np.float32)
    for i in range(0,64):
        matrix += np.loadtxt(fileList[i])
        print(i)
    matrix = matrix/Nmeas
    return matrix

def readmat64128(fileList):
    matrix = np.zeros((Nchain, Nchain), dtype=np.float32)
    for i in range(64,128):
        matrix += np.loadtxt(fileList[i])
    matrix = matrix/Nmeas
    return matrix

def readmat128192(fileList):
    matrix = np.zeros((Nchain, Nchain), dtype=np.float32)
    for i in range(128,192):
        matrix += np.loadtxt(fileList[i])
    matrix = matrix/Nmeas
    return matrix        

ContactProb064 = np.zeros((Nchain, Nchain), dtype=np.float32)
ContactProb64128 = np.zeros((Nchain, Nchain), dtype=np.float32)
#ContactProb128192 = np.zeros((Nchain, Nchain), dtype=np.float32)

ContactProb064 = readmat064(fileList_1)
ContactProb64128 = readmat64128(fileList_1)
#ContactProb128192 = readmat128192(fileList_1)

endtoendFile1  = os.path.join(dpath, "avgContactProb_Ext0_Topo0.res")
endtoendFile2  = os.path.join(dpath, "avgContactProb_Ext1000_Topo0.res")
#endtoendFile3  = os.path.join(dpath, "avgHomogeneity_%d.res" %i)

np.savetxt(endtoendFile1, ContactProb064)
np.savetxt(endtoendFile2, ContactProb64128)
#np.savetxt(endtoendFile3, ContactProb128192)

print("\033[1;32mChain1:Printed avg contact probability '%s'\033[0m" %endtoendFile1)


