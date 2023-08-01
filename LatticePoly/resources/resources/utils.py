##
##  utils.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import numpy as np
	

def getInputParam(key, paramFile):
	with open(paramFile, "r") as f:
		keyLines = [line for line in f if key in line]
		
		if len(keyLines) == 0:
			print("Could not locate parameter '%s' in file '%s'" % (key, paramFile))
			
			return None
		
		else:
			keyLine = keyLines[0]
			valLine = keyLine.split(';')[0]
			
			val = valLine.split('=')[1]
			
			return val.strip()


def fetchFiles(dirPath, fileName):
	fileList = []

	for root, dirs, files in os.walk(dirPath):
		if not dirs:
			if fileName in files:
				filePath = os.path.join(root, fileName)
				fileList.append(filePath)
	
	return fileList
	
		
# Adapted from Calandrini et al. (https://doi.org/10.1051/sfn/201112010)
def msdFFT(posHist):
	N = posHist.shape[0]
	
	sqDist = np.square(posHist).sum(axis=1)
	sqDist = np.append(sqDist, 0)
	
	S2 = sum([_autoCorrFFT(posHist[:,i]) for i in range(3)])
	Q = 2*sqDist.sum()
	
	S1 = np.zeros(N)
	
	for m in range(N):
		Q -= sqDist[m-1] + sqDist[N-m]
		S1[m] = Q / (N-m)
		
	return S1 - 2*S2


def _autoCorrFFT(x):
	N = x.shape[0]

	F = np.fft.fft(x, n=2*N)
	PSD = F * F.conjugate()

	autoCorr = np.fft.ifft(PSD)

	autoCorr = autoCorr[:N]
	autoCorr = autoCorr.real

	n = N*np.ones(N) - np.arange(N)

	return autoCorr / n
