##
##  LiqDensity.py
##  LatticePoly
##
##  Created by ppuel on 18/10/2024.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from hdf5Reader import hdf5Reader
import h5py

class LiqDensity():

        def __init__(self, inputDir, outputDir):
                self.reader = hdf5Reader(inputDir, 'traj.h5', -1, readLiq=True, readPoly=False, backInBox=False)
                self.processPath = os.path.join(outputDir, "process.h5")
                

        def Compute(self):
                self.ldensMean = np.zeros(self.reader.N, dtype=np.float32)
                self.ldensStd = np.zeros(self.reader.N, dtype=np.float32)
                self.ldensHist = np.zeros((self.reader.N,13), dtype=np.int32)

                for i in range(self.reader.N):
                        self.ProcessFrame(i)
                        
                       
                                        
        def ProcessFrame(self, i):                                                
                       
                data = next(self.reader)

                self.ldensMean[i] = data.liqDens.mean()

                self.ldensStd[i] = np.square(data.liqDens - self.ldensMean[i]).sum()

                for j in np.int8((data.liqDens+0.001)*12):
                        self.ldensHist[i][j] += 1

        def Print(self):
                self.reader.Close()
                with h5py.File(self.processPath, 'a') as processFile:
                        
                        self.PrintDataset(processFile, "liqMean", data = self.ldensMean)
                        self.PrintDataset(processFile, "liqSTD", data = np.sqrt(self.ldensStd / self.reader.nLiq))
                        self.PrintDataset(processFile, "liqHist", data = self.ldensHist)
                
        def PrintDataset(self, processFile, dataset_name, data):
                
                if dataset_name in processFile.keys():
                        tmp = processFile[dataset_name]
                        tmp[:] = data
                else:
                        processFile.create_dataset(dataset_name, data = data)
               

if __name__ == "__main__":
        if len(sys.argv) != 3:
                print("\033[1;31mUsage is %s inputDir outputDir\033[0m" % sys.argv[0])
                sys.exit()

        inputDir = sys.argv[1]
        outputDir = sys.argv[2]

        density = LiqDensity(inputDir, outputDir)

        density.Compute()
        density.Print()
        