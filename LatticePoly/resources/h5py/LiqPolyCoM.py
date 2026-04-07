##
##  PolyMSD.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys
import psutil

import numpy as np

from utils import msdFFT
from hdf5Reader import hdf5Reader
import h5py

class LiqPolyCoM():

        def __init__(self, inputDir : str):
                print(f"LiqPolyCoM : Init {inputDir}")
                self.reader = hdf5Reader(inputDir, 'traj.h5', -1, readLiq=True, readPoly=True)
                self.processPath = os.path.join(inputDir, "process.h5")
        
        
        def Compute(self):
                print("\n")
                
                vMem = psutil.virtual_memory()
                
                sizeTot = self.reader.N * (self.reader.polyPos.nbytes + self.reader.liqPos.nbytes)
                self.liqPosInit = self.reader.liqPos
 
                if sizeTot < vMem.available:
                
                        posHist = self.ReadHist()
                        self.LiqPolyCoM = msdFFT(np.mean(posHist,1))                
                else:
                        print("Memory overflow likely - reduce chosen number of frames")
                        sys.exit()


        def ReadHist(self):
                posHist = np.zeros((self.reader.N, self.reader.nTad+self.reader.nLiq, 3), dtype=np.float32)
                
                for i in range(self.reader.N):
                        data = next(self.reader)
                        posHist[i] = np.concatenate((data.polyPos,self.liqPosInit + data.liqDisp),axis=0)
                        
                        if (i+1) % 10 == 0:
                                        print("Processed %d out of %d frame" % (i+1, self.reader.N))
                        
                return posHist
        

        def Print(self):
                print("\n")
                
                self.reader.Close()
                
                with h5py.File(self.processPath, 'a') as processFile:
                        
                        self.PrintDataset(processFile, "LiqPolyCoM", data = self.LiqPolyCoM)
                       

        def PrintDataset(self, processFile : h5py.File, dataset_name : str, data : np.ndarray[tuple[int, ...], np.dtype[np.int32 | np.float32]]):
                
                if dataset_name in processFile.keys():
                        tmp = cast(h5py.Dataset, processFile[dataset_name])
                        tmp[:] = data
                else:
                        processFile.create_dataset(dataset_name, data = data)
                        
                print(f"Dataset {dataset_name} printed")
        
        
if __name__ == "__main__":
        if len(sys.argv) != 2:
                print("\033[1;31mUsage is %s inputDir\033[0m" % sys.argv[0])
                sys.exit()
                
        inputDir = sys.argv[1]

        CoM = LiqPolyCoM(inputDir)

        CoM.Compute()
        CoM.Print()
        
        print("\n")
        print("LiqPolyCoM : Done\n\n")        