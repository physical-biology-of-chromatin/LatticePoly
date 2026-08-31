##
##  PolyGyration.py
##  LatticePoly
##
##  Created by mtortora on 02/08/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from hdf5Reader import hdf5Reader
import h5py


class PolyTad():
        def __init__(self, inputDir, outputDir):
                print(f"PolyTad : Init {inputDir} {outputDir}")
                self.reader = hdf5Reader(inputDir, "traj.h5", -1, readLiq=False, readPoly=True)
                self.filePath = os.path.join(outputDir, "process.h5")
                
                self.NPcG = len(self.indata_PcG)
                self.NNoPcG = len(self.indata_NoPcG)
                
                self.PosPcG = np.zeros((self.reader.N, self.NPcG, 3), dtype=np.float32)
                self.PosNoPcG = np.zeros((self.reader.N, self.NNoPcG, 3), dtype=np.float32)


        def Compute(self):
                print("\n")
                self.PcGmeanDistance = np.zeros((self.reader.N, self.NPcG*(self.NPcG-1)//2), dtype=np.float32)
                self.PcGtoNoPcGmeanDistance = np.zeros((self.reader.N, self.NNoPcG*self.NPcG), dtype=np.float32)
                
                for i in range(self.reader.N):
                        self.ProcessFrame(i)
                        
                        if (i+1) % 10 == 0:
                                print("PolyTad : Processed %d out of %d configurations" % (i+1, self.reader.N))


        def ProcessFrame(self, i):
                data = next(self.reader)
                
                for ids, values in enumerate(self.indata_PcG.values()):
                        start, stop = values
                        self.PosPcG[i, ids] = np.mean(data.polyPos[start:stop], axis=0)
                for ids, values in enumerate(self.indata_NoPcG.values()):
                        start, stop = values
                        self.PosNoPcG[i, ids] = np.mean(data.polyPos[start:stop], axis=0)
                
                idsDistance = 0
                for PcG in range(self.NPcG):
                        for NoPcG in range(self.NNoPcG):
                                self.PcGtoNoPcGmeanDistance[i, idsDistance] = np.linalg.norm(self.PosPcG[i,PcG] - self.PosNoPcG[i,NoPcG])
                                idsDistance += 1
                
                idsDistance = 0
                for PcG in range(self.NPcG):
                        for PcG_bis in range(PcG+1, self.NPcG):
                                self.PcGmeanDistance[i, idsDistance] = np.linalg.norm(self.PosPcG[i,PcG] - self.PosPcG[i,PcG_bis])
                                idsDistance += 1


        def Print(self):
                print("\n")
                self.reader.Close()
                self.hfile = h5py.File(self.filePath, 'a')
                
                self.PrintDataset("PcGmeanDistanceHist", data = self.PcGmeanDistance)
                self.PrintDataset("PcGtoNoPcGmeanDistanceHist", data = self.PcGtoNoPcGmeanDistance)
                
                self.hfile.close()


        def PrintDataset(self, dataset_name, data):
                
                if dataset_name in self.hfile.keys():
                        tmp = self.hfile[dataset_name]
                        tmp[:] = data
                else:
                        self.hfile.create_dataset(dataset_name, data = data)
        
                print(f"Dataset {dataset_name} printed")


if __name__ == "__main__":
        if len(sys.argv) != 3:
                print("\033[1;31mUsage is %s inputDir outputDir\033[0m" % sys.argv[0])
                sys.exit()

        inputDir = sys.argv[1]
        outputDir = sys.argv[2]

        Tad = PolyTad(inputDir, outputDir)

        Tad.Compute()
        Tad.Print()
        
        print("\n")
        print("PolyTad : Done\n\n")
