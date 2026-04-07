##
##  PolyGyration.py
##  LatticePoly
##
##  Created by mtortora on 02/08/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

import os
import sys
import h5py

import numpy as np

from hdf5Reader import hdf5Reader

class PolyGyration():
        def __init__(self, inputDir : str):
                print(f"PolyGyration : Init {inputDir}")
                self.reader = hdf5Reader(inputDir, 'traj.h5', -1, readLiq=False, readPoly=True)
                self.processPath = os.path.join(inputDir, "process.h5")


        def Compute(self):
                print("\n")
                self.polyAniso = np.zeros(self.reader.N, dtype=np.float32)
                self.polyGyration = np.zeros(self.reader.N, dtype=np.float32)

                for i in range(self.reader.N):
                        self.ProcessFrame(i)
                        
                        if (i+1) % 10 == 0:
                                print("Processed %d out of %d frame" % (i+1, self.reader.N))


        def ProcessFrame(self, i : int):
                data = next(self.reader)
                norm = 0
                
                for d in self.reader.domains:
                        if d.size > 2:
                                pos = data.polyPos[d]
                                pos -= pos.mean(axis=0, keepdims=True)
                        
                                diag = np.linalg.svd(pos, compute_uv=False) * np.sqrt(12)/d.size
                                r2_gyr = np.square(diag).sum(axis=-1)
                        
                                r_gyr = np.sqrt(r2_gyr)
                                aniso = 3/2.*(diag**4).sum(axis=-1)/r2_gyr**2 - 1/2.
                        
                                norm += d.size
                                
                                self.polyAniso[i] += aniso * d.size
                                self.polyGyration[i] += r_gyr * d.size
                                                                        
                self.polyAniso[i] /= norm if norm > 0 else 1
                self.polyGyration[i] /= norm if norm > 0 else 1


        def Print(self):
                print("\n")
                
                self.reader.Close()
                
                with h5py.File(self.processPath, 'a') as processFile:
                        
                        self.PrintDataset(processFile, "polyAniso", data = self.polyAniso)
                        self.PrintDataset(processFile, "polyGyration", data = self.polyGyration)
                

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

        gyr = PolyGyration(inputDir)

        gyr.Compute()
        gyr.Print()
        
        print("\n")
        print("PolyGyration : Done\n\n")