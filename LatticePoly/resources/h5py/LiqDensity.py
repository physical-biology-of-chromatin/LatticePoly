##
##  LiqDensity.py
##  LatticePoly
##
##  Created by ppuel on 18/10/2024.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys
from typing import cast

import h5py
import numpy as np
from hdf5Reader import hdf5Reader


class LiqDensity():

class LiqDensity:
    def __init__(self, inputDir: str):
        print(f"LiqDensity : Init {inputDir}")
        self.reader = hdf5Reader(
            inputDir, "traj.h5", -1, readLiq=True, readPoly=False, backInBox=False
        )
        self.processPath = os.path.join(inputDir, "process.h5")

    def Compute(self):
        print("\n")
        self.ldensMean = np.zeros(self.reader.N, dtype=np.float32)
        self.ldensStd = np.zeros(self.reader.N, dtype=np.float32)
        self.ldensHist = np.zeros((self.reader.N, 13), dtype=np.int32)

        for i in range(self.reader.N):
            self.ProcessFrame(i)

            if (i + 1) % 10 == 0:
                print("Process %d out of %d trajectories" % (i + 1, self.reader.N))

    def ProcessFrame(self, i):

        data = next(self.reader)

        self.ldensMean[i] = data.liqDens.mean()

        self.ldensStd[i] = np.square(data.liqDens - self.ldensMean[i]).sum()

        for j in np.asarray((data.liqDens + 0.001) * 12, dtype=np.int32):
            self.ldensHist[i][j] += 1

    def Print(self):
        print("\n")
        self.reader.Close()
        with h5py.File(self.processPath, "a") as processFile:
            self.PrintDataset(processFile, "liqMean", data=self.ldensMean)
            self.PrintDataset(
                processFile, "liqSTD", data=np.sqrt(self.ldensStd / self.reader.nLiq)
            )
            self.PrintDataset(processFile, "liqHist", data=self.ldensHist)

    def PrintDataset(
        self,
        processFile: h5py.File,
        dataset_name: str,
        data: np.ndarray[tuple[int, ...], np.dtype[np.float32 | np.int32]],
    ):

        if dataset_name in processFile.keys():
            tmp = cast(h5py.Dataset, processFile[dataset_name])
            tmp[:] = data
        else:
            processFile.create_dataset(dataset_name, data=data)

        print(f"Dataset {dataset_name} printed")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("\033[1;31mUsage is %s inputDir\033[0m" % sys.argv[0])
        sys.exit()

    inputDir = sys.argv[1]

    density = LiqDensity(inputDir)

    density.Compute()
    density.Print()

    print("\n")
    print("LiqDensity : Done\n\n")
