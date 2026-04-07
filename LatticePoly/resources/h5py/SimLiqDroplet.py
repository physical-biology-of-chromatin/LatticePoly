##
##  LiqPolyCoincidence.py
##  LatticePoly
##Liq
##  Created by ppuel on ../../2023.
##  Copyright © 2023 ENS Lyon. All rights reserved.
##

import os
import pickle
import sys

import h5py
import numpy as np
import utils
from hdf5Reader import hdf5Reader
from LiqDroplet import Droplet, Event, LifeTime


class SimLiqDroplet:
    def __init__(self, inputDir, rMax):
        print(f"SimLiqDroplet : Init {inputDir} {rMax}")

        self.inputDir = inputDir
        self.rMax = rMax

        self.N = len(os.listdir(os.path.join(inputDir, "N/")))

        reader = hdf5Reader(os.path.join(inputDir, "N/0/"), "traj.h5", read_liq=True)
        self.Nmeas = reader.n_frame
        self.nLiq = reader.n_liq
        reader.close()

    def Process(self):
        print("\n")

        liqFraction = np.zeros(self.Nmeas)
        liqDropNum = np.zeros(self.Nmeas)
        liqDropSize = np.zeros(self.Nmeas)

        for n in range(self.N):
            if n != 6:
                with open(os.path.join(self.inputDir, f"N/{n}/liq_droplets.pickle"), "rb") as pfile:
                    dropletDict = pickle.load(pfile)
        
                for droplet in dropletDict.values():
                    rDroplet = ((droplet.tau - 1) ** 2 
                        + (np.max(droplet.sizes) - 2) ** 2
                        ) ** (1 / 2)
        
                    if rDroplet > self.rMax:
                        frame = np.fromiter(
                                map(lambda x: x[0], droplet.nodes), dtype=np.int32
                        )
        
                        liqFraction[frame] += droplet.sizes
                        liqDropNum[frame] += 1
        
                if (n + 1) % 10 == 0:
                        print("Process %d out of %d trajectories" % (n + 1, self.N))

        liqDropSize = (
            np.where(liqDropNum > 0, liqFraction / liqDropNum, np.nan) / self.N
        )
        liqFraction /= self.N * self.nLiq
        liqDropNum /= self.N

        with h5py.File(
            os.path.join(self.inputDir, "aggregated_process.h5"), "a"
        ) as processFile:
            utils.PrintDataset(processFile, "liqFraction", data=liqFraction)
            utils.PrintDataset(processFile, "liqDropNum", data=liqDropNum)
            utils.PrintDataset(processFile, "liqDropSize", data=liqDropSize)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\033[1;31mUsage is %s inputDir rMax\033[0m" % sys.argv[0])
        sys.exit()

    inputDir = sys.argv[1]

    rMax = float(sys.argv[2])

    SimDroplet = SimLiqDroplet(inputDir, rMax)

    SimDroplet.Process()

    print("\n")
    print("SimLiqDroplet : Done\n\n")
