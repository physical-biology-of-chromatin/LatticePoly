##
##  LiqPolyCoincidence.py
##  LatticePoly
##Liq
##  Created by ppuel on ../../2023.
##  Copyright © 2023 ENS Lyon. All rights reserved.
##

import os
import sys

import h5py
import numpy as np
import utils


class SimAggregateData:
    def __init__(self, inputDir, datasetNameList):
        print(f"SimAggregateData : Init {inputDir} {datasetNameList}")

        self.inputDir = inputDir
        self.datasetNameList = datasetNameList

        self.N = len(os.listdir(os.path.join(self.inputDir, "N/")))

        self.datasetShape = {}

        with h5py.File(os.path.join(self.inputDir, "N/0/process.h5"), "r") as tmpfile:
            for datasetName in datasetNameList:
                self.datasetShape[datasetName] = np.size(tmpfile[datasetName], axis=1)

    def Process(self):
        print("\n")

        self.datasetArray = {}

        for datasetName in datasetNameList:
            self.datasetArray[datasetName] = np.zeros(self.datasetShape[datasetName])

        for n in range(self.N):
            if n != 6:
                print(n)
                with h5py.File(
                    os.path.join(self.inputDir, f"N/{n}/process.h5"), "r"
                ) as dataFile:
                    for datasetName in datasetNameList:
                        self.datasetArray[datasetName] += np.mean(
                            dataFile[datasetName][-10:], axis=0
                        )

            if (n + 1) % 10 == 0:
                print("Process %d out of %d trajectories" % (n + 1, self.N))

    def Print(self):

        with h5py.File(
            os.path.join(self.inputDir, "aggregated_process.h5"), "a"
        ) as aggregatedFile:
            for datasetName in datasetNameList:
                utils.PrintDataset(
                    aggregatedFile,
                    datasetName + "_2de_axis",
                    data=self.datasetArray[datasetName] / (self.N - 1),
                )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(
            "\033[1;31mUsage is %s inputDir datasetName_1 ... datasetName_n\033[0m"
            % sys.argv[0]
        )
        sys.exit()

    inputDir = sys.argv[1]
    if len(sys.argv) == 3:
        datasetNameList = [sys.argv[2]]
    else:
        datasetNameList = sys.argv[2:]

    SimAggregate = SimAggregateData(inputDir, datasetNameList)

    SimAggregate.Process()
    SimAggregate.Print()

    print("\n")
    print("SimAggregateData : Done\n\n")
