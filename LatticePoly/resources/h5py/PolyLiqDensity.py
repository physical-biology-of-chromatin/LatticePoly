##
##  LiqDensity.py
##  LatticePoly
##
##  Created by ppuel on 18/10/2024.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##
##  MeM usage : ~ 250 MiB

import os
import sys

import h5py

# import networkx as nx
# import numba
import numpy as np
from hdf5Reader import hdf5Reader
from scipy.stats import pearsonr

# from memory_profiler import profile


class PolyLiqDensity:
    def __init__(self, inputDir: str, nCut: int):

        print(f"PolyLiqDensity : Init {inputDir} {nCut}")

        self.reader = hdf5Reader(
            inputDir,
            "traj.h5",
            init_frame=-1,
            read_liq=True,
            read_poly=True,
            back_in_box=True,
        )
        self.processPath = os.path.join(inputDir, "process.h5")

        self.nCut = nCut
        self.L = self.reader.box_dim[0]
        self.cut = self.L // self.nCut
        self.Ncube = self.nCut**3
        self.Nsite_per_cut = self.L**3 * 4 // self.Ncube

        self.Cube = np.zeros((self.Ncube, 6), dtype=np.int32)

        self.x_y_z = np.array([1, self.nCut, self.nCut**2], dtype=np.int32)

        self.list_threshold = [0, 1, 2, 3, 4, 5, 10, 20]
        self.Nthreshold = len(self.list_threshold)

        for vi in range(self.Ncube):
            iz = vi // self.nCut**2
            iy = vi // self.nCut - self.nCut * iz
            ix = vi - self.nCut * (self.nCut * iz + iy)

            self.Cube[vi, 0] = self.iCube((ix + 1) % self.nCut, iy, iz)
            self.Cube[vi, 1] = self.iCube((ix - 1) % self.nCut, iy, iz)
            self.Cube[vi, 2] = self.iCube(ix, (iy + 1) % self.nCut, iz)
            self.Cube[vi, 3] = self.iCube(ix, (iy - 1) % self.nCut, iz)
            self.Cube[vi, 4] = self.iCube(ix, iy, (iz + 1) % self.nCut)
            self.Cube[vi, 5] = self.iCube(ix, iy, (iz - 1) % self.nCut)

    def iCube(self, ix: int, iy: int, iz: int) -> np.int32:

        return np.sum(np.multiply(np.array([ix, iy, iz], dtype=np.int32), self.x_y_z))

    def Compute(self):
        print("\n")
        self.PolyVolume = np.zeros(
            (self.reader.n_frame, self.Nthreshold), dtype=np.float32
        )
        self.LiqVolume = np.zeros(
            (self.reader.n_frame, self.Nthreshold), dtype=np.float32
        )
        self.PixelLiqHist = np.zeros(
            (self.reader.n_frame, self.Nsite_per_cut + 1), dtype=np.float32
        )
        self.PixelHetHist = np.zeros(
            (self.reader.n_frame, self.Nsite_per_cut * 2 + 1), dtype=np.float32
        )
        self.PixelLiqHist_With_Het = np.zeros(
            (self.reader.n_frame, self.Nsite_per_cut + 1), dtype=np.float32
        )
        self.PixelLiqHist_Without_Het = np.zeros(
            (self.reader.n_frame, self.Nsite_per_cut + 1), dtype=np.float32
        )

        self.PolyLiqCoincidence = np.zeros(self.reader.n_frame, dtype=np.float32)
        self.PolyLiqCorrelation = np.zeros(self.reader.n_frame, dtype=np.float32)

        self.LiqDensityAroundPoly = np.zeros(self.reader.n_frame, dtype=np.float32)
        self.LiqDensityGazPhase = np.zeros(self.reader.n_frame, dtype=np.float32)
        self.RatioLiqDensityInAndOut = np.zeros(self.reader.n_frame, dtype=np.float32)

        # self.PunctaNumber = np.zeros((self.reader.n_frame, self.Nthreshold), dtype=np.int32)
        # self.PunctaVolume = np.zeros((self.reader.n_frame, self.Nthreshold), dtype=np.int32)

        # self.PunctaLiqNumber = np.zeros(
        #     (self.reader.n_frame, self.Nthreshold), dtype=np.int32
        # )
        # self.PunctaLiqFraction = np.zeros(
        #     (self.reader.n_frame, self.Nthreshold), dtype=np.float32
        # )
        # self.PunctaLiqDensity = np.zeros(
        #     (self.reader.n_frame, self.Nthreshold), dtype=np.float32
        # )

        # self.PunctaHetNumber = np.zeros(
        #     (self.reader.n_frame, self.Nthreshold), dtype=np.int32
        # )
        # self.PunctaHetFraction = np.zeros(
        #     (self.reader.n_frame, self.Nthreshold), dtype=np.float32
        # )
        # self.PunctaHetDensity = np.zeros(
        #     (self.reader.n_frame, self.Nthreshold), dtype=np.float32
        # )

        for i in range(self.reader.n_frame):
            self.ProcessFrame(i)

            if (i + 1) % 10 == 0:
                print(
                    "Processed %d out of %d configurations"
                    % (i + 1, self.reader.n_frame)
                )

    def ProcessFrame(self, i):
        _CubeHet = np.zeros(self.Ncube, dtype=np.int32)
        _CubeLiq = np.zeros(self.Ncube, dtype=np.int32)

        data = next(self.reader)

        Liqindex, Liqcount = np.unique(
            np.sum(
                np.multiply(
                    np.divide(data.liq_pos, self.cut).astype(np.int32), self.x_y_z
                ),
                axis=1,
            ),
            return_counts=True,
        )
        _CubeLiq[Liqindex] = Liqcount

        Hetindex, Hetcount = np.unique(
            np.sum(
                np.multiply(
                    np.divide(data.poly_pos[data.poly_type == 1], self.cut).astype(
                        np.int32
                    ),
                    self.x_y_z,
                ),
                axis=1,
            ),
            return_counts=True,
        )
        _CubeHet[Hetindex] = Hetcount

        LiqNumber, LiqProba = np.unique(_CubeLiq, return_counts=True)
        self.PixelLiqHist[i][LiqNumber] = LiqProba

        HetNumber, HetProba = np.unique(_CubeHet, return_counts=True)
        self.PixelHetHist[i][HetNumber] = HetProba

        LiqCount_With_Het, LiqProba_With_Het = np.unique(
            _CubeLiq[_CubeHet > 1], return_counts=True
        )
        self.PixelLiqHist_With_Het[i][LiqCount_With_Het] = LiqProba_With_Het

        LiqCount_Without_Het, LiqProba_Without_Het = np.unique(
            _CubeLiq[_CubeHet <= 1], return_counts=True
        )
        self.PixelLiqHist_Without_Het[i][LiqCount_Without_Het] = LiqProba_Without_Het

        self.LiqDensityAroundPoly[i] = np.sum(_CubeLiq[_CubeHet > 1]) / np.sum(
            np.astype(_CubeHet > 1, np.int32)
        )

        self.LiqDensityGazPhase[i] = np.sum(_CubeLiq[_CubeHet <= 1]) / np.sum(
            np.astype(_CubeHet <= 1, np.int32)
        )

        self.RatioLiqDensityInAndOut[i] = (
            self.LiqDensityAroundPoly[i] / self.LiqDensityGazPhase[i]
        )

        self.PolyLiqCoincidence[i] = np.sum(np.minimum(_CubeLiq, _CubeHet)) / np.sum(
            np.maximum(_CubeLiq, _CubeHet)
        )

        self.PolyLiqCorrelation[i] = pearsonr(x=_CubeLiq, y=_CubeHet).statistic

        for ids, threshold in enumerate(self.list_threshold):
            self.PolyVolume[i, ids] = (
                np.count_nonzero(Hetcount > threshold) / self.Ncube
            )
            self.LiqVolume[i, ids] = np.count_nonzero(Liqcount > threshold) / self.Ncube

            # LiqPolyindex = np.intersect1d(
            #     Liqindex[Liqcount > threshold], Hetindex[Hetcount > threshold]
            # )

            # connect = self._connectPBC(self.Cube, LiqPolyindex)
            # graph = nx.from_numpy_array(connect)
            # clusters = nx.connected_components(graph)
            # clusters_index = [LiqPolyindex[list(cluster)] for cluster in clusters]

            # for index in clusters_index:
            #     _LiqSum = np.sum(Liqcount[np.isin(Liqindex, index)])
            #     _HetSum = np.sum(Hetcount[np.isin(Hetindex, index)])

            #     self.PunctaNumber[i, ids] += 1
            #     self.PunctaVolume[i, ids] += len(index)

            #     self.PunctaLiqNumber[i, ids] += _LiqSum
            #     self.PunctaHetNumber[i, ids] += _HetSum

            # self.PunctaLiqDensity[i, ids] = (
            #     self.PunctaLiqNumber[i, ids] / self.PunctaVolume[i, ids]
            #     if self.PunctaVolume[i, ids] > 0
            #     else 0
            # )
            # self.PunctaLiqFraction[i, ids] = (
            #     self.PunctaLiqNumber[i, ids] / self.reader.n_liq
            # )
            # self.PunctaHetDensity[i, ids] = (
            #     self.PunctaHetNumber[i, ids] / self.PunctaVolume[i, ids]
            #     if self.PunctaVolume[i, ids] > 0
            #     else 0
            # )
            # self.PunctaHetFraction[i, ids] = (
            #     self.PunctaHetNumber[i, ids] / self.reader.n_het
            # )

    # @staticmethod
    # @numba.njit
    # def _connectPBC(
    #     Cube: np.ndarray[tuple[int, ...], np.dtype[np.int32]],
    #     LiqPolyindex: np.ndarray[tuple[int, ...], np.dtype[np.int32]],
    # ) -> np.ndarray[tuple[int, int], np.dtype[np.int32]]:

    #     nIndex = LiqPolyindex.shape[0]
    #     connect = np.zeros((nIndex, nIndex), dtype=np.int32)

    #     for i in range(nIndex):
    #         for j in range(i, nIndex):
    #                 if LiqPolyindex[i] in Cube[LiqPolyindex[j]]:
    #                         connect[i, j] = 1
    #                         connect[j, i] = 1

    #     return connect

    def Print(self):
        print("\n")

        self.reader.close()
        self.hfile = h5py.File(self.processPath, "a")

        self.PrintDataset("PolyVolume", data=self.PolyVolume)
        self.PrintDataset("LiqVolume", data=self.LiqVolume)

        self.PrintDataset("PixelHetHist", data=self.PixelHetHist)
        self.PrintDataset("PixelLiqHist", data=self.PixelLiqHist)
        self.PrintDataset("PixelLiqHist_With_Het", data=self.PixelLiqHist_With_Het)
        self.PrintDataset(
            "PixelLiqHist_Without_Het", data=self.PixelLiqHist_Without_Het
        )

        self.PrintDataset("PolyLiqCoincidence", data=self.PolyLiqCoincidence)
        self.PrintDataset("PolyLiqCorrelation", data=self.PolyLiqCorrelation)

        self.PrintDataset("LiqDensityAroundPoly", data=self.LiqDensityAroundPoly)
        self.PrintDataset("LiqDensityGazPhase", data=self.LiqDensityGazPhase)
        self.PrintDataset("RatioLiqDensityInAndOut", data=self.RatioLiqDensityInAndOut)

        # self.PrintDataset("PunctaNumber", data=self.PunctaNumber)
        # self.PrintDataset("PunctaVolume", data=self.PunctaVolume)

        # self.PrintDataset("PunctaLiqNumber", data=self.PunctaLiqNumber)
        # self.PrintDataset("PunctaLiqFraction", data=self.PunctaLiqFraction)
        # self.PrintDataset("PunctaLiqDensity", data=self.PunctaLiqDensity)

        # self.PrintDataset("PunctaHetNumber", data=self.PunctaHetNumber)
        # self.PrintDataset("PunctaHetFraction", data=self.PunctaHetFraction)
        # self.PrintDataset("PunctaHetDensity", data=self.PunctaHetDensity)

        self.hfile.close()

    def PrintDataset(self, dataset_name, data):

        if dataset_name in self.hfile.keys():
            del self.hfile[dataset_name]

        self.hfile.create_dataset(dataset_name, data=data)

        print(f"Dataset {dataset_name} printed")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\033[1;31mUsage is %s inputDir nCut\033[0m" % sys.argv[0])
        sys.exit()

    inputDir = sys.argv[1]
    nCut = int(sys.argv[2])

    density = PolyLiqDensity(inputDir, nCut)

    density.Compute()
    density.Print()

    print("\n")
    print("PolyLiqDensity : Done\n\n")
