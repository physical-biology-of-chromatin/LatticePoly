##
##  LiqCluster.py
##  LatticePoly
##
##  Created by ppuel on 18/10/2024.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

import math as m
import os
import sys

import h5py
import networkx as nx
import numpy as np
from hdf5Reader import hdf5Reader
from numba import njit, typed, types

# from itertools import chain
# from itertools import chain
from scipy.spatial import KDTree


class LiqCluster:
        def __init__(self, inputDir, Maxime=False, cutoff=1 / 2**0.5 + 1e-3):
                print(f"LiqCluster : Init {inputDir}")
                self.reader = hdf5Reader(
                        inputDir, "traj.h5", -1, readLiq=True, readPoly=True, backInBox=True
                )
                self.processPath = os.path.join(inputDir, "process.h5")
                self.Maxime = Maxime
                self.boxDim = self.reader.boxDim

                self.cutoff = cutoff

                self.nLiq = self.reader.nLiq

                if self.Maxime:
                        self.threshold = 0.5

                        self.dropNum = np.zeros(self.reader.N, dtype=np.int32)
                        self.dropVolume = np.zeros(self.reader.N, dtype=np.float32)
                        self.dropVolHist = np.zeros(
                                (self.reader.N, self.nLiq - 1), dtype=np.float32
                        )
                        self.liqFraction = np.zeros(self.reader.N, dtype=np.float32)
                else:
                        self.dropInfo = np.zeros(
                                (self.reader.N, self.nLiq // 2),
                                dtype=[
                                        ("size", "i4"),
                                        ("center_of_mass", "f4", (3,)),
                                        ("mean_number_of_neighbors", "f4"),
                                        ("r_gyr", "f4"),
                                        ("aniso", "f4"),
                                ],
                        )

                        self.liqInfo = np.zeros((self.reader.N, self.nLiq), dtype=np.int32) - 1

        def Compute(self):
                print("\n")
                for i in range(self.reader.N):
                        self.ProcessFrame(i)

                        if (i + 1) % 10 == 0:
                                print("Process %d out of %d trajectories" % (i + 1, self.reader.N))

        def ProcessFrame(self, i):

                data = next(self.reader)

                if self.Maxime:
                        liqTree = KDTree(data.liqPos, boxsize=data.boxDim)
                        PolyTree = KDTree(data.polyPos[data.polyType == 1], boxsize=data.boxDim)

                        liqPolyNeighbor = liqTree.query_ball_tree(PolyTree, self.cutoff)
                        liqPolyDensity = np.fromiter(
                                map(lambda x: len(x) / 12, liqPolyNeighbor), dtype=np.float32
                        )

                        dropPos = data.liqPos[(liqPolyDensity + data.liqDens) > self.threshold]
                        connect = self._connectPBC(data.boxDim, dropPos, self.cutoff * 2)

                        graph = nx.from_numpy_array(connect)
                        clusters = nx.connected_components(graph)

                        sizes = [len(cluster) for cluster in clusters]

                        n = len(sizes)

                        self.liqFraction[i] = np.sum(np.asarray(sizes)) / self.nLiq
                        self.dropNum[i] = n
                        self.dropVolume[i] = np.mean(np.asarray(sizes)) if n > 0 else 0.0

                        for s in sizes:
                                self.dropVolHist[i, s - 1] += 1

                else:
                        connect = self._connectPBC(data.boxDim, data.liqPos, self.cutoff)
                   
                        graph = nx.from_numpy_array(connect)
                        
                        print('ici')
                        
                        clusters = nx.connected_components(graph)
                        
                        print('ou là')
                
                        clusters = sorted(clusters, key=len, reverse=True)
                        
                        print('là')
                        
                        cluster_ids = [
                                np.asarray(list(cluster), dtype=np.int32) for cluster in clusters
                        ]      
                                
                        # cluster_ids = [
                        #         np.asarray(list(cluster), dtype = np.int32) for cluster in clusters
                        # ]

                        print('pas là')

                        c_id = 0

                        while (
                                len(cluster_ids) > 0
                                and c_id < len(cluster_ids)
                                and len(cluster_ids[c_id]) > 1.5
                        ):  # 0.001 > self.pValue_dropSize(len(cluster_ids[c_id]), q, self.nLiq): #*cluster_number**self.exponant:
                                
                                
                                ids = cluster_ids[c_id]
                                size = len(ids)
                                
                                print(size)

                                clusterPBC = data.liqPos[ids]

                                # Calculation of the mean number of neighbors of particles in the cluster

                                mean_number_of_neighbors = np.mean(
                                        np.floor(12 * data.liqDens[ids] + 0.001)
                                )

                                # Calculation of the Center of Mass by fixing PBC

                                succDict = typed.Dict.empty(types.int32, types.int32[:])  # numba magic

                                for key, value in nx.dfs_successors(
                                        graph.subgraph(clusters[c_id]), cluster_ids[c_id][0]
                                ).items():
                                        succDict[np.int32(key)] = np.array(value, dtype=np.int32)
                                        # numba need to reconstruct the dictionary

                                self._fixClusterPBC(data.boxDim, clusterPBC, succDict)

                                center_of_mass = clusterPBC.mean(axis=0, keepdims=True)

                                # Calculation of the Radius of Gyration and the Anysotropy

                                clusterPBC -= center_of_mass

                                diag = np.linalg.svd(clusterPBC, compute_uv=False) * np.sqrt(12) / size
                                r2_gyr = np.square(diag).sum(axis=-1)

                                r_gyr = np.sqrt(r2_gyr)
                                aniso = 3 / 2.0 * (diag**4).sum(axis=-1) / r2_gyr**2 - 1 / 2.0

                                # Saving all the information of the cluster

                                self.dropInfo[i][c_id] = (
                                        size,
                                        center_of_mass,
                                        mean_number_of_neighbors,
                                        r_gyr,
                                        aniso,
                                )  # pass as a tuple for structured ndarray

                                # Saving all the cluster of each particles

                                for particle in ids:
                                        self.liqInfo[i][particle] = c_id

                                c_id += 1

        def Print(self):
                print("\n")
                self.reader.Close()

                if self.Maxime:
                        with h5py.File(self.processPath, "a") as processFile:
                                self.PrintDataset(processFile, "maximeDropVolHist_PC", self.dropVolHist)
                                self.PrintDataset(processFile, "maximeDropNum_PC", self.dropNum)
                                self.PrintDataset(processFile, "maximeDropVolume_PC", self.dropVolume)
                                self.PrintDataset(processFile, "maximeLiqFraction_PC", self.liqFraction)

                else:
                        with h5py.File(self.processPath, "a") as processFile:
                                self.PrintDataset(processFile, "liqDropInfo", data=self.dropInfo)
                                self.PrintDataset(processFile, "liqInfo", data=self.liqInfo)

        def PrintDataset(self, processFile, dataset_name, data):
                if dataset_name in processFile.keys():
                        del processFile[dataset_name]

                processFile.create_dataset(dataset_name, data=data)

                print(f"Dataset {dataset_name} printed")

        @staticmethod
        @njit
        def _fixClusterPBC(dims, pts, succDict):

                for i, succList in succDict.items():
                        for succ in succList:
                                delta = pts[succ] - pts[i]
                                pts[succ] = np.where(
                                        np.abs(delta) > dims / 2.0,
                                        pts[succ] - np.copysign(dims, delta),
                                        pts[succ],
                                )

        @staticmethod
        @njit
        def _connectPBC(dims, liqPos, cutoff):
                nPoints = liqPos.shape[0]
                connect = np.zeros((nPoints, nPoints), dtype=np.int32)

                for i in range(nPoints):
                        for j in range(i, nPoints):
                                pDist = 0.0

                                for k in range(3):
                                        delta = liqPos[j, k] - liqPos[i, k]

                                        while abs(delta) > dims[k] / 2.0:
                                                shift = m.copysign(dims[k], delta)

                                                delta -= shift

                                        pDist += delta**2

                                if pDist < cutoff**2:
                                        connect[i, j] = 1
                                        connect[j, i] = 1

                return connect


if __name__ == "__main__":
        if len(sys.argv) not in [2, 3]:
                print("\033[1;31mUsage is %s inputDir [Maxime]\033[0m" % sys.argv[0])
                sys.exit()

        inputDir = sys.argv[1]
        maxime = len(sys.argv) == 3

        cluster = LiqCluster(inputDir, maxime)

        cluster.Compute()
        cluster.Print()

        print("\n")
        print("LiqCluster : Done\n\n")
