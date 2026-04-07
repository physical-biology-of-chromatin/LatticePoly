##
##  LiqCluster.py
##  LatticePoly
##
##  Created by ppuel on 18/10/2024.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

import os
import sys

import h5py
import networkx as nx
import numpy as np

from hdf5Reader import hdf5Reader


class ProcessFrame:
    def __init__(self, liq_pos, liq_dens, box_dim, cutoff=1 / np.sqrt(2) + 1e-3):
        self.liq_pos = liq_pos
        self.liq_dens = liq_dens
        self.box_dim = box_dim
        self.cutoff = cutoff

        self.clusters: list[set[int]] = []

        self.generate_graph()

        self.drop_info = np.zeros(
            len(self.clusters),
            dtype=[
                ("size", "i4"),
                ("center_of_mass", "f4", (3,)),
                ("mean_degree", "f4"),
                ("r_gyr", "f4"),
                ("aniso", "f4"),
            ],
        )
        self.liq_info = np.zeros(len(self.liq_pos), dtype=np.int32) - 1

    def generate_graph(self):
        delta = self.liq_pos[None, :] - self.liq_pos[:, None]
        deltaPCB = np.sum(
            np.minimum(np.abs(delta), np.abs(self.box_dim[None, None, :] - delta)) ** 2,
            axis=2,
        )
        self.connect = np.where(deltaPCB < self.cutoff**2, 1, 0)
        self.graph = nx.from_numpy_array(self.connect)
        self.clusters = sorted(
            nx.connected_components(self.graph),
            key=len,
            reverse=True,
        )

    def process(self):
        for cindex, cluster in enumerate(self.clusters):
            size = len(cluster)
            print(size)
            if size <= 1:
                break

            clusterPos = self.liq_pos[list(cluster)]
            clusterDens = self.liq_dens[list(cluster)]

            mean_degree = np.mean(np.floor(12 * clusterDens + 0.001))

            centroids = (
                np.mean(clusterPos)
                - np.linspace(0, 1, num=np.size(clusterPos))[:, None]
                * self.box_dim[None, :]
            )

            delta = clusterPos[None, :, :] - centroids[:, None, :]
            delta = np.sum(
                np.minimum(np.abs(delta), np.abs(self.box_dim[None, None, :] - delta))
                ** 2,
                axis=1,
            ).T

            centroid = np.diag(centroids[np.argmin(delta, axis=1)])

            delta = clusterPos - centroid[None, :]
            absolutePos = np.where(
                np.abs(delta) < np.abs(self.box_dim - delta),
                clusterPos,
                clusterPos - self.box_dim,
            )
            centeredPos = absolutePos - centroid[None, :]

            diag = np.linalg.svd(centeredPos, compute_uv=False) * np.sqrt(12) / size
            r2_gyr = np.sum(np.square(diag), axis=-1)
            r_gyr = np.sqrt(r2_gyr)
            aniso = 3 / 2.0 * np.sum(diag**4, axis=-1) / r2_gyr**2 - 1 / 2.0
            self.drop_info[cindex] = (
                size,
                centroid,
                mean_degree,
                r_gyr,
                aniso,
            )
            self.liq_info[list(cluster)] = cindex
        return self.drop_info, self.liq_info


class LiqCluster:
    def __init__(self, inputDir, cutoff=1 / 2**0.5 + 1e-3):
        print(f"LiqCluster : Init {inputDir}")
        self.reader = hdf5Reader(
            inputDir, "traj.h5", -1, read_liq=True, read_poly=True, back_in_box=True
        )
        self.processPath = os.path.join(inputDir, "process.h5")
        self.box_dim = self.reader.box_dim

        self.drop_info = np.zeros(
            (self.reader.n_frame, self.reader.n_liq // 2),
            dtype=[
                ("size", "i4"),
                ("centroid", "f4", (3,)),
                ("mean_degree", "f4"),
                ("r_gyr", "f4"),
                ("aniso", "f4"),
            ],
        )

        self.liq_info = np.zeros((self.reader.n_frame, self.reader.n_liq), dtype=np.int32) - 1
        next(self.reader)  # first frame is random noise

    def compute(self):
        print("+------------------------- compute -------------------------+")
        for findex, fdata in enumerate(self.reader):
            self.process_frame(findex, fdata)
            if findex % 10 == 0:
                print(f"Process {findex} out of {len(self.reader)} trajectories")

    def process_frame(self, findex, fdata):
        pframe = ProcessFrame(fdata.liqPos, fdata.liqDens, fdata.boxDim)
        drop_info, liq_info = pframe.process()
        self.drop_info[findex, : len(drop_info)] = drop_info
        self.liq_info[findex] = liq_info

    def print(self):
        print("\n")
        self.reader.close()

        with h5py.File(self.processPath, "a") as processFile:
            self.print_dataset(processFile, "liqDropInfo", data=self.drop_info)
            self.print_dataset(processFile, "liqInfo", data=self.liq_info)

    def print_dataset(self, processFile, dataset_name, data):
        if dataset_name in processFile.keys():
            del processFile[dataset_name]
        processFile.create_dataset(dataset_name, data=data)
        print(f"Dataset {dataset_name} printed")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("\033[1;31mUsage is %s inputDir\033[0m" % sys.argv[0])
        sys.exit()

    inputDir = sys.argv[1]

    cluster = LiqCluster(inputDir)

    cluster.compute()
    cluster.print()

    print("\n")
    print("LiqCluster : Done\n\n")
