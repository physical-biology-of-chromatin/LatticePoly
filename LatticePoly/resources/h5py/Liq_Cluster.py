##
##  Liq_Cluster.py
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

from Reader import Reader


class ProcessFrame:
    def __init__(self, liq_pos, liq_dens, box_dim, cutoff=1 / np.sqrt(2) + 1e-3):
        self.liq_pos = liq_pos
        self.liq_dens = liq_dens
        self.box_dim = box_dim
        self.cutoff = cutoff

        self.clusters: list[set[int]] = []

        self.generate_graph()

        self.drop_info = np.zeros(
            len(self.liq_pos) // 2,
            dtype=[
                ("size", "i4"),
                ("center_of_mass", "f4", (3,)),
                ("mean_number_of_neighbors", "f4"),
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
            if size <= 1:
                break

            clusterPos = self.liq_pos[list(cluster)]
            clusterDens = self.liq_dens[list(cluster)]

            mean_number_of_neighbors = np.mean(np.floor(12 * clusterDens + 0.001))

            array_center_of_mass = (
                np.mean(clusterPos, axis = 0)
                - np.linspace(0, 1, num=np.size(clusterPos))[:, None]
                * self.box_dim[None, :]
            )

            delta = clusterPos[None, :, :] - array_center_of_mass[:, None, :]
            delta = np.sum(
                np.minimum(np.abs(delta), np.abs(self.box_dim[None, None, :] - delta))
                ** 2,
                axis=1,
            ).T

            center_of_mass = np.diag(array_center_of_mass[np.argmin(delta, axis=1)])

            delta = clusterPos - center_of_mass[None, :]
            absolutePos = np.where(
                np.abs(delta) < np.abs(self.box_dim - delta),
                clusterPos,
                clusterPos - self.box_dim,
            )
            centeredPos = absolutePos - center_of_mass[None, :]

            diag = np.linalg.svd(centeredPos, compute_uv=False) * np.sqrt(12) / size
            r2_gyr = np.sum(np.square(diag), axis=-1)
            r_gyr = np.sqrt(r2_gyr)
            aniso = 3 / 2.0 * np.sum(diag**4, axis=-1) / r2_gyr**2 - 1 / 2.0

            center_of_mass_PBCs = np.where(center_of_mass < 0, center_of_mass + self.box_dim, center_of_mass) 
            # center_of_mass_PBCs = np.where(center_of_mass > self.box_dim, center_of_mass_PBCs - self.box_dim, center_of_mass_PBCs) 
            
            self.drop_info[cindex] = (
                size,
                center_of_mass_PBCs,
                mean_number_of_neighbors,
                r_gyr,
                aniso,
            )
            self.liq_info[list(cluster)] = cindex
        return self.drop_info, self.liq_info


class Liq_Cluster:
    def __init__(self, input_dir, cutoff=1 / 2**0.5 + 1e-3):
        self.input_dir = input_dir
        self.process_path = os.path.join(self.input_dir, "process.h5")

        self.reader = Reader(self.input_dir, read_liq=True, read_poly=False, back_in_box=True)
        
        self.box_dim = self.reader.box_dim

        self.drop_info = np.zeros(
            (self.reader.n_frame, self.reader.n_liq // 2),
            dtype=[
                ("size", "i4"),
                ("center_of_mass", "f4", (3,)),
                ("mean_number_of_neighbors", "f4"),
                ("r_gyr", "f4"),
                ("aniso", "f4"),
            ],
        )

        self.liq_info = np.zeros((self.reader.n_frame, self.reader.n_liq), dtype=np.int32) - 1
            
    def compute(self):
        with self.reader as iterator:
            next(iterator)  # first frame is random noise
            for findex, fdata in enumerate(iterator):
                self.process_frame(findex, fdata)
                if findex % 10 == 0:
                    print(f"Process {findex} out of {len(iterator)} trajectories")

    def process_frame(self, findex, fdata):
        pframe = ProcessFrame(fdata.liq_pos, fdata.liq_dens, fdata.box_dim)
        drop_info, liq_info = pframe.process()
        self.drop_info[findex, : len(drop_info)] = drop_info
        self.liq_info[findex] = liq_info

    def print(self):
        with h5py.File(self.process_path, "a") as process_file:
            self.print_dataset(process_file, "liq_drop_info", data=self.drop_info)
            self.print_dataset(process_file, "liq_info", data=self.liq_info)

    def print_dataset(self, process_file, dataset_name, data):
        if dataset_name in process_file.keys():
            del process_file[dataset_name]
        process_file.create_dataset(dataset_name, data=data)
        print(f"Dataset {dataset_name} printed")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage is {sys.argv[0]} input_dir")
        sys.exit()

    input_dir = sys.argv[1]

    cluster = Liq_Cluster(input_dir)

    cluster.compute()
    cluster.print()
