##
##  Liq_Density.py
##  LatticePoly
##
##  Created by ppuel on 18/10/2024.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys

import h5py
import numpy as np
from Reader import Reader

class Liq_Density:
    def __init__(self, input_dir: str):
        self.reader = Reader(
            input_dir, 
            read_liq=True,
            read_poly=False
        )
        self.process_path = os.path.join(input_dir, "process.h5")
        self.liq_dens_mean = np.zeros(self.reader.n_frame, dtype=np.float32)
        self.liq_dens_std = np.zeros(self.reader.n_frame, dtype=np.float32)
        self.liq_dens_hist = np.zeros((self.reader.n_frame, 13), dtype=np.int32)


    def compute(self):
        with self.reader as iterator:
            next(iterator)  # first frame is random noise
            for findex, fdata in enumerate(iterator):
                self.process_frame(findex, fdata)
                if findex % 10 == 0:
                    print(f"Process {findex} out of {len(iterator)} trajectories")

    def process_frame(self, findex, fdata):
        self.liq_dens_mean[findex] = fdata.liq_dens.mean()
        self.liq_dens_std[findex] = np.square(fdata.liq_dens - self.liq_dens_mean[findex]).sum()
        for j in np.asarray((fdata.liq_dens + 0.001) * 12, dtype=np.int32):
            self.liq_dens_hist[findex][j] += 1

    def print(self):
        with h5py.File(self.process_path, "a") as process_file:
            self.print_dataset(
                process_file,
                "liq_dens_mean",
                data=self.liq_dens_mean
            )
            self.print_dataset(
                process_file,
                "liq_dens_std",
                data=np.sqrt(self.liq_dens_std / self.reader.n_liq)
            )
            self.print_dataset(
                process_file,
                "liq_dens_hist",
                data=self.liq_dens_hist
            )

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

    density = Liq_Density(input_dir)

    density.compute()
    density.print()

    print("\n")
    print("Liq_Density : Done\n\n")
