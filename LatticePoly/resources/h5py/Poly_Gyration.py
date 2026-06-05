##
##  Poly_Gyration.py
##  LatticePoly
##
##  Created by mtortora on 02/08/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

import os
import sys

import h5py
import numpy as np
from Reader import Reader

class Poly_Gyration():
    def __init__(self, input_dir : str):
        self.reader = Reader(input_dir)
        self.process_path = os.path.join(input_dir, "process.h5")
        self.poly_aniso = np.zeros(self.reader.n_frame, dtype=np.float32)
        self.poly_gyration = np.zeros(self.reader.n_frame, dtype=np.float32)
        
    def compute(self):
        with self.reader as iterator:
            next(iterator)  # first frame is random noise
            for findex, fdata in enumerate(iterator):
                self.process_frame(findex, fdata)
                if findex % 10 == 0:
                    print(f"Process {findex} out of {len(iterator)} trajectories")

    def process_frame(self, findex, fdata):
        norm = 0
        
        for d in fdata.domains:
            if d.size > 2:
                pos = fdata.poly_pos[d]
                pos -= pos.mean(axis=0, keepdims=True)
                diag = np.linalg.svd(pos, compute_uv=False) * np.sqrt(12)/d.size
                r2_gyr = np.square(diag).sum(axis=-1)
                r_gyr = np.sqrt(r2_gyr)
                aniso = 3/2.*(diag**4).sum(axis=-1)/r2_gyr**2 - 1/2.
                norm += d.size
                self.poly_aniso[findex] += aniso * d.size
                self.poly_gyration[findex] += r_gyr * d.size
                                    
        self.poly_aniso[findex] /= norm if norm > 0 else 1
        self.poly_gyration[findex] /= norm if norm > 0 else 1

    def print(self):
        with h5py.File(self.process_path, 'a') as process_file:
            self.print_dataset(process_file, "poly_aniso", data = self.poly_aniso)
            self.print_dataset(process_file, "poly_gyration", data = self.poly_gyration)
        
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

    gyr = Poly_Gyration(input_dir)

    gyr.compute()
    gyr.print()