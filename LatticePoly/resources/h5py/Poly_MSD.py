##
##  Poly_MSD.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys
import psutil

import h5py
import numpy as np
from utils import msdFFT
from Reader import Reader

class Poly_MSD():
    def __init__(self, input_dir : str):
        self.reader = Reader(input_dir)
        self.process_path = os.path.join(input_dir, "process.h5")
        
    def compute(self):
        v_mem = psutil.virtual_memory()
        size_tot = self.reader.n_frame * self.reader.poly_pos.nbytes
        
        if size_tot < v_mem.available:
            self.cumul_dist_het = 0
            self.cumul_dist_hom = 0
            self.cumul_dist_PRE = 0
        
            pos_hist = self.read_hist()
            
            for tindex in range(self.reader.n_tad):
                if self.reader.poly_painter[tindex] == 1:
                    self.cumul_dist_PRE += msdFFT(pos_hist[:, tindex])
                elif self.reader.poly_type[tindex] == 1:
                    self.cumul_dist_het += msdFFT(pos_hist[:, tindex])
                else:
                    self.cumul_dist_hom += msdFFT(pos_hist[:, tindex])
                            
                if (tindex+1) % 1000 == 0:
                    print(f"Processed {tindex+1} out of {self.reader.n_tad} TADs")
                    
        else:
            print("Memory overflow likely - reduce chosen number of frames")
            sys.exit()


    def compute_tad(self, tindex):
        tad_pos_hist = np.zeros((self.reader.n_frame, 3), dtype=np.float32)

        with self.reader as iterator:
            next(iterator)  # first frame is random noise
            for findex, fdata in enumerate(iterator):
                tad_pos_hist[findex] = fdata.poly_pos[tindex]
                if findex % 10 == 0:
                    print(f"Process {findex} out of {len(iterator)} trajectories")
            
        self.dist_tad = msdFFT(tad_pos_hist)

    def read_hist(self):
        pos_hist = np.zeros((self.reader.n_frame, self.reader.n_tad, 3), dtype=np.float32)
        
        with self.reader as iterator:
            next(iterator)  # first frame is random noise
            for findex, fdata in enumerate(iterator):
                pos_hist[findex] = fdata.poly_pos
                if findex % 10 == 0:
                    print(f"Process {findex} out of {len(iterator)} trajectories")
            
        return pos_hist
    
    def print(self):
        with h5py.File(self.process_path, 'a') as process_file:
            if np.count_nonzero(self.reader.poly_painter == 1) > 0:
                poly_msd_PRE = self.cumul_dist_PRE / np.count_nonzero(self.reader.poly_painter == 1)
                self.print_dataset(process_file, "poly_msd_PRE", data = poly_msd_PRE)
            
            if self.reader.n_het > 0:
                poly_msd_het = self.cumul_dist_het /  self.reader.n_het
                self.print_dataset(process_file, "poly_msd_het", data = poly_msd_het)

            if self.reader.n_euc > 0:
                poly_msd_hom = self.cumul_dist_hom / self.reader.n_euc
                self.print_dataset(process_file, "poly_msd_hom", data = poly_msd_hom)    

    def print_tad(self, tindex):
        with h5py.File(self.process_path, 'a') as process_file:
            self.print_dataset(process_file, f"msd_tad_{tindex}", data = self.dist_tad)
    
    def print_dataset(self, process_file, dataset_name, data):
        if dataset_name in process_file.keys():
            del process_file[dataset_name]
        process_file.create_dataset(dataset_name, data=data)
        print(f"Dataset {dataset_name} printed")

if __name__ == "__main__":
    if len(sys.argv) not in [2, 3]:
        print(f"Usage is {sys.argv[0]} input_dir [tindex]")
        sys.exit()

    input_dir = sys.argv[1]
    msd = Poly_MSD(input_dir)

    if len(sys.argv) == 2:
        msd.compute()
        msd.print()        
    elif len(sys.argv) == 3:
        tindex = int(sys.argv[2])
        msd.compute_tad(tindex)
        msd.print_tad(tindex)