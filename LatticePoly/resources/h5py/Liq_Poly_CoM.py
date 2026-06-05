##
##  PolyMSD.py
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

class liq_poly_CoM():
    def __init__(self, input_dir : str):
        self.reader = Reader(input_dir, read_liq=True)
        self.process_path = os.path.join(input_dir, "process.h5")
    
    def compute(self):
        v_mem = psutil.virtual_memory()
        size_tot = self.reader.n_frame * (self.reader.poly_pos.nbytes + self.reader.liq_pos.nbytes)
        self.liq_pos_init = self.reader.liq_pos
 
        if size_tot < v_mem.available:
            pos_hist = self.read_hist()
            self.liq_poly_CoM = msdFFT(np.mean(pos_hist, 1))        
        else:
            print("Memory overflow likely - reduce chosen number of frames")
            sys.exit()


    def read_hist(self):
        pos_hist = np.zeros((self.reader.n_frame, self.reader.n_tad + self.reader.n_liq, 3), dtype=np.float32)
        
        with self.reader as iterator:
            next(iterator)  # first frame is liq_pos_init
            for findex, fdata in enumerate(iterator):
                pos_hist[findex] = np.concatenate((fdata.poly_pos, self.liq_pos_init + fdata.liq_disp),axis=0)
                if findex % 10 == 0:
                    print(f"Process {findex} out of {len(iterator)} trajectories")
        
        return pos_hist
        
    def print(self):
        with h5py.File(self.process_path, 'a') as process_file:
            self.print_dataset(process_file, "liq_poly_CoM", data = self.liq_poly_CoM)
          
    
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

    CoM = liq_poly_CoM(input_dir)

    CoM.compute()
    CoM.print()