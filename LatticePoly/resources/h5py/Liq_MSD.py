##
##  liq_MSD.py
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


class Liq_MSD():
    """
    Liq_MSD is a standard class for computing mean squared displacement of
    the liquid particles. It depends on Reader to open the h5 trajectory file
    and on msdFFT for the MSD computation. It computes the MSD of both the
    center of mass of the particles and the mean MSD of all the , the standard deviation and the histogram of 
    liquid local density and saves them in a process.h5 file on input_dir
    
    Attributes
    ----------
    input_dir : str
        Path to the directory where the trajectory file is located.
    
    
    Methods
    -------
    compute()
        Standard computing function that handle the open and iteration on
        the traj.h5 file.
    
    process_frame(findex: int, fdata: reader)
        Standard processing function that calculate statical data for the
        frame findex.
    
    print() / print_dataset()
        Standard output function write statical data in the process.h5 file.
        print_dataset uses del function to remove a dataset with the same name
        in the process.h5 file. This allows to recompute statical data in the
        same process.h5 file.
    """ 

    def __init__(self, input_dir : str):
        self.reader = Reader(input_dir, read_liq=True, read_poly=False)
        self.process_path = os.path.join(input_dir, "process.h5")
        


    def compute(self):
        print("\n")
        
        v_mem = psutil.virtual_memory()
        size_tot = self.reader.n_frame * self.reader.liq_pos.nbytes
        
        if size_tot < v_mem.available:
            self.cumul_dist = 0
            
            self.liq_pos_init = self.reader.liq_pos

            pos_hist = self.read_hist()
            
            print('\n')

            for lindex in range(self.reader.n_liq):
                self.cumul_dist += msdFFT(pos_hist[:, lindex])

                if (lindex+1) % 1000 == 0:
                    print("Processed %d out of %d protein" % (lindex+1, self.reader.n_liq))
            
            
            self.liq_MSD_CoM = msdFFT(np.mean(pos_hist,1))
            
            self.liq_MSD = self.cumul_dist / self.reader.n_liq
        
        else:
            print("Memory overflow likely - reduce chosen number of frames")
            sys.exit()


    def read_hist(self):
        pos_hist = np.zeros((self.reader.n_frame, self.reader.n_liq, 3), dtype=np.float32)

        with self.reader as iterator:
            next(iterator)  # first frame is liq_pos_init
            for findex, fdata in enumerate(iterator):
                pos_hist[findex] = self.liq_pos_init + fdata.liq_disp
                if findex % 10 == 0:
                    print(f"Process {findex} out of {len(iterator)} trajectories")
    
        return pos_hist


    def print(self):
        with h5py.File(self.process_path, 'a') as process_file:
            
            self.print_dataset(process_file, "liq_MSD", data = self.liq_MSD)
            self.print_dataset(process_file, "liq_MSD_CoM", data = self.liq_MSD_CoM)
        
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

    msd = Liq_MSD(input_dir)

    msd.compute()
    msd.print()