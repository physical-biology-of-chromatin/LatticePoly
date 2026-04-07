import os
import h5py
import numpy as np
from typing import cast

from hdf5Reader import hdf5Reader

class Process():
        def __init__(self, inputDir : str):
                self.reader = hdf5Reader(inputDir, "traj.h5", -1, readLiq=True, readPoly=True, backInBox=True)
                self.processPath = os.path.join(inputDir, "process.h5")
        
        def PrintDataset(self, processFile : h5py.File, dataset_name : str, data : np.ndarray[tuple[int, ...], np.dtype[np.float32 | np.int32]]):
                
                if dataset_name in processFile.keys():
                        tmp = cast(h5py.Dataset, processFile[dataset_name])
                        tmp[:] = data
                else:
                        processFile.create_dataset(dataset_name, data = data)
                        
                print(f"Dataset {dataset_name} printed")
        