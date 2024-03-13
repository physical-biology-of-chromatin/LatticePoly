import os
import sys
import pandas as pd
import numpy as np
import cooler


directory=["nagano_10kb_early-S_1.cool","nagano_10kb_early-S_2.cool","nagano_10kb_early-S_3.cool"]
for weight_columm in ["raw"]:
    for i,dire in enumerate(directory):
        clr=cooler.Cooler("///Users/mariachiara/Desktop/singlecell_mesc/"+dire)
        with clr.open("r+") as f:
            if(np.sum(clr.bins()[:].columns=="raw")==0):
                f["bins"].create_dataset("raw", data=np.ones(len(clr.bins()[:])), compression="gzip", compression_opts=6)
           