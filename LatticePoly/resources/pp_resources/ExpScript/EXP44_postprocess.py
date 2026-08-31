import subprocess
import os
from LiqCluster_lifeTime import LifeTime, Droplet, Event
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import h5py
import time

c = 0


font = fm.FontProperties(weight='bold',
                                style='normal', size=20)
fontlabel = {"labelsize" : 20}


ninter_list = ["10000","100000"]
ldens_list = ["0.00068","0.00135","0.00270"]

# jll_list = ["1.5"]#,"1.5","2.0","4.0"]
# jll_valency_list = ["4"]#["2","4","6","8","10","12"]
# ldens_list = ["0.010"]#["0.005","0.010","0.020","0.040"]



N = 960

plasma = mpl.colormaps["plasma"].resampled(8)

exp = "EXP44_Nazli_first_simulation_3R_null_model"

os.makedirs(f"/home/ppuel/data/{exp}/data/", exist_ok=True)



c_max = 2*3*960
c=0

for ldens in ldens_list:
        for ninter in ninter_list:
                local_density_mean = []
                sizes_max = []
                tau = []
                for n in range(N):
                        c += 1
                        print(f"{c/c_max*100:2.2f} %", end = '\r')
                        try : 
                                file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/LDENS/{ldens}/NINTER/{ninter}/N/{n}/liq_droplets.pickle", "rb")
                                droplet_dict = pickle.load(file)
                                file.close()
                                for droplet in droplet_dict.values():
                                        local_density_mean.append(np.mean(droplet.local_density))
                                        sizes_max.append(np.max(droplet.sizes))
                                        tau.append(droplet.tau)
                        except :
                                pass

                
                                                
                np.save(f"/home/ppuel/data/{exp}/data/local_density_mean_ldens_{ldens}_ninter_{ninter}.npy", local_density_mean)
                np.save(f"/home/ppuel/data/{exp}/data/sizes_max_ldens_{ldens}_ninter_{ninter}.npy", sizes_max)
                np.save(f"/home/ppuel/data/{exp}/data/tau_ldens_{ldens}_ninter_{ninter}.npy", tau)
