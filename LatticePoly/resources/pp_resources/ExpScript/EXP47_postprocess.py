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


ldens_list = ["0.00270","0.00405","0.00540"]
N = 400

plasma = mpl.colormaps["plasma"].resampled(8)

exp = "EXP47_Nazly_second_simulation_3R_ref"

os.makedirs(f"/home/ppuel/data/{exp}/data/", exist_ok=True)



# c_max = 3*400
# c=0

# for ldens in ldens_list:
#         local_density_mean = []
#         sizes_max = []
#         tau = []
#         for n in range(N):
#                 c += 1
#                 print(f"{c/c_max*100:2.2f} %", end = '\r')
#                 try : 
#                         file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/LDENS/{ldens}/N/{n}/liq_droplets.pickle", "rb")
#                         droplet_dict = pickle.load(file)
#                         file.close()
#                         for droplet in droplet_dict.values():
#                                 local_density_mean.append(np.mean(droplet.local_density))
#                                 sizes_max.append(np.max(droplet.sizes))
#                                 tau.append(droplet.tau)
#                 except :
#                         pass

        
                                        
#         np.save(f"/home/ppuel/data/{exp}/data/local_density_mean_ldens_{ldens}.npy", local_density_mean)
#         np.save(f"/home/ppuel/data/{exp}/data/sizes_max_ldens_{ldens}.npy", sizes_max)
#         np.save(f"/home/ppuel/data/{exp}/data/tau_ldens_{ldens}.npy", tau)


for ldens in ldens_list:

        # local_density_mean = np.load(f"/home/ppuel/data/{exp}/data/local_density_mean_ldens_{ldens}.npy")
        sizes_max = np.load(f"/home/ppuel/data/{exp}/data/sizes_max_ldens_{ldens}.npy")
        tau = np.load(f"/home/ppuel/data/{exp}/data/tau_ldens_{ldens}.npy")
        
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot()

        
        A, B = np.unique(np.array([sizes_max, tau]), return_counts=True, axis=1)

        A_min = np.min(A, axis=1)

        A_max = np.max(A, axis=1)


        r_max = 0

        for s, t in A.transpose():
                r_max = max(r_max, ((s-2)**2+(t-1)**2)**(1/2))

        ax.scatter(A[0], A[1], s=20*(np.log10(B)+1))
        
        u = np.linspace(0,np.pi/2,100)
        x = r_max*np.cos(u)+2
        y = r_max*np.sin(u)+1
        ax.plot(x, y, color="r")
        
        # ax.plot(1 + r_max*np.cos(np.linspace(0,np.pi/2,100)), 2 + r_max*np.sin(np.linspace(0,np.pi/2,100)),1)
        fontsize = 30
        ax.text( A_min[0]+.25,  A_max[1]+.25, f"r_max : {r_max:0.1f}", fontdict = {"fontsize" : fontsize})
        
        ax.set_xlim(2-.25, max(x))
        ax.set_ylim(1-.25, max(y))

        ax.set_xlabel("sizes max", fontdict = {"fontsize" : fontsize})
        ax.set_ylabel("life time", fontdict = {"fontsize" : fontsize})
        
        fig.suptitle(f"rho  = {ldens}")
        fig.savefig(f"/home/ppuel/data/{exp}/fig_r_max_3D_ldens_{ldens}.png")

        np.save(f"/home/ppuel/data/{exp}/r_max_3D_ldens_{ldens}.npy", r_max)
