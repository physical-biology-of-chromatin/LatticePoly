import subprocess
import os
from LiqCluster_lifeTime import LifeTime, Droplet, Event
import pickle
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import h5py
import scipy.stats as st
import itertools
import scipy.optimize as op
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression

font = fm.FontProperties(weight='bold',
                                style='normal', size=20)
fontlabel = {"labelsize" : 20}



exp = "EXP45_Lucy_first_attempt"

jll_list = ["0.5", "1.0", "2.0"]
jlp_list = ["0.2", "0.5", "1.0"]
jll_valency_list = ["2", "4", "6", "8"]
jlp_valency_list = ["2", "4", "6", "8"]
ldens_list = ["0.019", "0.038", "0.076"]


os.makedirs(f"/home/ppuel/data/{exp}/hist/", exist_ok=True)

hist_gyration = [[] for i in range(3)]
N = 10

# cmax = N*3*3*3*4*4

# c = 1
# for jll in jll_list:
#         for jlp in jlp_list:
#                 for jll_valency in jll_valency_list:
#                         for jlp_valency in jlp_valency_list:
#                                 for ldens in ldens_list:
#                                         i = ldens_list.index(ldens)
#                                         for n in range(N):
#                                                 if c%20 == 0:
#                                                         print(f"{c/cmax:0.2%}", end = "\r")
#                                                 if "process.h5" in os.listdir(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLL_VALENCY/{jll_valency}/JLP_VALENCY/{jlp_valency}/LDENS/{ldens}/N/{n}/"):
#                                                         file = h5py.File(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLL_VALENCY/{jll_valency}/JLP_VALENCY/{jlp_valency}/LDENS/{ldens}/N/{n}/process.h5")
#                                                         ratio = np.mean(file["polyGyration"][-10:])/file["polyGyration"][0]
#                                                         if ratio < 0.8:
#                                                                 print(f"45 {jll} {jlp} {jll_valency} {jlp_valency} {ldens} {n}")
#                                                         hist_gyration[i].append(ratio)
#                                                 else:
#                                                         pass
#                                                 c+=1


gyration_init = 0

cmax = N*3*3*3*4*4

c = 1
# for jll in jll_list:
#         i = jll_list.index(jll)
#         for jlp in jlp_list:
#                 j = jlp_list.index(jlp)
#                 for jll_valency in jll_valency_list:
#                         k = jll_valency_list.index(jll_valency)
#                         for jlp_valency in jlp_valency_list:
#                                 l = jlp_valency_list.index(jlp_valency)
#                                 for ldens in ldens_list:

#                                         m = ldens_list.index(ldens)
#                                         for n in range(N):
#                                                 if c%20 == 0:
#                                                         print(f"{c/cmax:0.2%}", end = "\r")
#                                                 if "process.h5" in os.listdir(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLL_VALENCY/{jll_valency}/JLP_VALENCY/{jlp_valency}/LDENS/{ldens}/N/{n}/"):
#                                                         file = h5py.File(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLL_VALENCY/{jll_valency}/JLP_VALENCY/{jlp_valency}/LDENS/{ldens}/N/{n}/process.h5")
#                                                         gyration_init += file["polyGyration"][0]
#                                                 else:
#                                                         pass
#                                                 c+=1
                                        
# gyration_init /= cmax
# file = open(f"/home/ppuel/data/{exp}/hist/gyr_init.float",'w')
# file.write(str(gyration_init))
# file.close()

file = open(f"/home/ppuel/data/{exp}/hist/gyr_init.float",'r')
gyration_init = float(file.readline().strip())



jll_list = ["0.5", "1.0", "2.0"]
jlp_list = ["0.2", "0.5", "1.0"]
jll_valency_list = ["4"]
jlp_valency_list = ["4"]
ldens_list = ["0.019"]#, "0.076"]


map_gyration = np.zeros((3,3))


for jll in jll_list:
        i = jll_list.index(jll)
        for jlp in jlp_list:
                j = jlp_list.index(jlp)
                for jll_valency in jll_valency_list:
                        k = jll_valency_list.index(jll_valency)
                        for jlp_valency in jlp_valency_list:
                                l = jlp_valency_list.index(jlp_valency)
                                for ldens in ldens_list:

                                        tmp_n = 0
                                        m = ldens_list.index(ldens)
                                        for n in range(N):
                                                if c%20 == 0:
                                                        print(f"{c/cmax:0.2%}", end = "\r")
                                                if "process.h5" in os.listdir(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLL_VALENCY/{jll_valency}/JLP_VALENCY/{jlp_valency}/LDENS/{ldens}/N/{n}/"):
                                                        tmp_n += 1
                                                        file = h5py.File(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLL_VALENCY/{jll_valency}/JLP_VALENCY/{jlp_valency}/LDENS/{ldens}/N/{n}/process.h5")
                                                        map_gyration[i,j] += np.mean(file["polyGyration"][-10:])
                                                else:
                                                        pass
                                                c+=1
                                        if tmp_n != 0:
                                                map_gyration[i,j] /= tmp_n*gyration_init
                                        if map_gyration[i,j] < 0.1:
                                                map_gyration[i,j] = np.nan
                                        # if map_gyration[i,j] < 0.95:
                                        #         print(f'45 {jll} {jlp} {jll_valency} {jlp_valency} {ldens}')



gyr_min = np.nanmin(map_gyration)
gyr_max = np.nanmax(map_gyration)
print(gyr_max, gyr_min)
fig, ax = plt.subplots(1,1)
fig.set_figheight(16)
fig.set_figwidth(24)
# for i in range(3):
im = ax.imshow(np.transpose(map_gyration), vmax = gyr_max, vmin = gyr_min, cmap = "plasma")
cbar = fig.colorbar(im, ax=ax, shrink=0.7)
cbar.set_ticks(ticks=cbar.get_ticks(), labels = [f"{tik:0.2f}" for tik in cbar.get_ticks()], font = font)
        

fig.savefig(f"/home/ppuel/data/{exp}/hist/PolyGyr_map_4_4")