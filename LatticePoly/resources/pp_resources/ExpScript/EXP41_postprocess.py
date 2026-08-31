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
import seaborn as sns

font = fm.FontProperties(weight='bold',
                                style='normal', size=20)
fontlabel = {"labelsize" : 20}


jll_list = ["1.0","1.5","2.0","4.0"]
jll_valency_list = ["2","4","6","8","10","12"]
ldens_list = ["0.005","0.010","0.020","0.040"]

# jll_list = ["1.5"]#,"1.5","2.0","4.0"]
# jll_valency_list = ["4"]#["2","4","6","8","10","12"]
# ldens_list = ["0.010"]#["0.005","0.010","0.020","0.040"]

N = 105

plasma = mpl.colormaps["plasma"].resampled(8)

exp = "EXP41_liqFraction_phase_diagram_droplet"

os.makedirs(f"/home/ppuel/data/{exp}/classifier/", exist_ok=True)

local_density_mean = []
sizes_max = []
tau = []

# c_max = 4*4*6*105
# c=0

# for ldens in ldens_list:
#         for jll in jll_list:
#                 for jll_valency in jll_valency_list:
#                         for n in range(N):
#                                 c += 1
#                                 print(f"{c/c_max*100:2.2f} %",end='\r')
#                                 file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLL_VALENCY/{jll_valency}/LDENS/{ldens}/N/{n}/liq_droplets.pickle", "rb")
#                                 droplet_dict = pickle.load(file)
#                                 file.close()
#                                 for droplet in droplet_dict.values():
#                                         local_density_mean.append(np.mean(droplet.local_density))
#                                         sizes_max.append(np.max(droplet.sizes))
#                                         tau.append(droplet.tau)

local_density_mean = np.load(f"/home/ppuel/data/{exp}/classifier/local_density_mean.npy")
sizes_max = np.load(f"/home/ppuel/data/{exp}/classifier/sizes_max.npy")
tau = np.load(f"/home/ppuel/data/{exp}/classifier/tau.npy")


tmp_fig, ax = plt.subplots()
fig, axs = plt.subplots(3,1)

fig.set_figheight(16)
fig.set_figwidth(24)
ldm_hist = ax.hist(local_density_mean)
smx_hist = ax.hist(sizes_max)
tau_hist = ax.hist(tau)

axs[0].plot((ldm_hist[1][:-1]+ldm_hist[1][1:])/2, np.log10(ldm_hist[0]))
axs[1].plot((smx_hist[1][:-1]+smx_hist[1][1:])/2, np.log10(smx_hist[0]))
axs[2].plot((tau_hist[1][:-1]+tau_hist[1][1:])/2, np.log10(tau_hist[0]))


fig.savefig(f"/home/ppuel/data/{exp}/classifier/histogram")
plt.close('all')