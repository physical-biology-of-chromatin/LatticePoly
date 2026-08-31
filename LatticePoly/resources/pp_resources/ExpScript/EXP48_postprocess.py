# import subprocess
import os
from LiqDroplet import LifeTime, Droplet, Event
import pickle
# import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import h5py
# import scipy.stats as st
# import itertools
# import scipy.optimize as op
# from sklearn.decomposition import PCA
# from sklearn.linear_model import LogisticRegression
# import seaborn as sns
import sys

plt.style.use('./resources/h5py/presentation.mplstyle')

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

font = fm.FontProperties(weight='bold',
                                style='normal', size=20)
fontlabel = {"labelsize" : 20}

jll_list = ["0.5","1.0","2.0"]
jlp_list = ["0.2","0.5","1.0"]
jll_valency_list = ["2", "6", "10"]
jlp_valency_list = ["2", "6", "10"]
jpl_valency_list = ["2", "6", "10"]
ldens_list = ["0.019","0.038"]


def conversion_kBT_to_kJ_per_mol(J):
        J = np.array([float(i) for i in J])
        T = 300
        kB = 1.380e-23
        Na = 6.022e23
        res = J * kB * Na * T * 1e-3
        return([f"{i:0.1f}" for i in res])

jll_list_conv = conversion_kBT_to_kJ_per_mol(jll_list)
jlp_list_conv = conversion_kBT_to_kJ_per_mol(jlp_list)
# jll_list = ["1.0"]
# jlp_list = ["0.5"]
# EV_list = ["0"]
# ninter_list = ["100000"]
# ldens_list = ["0.00135"]


N = 16

range_data = 10

# plasma = mpl.colormaps["plasma"].resampled(8)

exp = 48

init_path = "/Xnfs/physbiochrom/ppuel/data/" 

verif = False

for tmp_exp in os.listdir(init_path):
        if tmp_exp.split("_")[0] == f'EXP{exp}':
                exp = tmp_exp
                verif = not(verif)

if not(verif):
        print(f"The experience {exp} is not found or found multiple time")
        sys.exit()

output_path = f"/home/ppuel/data/{exp}/"

os.makedirs(output_path, exist_ok=True)

exp_path = os.path.join(init_path, exp)

c_max = 3**6*16
c=0

r_matrix_2D = [[],[]]
r_matrix_3D = [[],[]]
finish = np.zeros((3,3,3,3,3,3))

fig = plt.figure()
ax = fig.add_subplot()

# for jll in jll_list:
#         i = jll_list.index(jll)
#         for jlp in jlp_list:
#                 j = jlp_list.index(jlp)
#                 for jll_valency in jll_valency_list:
#                         k = jll_valency_list.index(jll_valency)
#                         for jlp_valency in jlp_valency_list:
#                                 l = jlp_valency_list.index(jlp_valency)
#                                 for jpl_valency in jpl_valency_list:
#                                         m = jpl_valency_list.index(jpl_valency)
#                                         for ldens in ldens_list:
#                                                 o = ldens_list.index(ldens)
#                                                 for n in range(N):
#                                                         print(f"{c/c_max*100:0.2f}", end = "\r")
#                                                         c += 1
#                                                         try : 
#                                                                 file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLL_VALENCY/{jll_valency}/JLP_VALENCY/{jlp_valency}/JPL_VALENCY/{jpl_valency}/LDENS/{ldens}/N/{n}/liq_droplets.pickle", "rb")
#                                                                 droplet_dict = pickle.load(file)
#                                                                 file.close()
#                                                                 for droplet in droplet_dict:
#                                                                         r_matrix_2D[o].append(((droplet.tau-1)**2+(np.max(droplet.size)-2)**2)**(1/2))
#                                                                         r_matrix_3D[o].append(((droplet.tau-1)**2+(np.max(droplet.size)-2)**2+(np.mean(droplet.local_density)-1)**2)**(1/2))
#                                                         except :
#                                                                 finish[i,j,k,l,m,o] += 1
# np.save(os.path.join(output_path, "r_matrix_2D_0.019.npy"), r_matrix_2D[0])
# np.save(os.path.join(output_path, "r_matrix_2D_0.038.npy"), r_matrix_2D[1])
# np.save(os.path.join(output_path, "r_matrix_3D_0.019.npy"), r_matrix_3D[0])
# np.save(os.path.join(output_path, "r_matrix_3D_0.038.npy"), r_matrix_3D[1])


# np.save(os.path.join(output_path, "finish.npy"), finish)

                                                

