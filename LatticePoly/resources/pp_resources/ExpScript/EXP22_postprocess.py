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
import itertools
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

ldens_list = [f"{i*0.003+0.001:0.3f}" for i in range(12)]

def conversion_density_to_yM(D, b=20e-9):
        D = np.array([float(i) for i in D])
        V_maille = b**3
        Na = 6.022e23
        n_particules_per_maille = 4*D
        yM = 1e3 * n_particules_per_maille/V_maille/Na
        return [f"{i:0.2f}" for i in yM]

def conversion_kBT_to_kJ_per_mol(J):
        J = np.array([float(i) for i in J])
        T = 300
        kB = 1.380e-23
        Na = 6.022e23
        res = J * kB * Na * T * 1e-3
        return([f"{i:0.2f}" for i in res])

ldens_list_conv = conversion_density_to_yM(ldens_list)



N = 1200

range_data = 10

exp = 22

init_path = "/Xnfs/physbiochrom/ppuel/data/"

verif = 0

for tmp_exp in os.listdir(init_path):
        if tmp_exp.split("_")[0] == f'EXP{exp}':
                exp = tmp_exp
                verif += 1

if verif != 1:
        print(f"The experience {exp} is not found or found multiple time")
        sys.exit()

output_path = f"/home/ppuel/data/{exp}/"

os.makedirs(output_path, exist_ok=True)

exp_path = os.path.join(init_path, exp)

c_max = 12*1200
c=0

# def add_len_2D(array, l):
#         l_tmp = len(array)
#         array = np.concatenate((array, np.zeros((l, l_tmp), dtype=np.int64)))
#         array = np.concatenate((array, np.zeros((l+l_tmp, l), dtype=np.int64)), axis=1)
#         return(array)


# for ldens in ldens_list:

#         matrix_data = np.astype(np.load(os.path.join(output_path, f"histogramme_goutte_{ldens}.npy")), np.uint64)

#         l = len(matrix_data)

#         for n in range(N):
#                 print(f"{n/N:0.3E}", end = '\r')
#                 file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/LDENS/{ldens}/N/{n+1200}/liq_droplets.pickle", "rb")
#                 droplet_dict = pickle.load(file)
#                 file.close()
#                 for droplet in droplet_dict.values():
#                         y = droplet.tau - 1
#                         x = np.max(droplet.sizes) - 2
#                         if max(x,y) > l-1:
#                                 matrix_data = add_len_2D(matrix_data, max(x,y) - l + 1)
#                                 l = max(x,y) + 1
#                         matrix_data[x,y] += 1

#         np.save(os.path.join(output_path, f"histogramme_goutte_{ldens}_2.npy"), matrix_data)
#         print(matrix_data)

fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot()


plasma = mpl.colormaps["plasma"].resampled(14)


for ldens in ldens_list:
        index = ldens_list.index(ldens)
        matrix_data = np.load(os.path.join(output_path, f"histogramme_goutte_{ldens}_2.npy"))
        print(matrix_data.sum())

        # set up the figure and Axes
        # fig = plt.figure(figsize=(12, 9))
        # ax1 = fig.add_subplot()
        # l_tmp = len(matrix_data)
        l = int(np.sqrt(len(matrix_data)))

        # for x in range(l_tmp):
        #         if np.isclose(matrix_data[-x-1,:].sum(),0,1e-5) and np.isclose(matrix_data[:,-x-1].sum(),0,1e-5):
        #                 l -= 1
        #         else:
        #                 break

        # matrix_data = matrix_data[:l,:l]

        _x = np.arange(l, dtype=np.int16)
        # print(np.shape(_x))
        # _xx, _yy = np.meshgrid(_x, _x)
        # x, y = _xx.ravel(), _yy.ravel()

        # print(np.shape(matrix_data))
        matrix_data = np.reshape(matrix_data, (l,l))

        list_diag = np.zeros(l*2)

        for i, j in itertools.product(_x, _x):
                list_diag[i+j] += matrix_data[i,j]

        list_diag /= list_diag.sum()

        list_diag /= (np.arange(l*2)+1)

        # print(np.shape(matrix_data))


        # top = # bottom = np.zeros_like(matrix_data)
        # width = depth = 1
        # print(top)

        ax.plot(np.where(list_diag != 0.0, np.log10(list_diag), None), color = plasma(1+index), label = f"{ldens_list_conv[index]}")

ax.set_xlabel('tau+size')
ax.set_ylabel('P_{gouttes} (log10)')
# ax.set_xlim(xmin=-0.5, xmax=5)

leg = fig.legend()
leg.set_title("densité de particules")
fig.savefig(os.path.join(output_path, f"hist/density_over_diag.png"))


