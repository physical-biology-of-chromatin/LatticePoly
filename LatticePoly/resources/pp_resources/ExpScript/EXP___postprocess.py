import os, sys, subprocess, pickle, h5py, itertools
from LiqDroplet import LifeTime, Droplet, Event
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
plt.style.use('./resources/h5py/presentation.mplstyle')
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

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

def find_exp(exp, init_path = "/Xnfs/physbiochrom/ppuel/data/"):
        verif = 0
        for tmp_exp in os.listdir(init_path):
                if tmp_exp.split("_")[0] == f'EXP{exp}':
                        exp_name = tmp_exp
                        verif += 1
        if verif != 1:
                print(f"The experience {exp} is not found or found multiple time")
                sys.exit()
        return(exp_name, os.path.join(init_path, exp_name))




ldens_list = [f"{i*0.003+0.001:0.3f}" for i in range(12)]
jll_list = [f"{i*0.2:0.1f}" for i in range(12)]
jll_valency_list = [f"{i+1:0d}" for i in range(12)]

ldens_list = ["0.031"]
jll_list = [f"{i*0.2:0.1f}" for i in range(12)]
jll_valency_list = [f"{i+1:0d}" for i in range(4)]

plasma = mpl.colormaps["plasma"].resampled(len(jll_list)+2)

N = 20

exp = 12
exp_name, exp_path = find_exp(exp)
output_path = f"/home/ppuel/data/{exp_name}/"
output_path_figure = f"/home/ppuel/data/{exp_name}/figure/"
os.makedirs(output_path, exist_ok=True)
os.makedirs(output_path_figure, exist_ok=True)

exp_tmp = 50
exp_name_tmp, exp_path_tmp = find_exp(exp_tmp)

c_max = len(jll_list)*len(jll_valency_list)*len(ldens_list)*N
c = 0

# file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp_name}/JLL/{jll}/JLL_VALENCY/{jll_valency}/LDENS/{ldens}/N/{n}/liq_droplets.pickle", "rb")
# droplet_dict = pickle.load(file)
# file.close()
