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




ldens_list = ["0.019", "0.038"]
ninter_list = ["10000", "100000"]

plasma = mpl.colormaps["plasma"].resampled(len(ldens_list)+2)

N = 48

exp = 49
exp_name, exp_path = find_exp(exp)
output_path = f"/home/ppuel/data/{exp_name}/"
output_path_figure = f"/home/ppuel/data/{exp_name}/figure/"
os.makedirs(output_path, exist_ok=True)
os.makedirs(output_path_figure, exist_ok=True)

# exp_tmp = 50
# exp_name_tmp, exp_path_tmp = find_exp(exp_tmp)



def extract_exp(exp):

        outputPathList = []
        
        _, exp_path = find_exp(exp)
        
        tmpList = list(map(lambda x: os.path.join(exp_path, x), os.listdir(exp_path)))
        
        bool_metadata = 0
        a = 1e5
        c = 0
        
        while len(tmpList) > 0 and c < a:
                c += 1
                tmpPath = tmpList[0]
                
                if os.path.isfile(tmpPath) or os.listdir(tmpPath) == []:
                        tmpList.remove(tmpPath)
                
                elif "N" in os.listdir(tmpPath):
                        outputPathList.append(tmpPath)
                        tmpList.remove(tmpPath)
                        if not(bool_metadata):
                                N = len(os.listdir(os.path.join(tmpPath,'N')))
                                bool_metadata = 1
                
                else:
                        for tmpFile in os.listdir(tmpPath):
                                if os.path.isdir(os.path.join(tmpPath,tmpFile)):
                                        tmpList.append(os.path.join(tmpPath,tmpFile))
                        tmpList.remove(tmpPath)

        return(outputPathList, N)

def extract_droplet(exp):

        exp_name, _ = find_exp(exp)
        output_path = f"/home/ppuel/data/{exp_name}/"
        os.makedirs(output_path, exist_ok=True)

        # outputPathList, N = extract_exp(exp)

        c_max = len(ninter_list)*len(ldens_list)*N
        c = 0

        for ldens, ninter in itertools.product(ldens_list, ninter_list):
                matrix_data = np.zeros(int(float(0.019)*48**3*4), dtype=np.int64)
                for n in range(N):
                        print(f"{c/c_max*100:0.1f}%", end="\r")
                        c += 1
                        try: 
                                file = open(f"{exp_path}/LDENS/{ldens}/NINTER/{ninter}/N/{n}/liq_droplets.pickle", "rb")
                                droplet_dict = pickle.load(file)
                                file.close()
                                for droplet in droplet_dict.values():
                                        indice = int(np.max(droplet.sizes))
                                        matrix_data[indice] += 1
                        except:
                                pass

                np.save(os.path.join(output_path, f'matrix_size_{ldens}_{ninter}.npy'), matrix_data)

# extract_droplet(exp)

for ninter in ninter_list:
        fig = plt.figure(figsize=(9,6))
        ax = fig.add_subplot()
        for ldens in ldens_list:
                matrix_data = np.load(os.path.join(output_path, f"matrix_size_{ldens}_{ninter}.npy"))
                ax.plot(np.where(matrix_data > 0, np.log10(matrix_data), np.nan), linewidth = 1.5, label = f"ldens = {ldens}", color = plasma(ldens_list.index(ldens)+1))
                i = 0
                while i < len(matrix_data) and matrix_data[len(matrix_data)-1-i] == 0:
                        i+=1
                # ax.plot([(len(matrix_data)-1-i), (len(matrix_data)-1-i)], [np.max(np.where(matrix_data > 0, np.log10(matrix_data), -1))/2, np.max(np.where(matrix_data > 0, np.log10(matrix_data), -1))],linewidth = 1, linestyle = 'dashed', color = plasma(ldens_list.index(ldens)+1))
        fig.legend()
        fig.savefig(os.path.join(output_path_figure, f"fig_size_ninter={ninter}.png"))
