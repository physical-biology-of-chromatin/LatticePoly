from matplotlib import axis, colormaps
import matplotlib.pyplot as plt

import numpy as np
import h5py, importlib, sys, os, pickle, subprocess, itertools

from hdf5Reader import hdf5Reader
from utils import *
from itertools import product, zip_longest
from scipy.spatial.distance import squareform
from matplotlib.colors import LogNorm

plt.style.use('./resources/h5py/presentation.mplstyle')

expNum = 52
XnfsDir = "/Xnfs/physbiochrom/ppuel/data/"
expName = find_exp(expNum, XnfsDir)
expPath = os.path.join(XnfsDir, expName)
figDir = f"/home/ppuel/figure/{expName}/Distance/"


os.makedirs(figDir, exist_ok=True)

dict_parameters, _, metaParameterN, metaParameterNmeas = exp_mapping(expPath)

pathList = exp_pathList(dict_parameters, -1, expPath) 

listDatasetName = ["distanceMap"]

# AggregateData(expPath, listDatasetName, pathList, metaParameterN)

# PrintAggregateData(expPath, figDir, listDatasetName, dict_parameters, metaParameterN, metaParameterNmeas, "map")

with h5py.File(os.path.join(pathList[0], "N/0/process.h5"), 'r') as hfile:
        print(hfile.keys())
















# def conversion_density_to_yM(D, b=20e-9):
#         D = np.array([float(i) for i in D])
#         V_maille = b**3
#         Na = 6.022e23
#         n_particules_per_maille = 4*D
#         yM = 1e3 * n_particules_per_maille/V_maille/Na
#         return [f"{i:0.2f}" for i in yM]

# def conversion_kBT_to_kJ_per_mol(J):
#         J = np.array([float(i) for i in J])
#         T = 300
#         kB = 1.380e-23
#         Na = 6.022e23
#         res = J * kB * Na * T * 1e-3
#         return([f"{i:0.2f}" for i in res])

# def find_exp(exp, init_path = "/scratch/Cascade/ppuel/data/"):
#         verif = 0
#         for tmp_exp in os.listdir(init_path):
#                 if tmp_exp.split("_")[0] == f'EXP{exp}':
#                         exp_name = tmp_exp
#                         verif += 1
#         if verif != 1:
#                 print(f"The experience {exp} is not found or found multiple time")
#                 sys.exit()
#         return(exp_name, os.path.join(init_path, exp_name))




# jll_list = [f"{1.0*2**(i-1):0.1f}" for i in range(3)]
# jll_valency_list = [f"{i+1}" for i in range(6)]
# jlp_list = [f"{0.49*2**(i-1):0.1f}" for i in range(3)]
# jpl_valency_list = [f"{i+1}" for i in range(6)]
# domain_list = ["Lucy_domain.in","Lucy_domain_50.in"]

# jll_list = jll_list
# jll_valency_list = jll_valency_list
# jlp_list = jlp_list
# jpl_valency_list = jpl_valency_list
# domain_list = domain_list[:1]

# plasma = mpl.colormaps["plasma"].resampled(len(jll_list)+2)

# N = 2

# exp = 52
# exp_2 = 50

# exp_name, exp_path = find_exp(exp)
# exp_name_2, exp_path_2 = find_exp(exp_2)

# jll_list_2 = [f"{i*0.2:0.1f}" for i in range(12)]

# output_path_classifier = f"/home/ppuel/data/{exp_name}/classifier/"
# output_path_figure_classifier = f"/home/ppuel/data/{exp_name}/figure/"

# os.makedirs(output_path_classifier, exist_ok=True)
# os.makedirs(output_path_figure_classifier, exist_ok=True)


# c_max = len(jll_list)*len(jll_valency_list)*len(jlp_list)*len(jpl_valency_list)*len(domain_list)*2
# c=0
# nliq = 0.019*48**3*4

# def f(x):
#         if x<1:
#                 return(x**4) 
#         elif x<3:
#                 return(1)
#         else:
#                 return(np.exp(3-x))                         

# def metric_WT(liqDropNum, liqDropSize, liqFraction):
#         return((liqFraction**4)*((np.tanh((liqDropSize/(nliq*liqFraction)-0.6)*6)+1)/2)*(f(liqDropNum)))

# def metric_TKO(liqDropNum, liqDropSize, liqFraction):
#         return((liqFraction**4)*((np.tanh((liqDropSize/(nliq*liqFraction)-0.6)*6)+1)/2)*(f(liqDropNum)))


# for domain in domain_list:
#         for jll in jll_list:
#                 i = jll_list.index(jll)
#                 for jll_valency in jll_valency_list:
#                         j = jll_valency_list.index(jll_valency)
#                         for jlp in jlp_list:
#                                 k = jlp_list.index(jlp)

#                                 for jpl_valency in jpl_valency_list:
#                                         l = jpl_valency_list.index(jpl_valency)
#                                         fig = plt.figure(figsize=(24,12))
#                                         ax_LiqHist = fig.add_subplot(1,2,1)
#                                         ax_HetHist = fig.add_subplot(1,2,2)
#                                         HetHist = np.zeros((100,108), dtype=np.float64)
#                                         LiqHist = np.zeros((100,108), dtype=np.float64)
                                                        
#                                         for n in range(N):
#                                                 print(f"{c/c_max*100:0.2f}%",end='\r')
#                                                 c+=1
                                                                
#                                                 try:
#                                                         with h5py.File(f"/scratch/Cascade/ppuel/data/{exp_name}/JLL/{jll}/JLP/{jlp}/JLL_VALENCY/{jll_valency}/JPL_VALENCY/{jpl_valency}/DOMAINPATH/{domain}/N/{n}/process.h5") as hfile:
#                                                                 LiqHist += hfile["LiqHist"][1:,:]
#                                                                 HetHist += hfile["HetHist"][1:,:]
                                                
#                                                 except:
#                                                         print(f"/scratch/Cascade/ppuel/data/{exp_name}/JLL/{jll}/JLP/{jlp}/JLL_VALENCY/{jll_valency}/JPL_VALENCY/{jpl_valency}/DOMAINPATH/{domain}/process.h5")
                        


#                                         im = ax_LiqHist.imshow(np.where(LiqHist>0, LiqHist, np.nan), cmap='plasma', origin='lower')
#                                         im = ax_HetHist.imshow(np.where(HetHist>0, HetHist, np.nan), cmap='plasma', origin='lower')
#                                         fig.suptitle(f"{jll}, {jlp}, {jll_valency}, {jpl_valency}")
#                                         fig.savefig(os.path.join(output_path_figure_classifier,f"figure_hist_{jll}_{jlp}_{jll_valency}_{jpl_valency}_{domain}.png"))
#                                         plt.close(fig=fig)
# # save matrix size
# if 0:
#         for jll in jll_list:
#                 for jll_valency in jll_valency_list:
#                         matrix_data_size = np.zeros((1800))
#                         matrix_data_rholoc = np.zeros((11*100))
                
#                         for n in range(N):
#                                 print(f"{c/c_max*100:0.1f}%", end = '\r')
#                                 c += 1
#                                 file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLL_VALENCY/{jll_valency}/N/{n}/liq_droplets.pickle", "rb")
#                                 droplet_dict = pickle.load(file)
#                                 file.close()
#                                 for droplet in droplet_dict.values():
#                                         y = np.max(droplet.sizes) - 2
#                                         matrix_data_size[y] += 1
#                                         rho_mean = int((np.mean(droplet.local_density)-1)*100)
#                                         matrix_data_rholoc[rho_mean] += 1

#                         np.save(os.path.join(output_path_classifier, f"data_{jll}_{jll_valency}_size.npy"), matrix_data_size)
#                         np.save(os.path.join(output_path_classifier, f"data_{jll}_{jll_valency}_rho.npy"), matrix_data_rholoc)
        


        
# if 0:
#         for jll_valency in jll_valency_list:
#                 fig = plt.figure(figsize=(9,6))
#                 ax = fig.add_subplot()
#                 for jll in jll_list:
#                         matrix_data = np.load(os.path.join(output_path_classifier, f"data_{jll}_{jll_valency}_rho.npy"))
#                         ax.plot(np.arange(len(matrix_data))/100+1, np.where(matrix_data > 0, np.log10(matrix_data), np.nan), linewidth = 1.5, label = f"jll = {jll}", color = plasma(jll_list.index(jll)+1))
#                         # ax.scatter(np.arange(len(matrix_data))/100+1, np.where(matrix_data > 0, np.log10(matrix_data), np.nan), s = 40, label = f"jll = {jll}", color = plasma(jll_list.index(jll)+1))
#                         i = 0
#                         while matrix_data[len(matrix_data)-1-i] == 0:
#                                 i+=1
#                         ax.plot([(len(matrix_data)-1-i)/100+1, (len(matrix_data)-1-i)/100+1], [np.max(np.where(matrix_data > 0, np.log10(matrix_data), -1))/2, np.max(np.where(matrix_data > 0, np.log10(matrix_data), -1))],linewidth = 1, linestyle = 'dashed', color = plasma(jll_list.index(jll)+1))
#                 fig.legend()
#                 fig.savefig(os.path.join(output_path_figure_classifier, f"fig_rho_valency={jll_valency}.png"))


# # save matrix droplet
# if 0:
#         def add_len_2D(array, l):
#                 l_tmp = len(array)
#                 array = np.concatenate((array, np.zeros((l, l_tmp), dtype=np.int64)))
#                 array = np.concatenate((array, np.zeros((l+l_tmp, l), dtype=np.int64)), axis=1)
#                 return(array)

#         for ldens in ldens_list:
#                 for jll in jll_list:
#                         for jll_valency in jll_valency_list:


#                                 matrix_data = np.zeros((3,3))
                                
#                                 l = len(matrix_data)

#                                 for n in range(N):
#                                         print(f"{c/c_max*100:0.1f}%", end = '\r')
#                                         c += 1
#                                         file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLL_VALENCY/{jll_valency}/LDENS/{ldens}/N/{n}/liq_droplets.pickle", "rb")
#                                         droplet_dict = pickle.load(file)
#                                         file.close()
#                                         for droplet in droplet_dict.values():
#                                                 y = droplet.tau - 1
#                                                 x = np.max(droplet.sizes) - 2
#                                                 if max(x,y) > l-1:
#                                                         matrix_data = add_len_2D(matrix_data, max(x,y) - l + 1)
#                                                         l = max(x,y) + 1
#                                                 matrix_data[x,y] += 1

#                                 np.save(os.path.join(output_path_classifier, f"data_{jll}_{jll_valency}_{ldens}.npy"), matrix_data)
#                                 print(matrix_data)

# # try to plot image + hist
# if 0:
#         plasma = mpl.colormaps["plasma"].resampled(14)

#         for ldens in ldens_list:
#                 for jll in jll_list:

#                         fig = plt.figure(figsize=(12, 9))
#                         ax = fig.add_subplot()

#                         for jll_valency in jll_valency_list:
#                                 print(jll, end = '\r')
#                                 index = jll_valency_list.index(jll_valency)
#                                 matrix_data = np.load(os.path.join(output_path_classifier, f"data_{jll}_{jll_valency}_{ldens}.npy"))
#                                 print(matrix_data.sum())
#                                 L = len(matrix_data)
#                                 # set up the figure and Axes
#                                 fig1 = plt.figure(figsize=(12, 9))
#                                 ax1 = fig1.add_subplot()
                                
#                                 l = int(np.sqrt(len(matrix_data)))
#                                 _x = np.arange(l, dtype=np.int16)
#                                 list_diag = np.zeros(l*2)

#                                 for i, j in itertools.product(_x, _x):
#                                         list_diag[i+j] += matrix_data[i,j]
#                                 # list_diag /= (np.arange(l*2)+1)
#                                 ax.plot(np.where(list_diag > 0, np.log10(list_diag), None), color = plasma(1+index), label = f"{jll_valency_list[index]}")
                                
#                                 xlim = 0
#                                 while np.sum(matrix_data[:,L-1-xlim]) < 1:
#                                         xlim += 1
                                
#                                 ylim = 0
#                                 while np.sum(matrix_data[L-1-ylim,:]) < 1:
#                                         ylim += 1


#                                 matrix_data = matrix_data[:L-ylim,:L-xlim]
#                                 ratio = (L-xlim)/(L-ylim)
                                
#                                 print('aspect', (L-xlim)/(L-ylim))

#                                 # the scatter plot:
#                                 im = ax1.imshow(np.where(matrix_data > 0, np.log10(matrix_data), np.nan), cmap = 'plasma', origin='lower')
                                

#                                 # create new Axes on the right and on the top of the current Axes
#                                 divider = make_axes_locatable(ax1)
#                                 # below height and pad are in inches
#                                 # ax1_histx = divider.append_axes("top", 1, pad=0.1, sharex=ax1) #, aspect = ratio)
#                                 # ax1_histy = divider.append_axes("right", ratio, pad=0.1, sharey=ax1, aspect = ratio)
#                                 # print(ax1_histx.get_aspect(), ax1_histy.get_aspect(), ax1.get_aspect())
#                                 # # make some labels invisible
#                                 # ax1_histx.xaxis.set_tick_params(labelbottom=False)
#                                 # ax1_histy.yaxis.set_tick_params(labelleft=False)


#                                 ax1.set_aspect('auto')
                                

#                                 # now determine nice limits by hand:
                                
#                                 # ax1_histx.plot(np.where(np.sum(matrix_data,axis=0) > 0, np.log10(np.sum(matrix_data,axis=0)), None))
#                                 # ax1_histy.plot(np.where(np.sum(matrix_data,axis=1) > 0, np.log10(np.sum(matrix_data,axis=1)), None), np.arange(len(matrix_data)))

#                                 # the xaxis of ax_histx and yaxis of ax_histy are shared with ax,
#                                 # thus there is no need to manually adjust the xlim and ylim of these
#                                 # axis.

#                                 ax1.set_xlabel('tau')
#                                 ax1.set_ylabel('size')
#                                 fig1.savefig(os.path.join(output_path_figure_classifier, f"density_{jll}_{jll_valency}_{ldens}.png"))

#                         ax.set_xlabel('tau+size')
#                         ax.set_ylabel('Pgouttes (log10)')

#                         leg = fig.legend()
#                         fig.savefig(os.path.join(output_path_figure_classifier, f"density_over_diag_{jll}_{ldens}.png"))


#                                 # for x in range(l_tmp):
#                                 #         if np.isclose(matrix_data[-x-1,:].sum(),0,1e-5) and np.isclose(matrix_data[:,-x-1].sum(),0,1e-5):
#                                 #                 l -= 1
#                                 #         else:
#                                 #                 break

#                                 # matrix_data = matrix_data[:l,:l]

                        
#                                 # print(np.shape(_x))
#                                 # _xx, _yy = np.meshgrid(_x, _x)
#                                 # x, y = _xx.ravel(), _yy.ravel()

#                                 # print(np.shape(matrix_data))

                                
#                                 # list_diag /= list_diag.sum()

                                
#                                 # print(list_diag)
#                                 # print(np.shape(matrix_data))


#                                 # top = # bottom = np.zeros_like(matrix_data)
#                                 # width = depth = 1
#                                 # print(top)  
                                
                        

