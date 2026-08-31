from matplotlib import axis, colormaps, layout_engine
import matplotlib.pyplot as plt
from networkx import rescale_layout
import pandas

import numpy as np
import h5py
import importlib
import sys
import os
import pickle
import subprocess
import itertools

from seaborn import reset_defaults

from hdf5Reader import hdf5Reader
from utils import find_exp, exp_mapping, exp_path_list

from scipy.spatial.distance import squareform
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit, Bounds
from tqdm import tqdm

plt.style.use('./resources/h5py/presentation.mplstyle')

expNum = 55
XnfsDir = "/Xnfs/physbiochrom/ppuel/data/"
expName = find_exp(expNum, XnfsDir)
expPath = os.path.join(XnfsDir, expName)
figDir = f"/home/ppuel/figure/{expName}/Distance/"


os.makedirs(figDir, exist_ok=True)

dict_parameters, _, metaParameterN, metaParameterNmeas = exp_mapping(expPath)

pathList = exp_path_list(dict_parameters, -1, expPath)

listDatasetName = ["liqMean"]

plotType = "plot"

fileName = "aggregated_process.h5"


# with h5py.File(os.path.join(pathList[0], "N/0/process.h5"), 'r') as hfile:
#         for dataset in hfile.keys():
#                 print(hfile[dataset].shape, dataset)


# AggregateData(listDatasetName, pathList, metaParameterN)

# PrintAggregateData(fileName, listDatasetName, plotType, expName)


def f(x, a, b, k, h):
        return(a + b * ((x / k)**h) / (1 + (x / k)**h))


def Fig1_mapping_jll_over_valency_with_liqMean():
    """
    This function is built to match for each self-interaction valency,
    the energy of phase transition. 
    The data are located in EXP12, 57 and 58
    To do that, we compute the average local density for each particle
    for each self-interaction energy. We fit the data by :
        
                            b x (x/k)^h
                f(x) = a + -------------
                            1 + (x/k)^h

    where f(0) = a, f(+∞) = a+b and f(k) = a + b/2 = (f(0) + f(+∞))/2
    so k is the energy of phase transition
    At h = 1, this function approximate a 2de order phase transition. 
    At h > 1, this function approximate a 1st order phase transition.

    For each valency, we restrict ourself to the increasing part of the 
    mean local density curve. We cut the data at index "c" to better fit
    them.
    """

    fig = plt.figure(figsize=[24, 14], layout = 'constrained')
    sfigs = fig.subfigures(2, 1, height_ratios=[4, 1])
    sfigs[0].suptitle("Initial density : 0.034")
    axs = sfigs[0].subplots(2,6)
    
    k_array = np.zeros(12)
    
    tqdm_55 = tqdm(total = 12*36)
    
    for val_id, jll_valency in enumerate([i+1 for i in range(12)]):
            ax = axs[val_id//6][val_id%6]
            jll_range = 36
            if jll_valency<9:
                    jll_array = [f"{k*0.2:0.1f}" for k in range(36)]
            else:
                    jll_array = ["0.0", "0.05", "0.10", "0.15", "0.2", "0.25", "0.30", "0.35", "0.4", "0.45", "0.50", "0.55", "0.6", "0.65", "0.70", "0.75", "0.8", "0.85", "0.90", "0.95", "1.0", "1.05", "1.10", "1.15", "1.2", "1.25", "1.30", "1.35", "1.4", "1.45", "1.50", "1.55", "1.6", "1.8", "2.0", "2.2"]
    
            dataset = np.zeros(jll_range)
    
            for jll_id, jll_str in enumerate(jll_array):
                    if (jll_valency < 9 and  jll_id < 12) or (jll_valency > 8 and len(jll_str) == 3) :
                            expNum = 12
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}/LDENS/0.034")
                            
                    elif (jll_valency < 9 and  jll_id > 11):
                            expNum = 57
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}")
                    else:
                            expNum = 58
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}")

                    with h5py.File(os.path.join(expPath, "aggregated_process.h5"), 'r') as hfile:
                            dataset[jll_id] = np.mean(hfile["liqMean"][-10:])
                    
                    tqdm_55.update()
    
            ax.set_title(f'valency = {jll_valency}')

            threshold = -np.inf
            c = 0
            if jll_valency > 4:
                    while c < 36 and threshold < dataset[c]:
                            threshold = dataset[c]
                            c += 1
            else:
                    c = 36
     
            jll_float_array = np.array(list(map(float,jll_array[:c])))
            dataset = dataset[:c]

            ax.scatter(jll_float_array, dataset)
            
            a, b, k, h = curve_fit(f, jll_float_array, dataset, p0 = [0.04,0.4,1,6], bounds = Bounds([0.03, 0.05, 0.5, 0], [0.05, 0.9, 2.8, 13]))[0]
            ax.text(0.7, 0.1, f"a = {a:0.3f}\nb = {b:0.3f}\nk = {k:0.3f}\nh = {k:0.3f}", transform=ax.transAxes)
            ax.set_xlabel("Energy of self-interaction")
            ax.set_ylabel("Mean local density")

            fit_linspace = np.linspace(0, jll_float_array[c-1], 100)
            ax.plot(fit_linspace, f(fit_linspace, a, b, k, h), '--',  color = 'tomato', alpha = 0.9)

            k_array[val_id] = k

            ax.plot([k, k], [min(dataset), max(dataset)], '--', color = 'black')
    
    ax = sfigs[1].subplots(1,1)
    ax.scatter(np.arange(1, 13), k_array, c = 'red')
    ax.tick_params(axis = 'y', labelcolor = 'red')
    ax.set_xlabel("Valency")
    ax.set_ylabel("Energy threshold")
    ax.set_box_aspect(1/3)
    fig.savefig(os.path.join(figDir, 'Fig1_mapping_jll_over_valency_with_liqMean.png'))
    
    plt.close(fig=fig)
    
    tqdm_55.close()


Fig1_mapping_jll_over_valency_with_liqMean()




def jll_jll_valency_matching_with_simTime():
    
    
    fig = plt.figure(figsize=[24, 12], layout = 'constrained')
    
    sfigs = fig.subfigures(2, 1, height_ratios=[4, 1])
    
    axs = sfigs[0].subplots(2,6)
    
    k_array = np.zeros(12)
    k_array2 = np.zeros(12)
    k_array3 = np.zeros(12)
    
    b_list =  [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    x0_list = [3,    3,3,3,   1.3,   1,   0.9,   0.8,   0.7, 0.7,  0.7,   0.7]
    l1_list = [0.01,    0.01,   0.01,   0.01,   3,   3,   3,   3,   3,   3,   3,   3]
    l2_list = [0.01,    0.01,   0.01,   0.01,   2,   2,   2,   2,   4,   4,   4,   7]
    
    tqdm_55 = tqdm(total = 12*36*4)
    
    for i in range(12):
    
            # b_tmp, x0_tmp, l1_tmp, l2_tmp = b_list[i], x0_list[i], l1_list[i], l2_list[i]
    
            ax = axs[i//6][i%6]
    
            jll_valency = i+1
    
            jll_range = 36
    
            if jll_valency<9:
                    jll_array = [f"{k*0.2:0.1f}" for k in range(36)]
            else:
                    jll_array = ["0.0", "0.05", "0.10", "0.15", "0.2", "0.25", "0.30", "0.35", "0.4", "0.45", "0.50", "0.55", "0.6", "0.65", "0.70", "0.75", "0.8", "0.85", "0.90", "0.95", "1.0", "1.05", "1.10", "1.15", "1.2", "1.25", "1.30", "1.35", "1.4", "1.45", "1.50", "1.55", "1.6", "1.8", "2.0", "2.2"]
    
    
            dataset = np.zeros(jll_range)
    
            color_list = []
    
            MCSmax = - np.inf
    
            for j, jll_str in enumerate(jll_array):
                    if (jll_valency < 9 and  j < 12) or (jll_valency > 8 and len(jll_str) == 3) :
                            expNum = 12
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}/LDENS/0.034")
    
                    elif (jll_valency < 9 and  j > 11):
                            expNum = 57
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}")
    
                    else:
                            expNum = 58
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}")
    
                    with h5py.File(os.path.join(expPath, "post_process.h5"), 'r') as hfile:
                            MCSmax = max(MCSmax, hfile["SimTime"][1])
                    tqdm_55.update()
            
            
    
            MCSmin = np.inf
    
            for j, jll_str in enumerate(jll_array):
                    if (jll_valency < 9 and  j < 12) or (jll_valency > 8 and len(jll_str) == 3) :
                            expNum = 12
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}/LDENS/0.034")
    
                    elif (jll_valency < 9 and  j > 11):
                            expNum = 57
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}")
    
                    else:
                            expNum = 58
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}")
    
                    with h5py.File(os.path.join(expPath, "post_process.h5"), 'r') as hfile:
                            if hfile["SimTime"][-1] >= MCSmax:
                                    MCSmin = min(MCSmin, hfile["SimTime"][-1]) 
                    tqdm_55.update()
    
            pathMCS = {}
    
            for j, jll_str in enumerate(jll_array):
                    if (jll_valency < 9 and  j < 12) or (jll_valency > 8 and len(jll_str) == 3) :
                            expNum = 12
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}/LDENS/0.034")
    
                    elif (jll_valency < 9 and  j > 11):
                            expNum = 57
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}")
    
                    else:
                            expNum = 58
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}")
    
                    with h5py.File(os.path.join(expPath, "post_process.h5"), 'r') as hfile:
                            if expNum < 13 : # hfile["SimTime"][-1] < MCSmax:
                                    pathMCS[expPath] = hfile["SimTime"][1]*10**(1/2) #np.nan
                            else:
                                    pathMCS[expPath] = hfile["SimTime"][1] #min(max(int(MCSmax // hfile["SimTime"][1])+1, int(MCSmin // hfile["SimTime"][1])), int(metaParameterNmeas))
                    tqdm_55.update()
    
            for j, jll_str in enumerate(jll_array):
                    if (jll_valency < 9 and  j < 12) or (jll_valency > 8 and len(jll_str) == 3) :
                            expNum = 12
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}/LDENS/0.034")
                            color_list.append("red")
    
                    elif (jll_valency < 9 and  j > 11):
                            expNum = 57
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}")
                            color_list.append("blue")
                    else:
                            expNum = 58
                            expName = find_exp(expNum, XnfsDir)
                            expPath = os.path.join(XnfsDir, expName, f"JLL/{jll_str}/JLL_VALENCY/{jll_valency:0d}")
                            color_list.append("green")
    
                    with h5py.File(os.path.join(expPath, "aggregated_process.h5"), 'r') as hfile:
                            dataset[j] = np.nan if np.isnan(pathMCS[expPath]) else pathMCS[expPath] #hfile["liqMean"][pathMCS[expPath]]
                    
                    
                    
                    
                    
                    tqdm_55.update()
    
            # ax_twinx = ax.twinx()
            # ax_twinx.plot(np.arange(jll_range-1)+0.5, np.diff(dataset)/2)
    
    
            ax.set_title(f'valency = {jll_valency}')
    
            # c_liqmean = -np.inf
            c = 36 # 0
            # error = 1
            # if jll_valency > 4:
            #         while c < 36 and c_liqmean < dataset[c] and error:
            #                 c_liqmean = dataset[c]
            #                 c += 1
            #                 error == 0
            # else:
            #         c = 36
    
            c -= 1
    
    
            jll_float_array = np.array(list(map(float,jll_array))[:c+1])
    
            ax.scatter(jll_float_array[:c+1], dataset[:c+1], c = color_list)
            ax.text(0.7, 0.1, f"MCSmax = {MCSmax:0.1f}\nMCSmin = {MCSmin:0.1f}", transform=ax.transAxes)
    
            # res = curve_fit(f, list(map(float,jll_array[:c+1])), dataset[:c+1], p0 = [0.04,0.4,1,6], bounds = Bounds([0.03, 0.05, 0.5, 0], [0.05, 0.9, 2.8, 13]))
            # ax.text(0.7, 0.1, f"a = {res[0][0]:0.3f}\nb = {res[0][1]:0.3f}\nk = {res[0][2]:0.3f}\nh = {res[0][3]:0.3f}\nf'(k) = {res[0][3]/res[0][2]*(res[0][1]-res[0][0])/4:0.3f}", transform=ax.transAxes)
    
            # res2 = curve_fit(lambda x, b, x0, lambda1, lambda2 : f2(x, b, np.max(dataset), x0, lambda1, lambda2), list(map(float,jll_array[:c+1])), dataset[:c+1], p0 = [b_tmp, x0_tmp, l1_tmp, l2_tmp], bounds = Bounds(0, 1))
            # ax.text(0.7, 0.4, f"b = {res2[0][0]:0.3f}\nx0 = {res2[0][1]:0.3f}\nk1 = {res2[0][2]:0.3f}\nk2 = {res2[0][3]:0.3f}", transform=ax.transAxes)
    
            # res3 = curve_fit(lambda x, alpha: alpha*f(x, *res[0]) + (1-alpha)*f2(x, res2[0][0], np.max(dataset), *res2[0][1:]), list(map(float,jll_array[:c+1])), dataset[:c+1], p0 = [0.5], bounds = bounds_alpha)
            # ax.text(0.7, 0.7, f"alpha = {res3[0][0]:0.3f}", transform=ax.transAxes)
    
    
            # ax.plot([res[0][2]*0.8, res[0][2]*1.2], [(res[0][0]+res[0][1])/2-res[0][2]*res[0][3]*0.2, (res[0][0]+res[0][1])/2+res[0][2]*res[0][3]*0.2], '.-', 'green')
    
            # def dev(x):
            #         return((res[0][0]+res[0][1])/2+(res[0][1]-res[0][0])/4/res[0][2]*res[0][3]*(x-res[0][2]))
    
    
    
    
            # tmp_array = np.linspace(res[0][2]*(1 - 2/res[0][3]), res[0][2]*(1 + 2/res[0][3]), 10)
            # ax.plot(tmp_array, dev(tmp_array), '--', color = 'green')
    
    
            # ax.plot(np.linspace(0, float(jll_array[c]), 100), f(np.linspace(0, float(jll_array[c]), 100), *res[0]), '--',  color = 'red', alpha = 0.7)
    
            # ax.plot(np.linspace(0, float(jll_array[c]), 100), [f2(x, res2[0][0], np.max(dataset), res2[0][1], res2[0][2], res2[0][3]) for x in np.linspace(0, float(jll_array[c]), 100)], '--',  color = 'green', alpha = 0.7)
    
            # ax.plot(np.linspace(0, float(jll_array[c]), 100), res3[0][0]*f(np.linspace(0, float(jll_array[c]), 100), *res[0]) + (1-res3[0][0]) * f2(np.linspace(0, float(jll_array[c]), 100), res2[0][0], np.max(dataset), res2[0][1], res2[0][2], res2[0][3]), color = 'purple')
    
            # alpha = res3[0][0]
    
            # k_array[i] = res[0][2]
    
            # print(np.sqrt(np.diag(res[1]))[2])
    
            # ax.plot([res[0][2], res[0][2]], [min(dataset), max(dataset)], '--', color = 'black')
    
            # k_array[i] = res[0][2] if jll_valency < 7 else res2[0][1]
            # k_array3[i] = alpha*res[0][2]+(1-alpha)*res2[0][1]
    
    
            # if i in [0, 3, 6, 9]:
            #         pandas.DataFrame({'jll' : list(map(float,jll_array)), 'liqMean' : dataset}).to_csv(os.path.join(figDir, f'geogebra_values_valency_{jll_valency}.csv'), float_format = '%.3f')
    
    ax = sfigs[1].subplots(1,2)[0]
    # ax_twinx = ax.twinx()
    
    # ax.scatter(np.arange(1, 13), k_array, c = 'red', label  = "k")
    # ax.scatter(np.arange(1, 13), k_array2, c = 'blue', label  = "k")
    # ax.scatter(np.arange(1, 13), k_array3, c = 'purple', label  = "k")
    # ax_twinx.scatter(np.arange(1, 13), h_array, c = 'blue', label  = "h")
    sfigs[1].legend()
    ax.tick_params(axis = 'y', labelcolor = 'red')
    # ax_twinx.tick_params(axis = 'y', labelcolor = 'blue')
    
    # np.savetxt(os.path.join(figDir, 'k_list.txt'), k_array)
    
    # res = curve_fit(lambda x, a, b : a*x+b, np.arange(1, 8), k_array[:7])[0]
    # res1 = curve_fit(lambda x, a, b : a*x+b, np.arange(7, 13), k_array[6:])[0]
    
    # ax.plot(np.arange(1, 8), res[0]*np.arange(1, 8)+res[1])
    # ax.plot(np.arange(7, 13), res1[0]*np.arange(7, 13)+res1[1])
    
    fig.savefig(os.path.join(figDir, 'concatenation_liqMean_SimTime.png'))
    
    plt.close(fig=fig)
    
    tqdm_55.close()
    
    
