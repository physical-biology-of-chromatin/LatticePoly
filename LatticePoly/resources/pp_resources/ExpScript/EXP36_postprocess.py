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

plasma = mpl.colormaps["plasma"].resampled(258)

exp = "EXP36_bithorax_light_sweep"
exp2 = "EXP37_bithorax_null_hypothesis"

os.makedirs(f"/home/ppuel/data/{exp}/droplet_tracking", exist_ok=True)

jll_list = ["1.0","1.5","2.0","4.0"]
jlp_list = ["1.0","1.5","2.0","4.0"]
jlpp_list = ["1.0","1.5","2.0","4.0"]


jll_valency_list = ["2","4","6"]
jlp_valency_list = ["2","4","6"]
jlpp_valency_list = ["2","4","6"]

ldens_list = ["0.0009","0.0018","0.0036"]
N = 1



# ldens_list = ["0.010","0.025"]
# val_list = [str(i*2) for i in range(1,7)]


# jll_list = ["0.0"]
# ldens_list = ["0.010"]
# val_list = ["2"]


def scatter_hist(x, y, corr_x, corr_y, color, ax, ax_histx, ax_histy):
        # no labels
        ax_histx.tick_params(axis="x", labelbottom=False)
        ax_histy.tick_params(axis="y", labelleft=False)

        # the scatter plot:
        ax.scatter(x, y, c = color)

        x_min = np.min(x)
        x_range = np.max(x)-np.min(x)
        x_bins = []
        for i in range(int(x_range)+2):
                x_bins += [i + np.min(x)-.5, i + np.min(x)-.5]
        x_hist = np.zeros(x_range*2+4)
        x_corr_hist = -np.zeros(x_range*2+4)
        for x_i in x:
                x_hist[(x_i-x_min)*2+1:(x_i-x_min)*2+3] += 1
        
        x_hist_log10 = np.where(x_hist < 1.5, 0.15, 0) + np.where(x_hist < .5, -0.15, 0) + np.where(x_hist > 1.5, np.log10(x_hist), 0)

        ax_histx.plot(x_bins, x_hist_log10)

        for x_i in corr_x:
                x_corr_hist[(x_i-x_min)*2+1:(x_i-x_min)*2+3] += 1/160*4**3

        x_corr_hist_log10 = np.where(x_corr_hist < 1.5, 0.15, 0) + np.where(x_corr_hist < .5, -0.15, 0) + np.where(x_corr_hist > 1.5, np.log10(x_corr_hist), 0)

        ax_histx.plot(x_bins, -x_corr_hist_log10)

        delta_x_hist = np.where((x_hist - x_corr_hist) > 1.5, np.log10(x_hist - x_corr_hist), 0) + np.where((x_hist - x_corr_hist) < -1.5, -np.log10(-(x_hist - x_corr_hist)), 0) + np.where(np.logical_and((x_hist - x_corr_hist) < 1.5, (x_hist - x_corr_hist) > .5), 0.15, 0) + np.where(np.logical_and((x_hist - x_corr_hist) > -1.5, (x_hist - x_corr_hist) < -.5), -0.15, 0)

        ax_histx.plot(x_bins, delta_x_hist)


        y_min = np.min(y)
        y_range = np.max(y)-np.min(y)
        y_bins = []
        for i in range(int(y_range)+2):
                y_bins += [i + np.min(y)-.5, i + np.min(y)-.5]
        y_hist = np.zeros(y_range*2+4)
        y_corr_hist = np.zeros(y_range*2+4)
        for y_i in y:
                y_hist[(y_i-y_min)*2+1:(y_i-y_min)*2+3] += 1

        y_hist_log10 = np.where(y_hist < 1.5, 0.15, 0) + np.where(y_hist < .5, -0.15, 0) + np.where(y_hist > 1.5, np.log10(y_hist), 0)

        ax_histy.plot(y_hist_log10, y_bins)
        

        for y_i in corr_y:
                y_corr_hist[(y_i-y_min)*2+1:(y_i-y_min)*2+3] += 1/160*4**3

        
        y_corr_hist_log10 = np.where(y_corr_hist < 1.5, 0.15, 0) + np.where(y_corr_hist < .5, -0.15, 0) + np.where(y_corr_hist > 1.5, np.log10(y_corr_hist), 0)

        ax_histy.plot(-y_corr_hist_log10, y_bins)

        delta_y_hist = np.where((y_hist - y_corr_hist) > 1.5, np.log10(y_hist - y_corr_hist), 0) + np.where((y_hist - y_corr_hist) < -1.5, -np.log10(-(y_hist - y_corr_hist)), 0) + np.where(np.logical_and((y_hist - y_corr_hist) < 1.5, (y_hist - y_corr_hist) > .5), 0.15, 0) + np.where(np.logical_and((y_hist - y_corr_hist) > -1.5, (y_hist - y_corr_hist) < -.5), -0.15, 0)

        ax_histy.plot(delta_y_hist, y_bins)


        ax_histx.set_xlim(np.min(x)-.75, np.max(x)+.75)
        ax_histy.set_ylim(np.min(y)-.75, np.max(y)+.75)
        ax.set_xlim(np.min(x)-.75, np.max(x)+.75)
        ax.set_ylim(np.min(y)-.75, np.max(y)+.75)



# liq_droplet_tau = [[] for ]
# liq_droplet_size = list(np.zeros((4,4,4,3,3,3)))

for ldens in ldens_list:
        o = ldens_list.index(ldens)

        liq_tau_array = np.load(f"/home/ppuel/data/{exp2}/tau_data_ldens{ldens}.npy")
        liq_size_array = np.load(f"/home/ppuel/data/{exp2}/size_data_ldens{ldens}.npy")
        

        # liq_MSD_array = np.load(f"/home/ppuel/data/{exp2}/liq_MSD_ldens{ldens}.npy")
        # poly_MSD_array = np.load(f"/home/ppuel/data/{exp2}/poly_MSD_ldens{ldens}.npy")
        # poly_gyr_array = np.load(f"/home/ppuel/data/{exp2}/gyr_data_ldens{ldens}.npy")
        # liq_MSD_null = np.mean([(np.log10(msd) - np.log10(liq_MSD_array[1]))/(np.log10(c+1)) for c, msd in enumerate(list(liq_MSD_array[1:50]))][1:])
        # poly_MSD_null = np.mean([(np.log10(msd) - np.log10(poly_MSD_array[1]))/(np.log10(c+1)) for c, msd in enumerate(list(poly_MSD_array[1:50]))][1:])
        # poly_gyr_null = np.mean(poly_gyr_array)

        
                        
        for jll_valency in jll_valency_list:

                print(jll_valency, end='\n')
                l = jll_valency_list.index(jll_valency)
                for jlp_valency in jlp_valency_list:
                        m = jlp_valency_list.index(jlp_valency)
                        print(jlp_valency, end='\r')

                        liq_droplet_tau = []
                        liq_droplet_size = []
                        color = []

                        for jll in jll_list:
                                i = jll_list.index(jll)
                                for jlp in jlp_list:

                                        j = jlp_list.index(jlp)
                                        for jlpp in jlpp_list:
                                                k = jlpp_list.index(jlpp)                              
                                                
                                                file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLPP/{jlpp}/JLL_VALENCY/{jll_valency}/JPL_VALENCY/{jlp_valency}/JLPP_VALENCY/6/LDENS/{ldens}/N/0/liq_droplets.pickle", "rb")
                                                droplet_dict = pickle.load(file)
                                                file.close()
                                                
                                                for droplet in droplet_dict.values():
                                                        liq_droplet_tau.append(droplet.tau)
                                                        liq_droplet_size.append(np.max(droplet.sizes))
                                                        color.append(3/4/(1/float(jll)+1/float(jlp)+1/float(jlpp))*256)


                        font = fm.FontProperties(weight='bold',
                                                        style='normal', size=12)
                        fontlabel = {"labelsize" : 12}


                        fig, axs = plt.subplot_mosaic([['histx', '.'],
                                                ['scatter', 'histy']],
                                                figsize=(10, 10),
                                                width_ratios=(4, 2), height_ratios=(2, 4),
                                                layout='constrained')
                        
                        fig.suptitle(f"JLL : {jll}, JLP : {jlp}, JLPP : {jlpp},\nJLL_VALENCY : {jll_valency}, JPL_VALENCY : {jlp_valency}, JLPP_VALENCY : 6,\nLDENS : {ldens}")
                        scatter_hist(liq_droplet_tau, liq_droplet_size, liq_tau_array, liq_size_array, color, axs['scatter'], axs['histx'], axs['histy'])
                        
                        axs['scatter'].set_ylabel("Size Droplet (Number of PRC1)", font = font)
                        axs['scatter'].set_xlabel("Time (kMCS)", font = font)
                        
                        axs['scatter'].tick_params(axis="both", **fontlabel)
                        axs['histx'].tick_params(axis="y", **fontlabel)
                        axs['histy'].tick_params(axis="x", **fontlabel)

                        fig.savefig(f"/home/ppuel/data/{exp}/Droplet_JLL_VALENCY_{jll_valency}_JPL_VALENCY_{jlp_valency}_JLPP_VALENCY_6_LDENS_{ldens}.png")
                        plt.close(fig=fig)
                        
                                                # file = h5py.File(os.path.join(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLPP/{jlpp}/JLL_VALENCY/{jll_valency}/JPL_VALENCY/{jlp_valency}/JLPP_VALENCY/6/LDENS/{ldens}/N/0/", "process.h5"),'r')
                                               

                                                # poly_gyr = np.mean(file["polyGyration"][-10:])
                                                # poly_gyr_matrix[j,k,l,m,o] = poly_gyr
                                                # liq_MSD = np.mean([(np.log10(msd) - np.log10(file["liqMSD"][1]))/(np.log10(c+1)) for c, msd in enumerate(list(file["liqMSD"][1:50]))][1:])
                                                # liq_MSD_matrix[j,k,l,m,o] = liq_MSD
                                                # poly_MSD = np.mean([(np.log10(msd) - np.log10(file["polyHetMSD"][1]))/(np.log10(c+1)) for c, msd in enumerate(list(file["polyHetMSD"][1:50]))][1:])
                                                # poly_MSD_matrix[j,k,l,m,o] = poly_MSD
                                                # if abs(poly_gyr-poly_gyr_null*0.8)/poly_gyr_null*0.8 < 0.1 and abs(liq_MSD-liq_MSD_null)/liq_MSD_null < 0.1 and abs(poly_MSD-poly_MSD_null)/poly_MSD_null < 0.1:
                                                #         print(jll,jlp, jlpp, jll_valency, jlp_valency, ldens)
                                                # file.close()

# np.save(f"/Xnfs/physbiochrom/ppuel/data/{exp}/poly_gyr_matrix.npy", poly_gyr_matrix.flatten())
# np.save(f"/Xnfs/physbiochrom/ppuel/data/{exp}/poly_MSD_matrix.npy", poly_MSD_matrix.flatten())
# np.save(f"/Xnfs/physbiochrom/ppuel/data/{exp}/liq_MSD_matrix.npy", liq_MSD_matrix.flatten())

