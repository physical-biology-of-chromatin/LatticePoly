import subprocess
import os
from LiqCluster_lifeTime import LifeTime, Droplet, Event
import pickle
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import h5py
import matplotlib.font_manager 

exp = "EXP37_bithorax_null_hypothesis"

os.makedirs(f"/home/ppuel/data/{exp}/droplet_tracking", exist_ok=True)

ldens_list = ["0.0009","0.0018","0.0036"]
N = 160

# ldens_list = ["0.010","0.025"]
# val_list = [str(i*2) for i in range(1,7)]


# jll_list = ["0.0"]
# ldens_list = ["0.010"]
# val_list = ["2"]

# for ldens in ldens_list:

        # liq_droplet_tau = []
        # liq_droplet_size = []

        # poly_gyr = []
        # liq_MSD = np.zeros(101)
        # poly_MSD = np.zeros(101)
        


        # for n in range(N):

                # if n%10 == 0:
                        # print(n, end='\r')

                # file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/LDENS/{ldens}/N/{n}/liq_droplets.pickle", "rb")
                # droplet_dict = pickle.load(file)
                # file.close()

                # for droplet in droplet_dict.values():
                        # liq_droplet_tau.append(droplet.tau)
                        # liq_droplet_size.append(np.max(droplet.sizes))

                # file = h5py.File(os.path.join(f"/Xnfs/physbiochrom/ppuel/data/{exp}/LDENS/{ldens}/N/{n}/", "process.h5"),'r')
                
                # poly_gyr.append(np.mean(file["polyGyration"][-10:]))
                # liq_MSD += file["liqMSD"]
                # poly_MSD += file["polyHetMSD"]
                
                # file.close()

        # np.save(f"/home/ppuel/data/{exp}/tau_data_ldens{ldens}.npy", np.array(liq_droplet_tau))
        # np.save(f"/home/ppuel/data/{exp}/size_data_ldens{ldens}.npy", np.array(liq_droplet_size))
        # np.save(f"/home/ppuel/data/{exp}/gyr_data_ldens{ldens}.npy", np.array(poly_gyr))
        # np.save(f"/home/ppuel/data/{exp}/liq_MSD_ldens{ldens}.npy", liq_MSD)
        # np.save(f"/home/ppuel/data/{exp}/poly_MSD_ldens{ldens}.npy", poly_MSD)

if 0:

        fig = plt.figure(figsize=(24,16))
        i = -1
        for ldens in ldens_list:
                i += 1

                liq_MSD = np.load(f"/home/ppuel/data/{exp}/liq_MSD_ldens{ldens}.npy")
                poly_MSD = np.load(f"/home/ppuel/data/{exp}/poly_MSD_ldens{ldens}.npy")

                ax = fig.add_subplot(1,3,i+1)
                ax.plot(np.log10([j for j in range(101)]), np.log10(liq_MSD))
                ax.plot(np.log10([j for j in range(101)]), np.log10(poly_MSD))

                
        fig.savefig(f"/home/ppuel/data/{exp}/MSD.png")


        fig = plt.figure(figsize=(24,16))
        i = -1
        for ldens in ldens_list:
                i += 1

                poly_gyr = np.load(f"/home/ppuel/data/{exp}/gyr_data_ldens{ldens}.npy")

                ax = fig.add_subplot(1,3,i+1)
                ax.plot(poly_gyr)
                ax.plot([0,len(poly_gyr)], [np.mean(poly_gyr),np.mean(poly_gyr)])
                print(np.mean(poly_gyr))
                
        fig.savefig(f"/home/ppuel/data/{exp}/Gyr.png")

def scatter_hist(x, y, ax, ax_histx, ax_histy):
        # no labels
        ax_histx.tick_params(axis="x", labelbottom=False)
        ax_histy.tick_params(axis="y", labelleft=False)

        # the scatter plot:
        ax.scatter(x, y)

        x_min = np.min(x)
        x_range = np.max(x)-np.min(x)
        x_bins = []
        for i in range(int(x_range)+2):
                x_bins += [i + np.min(x)-.5, i + np.min(x)-.5]
        x_hist = np.zeros(x_range*2+4)
        for x_i in x:
                x_hist[(x_i-x_min)*2+1:(x_i-x_min)*2+3] += 1

        print(x_bins, x_hist)
        x_hist[x_hist < 1] = 1

        y_min = np.min(y)
        y_range = np.max(y)-np.min(y)
        y_bins = []
        for i in range(int(y_range)+2):
                y_bins += [i + np.min(y)-.5, i + np.min(y)-.5]
        y_hist = np.zeros(y_range*2+4)
        for y_i in y:
                y_hist[(y_i-y_min)*2+1:(y_i-y_min)*2+3] += 1

        y_hist[y_hist < 1] = 1

        ax_histx.plot(x_bins, np.log10(x_hist))
        ax_histy.plot(np.log10(y_hist), y_bins)
        ax_histx.set_xlim(np.min(x)-.75, np.max(x)+.75)
        ax_histy.set_ylim(np.min(y)-.75, np.max(y)+.75)
        ax.set_xlim(np.min(x)-.75, np.max(x)+.75)
        ax.set_ylim(np.min(y)-.75, np.max(y)+.75)



for ldens in ldens_list:

        font = matplotlib.font_manager.FontProperties(family= 'Cambria', 
                                                        weight='bold',
                                                        style='normal', size=60)
        fontlabel = {"labelfontfamily" : 'Cambria', "labelsize" : 60}


        size_data = np.load(f"/home/ppuel/data/{exp}/size_data_ldens{ldens}.npy")
        tau_data = np.load(f"/home/ppuel/data/{exp}/tau_data_ldens{ldens}.npy")

        fig, axs = plt.subplot_mosaic([['histx', '.'],
                               ['scatter', 'histy']],
                              figsize=(6, 6),
                              width_ratios=(4, 1), height_ratios=(1, 4),
                              layout='constrained')
        scatter_hist(tau_data, size_data, axs['scatter'], axs['histx'], axs['histy'])

        fig.savefig(f"/home/ppuel/data/{exp}/Droplet_ldens_{ldens}.png")
