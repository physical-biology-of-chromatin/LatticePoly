import os
import sys
import subprocess
import pickle
import h5py
import itertools
from tqdm import tqdm
from Liq_Droplet import Droplet, Event, SDroplet
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
plt.style.use('./resources/h5py/presentation.mplstyle')
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import utils

import matplotlib.patches as mpatches


plasma = mpl.colormaps["plasma"].resampled(101)

exp = 60
exp_name = utils.find_exp(exp, "/Xnfs/physbiochrom/ppuel/data/")

figDir = f"/home/ppuel/figure/{exp_name}/"

os.makedirs(figDir, exist_ok=True)

exp_dir = os.path.join("/Xnfs/physbiochrom/ppuel/data/", exp_name)
output_path = f"/home/ppuel/data/{exp_name}/"
output_path_figure = f"/home/ppuel/data/{exp_name}/figure/"
os.makedirs(output_path, exist_ok=True)
os.makedirs(output_path_figure, exist_ok=True)

dict_parameters, is_poly, meta_parameter_N, meta_parameter_Nmeas = utils.exp_mapping(exp_dir)
path_list = utils.exp_path_list(dict_parameters, -1, exp_dir)



def FIG6_persistence_score_for_droplets():
    # calculate average residency time for a particle inside a droplet
    # |_ hard to do let's do mean persistance over size for each time frame for each droplet

    os.makedirs(os.path.join(figDir, f"droplet_analysis/FIG6_persistence_score/"), exist_ok=True)

    jll_valency_length = 12
    ninter_length = 4
    
    fig = plt.figure(layout='constrained', figsize = ((jll_valency_length+0.5)*4, ninter_length*4))

    ax = fig.subplots(ninter_length, jll_valency_length)

    for ninter_id, ninter in tqdm(enumerate([f"{int(10**(k+1))}" for k in range(ninter_length)]), total = ninter_length, leave = False):
     
        # i = 3 - i - 1
        
        jll_length = 12
        ldens_length = 12
        
        data = np.zeros((jll_valency_length, jll_length, ldens_length))

        for jll_valency_id, jll_valency in tqdm(enumerate([f"{k+1}" for k in range(jll_valency_length)]), total = jll_valency_length, leave = False):
            
            for jll_id, jll in tqdm(enumerate([f"{0.2*k+0.2:.1f}" for k in range(jll_length)]), total = jll_length, leave = False):
                for ldens_id, ldens in tqdm(enumerate([f"{0.003*k+0.001:.3f}" for k in range(ldens_length)]), total = ldens_length, leave = False):
            
                    input_dir = os.path.join(exp_dir, f"JLL/{jll}/JLL_VALENCY/{jll_valency}/LDENS/{ldens}/NINTER/{ninter}/N/0/liq_droplets.pickle")
                    
                    with open(input_dir, "rb") as pfile:
                        droplet_dict = pickle.load(pfile)
                
                    persistance_score_per_droplets = np.zeros(len(droplet_dict))
                    
                    tau_per_droplets = np.array(list(map(lambda x : x.tau, droplet_dict.values())))
                    # max_size_per_droplets = np.array(list(map(lambda x :  max(x.sizes), droplet_dict.values())))
                    
                    for droplet_id, droplet in droplet_dict.items():
                        if droplet.tau > 1 :
                            persistance_array = np.zeros(droplet.tau-1)
                            
                            for frame in range(droplet.tau - 1):
                                persistance_array[frame] = 2*droplet.out_event[frame][0].strength/(droplet.sizes[frame]+droplet.sizes[frame+1])
                    
                            persistance_score_per_droplets[droplet_id] = np.mean(persistance_array)
                
                        else:
                            persistance_score_per_droplets[droplet_id] = 0
                
                    data[jll_valency_id][ldens_id][jll_id] = np.average(persistance_score_per_droplets, weights = tau_per_droplets)
                
            
        for jll_valency_id, jll_valency in enumerate([f"{k+1}" for k in range(jll_valency_length)]):
            data_min = np.min(data)
            data_max = np.max(data)
            im = ax[ninter_id][jll_valency_id].imshow(data[jll_valency_id], cmap='plasma', vmin=data_min, vmax=data_max, origin='lower')

            ax[ninter_id][jll_valency_id].set_title(f"valency = {jll_valency}\nMCS = 10^{np.log10(int(ninter))+2}")
            ax[ninter_id][jll_valency_id].set_xticks(np.arange(jll_length), [f"{0.2*k+0.2:.1f}" for k in range(jll_length)], rotation = 90)
            ax[ninter_id][jll_valency_id].set_yticks(np.arange(ldens_length), [f"{0.003*k+0.001:.3f}" for k in range(ldens_length)])
            ax[ninter_id][jll_valency_id].set_xlabel("Self interaction")
            ax[ninter_id][jll_valency_id].set_ylabel("Particle density")
            
        fig.colorbar(im, ax = ax[ninter_id,:], orientation='vertical', label = "persistance score")
    
    fig.savefig(
        os.path.join(
            figDir,
            f"droplet_analysis/FIG6_persistence_score/FIG6.4_persistence_score_map.png",
        )       
    )

    plt.close(fig = fig)  
    
def print_droplet_proportion_events(ax, droplet_dict):
    
    start_event_dict = {'emergence' : 0, 'split' : 0, 'offload' : 0, 'splinter' : 0}
    stop_event_dict = {'amalgamate' : 0, 'blend' : 0, 'merge' : 0, 'evaporation' : 0 } 

    for droplet_id, droplet in droplet_dict.items():
        verif = False
        for event in droplet.in_event[0]:
            if event.name in start_event_dict.keys():
                start_event_dict[event.name] += 1
                if verif:
                    raise Exception()
                verif = True
            
        if droplet.frame_start + droplet.tau != 100:
            verif = False
            for event in droplet.out_event[-1]:
                if event.name in stop_event_dict.keys():
                    stop_event_dict[event.name] += 1
                    if verif:
                        raise Exception()
                    verif = True

    sum_start_event = sum(start_event_dict.values())
    sum_stop_event = sum(stop_event_dict.values())

    data = [count/sum_start_event for count in start_event_dict.values()] + [count/sum_stop_event for count in stop_event_dict.values()]
    label = [name for name in start_event_dict.keys()] + [name for name in stop_event_dict.keys()]
    shift = [0.8,0.4,-0.4,-0.8,-0.8,-0.4,0.4,0.8]
    color = ['darkred','firebrick','orangered','coral','deepskyblue', 'royalblue', 'slateblue', 'darkblue']
    pie = ax.pie(data, colors = color, wedgeprops=dict(width=0.5), startangle=90)
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
            bbox=bbox_props, zorder=0, va="center")
    
    for i, p in enumerate(pie.wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        # connectionstyle = f"angle,angleA=0,angleB={ang}"
        # kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(label[i], color = color[i], xy=(x, y), xytext=(1.35*np.sign(x), shift[i]),
                    horizontalalignment=horizontalalignment, **kw)

def FIG5_proportion_starting_stopping_events_as_valency_ninter():
    
    os.makedirs(os.path.join(figDir, f"droplet_analysis/FIG5_example/"), exist_ok=True)

    fig = plt.figure(layout='constrained', figsize=(13, 13))
    subfigs = fig.subfigures(2, 2, wspace=0.07, hspace=0.07, width_ratios=[1, 12], height_ratios=[12, 1])

    axs = subfigs[0, 1].subplots(3, 3)
    
    for i, jll_valency in tqdm(enumerate([f"{k+2}" for k in range(3)]), total = 3):
        i = 3 - i - 1
        for j, ninter in tqdm(enumerate([f"{int(10**(k+1))}" for k in range(3)]), total = 3, leave = False):
            input_dir = os.path.join(exp_dir, f"JLL/1.6/JLL_VALENCY/{jll_valency}/LDENS/0.016/NINTER/{ninter}/N/0/liq_droplets.pickle")
    
            with open(input_dir, "rb") as pfile:
                droplet_dict = pickle.load(pfile)

            print_droplet_proportion_events(axs[i][j], droplet_dict)

            axs[i][j].set_title(f"Valency = {jll_valency};\nMonte Carlo step = {int(ninter)*100:,}")

    fs = 20
    
    arrow = mpatches.FancyArrowPatch((0.3, 0), (0.3, 1),
                                     mutation_scale=50)
    ax = subfigs[0, 0].subplots(1,1)
    ax.add_patch(arrow)
    ax.text(0.6, 0.5, 'The ratio of offload/blend and of splinter/amalgamate over\nemergence/evaporation and split/merge increases as valency increases ', rotation = 90, 
        horizontalalignment='left',
        verticalalignment='center',
        fontsize = fs)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.spines[["left", "top", "right", "bottom"]].set_visible(False)

    arrow = mpatches.FancyArrowPatch((0, 0.3), (1, 0.3),
                                     mutation_scale=50)
    ax = subfigs[1, 1].subplots(1,1)
    ax.add_patch(arrow)
    ax.text(0.6, 0.5, 'The ratio of offload/blend and of splinter/amalgamate over\nemergence/evaporation and split/merge increases as MCS increases ', 
        horizontalalignment='center',
        verticalalignment='bottom',
        fontsize = fs)
    
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.spines[["left", "top", "right", "bottom"]].set_visible(False)



    fig.savefig(
        os.path.join(
            figDir,
            f"droplet_analysis/FIG5_example/FIG5_proportion_starting_stopping_events_as_valency_ninter.png",
        )       
    )

    plt.close(fig = fig)    

def FIG4_proportion_starting_stopping_events_as_jll_ldens():
    
    os.makedirs(os.path.join(figDir, f"droplet_analysis/FIG4_example/"), exist_ok=True)

    fig = plt.figure(layout='constrained', figsize=(13, 13))
    subfigs = fig.subfigures(2, 2, wspace=0.07, hspace=0.07, width_ratios=[1, 12], height_ratios=[12, 1])

    axs = subfigs[0, 1].subplots(3, 3)
    
    for i, jll in tqdm(enumerate([f"{k*0.2+1.4:0.1f}" for k in range(3)]), total = 3):
        i = 3 - i - 1
        for j, ldens in tqdm(enumerate([f"{k*0.003+0.013:0.3f}" for k in range(3)]), total = 3, leave = False):
            input_dir = os.path.join(exp_dir, f"JLL/{jll}/JLL_VALENCY/3/LDENS/{ldens}/NINTER/10/N/0/liq_droplets.pickle")
    
            with open(input_dir, "rb") as pfile:
                droplet_dict = pickle.load(pfile)

            print_droplet_proportion_events(axs[i][j], droplet_dict)

            axs[i][j].set_title(f"Self Interaction = {jll}k$_B$T;\nParticle density = {ldens}")

    fs = 20
    
    arrow = mpatches.FancyArrowPatch((0.3, 0), (0.3, 1),
                                     mutation_scale=50)
    ax = subfigs[0, 0].subplots(1,1)
    ax.add_patch(arrow)
    ax.text(0.6, 0.5, 'The ratio of merge/split over emergence/evaporation\n increases as self interaction increases ', rotation = 90, 
        horizontalalignment='left',
        verticalalignment='center',
        fontsize = fs)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.spines[["left", "top", "right", "bottom"]].set_visible(False)

    arrow = mpatches.FancyArrowPatch((0, 0.3), (1, 0.3),
                                     mutation_scale=50)
    ax = subfigs[1, 1].subplots(1,1)
    ax.add_patch(arrow)
    ax.text(0.5, 0.6, 'The ratio of merge/split over emergence/evaporation\n increases as particle density increases',
        horizontalalignment='center',
        verticalalignment='bottom',
        fontsize = fs)
    
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.spines[["left", "top", "right", "bottom"]].set_visible(False)



    fig.savefig(
        os.path.join(
            figDir,
            f"droplet_analysis/FIG4_example/FIG4_proportion_starting_stopping_events.png",
        )       
    )

    plt.close(fig = fig)    
    
def FIG3_Etude_de_cas():

    os.makedirs(os.path.join(figDir, f"droplet_analysis/FIG3_example/"), exist_ok=True)
    
    
    input_dir = os.path.join(exp_dir, "JLL/1.6/JLL_VALENCY/12/LDENS/0.016/NINTER/100/N/0/liq_droplets.pickle")
    
    with open(input_dir, "rb") as pfile:
        dropletDict = pickle.load(pfile)

    start_event_dict = {'emergence' : 0, 'split' : 0, 'offload' : 0, 'splinter' : 0}
    stop_event_dict = {'amalgamate' : 0, 'blend' : 0, 'merge' : 0, 'evaporation' : 0 } 
    for droplet_id, droplet in dropletDict.items():
        verif = False
        for event in droplet.in_event[0]:
            if event.name in start_event_dict.keys():
                start_event_dict[event.name] += 1
                if verif:
                    raise Exception()
                verif = True
            
        if droplet.frame_start + droplet.tau != 100:
            verif = False
            for event in droplet.out_event[-1]:
                if event.name in stop_event_dict.keys():
                    stop_event_dict[event.name] += 1
                    if verif:
                        raise Exception()
                    verif = True

    sum_start_event = sum(start_event_dict.values())
    # print("\n".join([f"{name} : {count/sum_start_event*100:0.1f}%" for name, count in start_event_dict.items()]))
    sum_stop_event = sum(stop_event_dict.values())
    # print("\n".join([f"{name} : {count/sum_stop_event*100:0.1f}%" for name, count in stop_event_dict.items()]))

    fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    
    data = [count/sum_start_event for count in start_event_dict.values()] + [count/sum_stop_event for count in stop_event_dict.values()]
    label = [name for name in start_event_dict.keys()] + [name for name in stop_event_dict.keys()]
    shift = [0.8,0.4,-0.4,-0.8,-0.8,-0.4,0.4,0.8]
    color = ['darkred','firebrick','orangered','coral','deepskyblue', 'royalblue', 'slateblue', 'darkblue']
    pie = ax.pie(data, colors = color, wedgeprops=dict(width=0.5), startangle=90)
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")
    
    for i, p in enumerate(pie.wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        # connectionstyle = f"angle,angleA=0,angleB={ang}"
        # kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(label[i], color = color[i], xy=(x, y), xytext=(1.35*np.sign(x), shift[i]),
                    horizontalalignment=horizontalalignment, **kw)

    fig.savefig(
        os.path.join(
            figDir,
            f"droplet_analysis/FIG3_example/FIG3_pie_chart_start_and_stop.png",
        )       
    )

    plt.close(fig = fig)    
    
    
                # print(f"{droplet_id} {droplet.tau} -> "
                #     + " ".join([f"{event.name}:{event.strength}" for event in droplet.in_event[i]])
                #     + ";  "
                #     + " ".join([f"{event.name}:{event.strength}" for event in droplet.out_event[i-1]]))

def FIG_1_and_2():
    
    os.makedirs(os.path.join(figDir, f"droplet_analysis/FIG2_tau_max/"), exist_ok=True)
    
    
    for fig_id, valency in tqdm(enumerate([f"{i+1}" for i in range(12)]), leave=True, total=12):
    
        fig = plt.figure(figsize=(24, 6))
        ax = fig.subplots(1, 4, sharex = True, sharey = True)
            
        for ax_k, ninter in tqdm(enumerate([f"{10**(i+1):d}" for i in range(4)]), leave=False, total=4):
        
            # for ax_j, jll in tqdm(enumerate([f"{0.2 * i + 0.2:0.1f}" for i in range(12)]), leave=False, total=12):
            #     ax[0][ax_j].text(0.5, 1.3, f"Ep-p = {jll}",
            #     horizontalalignment='center',
            #     verticalalignment='top',
            #     transform=ax[0][ax_j].transAxes)
        
            # for ax_i, ldens in tqdm(enumerate([f"{0.001+i*0.003:.3f}" for i in range(12)]), leave=False, total=12):
            #     ax[ax_i][11].text(1.3, 0.5, f"Density = {ldens}",
            #     horizontalalignment='right',
            #     verticalalignment='center',
            #     rotation='vertical',
            #     transform=ax[ax_i][11].transAxes)
        
            data = np.zeros((12,12))
            
            for row_i, ldens in tqdm(enumerate([f"{0.001+i*0.003:.3f}" for i in range(12)]), leave=False, total=12):
                for col_j, jll in tqdm(enumerate([f"{0.2 * i + 0.2:0.1f}" for i in range(12)]), leave=False, total=12):
                    with h5py.File(
                        os.path.join(exp_dir,f"JLL/{jll}/JLL_VALENCY/{valency}/LDENS/{ldens}/NINTER/{ninter}/aggregated_process.h5"),
                        "r",
                    ) as hfile:
                        data[row_i][col_j] = np.log2(max(hfile["list_tau"])) if len(hfile["list_tau"]) > 0 else -1
        
            ax[ax_k].imshow(data, cmap=plasma, origin="lower", vmin=-1, vmax=np.log2(100))
        
            ax[ax_k].set_xlabel("Self interaction")
            ax[ax_k].set_ylabel("Particle density")
        
            fig.savefig(
                os.path.join(
                    figDir,
                    f"droplet_analysis/FIG2_tau_max/FIG2.{fig_id}_valency={valency}.png",
                )
            )
            plt.close(fig)
        
        # ax[1][ax_index].scatter(
        #     np.array(hfile["list_local_time"]),
        #     np.array(hfile["list_CoM_displacement"]),
        #     c=[plasma(tau) for tau in np.array(hfile["list_tau"])],
        # )
        # ax[2][ax_index].scatter(
        #     np.array(hfile["list_CoM_displacement"]),
        #     np.array(hfile["list_size"]),
        #     c=[plasma(tau) for tau in np.array(hfile["list_tau"])],
        # )
        # ax[1][ax_index].set_xlabel("list_local_time")
        # ax[1][ax_index].set_ylabel("list_CoM_displacement")
        # ax[2][ax_index].set_xlabel("list_CoM_displacement")
        # ax[2][ax_index].set_ylabel("list_size")
        
        
        
        # for path in tqdm(path_list[1:], total = len(path_list)):
        #     list_file = os.listdir(os.path.join(path, "N/"))
        #     os.makedirs(os.path.join(path, "N/0/"), exist_ok=True)
        #     for file in list_file:
        #         os.rename(os.path.join(path, "N/", file), os.path.join(path, "N/0/", file))
        
        
        # def is_non_zero_file(fpath):
        #     return os.path.isfile(fpath) and os.path.getsize(fpath) > 0
        
        
        # for i in range(576):
        #     if is_non_zero_file(os.path.join(exp_dir,f"/tmp/SimProcess_14529469_{i}.err")):
        #         print(i)
        
        # file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp_name}/JLL/{jll}/JLL_VALENCY/{jll_valency}/LDENS/{ldens}/N/{n}/liq_droplets.pickle", "rb")
        # droplet_dict = pickle.load(file)
        # file.close()
