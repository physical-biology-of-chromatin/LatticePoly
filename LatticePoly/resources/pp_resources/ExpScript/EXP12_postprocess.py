# import math as m
import os
# import pickle
# from itertools import chain, product  # , zip_longest
from EXP60_postprocess import print_droplet_proportion_events
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from Liq_Droplet import Droplet, Event
from tqdm import tqdm
from utils import exp_mapping, exp_path_list, find_exp
from sklearn.decomposition import PCA
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import StandardScaler, QuantileTransformer

from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import pickle

plt.style.use("./resources/h5py/presentation.mplstyle")

exp_num = 12
Xnfs_dir = "/Xnfs/physbiochrom/ppuel/data/"
exp_name = find_exp(exp_num, Xnfs_dir)
exp_dir = os.path.join(Xnfs_dir, exp_name)
fig_dir = f"/home/ppuel/figure/{exp_name}/"

scaler_dir = os.path.join(exp_dir, "scaler")
pca_dir = os.path.join(exp_dir, "pca")

os.makedirs(scaler_dir, exist_ok=True)
os.makedirs(pca_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

dict_parameters, _, meta_parameter_Nstat, meta_parameter_Nframe = exp_mapping(exp_dir)

path_list = exp_path_list(dict_parameters, "-1", exp_dir)

def exploration_persistance_over_jll():
    for jll in [f"{1.0+0.2*i:.1f}" for i in range(5)]:
        input_dir = os.path.join(exp_dir, f"JLL/{jll}/JLL_VALENCY/12/LDENS/0.016/N/0/liq_droplets.pickle")

        with open(input_dir, "rb") as pfile:
            droplet_dict = pickle.load(pfile)
            print(len(droplet_dict))
            droplet_max_list = sorted(list(droplet_dict.values()), key = lambda x: x.tau, reverse=True)[:5]

        with open(os.path.join(fig_dir, "example_droplet_persistance.txt"), 'a') as txt_file:    
            for droplet in droplet_max_list:
                txt_file.write(f"{jll} : " + "".join([f"{droplet.sizes[i]: 4d}--{sum(list(map(lambda x: (x.name == 'persistance')*x.strength, droplet.out_event[i]))): 4d}-->" for i in range(droplet.tau - 1)])+f"{droplet.sizes[droplet.tau - 1]: 4d}"+'\n')



def ACP_on_SI_droplets(jll, ldens, scaler = None, pca = None, data_SI = np.empty((0,))):
    dataset_list = ["list_droplet_mean_size", "list_droplet_tau", "list_droplet_mean_CoM_displacement", "list_droplet_mean_number_of_neighbors"]
    key_list = ["Size", "Tau", "CoM Displacement", "Local density"]
    
    valency_length = 11
    valency_list = [f"{i+2}" for i in range(valency_length)]
    valency = 4
    
    with h5py.File(os.path.join(exp_dir, f"JLL/{jll}/JLL_VALENCY/{valency}/LDENS/{ldens}/aggregated_process.h5"), 'r') as hfile:
        data = np.zeros((len(hfile[dataset_list[0]]), 4)  if len(hfile[dataset_list[0]]) < 1e3 else (int(1e3), 4))
        
        for i in range(4):
            data[:,i] = np.array(hfile[dataset_list[i]])[:int(1e3)] #[:int(1e3)]

    # if jll == '0.0':
    #     for jll_valency in tqdm(valency_list, total = valency_length, leave = False):
    #         with h5py.File(os.path.join(exp_dir, f"JLL/{jll}/JLL_VALENCY/{jll_valency}/LDENS/{ldens}/aggregated_process.h5"), 'r') as hfile:
    #             data_tmp = np.zeros((len(hfile[dataset_list[0]]), 4))
                
    #             for i in range(4):
    #                 data_tmp[:,i] = np.array(hfile[dataset_list[i]])
    
    #         data = np.concatenate((data, data_tmp))
            

    if scaler == None:
        scaler = QuantileTransformer(output_distribution = 'normal')
        data = scaler.fit_transform(data)
        if len(data_SI) == 0:
            with open(os.path.join(scaler_dir, f'QuantileTransformer_ldens_{ldens}.bin'),'wb') as f:
                pickle.dump(scaler, f)
        
    else:
        data = scaler.transform(data)
        
    # print(np.mean(data, axis = 0),np.std(data, axis = 0))
    # print(scaler.mean_, scaler.var_)

    # df = pd.DataFrame(data, columns = key_list)

    if pca == None:
        pca = PCA(n_components=4)
        data_pca = pca.fit_transform(data)
        df = pd.DataFrame(data_pca, columns = [f"composante {i+1}" for i in range(4)])
        
        if len(data_SI) == 0:
            with open(os.path.join(pca_dir, f'pca_ldens_{ldens}.bin'),'wb') as f:
                pickle.dump(pca, f)
            
    else:
        if len(data_SI) == 0:
            data_pca = pca.transform(data)
            df = pd.DataFrame(data_pca, columns = [f"composante {i+1}" for i in range(4)])
            
        else:
            data_pca = np.concatenate((pca.transform(data), data_SI), axis = 0) #
            data_type = np.column_stack((data_pca, np.array([0 for i in range(int(1e3))] + [1 for i in range(int(1e3))])))
            
            df = pd.DataFrame(data_type, columns = [f"composante {i+1}" for i in range(4)] + ["type"])
            # df_type = pd.DataFrame({'type' : ['AI' for i in range(len(data)-len(data_SI))] + ['SI' for i in range(len(data_SI))]})
            # print(df_type.head())
            
    # if df_reduced_SI != None:
    #     df_reduced.join(pd.DataFrame(['AI' for i in range(len(df_reduced))], columns = 'type'))
    #     df_reduced_SI.join(pd.DataFrame(['SI' for i in range(len(df_reduced))], columns = 'type'))


    if len(data_SI) == 0:
        pp = sns.pairplot(df, diag_kind='kde')
    else:
        # print(df_reduced.head())
        
        # df_concat = pd.concat([df, df_type], axis=1)
        # print(df_reduced.head())
        # print(df_reduced.shape)
        
        pp = sns.pairplot(df, diag_kind='kde', vars=[f"composante {i+1}" for i in range(4)], hue= 'type', palette = {0 : 'orange', 1 : 'blue'})
    
    for i in range(4):
        pp.figure.axes[i].set_title(f"var_ratio = {pca.explained_variance_ratio_[i]*100:0.1f}%\n{'\n'.join(list(map(lambda x: f'{x:0.3f}' , pca.components_[i])))}")
    pp.figure.suptitle(f'Particle density : {ldens}')
    pp.add_legend()
    pp.figure.set_layout_engine('constrained')
    pp.figure.savefig(os.path.join(fig_dir, f"droplet_analysis/quantile_transform_uniform_acp_valency_1_on_data_jll_{jll}_ldens_{ldens}.png"))

    return(scaler, pca, data_pca)
    
jll = "0.0"
# for ldens in tqdm([f"{0.001+i*0.003:0.3f}" for i in range(12)], total = 12):
#     scaler, pca, data_SI = ACP_on_SI_droplets(jll, ldens)

print("----SI-done------")

ldens = "0.016"
print(jll, ldens)
scaler, pca, data_SI = ACP_on_SI_droplets(jll, ldens)

for jll in tqdm([f"{0.2*(i+1):0.1f}" for i in range(11)], total = 11):
    
    print(jll, ldens)
    ACP_on_SI_droplets(jll, ldens, scaler, pca, data_SI)

# print(pca.feature_names_in_)
# print(pca.n_samples_) 
# print(pca.n_components_) 
# print(pca.noise_variance_) 
# print(pca.explained_variance_ratio_) 
# print(pca.components_)

# pp = sns.pairplot(df, diag_kind='kde')
# for i in range(3):
#     pp.figure.axes[i].set_title(f"var_ratio = {pca.explained_variance_ratio_[i]*100:0.1f}%\n{'\n'.join(list(map(lambda x: f'{x:0.3f}' , pca.components_[i])))}")
# pp.figure.suptitle(f'Particle density : {ldens}')
    
# pp.figure.set_layout_engine('constrained')
# pp.figure.savefig(os.path.join(fig_dir, f"droplet_analysis/acp_ldens_{ldens}.png"))

# fig, ax = plt.subplots(subplot_kw = {"projection" : "3d"})
# ax.scatter(X_reduced[:,0], X_reduced[:,1], X_reduced[:,2])


# n_std = 5

# cov = np.cov(data.T)
# pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
# ell_radius_x = np.sqrt(1 + pearson)
# ell_radius_y = np.sqrt(1 - pearson)
# ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
#                     facecolor='none', edgecolor='red')

# Calculating the standard deviation of x from
# the squareroot of the variance and multiplying
# with the given number of standard deviations.
# scale_x = np.sqrt(cov[0, 0]) * n_std
# mean_x = np.mean(data[:, 0])

# calculating the standard deviation of y ...
# scale_y = np.sqrt(cov[1, 1]) * n_std
# mean_y = np.mean(data[:, 1])

# transf = transforms.Affine2D() \
#     .rotate_deg(45) \
#     .scale(scale_x, scale_y) \
#     .translate(mean_x, mean_y)

# ellipse.set_transform(transf + ax.transData)
# ax.add_patch(ellipse)
# fig.savefig(os.path.join(fig_dir, f"droplet_analysis/transform_ellipse.png"))


def FIG3_proportion_starting_stopping_events():
    
    os.makedirs(os.path.join(fig_dir, f"droplet_analysis/FIG3_events_SI/"), exist_ok=True)

    fig = plt.figure(layout='constrained', figsize=(12*4, 3*4))
    subfigs = fig.subfigures(2, 2, wspace=0.07, hspace=0.07, width_ratios=[1, 12], height_ratios=[12, 1])

    axs = subfigs[0, 1].subplots(3, 12)
    
    for i, jll in tqdm(enumerate([f"{k*0.2:0.1f}" for k in range(3)]), total = 3):
        i = 3 - i - 1
        for j, ldens in tqdm(enumerate([f"{k*0.003+0.001:0.3f}" for k in range(12)]), total = 12, leave = False):
            input_dir = os.path.join(exp_dir, f"JLL/{jll}/JLL_VALENCY/1/LDENS/{ldens}/N/0/liq_droplets.pickle")
    
            with open(input_dir, "rb") as pfile:
                droplet_dict = pickle.load(pfile)

            print_droplet_proportion_events(axs[i][j], droplet_dict)

            axs[i][j].set_title(f"Self Interaction = {jll}kBT;\nParticle density = {ldens}")

    # fs = 20
    
    # arrow = mpatches.FancyArrowPatch((0.3, 0), (0.3, 1),
    #                                  mutation_scale=50)
    # ax = subfigs[0, 0].subplots(1,1)
    # ax.add_patch(arrow)
    # ax.text(0.6, 0.5, 'The ratio of merge/split over emergence/evaporation\n increases as self interaction increases ', rotation = 90, 
    #     horizontalalignment='left',
    #     verticalalignment='center',
    #     fontsize = fs)
    # ax.yaxis.set_visible(False)
    # ax.xaxis.set_visible(False)
    # ax.spines[["left", "top", "right", "bottom"]].set_visible(False)

    # arrow = mpatches.FancyArrowPatch((0, 0.3), (1, 0.3),
    #                                  mutation_scale=50)
    # ax = subfigs[1, 1].subplots(1,1)
    # ax.add_patch(arrow)
    # ax.text(0.5, 0.6, 'The ratio of merge/split over emergence/evaporation\n increases as particle density increases',
    #     horizontalalignment='center',
    #     verticalalignment='bottom',
    #     fontsize = fs)
    
    # ax.yaxis.set_visible(False)
    # ax.xaxis.set_visible(False)
    # ax.spines[["left", "top", "right", "bottom"]].set_visible(False)



    fig.savefig(
        os.path.join(
            fig_dir,
            f"droplet_analysis/FIG3_events_SI/FIG3_proportion_starting_stopping_events.png",
        )       
    )

    plt.close(fig = fig)    


def FIG2_persistence_score_for_droplets_SI():
    # calculate average residency time for a particle inside a droplet
    # |_ hard to do let's do mean persistance over size for each time frame for each droplet

    os.makedirs(os.path.join(fig_dir, f"droplet_analysis/FIG2_persistence_score_SI/"), exist_ok=True)

    
    fig = plt.figure(layout='constrained', figsize = (12,12))

    ax = fig.subplots()

    jll_valency_length = 12
    ldens_length = 12
    
    data = np.zeros((jll_valency_length, ldens_length))

    for jll_valency_id, jll_valency in tqdm(enumerate([f"{k+1}" for k in range(jll_valency_length)]), total = jll_valency_length, leave = False):
        
        for ldens_id, ldens in tqdm(enumerate([f"{0.003*k+0.001:.3f}" for k in range(ldens_length)]), total = ldens_length, leave = False):
            
            input_dir = os.path.join(exp_dir, f"JLL/{0.0}/JLL_VALENCY/{jll_valency}/LDENS/{ldens}/N/0/liq_droplets.pickle")
            
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
        
            data[ldens_id][jll_valency_id] = np.average(persistance_score_per_droplets, weights = tau_per_droplets)
            
        
    data_min = np.min(data)
    data_max = np.max(data)
    im = ax.imshow(data, cmap='plasma', vmin=data_min, vmax=data_max, origin='lower')

    ax.set_title(f"jll = {0.0}k$_B$T")
    ax.set_xticks(np.arange(jll_valency_length), [f"{k+1}" for k in range(jll_valency_length)])
    ax.set_yticks(np.arange(ldens_length), [f"{0.003*k+0.001:.3f}" for k in range(ldens_length)])
    ax.set_xlabel("Valency")
    ax.set_ylabel("Particle density")
    
    fig.colorbar(im, ax = ax, orientation='vertical', label = "persistance score")
    
    fig.savefig(
        os.path.join(
            fig_dir,
            f"droplet_analysis/FIG2_persistence_score_SI/FIG2_persistence_score_SI.png",
        )       
    )

    plt.close(fig = fig)  


def FIG1_persistence_score_for_droplets():
    # calculate average residency time for a particle inside a droplet
    # |_ hard to do let's do mean persistance over size for each time frame for each droplet

    os.makedirs(os.path.join(fig_dir, f"droplet_analysis/FIG1_persistence_score/"), exist_ok=True)

    jll_valency_length = 12
    
    fig = plt.figure(layout='constrained', figsize = (jll_valency_length*4, 1*4))

    ax = fig.subplots(1, jll_valency_length)

    jll_length = 12
    ldens_length = 12
    
    data = np.zeros((jll_valency_length, jll_length, ldens_length))

    for jll_valency_id, jll_valency in tqdm(enumerate([f"{k+1}" for k in range(jll_valency_length)]), total = jll_valency_length, leave = False):
        
        for jll_id, jll in tqdm(enumerate([f"{0.2*k:.1f}" for k in range(jll_length)]), total = jll_length, leave = False):
            for ldens_id, ldens in tqdm(enumerate([f"{0.003*k+0.001:.3f}" for k in range(ldens_length)]), total = ldens_length, leave = False):
        
                input_dir = os.path.join(exp_dir, f"JLL/{jll}/JLL_VALENCY/{jll_valency}/LDENS/{ldens}/N/0/liq_droplets.pickle")
                
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
        im = ax[jll_valency_id].imshow(data[jll_valency_id], cmap='plasma', vmin=data_min, vmax=data_max, origin='lower')

        ax[jll_valency_id].set_title(f"valency = {jll_valency}")
        ax[jll_valency_id].set_xticks(np.arange(jll_length), [f"{0.2*k+0.2:.1f}" for k in range(jll_length)], rotation = 90)
        ax[jll_valency_id].set_yticks(np.arange(ldens_length), [f"{0.003*k+0.001:.3f}" for k in range(ldens_length)])
        ax[jll_valency_id].set_xlabel("Self interaction")
        ax[jll_valency_id].set_ylabel("Particle density")
        
    fig.colorbar(im, ax = ax[:], orientation='vertical', label = "persistance score")

    fig.savefig(
        os.path.join(
            fig_dir,
            f"droplet_analysis/FIG1_persistence_score/FIG1_persistence_score_map.png",
        )       
    )

    plt.close(fig = fig)  
    

def distance(list_pts):
    list_pts_dist = []
    for i in range(len(list_pts) - 1):
        dist = 0
        for j in range(3):
            delta = abs(list_pts[i][j] - list_pts[i + 1][j])
            dist += min(delta, 24 - delta) ** 2
        list_pts_dist.append(float(dist ** (1 / 2)))
    return list_pts_dist


# def distance(list_pts):
#     delta = np.abs(list_pts[:-1] - list_pts[1:])
#     pts_dist = np.sum(np.minimum(delta, 24 - delta) ** 2, axis=1) ** (1 / 2)
#     return pts_dist


plasma = mpl.colormaps["plasma"].resampled(101)


def aggregate_distance(path, ax):

    for n in tqdm(range(10), leave=True, total=10):
        with open(os.path.join(path, f"N/{n}/liq_droplets.pickle"), "rb") as pfile:
            droplet_dict = pickle.load(pfile)

            for droplet in droplet_dict.values():
                if droplet.tau > 2:
                    array = distance(np.array(droplet.center_of_mass))
                    sizes = np.array(
                        [droplet.frame_start + i for i in range(droplet.tau - 1)]
                    )
                    ax.plot(sizes, array, color=plasma(droplet.tau))

                    # print(droplet.tau, distance(np.array(droplet.center_of_mass)))
                    # , c = [plasma(size) for size in droplet.sizes[1:]])
    # except Exception:
    #     print(n)

def testDroplet():
    
    ldens = "0.034"
    
    os.makedirs(os.path.join(fig_dir, f"testDroplet/LDENS_{ldens}_scatter"), exist_ok=True)
    
    for valency in tqdm([f"{i + 3}" for i in range(6)], leave=True, total=6):
        for jll in tqdm([f"{0.2 * i:0.1f}" for i in range(12)], leave=False, total=12):
            fig = plt.figure(figsize=(30, 10))
            ax = fig.subplots(1, 3)
            with h5py.File(
                f"/Xnfs/physbiochrom/ppuel/data/EXP12_liqFraction_starting_of_the_phase_diagram/JLL/{jll}/JLL_VALENCY/{valency}/LDENS/{ldens}/aggregated_process.h5",
                "r",
            ) as hfile:
                ax[0].scatter(
                    np.array(hfile["list_local_time"]),
                    np.array(hfile["list_size"]),
                    c=[plasma(tau) for tau in np.array(hfile["list_tau"])],
                )
                ax[1].scatter(
                    np.array(hfile["list_local_time"]),
                    np.array(hfile["list_CoM_displacement"]),
                    c=[plasma(tau) for tau in np.array(hfile["list_tau"])],
                )
                ax[2].scatter(
                    np.array(hfile["list_CoM_displacement"]),
                    np.array(hfile["list_size"]),
                    c=[plasma(tau) for tau in np.array(hfile["list_tau"])],
                )
                ax[0].set_xlabel("list_local_time")
                ax[0].set_ylabel("list_size")
                ax[0].set_yscale("log", base=2)
                ax[1].set_xlabel("list_local_time")
                ax[1].set_ylabel("list_CoM_displacement")
                ax[2].set_xlabel("list_CoM_displacement")
                ax[2].set_ylabel("list_size")
                ax[0].set_yscale("log", base=2)
    
            fig.savefig(
                os.path.join(
                    fig_dir,
                    f"testDroplet/LDENS_{ldens}_scatter/scatter_JLL_{jll}_JLL_VALENCY_{valency}_LDENS_{ldens}.png",
                )
            )
            plt.close(fig)

# with h5py.File("/Xnfs/physbiochrom/ppuel/data/EXP12_liqFraction_starting_of_the_phase_diagram/JLL/1.6/JLL_VALENCY/12/LDENS/0.034/N/0/process.h5") as hfile:
#      print(hfile["liq_drop_info"][3:5,:5])

# with open("/Xnfs/physbiochrom/ppuel/data/EXP12_liqFraction_starting_of_the_phase_diagram/JLL/1.6/JLL_VALENCY/12/LDENS/0.034/N/0/liq_droplets.pickle", "rb"
# ) as pfile:
#     droplet_dict = pickle.load(pfile)
#     for droplet in droplet_dict.values():
#         if droplet.tau > 10:
#             array = np.array(droplet.center_of_mass)
#             print(array)

# fig.savefig(os.path.join(fig_dir, "testDroplet/test_comparaison_all_tau.png"))


# def is_non_zero_file(fpath):
#     return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


# for i in range(900):
#     if is_non_zero_file(
#         f"/Xnfs/physbiochrom/ppuel/data/EXP12_liqFraction_starting_of_the_phase_diagram/tmp/SimProcess_14489522_{i}.err"
#     ):
#         print(i)


# hist, bin_egde = np.histogram(
#     aggregate_distance("/Xnfs/physbiochrom/ppuel/data/EXP12_liqFraction_starting_of_the_phase_diagram/JLL/1.6/JLL_VALENCY/12/LDENS/0.034/")
# )ax.plot(aggregate_distance("/Xnfs/physbiochrom/ppuel/data/EXP12_liqFraction_starting_of_the_phase_diagram/JLL/1.6/JLL_VALENCY/12/LDENS/0.034/"))
# (bin_egde[1:] + bin_egde[:-1]) / 2, hist)


# hist, bin_egde = np.histogram(
#     list(

#             aggregate_distance(
#                 "/Xnfs/physbiochrom/ppuel/data/EXP12_liqFraction_starting_of_the_phase_diagram/JLL/1.8/JLL_VALENCY/6/LDENS/0.034/"
#             )

#     )
# )
# ax.semilogy((bin_egde[1:] + bin_egde[:-1]) / 2, hist)


# hist, bin_egde = np.histogram(
#     list(
#             aggregate_distance(
#                 "/Xnfs/physbiochrom/ppuel/data/EXP12_liqFraction_starting_of_the_phase_diagram/JLL/2.2/JLL_VALENCY/12/LDENS/0.034/"

#         )
#     )
# )
# ax.semilogy((bin_egde[1:] + bin_egde[:-1]) / 2, hist)


def Analyse() -> None:
    listDatasetName = ["liqMean"]

    plotType = "image"

    fileName = "aggregated_process.h5"

    # with h5py.File(os.path.join(path_list[0], "N/0/process.h5"), 'r') as hfile:
    #         for dataset in hfile.keys():
    #                 print(hfile[dataset].shape, dataset)

    # AggregateData(listDatasetName, path_list, metaParameterN)

    PrintAggregateData(fileName, listDatasetName, plotType, expName, SimTime=False)


def exp_diffusion():

    # DLiq = 1e6 # nm²/s
    N = int(metaParameterN)
    Nmeas_MSD = int(metaParameterNmeas) - 1
    nInter = 10000

    diffusion_CoM_matrix_polyfit_mean = np.zeros(
        (
            len(dict_parameters["JLL"]),
            len(dict_parameters["JLL_VALENCY"]),
            len(dict_parameters["LDENS"]),
        )
    )
    diffusion_CoM_matrix_polyfit_var = np.zeros(
        (
            len(dict_parameters["JLL"]),
            len(dict_parameters["JLL_VALENCY"]),
            len(dict_parameters["LDENS"]),
        )
    )
    diffusion_CoM_matrix_sum_mean = np.zeros(
        (
            len(dict_parameters["JLL"]),
            len(dict_parameters["JLL_VALENCY"]),
            len(dict_parameters["LDENS"]),
        )
    )
    diffusion_CoM_matrix_sum_var = np.zeros(
        (
            len(dict_parameters["JLL"]),
            len(dict_parameters["JLL_VALENCY"]),
            len(dict_parameters["LDENS"]),
        )
    )
    # time_matrix = np.zeros((len(dict_parameters["JLL"]), len(dict_parameters["JLL_VALENCY"]), len(dict_parameters["LDENS"])))

    for ids, (jll, jll_valency, ldens) in tqdm(
        enumerate(
            product(
                dict_parameters["JLL"],
                dict_parameters["JLL_VALENCY"],
                dict_parameters["LDENS"],
            )
        ),
        total=len(dict_parameters["JLL"])
        * len(dict_parameters["JLL_VALENCY"])
        * len(dict_parameters["LDENS"]),
    ):
        miniFig = plt.figure(figsize=(8, 4))
        axs: list[Axes] = miniFig.subplots(1, 2)

        liqMSD_CoM_matrix = np.zeros((N, Nmeas_MSD))
        for n in range(N):
            with h5py.File(
                os.path.join(
                    exp_dir,
                    f"JLL/{jll}/JLL_VALENCY/{jll_valency}/LDENS/{ldens}/N/{n}/process.h5",
                ),
                "r",
            ) as nFile:
                dset: h5py.Dataset = nFile.require_dataset(
                    "liqMSD_CoM", Nmeas_MSD + 1, np.float64
                )
                liqMSD_CoM_matrix[n] = dset[1 : Nmeas_MSD + 1]

        np.save(
            os.path.join(
                fig_dir,
                f"npyData/JLL_{jll}_JLL_VALENCY_{jll_valency}_LDENS_{ldens}_liqMSD_CoM_matrix.npy",
            ),
            liqMSD_CoM_matrix,
        )

        liqMSD_CoM_Mean = np.mean(liqMSD_CoM_matrix, axis=0)
        liqMSD_CoM_Var = np.var(liqMSD_CoM_matrix, axis=0)

        axs[0].loglog(np.arange(1, Nmeas_MSD + 1) * nInter, liqMSD_CoM_Mean, "darkblue")
        axs[0].tick_params(axis="y", labelcolor="darkblue")

        axs[0].set_xlabel("MCS")
        axs[0].set_ylabel("MSD of CoM", color="darkblue")

        axs[1].plot(
            np.arange(1, Nmeas_MSD + 1), np.sqrt(liqMSD_CoM_Var / (12 * N)), "blue"
        )
        axs[1].tick_params(axis="y", labelcolor="blue")

        axs[1].set_xlabel("MCS")
        axs[1].set_ylabel("Variation of MSD of CoM", color="blue")

        miniFig.suptitle(
            f"JLL = {jll}\nJLL_VALENCY = {jll_valency}\nLDENS = {ldens}", color="black"
        )

        miniFig.savefig(
            os.path.join(
                fig_dir,
                f"miniFig/JLL_{jll}_JLL_VALENCY_{jll_valency}_LDENS_{ldens}_MSD_CoM.png",
            )
        )
        plt.close(fig=miniFig)

        poly_CoM = np.polyfit(
            np.arange(1, Nmeas_MSD + 1) * nInter,
            liqMSD_CoM_Mean * 20**2 * 2,
            deg=1,
            w=1 / np.sqrt(liqMSD_CoM_Var),
            cov="unscaled",
        )

        diffusion_CoM_matrix_polyfit_mean[ids // (12 * 12)][(ids // 12) % 12][
            ids % 12
        ] = poly_CoM[0][0]
        diffusion_CoM_matrix_polyfit_var[ids // (12 * 12)][(ids // 12) % 12][
            ids % 12
        ] = np.sqrt(poly_CoM[1][0][0])

        D_CoM_ref_matrix = (
            liqMSD_CoM_matrix / np.arange(1, Nmeas_MSD + 1) / nInter * 20**2 * 2
        )  # nm²/MCS

        diffusion_CoM_matrix_sum_mean[ids // (12 * 12)][(ids // 12) % 12][ids % 12] = (
            np.mean(D_CoM_ref_matrix)
        )
        diffusion_CoM_matrix_sum_var[ids // (12 * 12)][(ids // 12) % 12][ids % 12] = (
            np.var(D_CoM_ref_matrix)
        )

    np.save(
        os.path.join(fig_dir, "diffusion_CoM_matrix_polyfit_mean.npy"),
        diffusion_CoM_matrix_polyfit_mean,
    )
    np.save(
        os.path.join(fig_dir, "diffusion_CoM_matrix_polyfit_var.npy"),
        diffusion_CoM_matrix_polyfit_var,
    )
    np.save(
        os.path.join(fig_dir, "diffusion_CoM_matrix_sum_mean.npy"),
        diffusion_CoM_matrix_sum_mean,
    )
    np.save(
        os.path.join(fig_dir, "diffusion_CoM_matrix_sum_var.npy"),
        diffusion_CoM_matrix_sum_var,
    )

    fig_polyfit = plt.figure(figsize=(4 * 5 + 2, 6 * 5), layout="constrained")
    axs_polyfit: list[list[Axes]] = fig_polyfit.subplots(5, 6)

    fig_sum = plt.figure(figsize=(4 * 5 + 2, 6 * 5), layout="constrained")
    axs_sum: list[list[Axes]] = fig_sum.subplots(5, 6)

    min_diffusion_CoM_matrix_mean, max_diffusion_CoM_matrix_mean = (
        min(
            np.min(diffusion_CoM_matrix_polyfit_mean),
            np.min(diffusion_CoM_matrix_sum_mean),
        ),
        max(
            diffusion_CoM_matrix_polyfit_mean.max(), diffusion_CoM_matrix_sum_mean.max()
        ),
    )
    min_diffusion_CoM_matrix_var, max_diffusion_CoM_matrix_var = (
        min(
            np.min(diffusion_CoM_matrix_polyfit_var),
            np.min(diffusion_CoM_matrix_sum_var),
        ),
        max(diffusion_CoM_matrix_polyfit_var.max(), diffusion_CoM_matrix_sum_var.max()),
    )

    for ids, jll_valency in enumerate(dict_parameters["JLL_VALENCY"]):
        row_ids = ids // 6
        col_ids = ids % 6

        ax_polyfit_mean = axs_polyfit[row_ids][col_ids]
        ax_polyfit_var = axs_polyfit[row_ids + 1][col_ids]

        ax_sum_mean = axs_sum[row_ids][col_ids]
        ax_sum_var = axs_sum[row_ids + 1][col_ids]

        im_poly_mean: AxesImage = ax_polyfit_mean.imshow(
            diffusion_CoM_matrix_polyfit_mean[:, ids, :],
            cmap="Purples",
            vmin=min_diffusion_CoM_matrix_mean,
            vmax=max_diffusion_CoM_matrix_mean,
            origin="lower",
        )
        im_polyfit_var: AxesImage = ax_polyfit_var.imshow(
            diffusion_CoM_matrix_polyfit_var[:, ids, :],
            cmap="Reds",
            vmin=min_diffusion_CoM_matrix_var,
            vmax=max_diffusion_CoM_matrix_var,
            origin="lower",
        )
        im_sum_mean: AxesImage = ax_sum_mean.imshow(
            diffusion_CoM_matrix_sum_mean[:, ids, :],
            cmap="Purples",
            vmin=min_diffusion_CoM_matrix_mean,
            vmax=max_diffusion_CoM_matrix_mean,
            origin="lower",
        )
        im_sum_var: AxesImage = ax_sum_var.imshow(
            diffusion_CoM_matrix_sum_var[:, ids, :],
            cmap="Reds",
            vmin=min_diffusion_CoM_matrix_var,
            vmax=max_diffusion_CoM_matrix_var,
            origin="lower",
        )

        ax_polyfit_mean.set_xlabel("Jll")
        ax_polyfit_mean.set_ylabel("Ldens")
        ax_polyfit_mean.set_xticks(
            np.arange(len(dict_parameters["JLL"])), dict_parameters["JLL"], rotation=90
        )
        ax_polyfit_mean.set_yticks(
            np.arange(len(dict_parameters["LDENS"])), dict_parameters["LDENS"]
        )

        if ids == 0:
            fig_polyfit.colorbar(
                im_poly_mean,
                cax=axs_polyfit[4][1],
                orientation="horizontal",
                label="Diffusion CoM",
            )
            fig_polyfit.colorbar(
                im_polyfit_var,
                cax=axs_polyfit[4][4],
                orientation="horizontal",
                label="Diffusion CoM variation",
            )
            fig_sum.colorbar(
                im_sum_mean,
                cax=axs_sum[4][1],
                orientation="horizontal",
                label="Diffusion CoM",
            )
            fig_sum.colorbar(
                im_sum_var,
                cax=axs_sum[4][4],
                orientation="horizontal",
                label="Diffusion CoM variation",
            )

    fig_polyfit.savefig(os.path.join(fig_dir, "Diffusion_CoM_polyfit_image.png"))
    fig_sum.savefig(os.path.join(fig_dir, "Diffusion_CoM_sum_image.png"))


def ref_diffusion(dict_parameters, metaParameterN, metaParameterNmeas):

    metaParameterN = int(metaParameterN)
    Nmeas_MSD = int(metaParameterNmeas) - 1
    nInter = 10000

    DLiq = 1e6

    fig = plt.figure(figsize=(32, 16), layout="constrained")
    axList = []
    subfiguresList = fig.subfigures(2, 2, height_ratios=(3, 1), width_ratios=(1, 1))

    axList.append(subfiguresList[0][0].subplots(3, 4))
    axList.append(subfiguresList[1][0].subplots(1, 2))
    axList.append(subfiguresList[0][1].subplots(3, 4))
    axList.append(subfiguresList[1][1].subplots(1, 2))

    D_liq_ref_mean = np.zeros(len(dict_parameters["LDENS"]))
    D_CoM_ref_mean = np.zeros(len(dict_parameters["LDENS"]))
    D_liq_ref_var = np.zeros(len(dict_parameters["LDENS"]))
    D_CoM_ref_var = np.zeros(len(dict_parameters["LDENS"]))

    D_liq_ref_polyfit_mean = np.zeros(len(dict_parameters["LDENS"]))
    D_liq_ref_polyfit_var = np.zeros(len(dict_parameters["LDENS"]))
    D_CoM_ref_polyfit_mean = np.zeros(len(dict_parameters["LDENS"]))
    D_CoM_ref_polyfit_var = np.zeros(len(dict_parameters["LDENS"]))

    listFile = os.listdir(fig_dir)

    for i, ldens in tqdm(
        enumerate(dict_parameters["LDENS"]), total=len(dict_parameters["LDENS"])
    ):
        if f"liqMSD_matrix_Ldens_{ldens}.npy" in listFile:
            liqMSD_matrix = np.load(
                os.path.join(fig_dir, f"liqMSD_matrix_Ldens_{ldens}.npy")
            )
        else:
            liqMSD_matrix = np.zeros((12 * metaParameterN, Nmeas_MSD))

        if f"liqMSD_matrix_CoM_Ldens_{ldens}.npy" in listFile:
            liqMSD_CoM_matrix = np.load(
                os.path.join(fig_dir, f"liqMSD_matrix_CoM_Ldens_{ldens}.npy")
            )
        else:
            liqMSD_CoM_matrix = np.zeros((12 * metaParameterN, Nmeas_MSD))

        if (
            f"liqMSD_matrix_Ldens_{ldens}.npy" not in listFile
            or f"liqMSD_matrix_CoM_Ldens_{ldens}.npy" not in listFile
        ):
            for j, jll_valency in tqdm(
                enumerate(dict_parameters["JLL_VALENCY"]),
                total=len(dict_parameters["JLL_VALENCY"]),
                leave=False,
            ):
                for n in range(metaParameterN):
                    with h5py.File(
                        os.path.join(
                            exp_dir,
                            f"JLL/0.0/JLL_VALENCY/{jll_valency}/LDENS/{ldens}/N/{n}/process.h5",
                        ),
                        "r",
                    ) as nFile:
                        if f"liqMSD_matrix_Ldens_{ldens}.npy" not in listFile:
                            liqMSD_matrix[j * metaParameterN + n] = (
                                nFile.require_dataset(
                                    "liqMSD", Nmeas_MSD + 1, np.float64
                                )[1 : Nmeas_MSD + 1]
                            )
                        if f"liqMSD_matrix_CoM_Ldens_{ldens}.npy" not in listFile:
                            liqMSD_CoM_matrix[j * metaParameterN + n] = (
                                nFile.require_dataset(
                                    "liqMSD_CoM", Nmeas_MSD + 1, np.float64
                                )[1 : Nmeas_MSD + 1]
                            )

        if f"liqMSD_matrix_Ldens_{ldens}.npy" not in listFile:
            np.save(
                os.path.join(fig_dir, f"liqMSD_matrix_Ldens_{ldens}.npy"), liqMSD_matrix
            )
        if f"liqMSD_matrix_CoM_Ldens_{ldens}.npy" not in listFile:
            np.save(
                os.path.join(fig_dir, f"liqMSD_matrix_CoM_Ldens_{ldens}.npy"),
                liqMSD_CoM_matrix,
            )

        ax = axList[0][i // 4][i % 4]
        ax_var = axList[2][i // 4][i % 4]

        liqMSD_Mean = np.mean(liqMSD_matrix, axis=0)
        liqMSD_Var = np.var(liqMSD_matrix, axis=0)

        liqMSD_CoM_Mean = np.mean(liqMSD_CoM_matrix, axis=0)
        liqMSD_CoM_Var = np.var(liqMSD_CoM_matrix, axis=0)

        ax.loglog(np.arange(1, Nmeas_MSD + 1), liqMSD_Mean, "darkred")
        ax.tick_params(axis="y", labelcolor="darkred")

        ax_twinx = ax.twinx()
        ax_twinx.loglog(np.arange(1, Nmeas_MSD + 1), liqMSD_CoM_Mean, "darkblue")
        ax_twinx.tick_params(axis="y", labelcolor="darkblue")

        ax.set_xlabel("MCS")
        ax.set_ylabel("MSD of particles", color="darkred")
        ax_twinx.set_ylabel("MSD of CoM", color="darkblue")

        ax_var.plot(
            np.arange(1, Nmeas_MSD + 1),
            np.sqrt(liqMSD_Var / (12 * metaParameterN)),
            "red",
        )
        ax_var.tick_params(axis="y", labelcolor="red")

        ax_var_twinx = ax_var.twinx()
        ax_var_twinx.plot(
            np.arange(1, Nmeas_MSD + 1),
            np.sqrt(liqMSD_CoM_Var / (12 * metaParameterN)),
            "blue",
        )
        ax_var_twinx.tick_params(axis="y", labelcolor="blue")

        ax_var.set_xlabel("MCS")
        ax_var.set_ylabel("Variation of MSD", color="red")
        ax_var_twinx.set_ylabel("Variation of MSD of CoM", color="blue")

        ax.set_title(f"Particles density : {ldens}", color="darkgreen")
        ax_var.set_title(f"Particles density : {ldens}", color="darkgreen")

        poly_liq = np.polyfit(
            np.arange(1, Nmeas_MSD + 1) * nInter,
            liqMSD_Mean * 20**2 * 2,
            deg=1,
            w=1 / np.sqrt(liqMSD_Var),
            cov="unscaled",
        )
        poly_CoM = np.polyfit(
            np.arange(1, Nmeas_MSD + 1) * nInter,
            liqMSD_CoM_Mean * 20**2 * 2,
            deg=1,
            w=1 / np.sqrt(liqMSD_CoM_Var),
            cov="unscaled",
        )

        D_liq_ref_polyfit_mean[i] = poly_liq[0][0]
        D_liq_ref_polyfit_var[i] = np.sqrt(poly_liq[1][0][0])
        D_CoM_ref_polyfit_mean[i] = poly_CoM[0][0]
        D_CoM_ref_polyfit_var[i] = np.sqrt(poly_CoM[1][0][0])

        D_liq_ref_matrix = (
            liqMSD_matrix / np.arange(1, Nmeas_MSD + 1) / nInter * 20**2 * 2
        )  # nm²/MCS
        D_CoM_ref_matrix = (
            liqMSD_CoM_matrix / np.arange(1, Nmeas_MSD + 1) / nInter * 20**2 * 2
        )  # nm²/MCS

        D_liq_ref_mean[i] = np.mean(D_liq_ref_matrix)
        D_CoM_ref_mean[i] = np.mean(D_CoM_ref_matrix)
        D_liq_ref_var[i] = np.var(D_liq_ref_matrix)
        D_CoM_ref_var[i] = np.var(D_CoM_ref_matrix)

    ax = axList[1][0]
    # ax.plot(np.arange(12)*0.003+0.001, D_liq_ref_mean, color = 'sandybrown')
    ax.plot(np.arange(12) * 0.003 + 0.001, D_liq_ref_polyfit_mean, color="lightcoral")
    ax.set_xticks(
        ticks=np.arange(12) * 0.003 + 0.001,
        labels=dict_parameters["LDENS"],
        color="darkgreen",
    )
    ax.tick_params(axis="x", labelrotation=90)
    ax_twinx = ax.twinx()
    ax_twinx.plot(
        np.arange(12) * 0.003 + 0.001, D_liq_ref_polyfit_mean / DLiq * nInter, alpha=0
    )
    ax_twinx.set_ylabel("Simulated Time (s)")

    ax.text(
        0.9,
        0.1,
        "The increase is due to\nmore swapping particles",
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax.transAxes,
    )

    ax.set_ylabel("Particles diffusion", color="lightcoral")
    ax.set_xlabel("Particles density", color="darkgreen")

    ax = axList[1][1]
    # ax.plot(np.arange(12)*0.003+0.001, D_CoM_ref_mean, color = 'lightseagreen')
    ax.plot(np.arange(12) * 0.003 + 0.001, D_CoM_ref_polyfit_mean, color="steelblue")
    ax.set_xticks(
        ticks=np.arange(12) * 0.003 + 0.001,
        labels=dict_parameters["LDENS"],
        color="darkgreen",
    )
    ax.tick_params(axis="x", labelrotation=90)

    ax.set_ylabel("CoM Diffusion", color="steelblue")
    ax.set_xlabel("Particles density", color="darkgreen")

    ax = axList[3][0]
    # ax.plot(np.arange(12)*0.003+0.001, np.sqrt(D_liq_ref_var/(12*metaParameterN*Nmeas_MSD)), color = 'lightcoral')
    ax.plot(np.arange(12) * 0.003 + 0.001, D_liq_ref_polyfit_var, color="sandybrown")
    ax.set_xticks(
        ticks=np.arange(12) * 0.003 + 0.001,
        labels=dict_parameters["LDENS"],
        color="darkgreen",
    )
    ax.tick_params(axis="x", labelrotation=90)

    ax.set_ylabel("Variation of Particles Diffusion", color="sandybrown")
    ax.set_xlabel("Particles density", color="darkgreen")

    ax = axList[3][1]
    # ax.plot(np.arange(12)*0.003+0.001, np.sqrt(D_CoM_ref_var/(12*metaParameterN*Nmeas_MSD)), color = 'steelblue')
    ax.plot(np.arange(12) * 0.003 + 0.001, D_CoM_ref_polyfit_var, color="lightseagreen")
    ax.set_xticks(
        ticks=np.arange(12) * 0.003 + 0.001,
        labels=dict_parameters["LDENS"],
        color="darkgreen",
    )
    ax.tick_params(axis="x", labelrotation=90)

    ax.set_ylabel("Variation of CoM Diffusion", color="lightseagreen")
    ax.set_xlabel("Particles density", color="darkgreen")

    poly_D_liq = np.polyfit(np.arange(12) * 0.003 + 0.001, D_liq_ref_mean, 1)
    poly_D_CoM = np.polyfit(1 / (np.arange(12) * 0.003 + 0.001), D_CoM_ref_mean, 1)

    ax = axList[1][0]
    ax.plot(
        np.arange(12) * 0.003 + 0.001,
        np.polyval(poly_D_liq, np.arange(12) * 0.003 + 0.001),
        "--",
        color="black",
    )
    ax = axList[1][1]
    ax.plot(
        np.arange(12) * 0.003 + 0.001,
        1 / np.polyval(poly_D_CoM, np.arange(12) * 0.003 + 0.001),
        "--",
        color="black",
    )

    print(poly_D_liq)
    print(poly_D_CoM)

    fig.savefig(os.path.join(fig_dir, "DiffusionCoefficientRef_fullMSD.png"))


# ref_diffusion(dict_parameters, metaParameterN, metaParameterNmeas)
