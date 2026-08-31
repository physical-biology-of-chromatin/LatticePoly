import subprocess
import os
# from LiqCluster_lifeTime import LifeTime, Droplet, Event
import pickle
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import h5py
import scipy.stats as st
import itertools
import scipy.optimize as op
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
import seaborn as sns


jll_list = ["1.0","1.5","2.0","4.0"]
jlp_list = ["1.0","1.5","2.0","4.0"]
jlpp_list = ["1.0","1.5","2.0","4.0"]


jll_valency_list = ["2","4","6"]
jlp_valency_list = ["2","4","6"]
jlpp_valency_list = ["2","4","6"]


font = fm.FontProperties(weight='bold',
                                style='normal', size=20)
fontlabel = {"labelsize" : 20}


ldens_list = ["0.0009","0.0018","0.0036"]

plasma = mpl.colormaps["plasma"].resampled(258)

exp = "EXP38_bithorax_light_sweep_corrigé"

# # exp2 = "EXP37_bithorax_null_hypothesis"

os.makedirs(f"/home/ppuel/data/{exp}/gyr_fit/", exist_ok=True)


# gyr_init = np.zeros((3,16,9))
# gyr_final = np.zeros((3,16,9))

# for m in range(3):
#         ldens = ldens_list[m]
#         for i in range(4):
#                 print(i)
#                 jll = jll_list[i]
#                 for j in range(4):
#                         jlp = jlp_list[j]
#                         for k in range(3):
#                                 print(k,end='\r')
#                                 jll_valency = jll_valency_list[k]
#                                 for l in range(3):
#                                         jlp_valency = jlp_valency_list[l]
#                                         gyr = np.zeros(101)
#                                         for n in range(5,10):
#                                                 file = h5py.File(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLPP/{1.0}/JLL_VALENCY/{jll_valency}/JLP_VALENCY/{jlp_valency}/JLPP_VALENCY/{2}/LDENS/{ldens}/N/{n}/process.h5", "r")
#                                                 gyr += np.array(file["polyGyration"])
                                
#                                         gyr = gyr/5
#                                         # f = lambda x, tau: np.exp(-x/tau+np.log(gyr[0]-gyr[-1])) +  gyr[-1]
#                                         # res = op.curve_fit(f, np.arange(101), gyr, p0=10)
                                        
#                                         gyr_init[m,i*4+j,k*3+l] = np.mean(gyr[:3])
#                                         gyr_final[m,i*4+j,k*3+l] = np.mean(gyr[-3:])

# np.save(f"/home/ppuel/data/{exp}/gyr_fit/gyr_init_relax.npy", gyr_init)
# np.save(f"/home/ppuel/data/{exp}/gyr_fit/gyr_final_relax.npy", gyr_final)

gyr_init = np.load(f"/home/ppuel/data/{exp}/gyr_fit/gyr_init.npy")
gyr_final = np.load(f"/home/ppuel/data/{exp}/gyr_fit/gyr_final.npy")
gyr_init_relax = np.load(f"/home/ppuel/data/{exp}/gyr_fit/gyr_init_relax.npy")
gyr_final_relax = np.load(f"/home/ppuel/data/{exp}/gyr_fit/gyr_final_relax.npy")

fig, axs = plt.subplots(1,2)

axs[0].hist(gyr_init.flatten(), density = True, histtype = 'step')
axs[0].hist(gyr_init_relax.flatten(), density = True, histtype = 'step')


axs[1].hist(gyr_final.flatten(), density = True, histtype = 'step')
axs[1].hist(gyr_final_relax.flatten(), density = True, histtype = 'step')

fig.savefig(f"/home/ppuel/data/{exp}/gyr_fit/hist_relax.png")

# tcov_imshow = np.load(f"/home/ppuel/data/{exp}/gyr_fit/tcov_imshow.npy")
# tau_imshow = np.load(f"/home/ppuel/data/{exp}/gyr_fit/tau_imshow.npy")

# fig, ax = plt.subplots()

# ax.hist((np.delete(gyr_final[0],[1,5,9,13,2,6,10,14,3,7,11,15],0)).flatten(), bins = 20, density=True, histtype = "step")
# ax.hist((np.delete(gyr_final[0],[0,4,8,12,2,6,10,14,3,7,11,15],0)).flatten(), bins = 20, density=True, histtype = "step")
# ax.hist((np.delete(gyr_final[0],[1,5,9,13,0,4,8,12,3,7,11,15],0)).flatten(), bins = 20, density=True, histtype = "step")
# ax.hist((np.delete(gyr_final[0],[1,5,9,13,0,4,8,12,2,6,10,14],0)).flatten(), bins = 20, density=True, histtype = "step")

# fig.savefig(f"/home/ppuel/data/{exp}/gyr_fit/hist_gyr_final_by_jlp.png")


# fig, ax = plt.subplots()

# ax.hist((gyr_final[0]).flatten(), bins = 20, histtype = "step")
# ax.hist((gyr_final[1]).flatten(), bins = 20, histtype = "step")
# ax.hist((gyr_final[2]).flatten(), bins = 20, histtype = "step")

# fig.savefig(f"/home/ppuel/data/{exp}/gyr_fit/hist_gyr_final.png")


# fig, ax = plt.subplots()

# ax.hist((gyr_init[0]).flatten(), bins = 20, histtype = "step")
# ax.hist((gyr_init[1]).flatten(), bins = 20, histtype = "step")
# ax.hist((gyr_init[2]).flatten(), bins = 20, histtype = "step")

# fig.savefig(f"/home/ppuel/data/{exp}/gyr_fit/hist_gyr_initial.png")



# fig, axs = plt.subplots(3,3, figsize = (36,24))

# for m in range(3):
#         ax = axs[0][m]
#         im = ax.imshow(np.flip(np.transpose(gyr_final[m]),0), cmap = 'plasma')
#         cbar = fig.colorbar(im, ax=ax, shrink=0.7)
#         cbar.set_ticks(ticks=cbar.get_ticks(), labels = [f"{tik:0.2f}" for tik in cbar.get_ticks()], font = font)
#         y_ticks = [f"{jll_valency_list[i]} {jlp_valency_list[j]}" for i in range(3) for j in range(3)]
#         x_ticks = [f"{jll_list[i]} {jlp_list[j]}" for i in range(4) for j in range(4)]
#         ax.set_yticks([8-i for i in range(9)], y_ticks, rotation = 0, font = font)
#         ax.set_xticks([i for i in range(16)], x_ticks, rotation = 90, font = font)
   
# for m in range(3):
#         ax = axs[1][m]
#         im = ax.imshow(np.flip(np.log10(np.transpose(tcov_imshow[m])),0), cmap = 'plasma')
#         cbar = fig.colorbar(im, ax=ax, shrink=0.7)
#         cbar.set_ticks(ticks=cbar.get_ticks(), labels = [f"{tik:0.2f}" for tik in cbar.get_ticks()], font = font)
#         y_ticks = [f"{jll_valency_list[i]} {jlp_valency_list[j]}" for i in range(3) for j in range(3)]
#         x_ticks = [f"{jll_list[i]} {jlp_list[j]}" for i in range(4) for j in range(4)]
#         ax.set_yticks([8-i for i in range(9)], y_ticks, rotation = 0, font = font)
#         ax.set_xticks([i for i in range(16)], x_ticks, rotation = 90, font = font)

  
# for m in range(3):
#         ax = axs[2][m]
#         im = ax.imshow(np.flip(np.transpose(tau_imshow[m])      ,0), cmap = 'plasma')
#         cbar = fig.colorbar(im, ax=ax, shrink=0.7)
#         cbar.set_ticks(ticks=cbar.get_ticks(), labels = [f"{tik:0.2f}" for tik in cbar.get_ticks()], font = font)
#         y_ticks = [f"{jll_valency_list[i]} {jlp_valency_list[j]}" for i in range(3) for j in range(3)]
#         x_ticks = [f"{jll_list[i]} {jlp_list[j]}" for i in range(4) for j in range(4)]
#         ax.set_yticks([8-i for i in range(9)], y_ticks, rotation = 0, font = font)
#         ax.set_xticks([i for i in range(16)], x_ticks, rotation = 90, font = font)
   
# fig.savefig(f"/home/ppuel/data/{exp}/gyr_fit/Comp.png")    

# fig, ax = plt.subplots()

# gyr = np.zeros(101)
# for n in range(5):
#         file = h5py.File(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{2.0}/JLP/{1.0}/JLPP/{1.0}/JLL_VALENCY/{2}/JLP_VALENCY/{2}/JLPP_VALENCY/{2}/LDENS/0.0009/N/{n}/process.h5", "r")
#         ax.plot(np.arange(101),np.array(file["polyGyration"]))
#         gyr += np.array(file["polyGyration"])
# ax.plot(np.arange(101),gyr/5, color = 'black')

# f = lambda x, tau: np.exp(-x/tau+np.log(gyr[0]-gyr[-1])) +  gyr[-1]                                         
# res = op.curve_fit(f, np.arange(101), gyr, p0=10)
# print(res)
# ax.plot(np.arange(101), f(np.arange(101), res[0][0]))


# fig.savefig(f"/home/ppuel/data/{exp}/gyr_fit/gyrtest.png")

plt.close("all")
# os.makedirs(f"/home/ppuel/data/{exp}/stat/logreg/", exist_ok=True)



# # jll_list = ["1.0"]#,"1.5","2.0","4.0"]
# # jlp_list = ["1.0"]#,"1.5","2.0","4.0"]
# # jlpp_list = ["1.0"]#,"1.5","2.0","4.0"]


# # jll_valency_list = ["2"]#,"4","6"]
# # jlp_valency_list = ["2"]#,"4","6"]
# # jlpp_valency_list = ["2"]#,"4","6"]

# # ldens_list = ["0.0009"]#,"0.0018","0.0036"]
# N = 5



# # ldens_list = ["0.010","0.025"]
# # val_list = [str(i*2) for i in range(1,7)]


# # jll_list = ["0.0"]
# # ldens_list = ["0.010"]
# # val_list = ["2"]


# def scatter_hist(x, y, corr_x, corr_y, color, ax, ax_histx, ax_histy):
#         # no labels
#         ax_histx.tick_params(axis="x", labelbottom=False)
#         ax_histy.tick_params(axis="y", labelleft=False)

#         # the scatter plot:
#         ax.scatter(x, y, c = color)

#         x_min = np.min(x)
#         x_range = np.max(x)-np.min(x)
#         x_bins = []
#         for i in range(int(x_range)+2):
#                 x_bins += [i + np.min(x)-.5, i + np.min(x)-.5]
#         x_hist = np.zeros(x_range*2+4)
#         x_corr_hist = -np.zeros(x_range*2+4)
#         for x_i in x:
#                 x_hist[(x_i-x_min)*2+1:(x_i-x_min)*2+3] += 1
        
#         x_hist_log10 = np.where(x_hist < 1.5, 0.15, 0) + np.where(x_hist < .5, -0.15, 0) + np.where(x_hist > 1.5, np.log10(x_hist), 0)

#         ax_histx.plot(x_bins, x_hist_log10)

#         for x_i in corr_x:
#                 x_corr_hist[(x_i-x_min)*2+1:(x_i-x_min)*2+3] += 1/160*4**3

#         x_corr_hist_log10 = np.where(x_corr_hist < 1.5, 0.15, 0) + np.where(x_corr_hist < .5, -0.15, 0) + np.where(x_corr_hist > 1.5, np.log10(x_corr_hist), 0)

#         ax_histx.plot(x_bins, -x_corr_hist_log10)

#         delta_x_hist = np.where((x_hist - x_corr_hist) > 1.5, np.log10(x_hist - x_corr_hist), 0) + np.where((x_hist - x_corr_hist) < -1.5, -np.log10(-(x_hist - x_corr_hist)), 0) + np.where(np.logical_and((x_hist - x_corr_hist) < 1.5, (x_hist - x_corr_hist) > .5), 0.15, 0) + np.where(np.logical_and((x_hist - x_corr_hist) > -1.5, (x_hist - x_corr_hist) < -.5), -0.15, 0)

#         ax_histx.plot(x_bins, delta_x_hist)


#         y_min = np.min(y)
#         y_range = np.max(y)-np.min(y)
#         y_bins = []
#         for i in range(int(y_range)+2):
#                 y_bins += [i + np.min(y)-.5, i + np.min(y)-.5]
#         y_hist = np.zeros(y_range*2+4)
#         y_corr_hist = np.zeros(y_range*2+4)
#         for y_i in y:
#                 y_hist[(y_i-y_min)*2+1:(y_i-y_min)*2+3] += 1

#         y_hist_log10 = np.where(y_hist < 1.5, 0.15, 0) + np.where(y_hist < .5, -0.15, 0) + np.where(y_hist > 1.5, np.log10(y_hist), 0)

#         ax_histy.plot(y_hist_log10, y_bins)
        

#         for y_i in corr_y:
#                 y_corr_hist[(y_i-y_min)*2+1:(y_i-y_min)*2+3] += 1/160*4**3

        
#         y_corr_hist_log10 = np.where(y_corr_hist < 1.5, 0.15, 0) + np.where(y_corr_hist < .5, -0.15, 0) + np.where(y_corr_hist > 1.5, np.log10(y_corr_hist), 0)

#         ax_histy.plot(-y_corr_hist_log10, y_bins)

#         delta_y_hist = np.where((y_hist - y_corr_hist) > 1.5, np.log10(y_hist - y_corr_hist), 0) + np.where((y_hist - y_corr_hist) < -1.5, -np.log10(-(y_hist - y_corr_hist)), 0) + np.where(np.logical_and((y_hist - y_corr_hist) < 1.5, (y_hist - y_corr_hist) > .5), 0.15, 0) + np.where(np.logical_and((y_hist - y_corr_hist) > -1.5, (y_hist - y_corr_hist) < -.5), -0.15, 0)

#         ax_histy.plot(delta_y_hist, y_bins)


#         ax_histx.set_xlim(np.min(x)-.75, np.max(x)+.75)
#         ax_histy.set_ylim(np.min(y)-.75, np.max(y)+.75)
#         ax.set_xlim(np.min(x)-.75, np.max(x)+.75)
#         ax.set_ylim(np.min(y)-.75, np.max(y)+.75)



# def ecdf(a):
#     x, counts = np.unique(a, return_counts=True)
#     cusum = np.cumsum(counts)
#     return x, cusum / cusum[-1]

# # liq_droplet_tau = [[] for ]
# # liq_droplet_size = list(np.zeros((4,4,4,3,3,3)))

# # liq_droplet_tau = [[[[[[[[[] for h in range(5)] for i in range(3)] for j in range(3)] for k in range(3)] for m in range(3)] for n in range(4)] for o in range(4)] for p in range(4)]
# # liq_droplet_size = [[[[[[[[[] for h in range(5)] for i in range(3)] for j in range(3)] for k in range(3)] for m in range(3)] for n in range(4)] for o in range(4)] for p in range(4)]
# # color = []


# # liqMSD = np.zeros((4,4,4,3,3,3,3,5))
# # polyHetMSD = np.zeros((4,4,4,3,3,3,3,5))
# # polyGyration = np.zeros((4,4,4,3,3,3,3,5))


# # liqMSD = np.load(f"/home/ppuel/data/{exp}/stat/liqMSD.npy") #, allow_pickle=True)
# # polyHetMSD = np.load(f"/home/ppuel/data/{exp}/stat/polyHetMSD.npy") #, allow_pickle=True)
# # polyGyration = np.load(f"/home/ppuel/data/{exp}/stat/MSD_gyr/polyGyration.npy") #, allow_pickle=True)

# # liqMSD = np.reshape(liqMSD,(4,4,4,3,3,3,3,5))
# # polyHetMSD = np.reshape(polyHetMSD,(4,4,4,3,3,3,3,5))
# # polyGyration = np.reshape(polyGyration,(4,4,4,3,3,3,3,5))

# # print(polyGyration[0][0][0][0][0][0][0])

# pca = PCA()#n_components=7)


# # liq_droplet_tau = np.load(f"/home/ppuel/data/{exp}/stat/liq_droplet_tau.npy", allow_pickle=True)
# # liq_droplet_size = np.load(f"/home/ppuel/data/{exp}/stat/liq_droplet_size.npy", allow_pickle=True)

# # liq_droplet_tau = np.reshape(liq_droplet_tau,(4,4,4,3,3,3,3,5))
# # liq_droplet_size = np.reshape(liq_droplet_size,(4,4,4,3,3,3,3,5))

# log_reg_coef = np.zeros((4,4,4,3,3,3,3,2))
# log_reg_score = np.zeros((4,4,4,3,3,3,3))


# log_reg_coef = np.reshape(np.load(f"/home/ppuel/data/{exp}/stat/logreg/log_reg_coef.npy"),(4,4,4,3,3,3,3,2))
# log_reg_score = np.reshape(np.load(f"/home/ppuel/data/{exp}/stat/logreg/log_reg_score.npy"),(4,4,4,3,3,3,3))

# min_coef1 = np.min(log_reg_coef[:,:,:,:,:,:,:,0])
# max_coef1 = np.max(log_reg_coef[:,:,:,:,:,:,:,0])
# min_coef2 = np.min(log_reg_coef[:,:,:,:,:,:,:,1])
# max_coef2 = np.max(log_reg_coef[:,:,:,:,:,:,:,1])
# min_score = np.min(log_reg_score[:,:,:,:,:,:,:])
# max_score = np.max(log_reg_score[:,:,:,:,:,:,:])

# dim = [4,4,4,3,3,3,3]
# label = ['jll', 'jlp', 'jlpp', 'jll_valency', 'jlp_valency', 'jlpp_valency', 'ldens']
# condition = ''#'ldens_0.0018_'

# for i in range(len(dim)):
#         print(label[i])
#         coef = np.moveaxis(log_reg_coef,i,-1)

#         log_reg_coef1_by_i = [coef[:,:,:,:,:,:,0,c].flatten() for c in range(dim[i])]
#         log_reg_coef2_by_i = [coef[:,:,:,:,:,:,1,c].flatten() for c in range(dim[i])]
#         log_reg_score_by_i = [np.moveaxis(log_reg_score,i,-1)[:,:,:,:,:,:,c].flatten() for c in range(dim[i])]
#         log_reg_coef1_by_i_cumsum = [None for c in range(dim[i])]
#         log_reg_coef2_by_i_cumsum = [None for c in range(dim[i])]
#         log_reg_score_by_i_cumsum = [None for c in range(dim[i])]

#         fig = plt.figure(figsize=(24,8))
#         ax1 = fig.add_subplot(1,3,1)
#         ax2 = fig.add_subplot(1,3,2)
#         ax3 = fig.add_subplot(1,3,3)

#         for c in range(dim[i]):
#                 # print(len(liq_droplet_tau_by_ldens[o]), liq_droplet_tau_by_ldens[o][0])
#                 x, y = ecdf(log_reg_coef1_by_i[c])
#                 x = np.insert(x, 0, x[0])
#                 y = np.insert(y, 0, 0.)
#                 x = np.insert(x, 0, min_coef1)
#                 y = np.insert(y, 0, 0.)
#                 x = np.insert(x, -1, max_coef1)
#                 y = np.insert(y, -1, 1.)
#                 log_reg_coef1_by_i_cumsum[c] = y
#                 ax1.plot(x, y)
#                 ax1.set_xlim(min_coef1-5/100*(max_coef1-min_coef1),max_coef1+5/100*(max_coef1-min_coef1))
#                 ax1.set_title("Coef 1 : Tau", font = font)

#                 x, y = ecdf(log_reg_coef2_by_i[c])
#                 x = np.insert(x, 0, x[0])
#                 y = np.insert(y, 0, 0.)
#                 x = np.insert(x, 0, min_coef2)
#                 y = np.insert(y, 0, 0.)
#                 x = np.insert(x, -1, max_coef2)
#                 y = np.insert(y, -1, 1.)
#                 log_reg_coef2_by_i_cumsum[c] = x, y
#                 ax2.plot(x, y)
#                 ax2.set_xlim(min_coef2-5/100*(max_coef2-min_coef2),max_coef2+5/100*(max_coef2-min_coef2))
#                 ax2.set_title("Coef 1 : Size", font = font)

#                 x, y = ecdf(log_reg_score_by_i[c])
#                 x = np.insert(x, 0, x[0])
#                 y = np.insert(y, 0, 0.)
#                 x = np.insert(x, 0, min_score)
#                 y = np.insert(y, 0, 0.)
#                 x = np.insert(x, -1, max_score)
#                 y = np.insert(y, -1, 1.)
#                 log_reg_score_by_i_cumsum[c] = x, y
#                 ax3.plot(x, y)
#                 ax3.set_xlim(min_score-5/100*(max_score-min_score),max_score+5/100*(max_score-min_score))
#                 ax3.set_title("Accuracy", font = font)

#         fig.savefig(f"/home/ppuel/data/{exp}/stat/logreg/logreg_well_fit_{condition}by_{label[i]}.png")

# # fig, axs = plt.subplots(figsize=(24,8), nrows=1, ncols=3)

# # c_max = 4*4*4*3*3*3*3*5
# # c = 0
# # for ldens in ldens_list:
# #         o = ldens_list.index(ldens)

# #         liq_tau_array_null = np.load(f"/home/ppuel/data/{exp2}/tau_data_ldens{ldens}.npy")
# #         liq_size_array_null = np.load(f"/home/ppuel/data/{exp2}/size_data_ldens{ldens}.npy")
# #         X_null = [[liq_tau_array_null[null],liq_size_array_null[null]] for null in range(len(liq_tau_array_null))]
# #         y_null = [0 for null in range(len(liq_tau_array_null))]
        
# #         for jll_valency in jll_valency_list:
# #                 l = jll_valency_list.index(jll_valency)

# #                 for jlp_valency in jlp_valency_list:
# #                         m = jlp_valency_list.index(jlp_valency)

# #                         for jlpp_valency in jlpp_valency_list:
# #                                 n = jlpp_valency_list.index(jlpp_valency)

# #                                 for jll in jll_list:
# #                                         i = jll_list.index(jll)

# #                                         for jlp in jlp_list:
# #                                                 j = jlp_list.index(jlp)

# #                                                 for jlpp in jlpp_list:
# #                                                         k = jlpp_list.index(jlpp)  
                                                        
# #                                                         log_reg = LogisticRegression(penalty=None, dual=False)

# #                                                         X = X_null
# #                                                         y = y_null

# #                                                         for p in range(N):                           
# #                                                                 print(f"{c/c_max:0.2f}", end='\r')
                                                                
# #                                                                 for interact in range(len(liq_droplet_tau[i,j,k,l,m,n,o,p])):
                                                                        
# #                                                                         X.append([liq_droplet_tau[i,j,k,l,m,n,o,p][interact], liq_droplet_size[i,j,k,l,m,n,o,p][interact]])
# #                                                                         y.append(1)
                                                                
# #                                                                 c+=1
                                                        
# #                                                         X = np.array(X)
# #                                                         y = np.array(y)

# #                                                         log_reg.fit(X, y)
                                                        
# #                                                         # X_test = np.array([[i,j] for i in range(1,100) for j in range(2,50)])
# #                                                         # c = log_reg.decision_function(X_test)
# #                                                         # c_max = np.max(np.abs(c))

# #                                                         # fig = plt.figure(figsize=(8, 6))
# #                                                         # ax = fig.subplots()
# #                                                         # ax.scatter(x=X_test[:, 0],
# #                                                         #         y=X_test[:, 1], 
# #                                                         #         c=c, 
# #                                                         #         cmap='coolwarm',
# #                                                         #         vmin = -c_max,
# #                                                         #         vmax = c_max,
# #                                                         #         marker='o')
# #                                                         # ax.set_xlabel("Tau", font = font)
# #                                                         # ax.set_ylabel("Size", font = font)
# #                                                         # ax.tick_params(axis="both", **fontlabel)
# #                                                         # fig.suptitle("Logistic Regression Decision Boundary", font = font)
# #                                                         # plt.savefig(f"/home/ppuel/data/{exp}/log_reg_test.png")

# #                                                         log_reg_coef[i,j,k,l,m,n,o,:] = log_reg.coef_
# #                                                         log_reg_score[i,j,k,l,m,n,o] = log_reg.score(X, y)
                                                      
# #                                                         # fig = plt.figure(figsize=(8, 6))
# #                                                         # ax = fig.subplots()
# #                                                         # ax.hist(log_reg.decision_function(X[y == 0]), density=True, histtype='step', color='blue')
# #                                                         # ax.hist(log_reg.decision_function(X[y == 1]), density=True, histtype='step', color='red')
# #                                                         # ax.set_xlabel("Confidence Score ", font = font)
# #                                                         # ax.set_ylabel("Density", font = font)
# #                                                         # ax.tick_params(axis="both", **fontlabel)
# #                                                         # fig.suptitle("Decision Function", font = font)
# #                                                         # plt.savefig(f"/home/ppuel/data/{exp}/Confidence_Score.png")


# # axs[0].hist(log_reg_coef[:,:,:,:,:,:,:,0].flatten(), density = True, histtype = 'step')

# # axs[1].hist(log_reg_coef[:,:,:,:,:,:,:,1].flatten(), density = True, histtype = 'step')

# # axs[2].hist(log_reg_score[:,:,:,:,:,:,:].flatten(), density = True, histtype = 'step')

# # fig.savefig(f"/home/ppuel/data/{exp}/log_reg_stat_total.png")

# # np.save(f"/home/ppuel/data/{exp}/stat/log_reg_coef.npy", log_reg_coef.flatten())
# # np.save(f"/home/ppuel/data/{exp}/stat/log_reg_score.npy", log_reg_score.flatten())


# # for i in range(len(liq_droplet_tau)):
# #        for j in range(len(liq_droplet_tau[i])):
# #               data.append([liq_droplet_tau[i][j],liq_droplet_size[i][j]])

# # pca.fit(data)

# # x=pca.transform(data)
# # print(pca.components_)
# # print(pca.explained_variance_ratio_)

# # plt.figure(figsize=(10,10))
# # plt.scatter(x[:,0],x[:,1])
# # plt.xlabel('pc1')
# # plt.ylabel('pc2')
# # plt.savefig(f'/home/ppuel/data/{exp}/PCA_total.png')

# # dim1 = 5
# # dim2 = 3

# # liqMSD_by_1 = [liqMSD[:,:,:,:,:,:,:,c].flatten() for c in range(dim1)]
# # polyHetMSD_by_1 = [polyHetMSD[:,:,:,:,:,:,:,c].flatten() for c in range(dim1)]
# # polyGyration_by_1 = [polyGyration[:,:,:,:,:,:,:,c].flatten() for c in range(dim1)]
# # liqMSD_by_1_cumsum = [None for c in range(dim1)]
# # polyHetMSD_by_1_cumsum = [None for c in range(dim1)]
# # polyGyration_by_1_cumsum = [None for c in range(dim1)]

# # liqMSD_by_2 = [liqMSD[:,:,:,:,:,:,c,:].flatten() for c in range(dim2)]
# # polyHetMSD_by_2 = [polyHetMSD[:,:,:,:,:,:,c,:].flatten() for c in range(dim2)]
# # polyGyration_by_2 = [polyGyration[:,:,:,:,:,:,c,:].flatten() for c in range(dim2)]
# # liqMSD_by_2_cumsum = [None for c in range(dim2)]
# # polyHetMSD_by_2_cumsum = [None for c in range(dim2)]
# # polyGyration_by_2_cumsum = [None for c in range(dim2)]


# # fig = plt.figure(figsize=(36,12))
# # ax1 = fig.add_subplot(1,3,1)
# # ax2 = fig.add_subplot(1,3,2)
# # ax3 = fig.add_subplot(1,3,3)

# # for c in range(dim1):
# #         # print(len(liq_droplet_tau_by_ldens[o]), liq_droplet_tau_by_ldens[o][0])
# #         x, y = ecdf(liqMSD_by_1[c])
# #         x = np.insert(x, 0, x[0])
# #         y = np.insert(y, 0, 0.)
# #         liqMSD_by_1_cumsum[c] = y
# #         ax1.plot(x, y)

# #         x, y = ecdf(polyHetMSD_by_1[c])
# #         x = np.insert(x, 0, x[0])
# #         y = np.insert(y, 0, 0.)
# #         polyHetMSD_by_1_cumsum[c] = x, y
# #         ax2.plot(x, y)

# #         x, y = ecdf(polyGyration_by_1[c])
# #         x = np.insert(x, 0, x[0])
# #         y = np.insert(y, 0, 0.)
# #         polyGyration_by_1_cumsum[c] = x, y
# #         ax3.plot(x, y)

# # fig.savefig(f"/home/ppuel/data/{exp}/stat/MSD_gyr_by_N.png")


# # fig = plt.figure(figsize=(36,12))
# # ax1 = fig.add_subplot(1,3,1)
# # ax2 = fig.add_subplot(1,3,2)
# # ax3 = fig.add_subplot(1,3,3)

# # for c in range(dim2):
# #         # print(len(liq_droplet_tau_by_ldens[o]), liq_droplet_tau_by_ldens[o][0])
# #         x, y = ecdf(liqMSD_by_2[c])
# #         x = np.insert(x, 0, x[0])
# #         y = np.insert(y, 0, 0.)
# #         liqMSD_by_2_cumsum[c] = y
# #         ax1.plot(x, y)

# #         x, y = ecdf(polyHetMSD_by_2[c])
# #         x = np.insert(x, 0, x[0])
# #         y = np.insert(y, 0, 0.)
# #         polyHetMSD_by_2_cumsum[c] = x, y
# #         ax2.plot(x, y)

# #         x, y = ecdf(polyGyration_by_2[c])
# #         x = np.insert(x, 0, x[0])
# #         y = np.insert(y, 0, 0.)
# #         polyGyration_by_2_cumsum[c] = x, y
# #         ax3.plot(x, y)

# # fig.savefig(f"/home/ppuel/data/{exp}/stat/MSD_gyr_by_density.png")


# # diff = 0

# # for c in range(dim1-1):
# #         diff += np.sum(np.abs(tau_by_jll_valency_cumsum[c] - tau_by_jll_valency_cumsum[c+1]))/100
# # print(diff/(dim1-1))


# # liq_droplet_tau_by_jll_valency = [list(itertools.chain.from_iterable(list(liq_droplet_tau[:,:,:,:,:,:,:,c].flatten()))) for c in range(dim1)]
# # liq_droplet_size_by_jll_valency = [list(itertools.chain.from_iterable(list(liq_droplet_size[:,:,:,:,:,:,:,c].flatten()))) for c in range(dim1)]
# # tau_by_jll_valency_cumsum = [None for c in range(dim1)]
# # size_by_jll_valency_cumsum = [None for c in range(dim1)]

# # fig = plt.figure(figsize=(24,12))
# # ax = fig.add_subplot(1,2,1)
# # ax2 = fig.add_subplot(1,2,2)

# # for c in range(dim1):
# #         # print(len(liq_droplet_tau_by_ldens[o]), liq_droplet_tau_by_ldens[o][0])
# #         x, y = ecdf(liq_droplet_tau_by_jll_valency[c])
# #         x = np.insert(x, 0, x[0])
# #         y = np.insert(y, 0, 0.)
# #         tau_by_jll_valency_cumsum[c] = y
# #         ax.plot(x, y)

# #         x, y = ecdf(liq_droplet_size_by_jll_valency[c])
# #         x = np.insert(x, 0, x[0])
# #         y = np.insert(y, 0, 0.)
# #         size_by_jll_valency_cumsum[c] = x, y
# #         ax2.plot(x, y)

# # fig.savefig(f"/home/ppuel/data/{exp}/stat/droplet_by_N.png")

# # diff = 0

# # for c in range(dim1-1):
# #         diff += np.sum(np.abs(tau_by_jll_valency_cumsum[c] - tau_by_jll_valency_cumsum[c+1]))/100
# # print(diff/(dim1-1))

# # liq_droplet_tau_by_jll = [list(itertools.chain.from_iterable(list(liq_droplet_tau[:,:,:,:,:,:,c,:].flatten()))) for c in range(dim2)]
# # liq_droplet_size_by_jll = [list(itertools.chain.from_iterable(list(liq_droplet_size[:,:,:,:,:,:,c,:].flatten()))) for c in range(dim2)]
# # tau_by_jll_cumsum = [None for c in range(dim2)]
# # size_by_jll_cumsum = [None for c in range(dim2)]


# # fig = plt.figure(figsize=(24,12))
# # ax = fig.add_subplot(1,2,1)
# # ax2 = fig.add_subplot(1,2,2)


# # for c in range(dim2):
# #         x, y = ecdf(liq_droplet_tau_by_jll[c])
# #         x = np.insert(x, 0, x[0])
# #         y = np.insert(y, 0, 0.)
# #         tau_by_jll_cumsum[c] = y
# #         ax.plot(x, y)

# #         x, y = ecdf(liq_droplet_size_by_jll[c])
# #         x = np.insert(x, 0, x[0])
# #         y = np.insert(y, 0, 0.)
# #         size_by_jll_cumsum[c] = y
# #         ax2.plot(x, y)


# # fig.savefig(f"/home/ppuel/data/{exp}/stat/droplet_by_ldens_valency.png")


# # diff = 0

# # for c in range(dim2-1):
# #         diff += np.sum(np.abs(tau_by_jll_cumsum[c] - tau_by_jll_cumsum[c+1]))/100
# # print(diff/(dim2-1))

#         # print(np.shape(liq_droplet_size[:,:,:,:,:,:,o,:].flatten()))
#         # print(np.shape(liq_droplet_size[:,:,:,:,:,:,o+1,:].flatten()))
       
#         # # liq_tau_array = np.load(f"/home/ppuel/data/{exp2}/tau_data_ldens{ldens}.npy")
#         # # liq_size_array = np.load(f"/home/ppuel/data/{exp2}/size_data_ldens{ldens}.npy")
        

#         # # liq_MSD_array = np.load(f"/home/ppuel/data/{exp2}/liq_MSD_ldens{ldens}.npy")
#         # # poly_MSD_array = np.load(f"/home/ppuel/data/{exp2}/poly_MSD_ldens{ldens}.npy")
#         # # poly_gyr_array = np.load(f"/home/ppuel/data/{exp2}/gyr_data_ldens{ldens}.npy")
#         # # liq_MSD_null = np.mean([(np.log10(msd) - np.log10(liq_MSD_array[1]))/(np.log10(c+1)) for c, msd in enumerate(list(liq_MSD_array[1:50]))][1:])
#         # # poly_MSD_null = np.mean([(np.log10(msd) - np.log10(poly_MSD_array[1]))/(np.log10(c+1)) for c, msd in enumerate(list(poly_MSD_array[1:50]))][1:])
#         # # poly_gyr_null = np.mean(poly_gyr_array)


# # liq_droplet_tau = [[[[[[[[[] for h in range(5)] for i in range(3)] for j in range(3)] for k in range(3)] for m in range(3)] for n in range(4)] for o in range(4)] for p in range(4)]
# # liq_droplet_size = [[[[[[[[[] for h in range(5)] for i in range(3)] for j in range(3)] for k in range(3)] for m in range(3)] for n in range(4)] for o in range(4)] for p in range(4)]

# # liqMSD = np.zeros((4,4,4,3,3,3,3,5))
# # polyHetMSD = np.zeros((4,4,4,3,3,3,3,5))
# # polyGyration = np.zeros((4,4,4,3,3,3,3,5))

# # for ldens in ldens_list:
# #         o = ldens_list.index(ldens)
      
                        
# #         for jll_valency in jll_valency_list:
# #                 print(jll_valency, end='\n')
# #                 l = jll_valency_list.index(jll_valency)

# #                 for jlp_valency in jlp_valency_list:
# #                         m = jlp_valency_list.index(jlp_valency)

# #                         for jlpp_valency in jlpp_valency_list:
# #                                 n = jlpp_valency_list.index(jlpp_valency)

# #                                 for jll in jll_list:
# #                                         i = jll_list.index(jll)

# #                                         for jlp in jlp_list:
# #                                                 print(jlp, end='\r')
# #                                                 j = jlp_list.index(jlp)

# #                                                 for jlpp in jlpp_list:
# #                                                         k = jlpp_list.index(jlpp)  

# #                                                         for p in range(N):                           
                                                                
# #                                                                 file = h5py.File(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLPP/{jlpp}/JLL_VALENCY/{jll_valency}/JLP_VALENCY/{jlp_valency}/JLPP_VALENCY/{jlpp_valency}/LDENS/{ldens}/N/{p}/process.h5", "r")
# #                                                                 liqMSD[i,j,k,l,m,n,o,p] = np.mean([(np.log10(msd) - np.log10(file["liqMSD"][1]))/(np.log10(c+1)) for c, msd in enumerate(list(file["liqMSD"][1:50]))][1:])
# #                                                                 polyHetMSD[i,j,k,l,m,n,o,p] = np.mean([(np.log10(msd) - np.log10(file["polyHetMSD"][1]))/(np.log10(c+1)) for c, msd in enumerate(list(file["polyHetMSD"][1:50]))][1:])
# #                                                                 polyGyration[i,j,k,l,m,n,o,p] = np.mean(file["polyGyration"][-10:])

# #                                                                 # droplet_dict = pickle.load(file)
# #                                                                 file.close()
                                                                
#                                                                 # for droplet in droplet_dict.values():
#                                                                 #         liq_droplet_tau[i][j][k][l][m][n][o][p].append(droplet.tau)
#                                                                 #         liq_droplet_size[i][j][k][l][m][n][o][p].append(np.max(droplet.sizes))
#                                                                 #         # color.append(3/4/(1/float(jll)+1/float(jlp)+1/float(jlpp))*256)


# # np.save(f"/home/ppuel/data/{exp}/stat/liqMSD.npy", liqMSD.flatten())
# # np.save(f"/home/ppuel/data/{exp}/stat/polyHetMSD.npy", polyHetMSD.flatten())
# # np.save(f"/home/ppuel/data/{exp}/stat/polyGyration.npy", polyGyration.flatten())


#                         # font = fm.FontProperties(weight='bold',
#                         #                                 style='normal', size=12)
#                         # fontlabel = {"labelsize" : 12}


#                         # fig, axs = plt.subplot_mosaic([['histx', '.'],
#                         #                         ['scatter', 'histy']],
#                         #                         figsize=(10, 10),
#                         #                         width_ratios=(4, 2), height_ratios=(2, 4),
#                         #                         layout='constrained')
                        
#                         # fig.suptitle(f"JLL : {jll}, JLP : {jlp}, JLPP : {jlpp},\nJLL_VALENCY : {jll_valency}, JPL_VALENCY : {jlp_valency}, JLPP_VALENCY : 6,\nLDENS : {ldens}")
#                         # scatter_hist(liq_droplet_tau, liq_droplet_size, liq_tau_array, liq_size_array, color, axs['scatter'], axs['histx'], axs['histy'])
                        
#                         # axs['scatter'].set_ylabel("Size Droplet (Number of PRC1)", font = font)
#                         # axs['scatter'].set_xlabel("Time (kMCS)", font = font)
                        
#                         # axs['scatter'].tick_params(axis="both", **fontlabel)
#                         # axs['histx'].tick_params(axis="y", **fontlabel)
#                         # axs['histy'].tick_params(axis="x", **fontlabel)

#                         # fig.savefig(f"/home/ppuel/data/{exp}/Droplet_JLL_VALENCY_{jll_valency}_JPL_VALENCY_{jlp_valency}_JLPP_VALENCY_6_LDENS_{ldens}.png")
#                         # plt.close(fig=fig)
                        
#                                                 # file = h5py.File(os.path.join(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLPP/{jlpp}/JLL_VALENCY/{jll_valency}/JPL_VALENCY/{jlp_valency}/JLPP_VALENCY/6/LDENS/{ldens}/N/0/", "process.h5"),'r')
                                               

#                                                 # poly_gyr = np.mean(file["polyGyration"][-10:])
#                                                 # poly_gyr_matrix[j,k,l,m,o] = poly_gyr
#                                                 # liq_MSD = np.mean([(np.log10(msd) - np.log10(file["liqMSD"][1]))/(np.log10(c+1)) for c, msd in enumerate(list(file["liqMSD"][1:50]))][1:])
#                                                 # liq_MSD_matrix[j,k,l,m,o] = liq_MSD
#                                                 # poly_MSD = np.mean([(np.log10(msd) - np.log10(file["polyHetMSD"][1]))/(np.log10(c+1)) for c, msd in enumerate(list(file["polyHetMSD"][1:50]))][1:])
#                                                 # poly_MSD_matrix[j,k,l,m,o] = poly_MSD
#                                                 # if abs(poly_gyr-poly_gyr_null*0.8)/poly_gyr_null*0.8 < 0.1 and abs(liq_MSD-liq_MSD_null)/liq_MSD_null < 0.1 and abs(poly_MSD-poly_MSD_null)/poly_MSD_null < 0.1:
#                                                 #         print(jll,jlp, jlpp, jll_valency, jlp_valency, ldens)
#                                                 # file.close()

# # np.save(f"/Xnfs/physbiochrom/ppuel/data/{exp}/poly_gyr_matrix.npy", poly_gyr_matrix.flatten())
# # np.save(f"/Xnfs/physbiochrom/ppuel/data/{exp}/poly_MSD_matrix.npy", poly_MSD_matrix.flatten())
# # np.save(f"/Xnfs/physbiochrom/ppuel/data/{exp}/liq_MSD_matrix.npy", liq_MSD_matrix.flatten())