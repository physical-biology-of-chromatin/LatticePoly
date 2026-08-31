# import subprocess
import os
# from LiqCluster_lifeTime import LifeTime, Droplet, Event
# import pickle
# import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import h5py
# import scipy.stats as st
# import itertools
# import scipy.optimize as op
# from sklearn.decomposition import PCA
# from sklearn.linear_model import LogisticRegression
# import seaborn as sns


from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

font = fm.FontProperties(weight='bold',
                                style='normal', size=20)
fontlabel = {"labelsize" : 20}

jll_list = ["0.5","1.0","2.0"]
jlp_list = ["0.2","0.5","1.0"]
EV_list = ["0","1","10"]
ninter_list = ["100000"]
ldens_list = ["0.00068","0.00135","0.00270"]
percentage_list = ["12.5%", "25%", "50%"]


def conversion_kJ_per_mol_to_kBT(E):
        E = np.array([float(i) for i in E])
        T = 300
        kB = 1.380e-23
        Na = 6.022e23
        E = E / ( kB * Na * T * 1e-3 )
        return([f"{i:0.1f}" for i in E])


EV_list_conv = conversion_kJ_per_mol_to_kBT(np.array(EV_list))
# jll_list = ["1.0"]
# jlp_list = ["0.5"]
# EV_list = ["0"]
# ninter_list = ["100000"]
# ldens_list = ["0.00135"]


N = 40

plasma = mpl.colormaps["plasma"].resampled(8)

exp = "EXP43_Nazli_first_simulation_3R"

os.makedirs(f"/home/ppuel/data/{exp}/mapping/", exist_ok=True)



c_max = 3*3*3*3
c=0



fig, axs = plt.subplots(3,3)

fig.set_figheight(22)
fig.set_figwidth(24)

# min_X = 1e22
# max_X = -1
# min_Y = 1e22
# max_Y = -1





liqFraction_matrix = np.load(f"/home/ppuel/data/{exp}/liqFraction_matrix.npy")

liqDropNum_matrix = np.load(f"/home/ppuel/data/{exp}/liqDropNum_matrix.npy")

liqDropSize_matrix = np.load(f"/home/ppuel/data/{exp}/liqDropSize_matrix.npy")

max_size = [np.max(liqDropSize_matrix[i]) for i in range(3)]
print(max_size)
for ldens in ldens_list:
        i = ldens_list.index(ldens)

        for EV in EV_list:

                j = EV_list.index(EV)

                
                # for jlp in jlp_list:
                #         k = jlp_list.index(jlp)


                #         for jll in jll_list:
                #                 l = jll_list.index(jll)

                ax = axs[i, j]
                im = ax.imshow(liqDropNum_matrix[i, j],origin = 'lower', cmap = 'plasma', vmin = 0, vmax = 8)#max_size[i])
                ax.set_xlabel("P <-> P", font = font)
                ax.set_ylabel("P <-> H", font = font)
                ax.set_xticks([i for i in range(3)], conversion_kJ_per_mol_to_kBT(np.array(jll_list)), font = font)
                ax.set_yticks([i for i in range(3)], conversion_kJ_per_mol_to_kBT(np.array(jlp_list)), font = font)
                ax.set_title(f"P% = {percentage_list[i]}, EV = {EV_list_conv[j]}", font = font)


        ax_divider = make_axes_locatable(axs[i,2])
        # Add an Axes to the right of the main Axes.
        cax1 = ax_divider.append_axes("right", size="7%", pad="2%")
        cb1 = fig.colorbar(im, cax=cax1)
fig.savefig(f"/home/ppuel/data/{exp}/liqDropNum_matrix.png")




#                 # local_density_mean = []
#                 # sizes_max = []
#                 # tau = []

#                                 # liqDropNum_dict = {}
#                                 # liqFraction_dict = {}
#                                 # liqDropSize_dict = {}
#                                 # count = {}

                                # for ninter in ninter_list:


                                #         c += 1
                                #         print(f"{c/c_max*100:2.2f} %", end = '\r')

                                #         file = h5py.File(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/EV/{EV}/LDENS/{ldens}/NINTER/{ninter}/post_process.h5")


                                #         liqDropNum = np.reshape(file['liqDropNum'],101)
                                #         liqFraction = np.reshape(file['liqFraction'],101)
                                #         liqDropSize = np.reshape(file['liqDropSize'],101)
                                        # polyHomSimTime = np.reshape(file['polyHetSimTime'],101)[-10:]


#                                         liqFraction_matrix[i,j,k,l] = np.mean(np.reshape(file['liqFraction'],101)[-10:])
#                                         liqDropNum_matrix[i,j,k,l] = np.mean(np.reshape(file['liqDropNum'],101)[-10:])
#                                         liqDropSize_matrix[i,j,k,l] = np.mean(np.reshape(file['liqDropSize'],101)[-10:])

#                                         file.close()


                                        # min_X = min(min_X, np.min(polyHomSimTime))
                                        # max_X = max(max_X, np.max(polyHomSimTime))
                                        # min_Y = min(min_Y, np.min(liqFraction))
                                        # max_Y = max(max_Y, np.max(liqFraction))


                                        # for i in range(1,101):
                                        #         if i*ninter in liqFraction_dict.keys():
                                        #                 liqDropNum_dict[polyHomSimTime[i]] += liqDropNum[i]
                                        #                 liqFraction_dict[polyHomSimTime[i]] += liqFraction[i]
                                        #                 liqDropSize_dict[polyHomSimTime[i]] += liqDropSize[i]
                                        #                 count[polyHomSimTime[i]] += 1

                                        #         else:
                                        #                 liqDropNum_dict[polyHomSimTime[i]] = liqDropNum[i]
                                        #                 liqFraction_dict[polyHomSimTime[i]] = liqFraction[i]
                                        #                 liqDropSize_dict[polyHomSimTime[i]] = liqDropSize[i]
                                        #                 count[polyHomSimTime[i]] = 1
                                # if ldens_list.index(ldens) == 0 and EV_list.index(EV) == 0:
                                #         ax.scatter(liqFraction_dict.keys(), [liqFraction_dict[frame]/count[frame] if liqFraction_dict[frame]/count[frame] > 1 else -1 for frame in liqFraction_dict.keys()], label = f"jll : {jll}; jlp : {jlp}")#, color = 'blue')
                                # else:
                                #         ax.scatter(liqFraction_dict.keys(), [liqFraction_dict[frame]/count[frame] if liqFraction_dict[frame]/count[frame] > 1 else -1 for frame in liqFraction_dict.keys()])


                                        # from mpl_toolkits import axisartist
                                        # from mpl_toolkits.axes_grid1 import host_subplot
                                        # fig = plt.figure(figsize=(8,6))
                                        # host = host_subplot(111, axes_class=axisartist.Axes, figure=fig)
                                        # plt.subplots_adjust(right=0.75)

                                        # par1 = host.twinx()
                                        # par2 = host.twinx()

                                        # par2.axis["right"] = par2.new_fixed_axis(loc="right", offset=(60, 0))

                                        # par1.axis["right"].toggle(all=True)
                                        # par2.axis["right"].toggle(all=True)

                                        # p1, = host.plot(liqFraction, label="Liquid fraction")
                                        # p2, = par1.plot(liqDropNum, label="Number of droplets")
                                        # par1.plot([0, 100], [1, 1], "--", color = 'grey')
                                        # p3, = par2.plot(liqDropSize, label="Sizes of the droplets")

                                        # host.set(xlim=(-.5, 100), ylim=(0, 1), xlabel="MCS", ylabel="Liquid fraction")
                                        # par1.set(ylabel="Number of droplets")
                                        # par2.set(ylim=(1, float(ldens)*51**3*4), ylabel="Sizes of the droplets")

                                        # host.legend()

                                        # host.axis["left"].label.set_color(p1.get_color())
                                        # par1.axis["right"].label.set_color(p2.get_color())
                                        # par2.axis["right"].label.set_color(p3.get_color())

                                        # fig.savefig(f"/home/ppuel/data/{exp}/mapping/{ldens}_{EV}_{jll}_{jlp}.png")

                                        # plt.close(fig=fig)

                                        # ax1.plot(liqFraction_dict.keys(), [liqFraction_dict[frame]/count[frame] for frame in liqFraction_dict.keys()], color = 'blue')
                                        # ax2.plot(liqDropNum_dict.keys(), [liqDropNum_dict[frame]/count[frame] for frame in liqDropNum_dict.keys()], color = 'orange')
                                        # ax3.plot(liqDropSize_dict.keys(), [liqDropSize_dict[frame]/count[frame] for frame in liqDropSize_dict.keys()], color = 'green')

                                        # ax1.set_yticklabels(ax1.get_xticklabels(), color = 'blue', font = font)
                                        # ax2.set_yticklabels(ax2.get_yticklabels(), color = 'orange', font = font)
                                        # ax3.set_yticklabels(ax3.get_yticklabels(), color = 'green', font = font)
                # ax.set_xscale('log')
                # ax.set_yscale('log')
                # if ldens_list.index(ldens) == 0:
                #         ax.set_ylabel("liqFraction_dict", color = 'blue')
                                # ax1.set_ylabel("liqFraction", color = 'blue', font = font)
                                # ax2.set_ylabel("liqDropNum", color = 'orange', font = font)
                                # ax3.set_ylabel("liqDropSize", color = 'green', font = font)

        # fig.suptitle(f'Jlp : {jlp}, Jll : {jll}, Ldens : {ldens}')

# np.save(f"/home/ppuel/data/{exp}/liqFraction_matrix.npy",liqFraction_matrix)

# np.save(f"/home/ppuel/data/{exp}/liqDropNum_matrix.npy",liqDropNum_matrix)

# np.save(f"/home/ppuel/data/{exp}/liqDropSize_matrix.npy",liqDropSize_matrix)

# fig.legend()

# fig.savefig(f"/home/ppuel/data/{exp}/figure/liqFraction_fix.png")

# print(min_X, max_X, min_Y, max_Y)

                                        # file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/EV/{EV}/LDENS/{ldens}/NINTER/{ninter}/N/{n}/liq_droplets.pickle", "rb")
                                        # droplet_dict = pickle.load(file)
                                        # file.close()
                                        # for droplet in droplet_dict.values():
                                        #         local_density_mean.append(np.mean(droplet.local_density))
                                        #         sizes_max.append(np.max(droplet.sizes))
                                        #         tau.append(droplet.tau)



                # np.save(f"/home/ppuel/data/{exp}/data/local_density_mean_ldens_{ldens}_ninter_{ninter}.npy", local_density_mean)
                # np.save(f"/home/ppuel/data/{exp}/data/sizes_max_ldens_{ldens}_ninter_{ninter}.npy", sizes_max)
                # np.save(f"/home/ppuel/data/{exp}/data/tau_ldens_{ldens}_ninter_{ninter}.npy", tau)
