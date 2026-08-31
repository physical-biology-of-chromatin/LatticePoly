import subprocess
import os
from LiqCluster_lifeTime import LifeTime, Droplet, Event
import pickle
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import h5py

exp = "EXP34_bitorax_first_trial"

os.makedirs(f"/home/ppuel/data/{exp}/droplet_tracking", exist_ok=True)

jll_list = ["0.0", "2.2", "4.4"]
jlp_list = ["2.2", "4.4"]
jlpp_list = ["2.2", "4.4"]
jll_valency_list = [str(i*2) for i in range(1,7)]
N = 1

# ldens_list = ["0.010","0.025"]
# val_list = [str(i*2) for i in range(1,7)]


# jll_list = ["0.0"]
# ldens_list = ["0.010"]
# val_list = ["2"]



list_gyr = []

for jll in jll_list:
        for jlp in jlp_list:
                for jlpp in jlpp_list:
                        for jll_valency in jll_valency_list:
                                        for n in range(N):

                                                print(f"jll {jll} jlp {jlp} jlpp {jlpp} jll_valency {jll_valency}")

                                                # subprocess.run(f".venv/bin/python3 resources/h5py/LiqCluster_hist.py /home/ppuel/data/{exp}/JLL/{jll}/JLL_VALENCY/{val}/LDENS/{ldens}/N/{n}/ -1 -1", shell=True, executable="/bin/bash")
                                                # subprocess.run(f".venv/bin/python3 resources/h5py/LiqCluster_lifeTime.py /Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLL_VALENCY/{val}/LDENS/{ldens}/N/{n}/ -1", shell=True, executable="/bin/bash")

                                                try:
                                                        # file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLPP/{jlpp}/JLL_VALENCY/{jll_valency}/N/{n}/liq_droplets.pickle", "rb")
                                                        # droplet_dict = pickle.load(file)
                                                        # print(len(droplet_dict))
                                                        # file.close()

                                                        file = h5py.File(os.path.join(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLP/{jlp}/JLPP/{jlpp}/JLL_VALENCY/{jll_valency}/N/{n}/", "process.h5"),'r')
                                                        list_gyr.append(np.mean(file["polyGyration"][-10:]))
                                                        file.close()
                                                except Exception as ex:
                                                        droplet_dict = {}
                                                        print("error")
                                                        pass#print("Error during unpickling object (Possibly unsupported):", ex)
                                                

                                                # if droplet_dict != {}:
                                                #         fig = plt.figure()
                                                #         ax = fig.add_subplot()

                                                #         max_tau = -1
                                                #         for droplet in droplet_dict.values():
                                                #                 max_tau = max(max_tau, droplet.tau)
                                                #         rainbow = mpl.colormaps['rainbow'].resampled(max_tau)

                                                #         for droplet in droplet_dict.values():
                                                #                 if droplet.tau > 1 and np.max(droplet.sizes) > 4:
                                                #                         ax.scatter([droplet.nodes[0][0]+ i for i in range(len(droplet.sizes))], droplet.sizes, color = rainbow(droplet.tau), s = 100)
                                                #                         ax.plot([droplet.nodes[0][0]+ i for i in range(len(droplet.sizes))], droplet.sizes, color = rainbow(droplet.tau), linewidth = 10)


                                                #         fig.savefig(f"/home/ppuel/data/{exp}/droplet_tracking/fig_droplet_size_jll_{jll}_jlp_{jlp}_jlpp_{jlpp}_jll_valency_{jll_valency}.png")
                                                #         plt.close(fig=fig)

list_gyr = np.array(list_gyr)
list_gyr = list_gyr/np.max(list_gyr)

fig = plt.figure()
ax = fig.add_subplot()
ax.scatter(range(len(list_gyr)), list_gyr)
plt.savefig(f"/home/ppuel/data/{exp}/gyration_data")
plt.close(fig=fig)
