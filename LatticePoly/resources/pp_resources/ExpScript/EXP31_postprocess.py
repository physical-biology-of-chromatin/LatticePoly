import subprocess
import os
from LiqCluster_lifeTime import LifeTime, Droplet, Event
import pickle
import networkx as nx

exp = "EXP31_liqFraction_liqlifeTime_test"

jll_list = ["0.0","1.0","1.6","2.2"]
ldens_list = ["0.010","0.025"]
val_list = [str(i*2) for i in range(1,7)]


# jll_list = ["0.0"]
# ldens_list = ["0.010"]
# val_list = ["2"]
N = 1

for jll in jll_list:
        for ldens in ldens_list:
                for val in val_list:
                        for n in range(N):

                                print(f"jll {jll} ldens {ldens} val {val}")

                                # subprocess.run(f".venv/bin/python3 resources/h5py/LiqCluster_hist.py /home/ppuel/data/{exp}/JLL/{jll}/JLL_VALENCY/{val}/LDENS/{ldens}/N/{n}/ -1 -1", shell=True, executable="/bin/bash")
                                subprocess.run(f".venv/bin/python3 resources/h5py/LiqCluster_lifeTime.py /Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLL_VALENCY/{val}/LDENS/{ldens}/N/{n}/ -1", shell=True, executable="/bin/bash")


                                try:
                                        file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLL_VALENCY/{val}/LDENS/{ldens}/N/0/liq_droplets.pickle", "rb")
                                        droplet_dict = pickle.load(file)
                                        print(type(droplet_dict[0].in_event[0][0]), droplet_dict[0].in_event[0][0])
                                        file.close()
                                except Exception as ex:
                                        pass#print("Error during unpickling object (Possibly unsupported):", ex)
                                try:
                                        G = nx.read_gml(f"/Xnfs/physbiochrom/ppuel/data/{exp}/JLL/{jll}/JLL_VALENCY/{val}/LDENS/{ldens}/N/0/graph.gml")
                                        print(type(G.nodes(data = True)), G.nodes(data = True))
                                except Exception as ex:
                                        pass#print("Error during unpickling object (Possibly unsupported):", ex)
