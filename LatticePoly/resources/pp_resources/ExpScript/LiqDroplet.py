##
##  LiqCluster.py
##  LatticePoly
##
##  Created by ppuel on 19/02/2025.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np
import networkx as nx

import h5py

import matplotlib.pyplot as plt

import matplotlib as mpl

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import matplotlib.font_manager as font_manager

import pickle

from hdf5Reader import hdf5Reader

from vtk import vtkXMLPolyDataWriter
from vtk import vtkPoints
from vtk import vtkFloatArray
from vtk import vtkPolyData
from vtk import vtkCellArray
from vtk import vtkIntArray
from vtk import vtkLine
from vtk import vtkCubeSource

import time

# BrBG = PRGn_10.mpl_colormap
viridis = plt.cm.Greens

class Droplet():
        def __init__(self, droplet_id, node, size, in_event, center_of_mass, local_density):
                self.droplet_id = droplet_id
                self.nodes = [node]
                self.tau = 1
                self.sizes = [size]
                self.in_event = [in_event]
                self.out_event = [[]]
                self.center_of_mass = [center_of_mass]
                self.local_density = [local_density]
                self.frame_start = node[0]

        def append(self, node, size, in_event, center_of_mass, local_density):
                self.nodes.append(node)
                self.tau += 1
                self.sizes.append(size)
                self.in_event.append(in_event)
                self.out_event.append([])
                self.center_of_mass.append(center_of_mass)
                self.local_density.append(local_density)
                
                

        def append_in_event(self, node, in_event):
                self.in_event[self.nodes.index(node)].append(in_event)

        def append_out_event(self, node, out_event):
                self.out_event[self.nodes.index(node)].append(out_event)

        def get_out_event(self, node):
                return(self.out_event[self.nodes.index(node)])

        def is_persistante_at_node(self, node):
                for event in self.out_event[self.nodes.index(node)]:
                        if event.name == 'persistance':
                                return(True)
                return(False)

        def was_persistante_at_node(self, node):
                for event in self.in_event[self.nodes.index(node)]:
                        if event.name == 'persistance':
                                return(True)
                return(False)

class Event():
        def __init__(self, name, strength, node):
                self.name = name
                self.strength = strength
                self.node = node

        def __str__(self):
                return(str(self.name) + " "  + str(self.strength) + " "  + str(self.node))

        def __repr__(self):
                return(str(self.name) + " " + str(self.strength) + " "  + str(self.node))

class LifeTime():
        def __init__(self, outputDir):
                
                self.outputDir = outputDir
                                
                self.filePath = os.path.join(outputDir, "process.h5")
                self.file = h5py.File(self.filePath, 'r')
                self.liqInfo = np.array(self.file["liqInfo"][:,:,0], dtype=int)
                self.liqCarac = np.array(self.file["liqInfo"][:,:,1:])
                self.liqtest = np.array(self.file["liqInfo"])
                self.dropInfo = np.array(self.file["liqDropInfo"])
                
                self.N, self.nLiq = np.shape(self.liqInfo)
                self.dropSizeMax = np.max(self.dropInfo[:,0])
                self.dropSizeMin = np.min(np.where(self.dropInfo[:,0] == 0, 2*self.nLiq, self.dropInfo[:,0]))
                
                self.G = nx.Graph()
                
                self.droplet_dict = {}
                self.droplet_id = 0
                self.node_neighbors_dict = {}


        def Compute_graph(self):
                self.G.add_node((0, int(-1)), size = self.nLiq, droplet = -1)

                for frame in range(1, self.N):

                        nLiq_gaz = int((self.liqInfo[frame] < -.5).sum())
                        
                        if nLiq_gaz > .5:
                                cm_x, cm_y, cm_z, rho_loc_mean = np.mean(self.liqCarac[frame][self.liqInfo[frame] < -.5], axis = 0)
                                self.G.add_node((frame, int(-1)), size = nLiq_gaz, droplet = -1, center_of_mass_x = float(cm_x), center_of_mass_y = float(cm_y), center_of_mass_z = float(cm_z), local_density = float(rho_loc_mean))
                                
                        cluster_id = 0
                        
                        while self.dropInfo[frame][cluster_id][0] > 0:
                                # print(self.liqInfo[frame])
                                # print(cluster_id)
                                # print(self.liqInfo[frame] == cluster_id)
                                # print(self.liqCarac[frame])
                                # print(self.liqCarac[frame][self.liqInfo[frame] == cluster_id])
                                # print(np.mean(self.liqCarac[frame][self.liqInfo[frame] == cluster_id], axis = 0))
                                # raise TypeError("test")
                                cm_x, cm_y, cm_z, rho_loc_mean = np.mean(self.liqCarac[frame][self.liqInfo[frame] == cluster_id], axis = 0)
                                self.G.add_node((frame, int(cluster_id)), size = int(self.dropInfo[frame][cluster_id][0]), droplet = -1, in_flux = 0, out_flux = 0, r_gyr = self.dropInfo[frame][cluster_id][1], aniso = self.dropInfo[frame][cluster_id][2], center_of_mass = (cm_x, cm_y, cm_z), local_density = rho_loc_mean)
                                cluster_id += 1
                        
                        if frame > 1:
                                for particule in range(self.nLiq):
                                        try:
                                                self.G.edges[(frame-1, int(self.liqInfo[frame-1][particule])),(frame, int(self.liqInfo[frame][particule]))]["weight"] += 1
                                        except:
                                                self.G.add_edge((frame-1, int(self.liqInfo[frame-1][particule])),(frame, int(self.liqInfo[frame][particule])), event = "")
                                                self.G.edges[(frame-1, int(self.liqInfo[frame-1][particule])),(frame, int(self.liqInfo[frame][particule]))]["weight"] = 1
                                        if int(self.liqInfo[frame-1][particule]) < -.5 and int(self.liqInfo[frame][particule]) > -.5:
                                                self.G.nodes[(frame,int(self.liqInfo[frame][particule]))]["in_flux"] += 1
                                        if int(self.liqInfo[frame-1][particule]) > -.5 and int(self.liqInfo[frame][particule]) < -.5:
                                                self.G.nodes[(frame-1,int(self.liqInfo[frame-1][particule]))]["out_flux"] += 1

                        elif frame > 0:
                                for particule in range(self.nLiq):
                                        try:
                                                self.G.edges[(0, -1),(frame, int(self.liqInfo[frame][particule]))]["weight"] += 1
                                        except:
                                                self.G.add_edge((0, -1),(frame, int(self.liqInfo[frame][particule])))
                                                self.G.edges[(0, -1),(frame, int(self.liqInfo[frame][particule]))]["weight"] = 1
                                        if int(self.liqInfo[frame][particule]) > -.5:
                                                self.G.nodes[(frame,int(self.liqInfo[frame][particule]))]["in_flux"] += 1

                self.dict_max_give = {(frame, cluster) : max({(neigh_frame, neigh_cluster) : weight if neigh_frame > frame else -1 for _, (neigh_frame, neigh_cluster), weight in self.G.edges((frame, cluster), data = "weight")}.items(), key=lambda k: k[1])[0] if frame > -.5 else (frame, cluster) for frame, cluster in self.G.nodes()}
                self.dict_max_receive = {(frame, cluster) : max({(neigh_frame, neigh_cluster) : weight if neigh_frame < frame else -1 for _, (neigh_frame, neigh_cluster), weight in self.G.edges((frame, cluster), data = "weight")}.items(), key=lambda k: k[1])[0] if frame < self.N else (frame, cluster) for frame, cluster in self.G.nodes()}

                self.dict_edges_per_frame_and_type = {(frame, frame+1) : {0 : [], 'A' : [], 'B' : [], 'C' : []} for frame in range(self.N)}
                
                for (u,v) in self.G.edges:
                        if v[0] < u[0]:
                                u, v = v, u
                        
                        if (u[1] > -.5 and v[1] > -.5 and self.dict_max_give[u] == v and self.dict_max_receive[v] == u) or (self.dict_max_receive[v] == u and u[1] < -.5 and v[1] > -.5):
                                self.dict_edges_per_frame_and_type[(u[0],v[0])][0].append((u,v))
                                if (self.dict_max_receive[v] == u and u[1] < -.5 and v[1] > -.5):
                                        self.dict_edges_per_frame_and_type[(u[0],v[0])]["A"].append((u,v))

                        elif self.dict_max_receive[v] == u and v[1] > -.5:
                                self.dict_edges_per_frame_and_type[(u[0],v[0])]["A"].append((u,v))
                                
                        elif self.dict_max_give[u] == v and u[1] > -.5:
                                self.dict_edges_per_frame_and_type[(u[0],v[0])]["B"].append((u,v))

                        elif u[1] > -.5 and v[1] > -.5:
                                self.dict_edges_per_frame_and_type[(u[0],v[0])]["C"].append((u,v))
                                                       
                        else:
                                pass

                # print('0', self.dict_edges_per_frame_and_type[(3,4)][0], '\n')
                # print('A', self.dict_edges_per_frame_and_type[(3,4)]['A'], '\n')
                # print('B', self.dict_edges_per_frame_and_type[(3,4)]['B'], '\n')
                # print('C', self.dict_edges_per_frame_and_type[(3,4)]['C'], '\n')

                # if self.dropSizeMax-self.dropSizeMin < 2:
                #         self.node_size = np.array([self.G.nodes[(frame,cluster)]["size"]/self.nLiq for frame, cluster in self.G.nodes()])
                # else:
                #         self.node_size = np.array([(self.G.nodes[(frame,cluster)]["size"]-self.dropSizeMin)/(self.dropSizeMax-self.dropSizeMin) if cluster != -1 else self.G.nodes[(frame,cluster)]["size"]/self.nLiq for frame, cluster in self.G.nodes()])
                # self.edge_size = np.array([self.G.edges[(u, v)]["weight"]/self.G.nodes[u]["size"] if u[1] != -1 else self.G.edges[(u, v)]["weight"]/self.G.nodes[v]["size"] for u, v in self.G.edges()])
                # self.edge_size = (self.edge_size - np.min(self.edge_size))/(np.max(self.edge_size) - np.min(self.edge_size))

        def Compute_droplet(self):
                # print(self.G.edges((4,0), data=True))
                for node in self.G.nodes():
                        if node[1] > -.5:
                                self.Compute_node_neighbors(node)
                


                for frame in range(self.N):
                        for edges in self.dict_edges_per_frame_and_type[(frame, frame+1)][0]:
                                self.Compute_edge_0(edges)
                                
                        for edges in self.dict_edges_per_frame_and_type[(frame, frame+1)]['A']:
                                self.Compute_edge_A(edges)
                                
                        for edges in self.dict_edges_per_frame_and_type[(frame, frame+1)]['B']:
                                self.Compute_edge_B(edges)
                                
                        for edges in self.dict_edges_per_frame_and_type[(frame, frame+1)]['C']:
                                self.Compute_edge_C(edges)

                self.Check_droplet()

        def Check_droplet(self):
                for node in self.G.nodes():
                        if  node[1] > -.5 and node[0] < self.N - 1 and self.G.nodes[node]["droplet"] == -1:
                                raise NameError(f"node {node} as no droplet associated")
                for edge in self.G.edges():
                        if (edge[0][1] > -0.5 or edge[1][1] > -0.5) and self.G.edges[edge]["event"] == "":
                                raise NameError(f"edge {edge} as no event associated")
                for droplet in self.droplet_dict.values():
                        for i in range(droplet.tau):
                                if droplet.nodes[i][0] != droplet.frame_start + i:
                                        raise NameError(f"node {droplet.nodes[i]}'s frame in droplet {droplet.droplet_id} is not consistant")

                

        def Compute_node_neighbors(self, node):
                list_in_node_A = []
                list_in_node_B = []
                list_out_node_A = []
                list_out_node_B = []

                for _,v in self.G.edges(node):
                        if node[0] > v[0] and v[1] > -.5:
                                if self.dict_max_give[v] == node:
                                        list_in_node_B.append(v)
                                else:
                                        list_in_node_A.append(v)
                        if node[0] < v[0] and v[1] > -.5:
                                if self.dict_max_receive[v] == node:
                                        list_out_node_B.append(v)
                                else:
                                        list_out_node_A.append(v)

                self.node_neighbors_dict[node] = {"in_A" : list_in_node_A, "in_B" : list_in_node_B, "out_A" : list_out_node_A, "out_B" : list_out_node_B}

               
        def Compute_edge_0(self, edge):

                u, v = edge
                list_in_v_B = self.node_neighbors_dict[v]["in_B"]

                if u[1] < -.5 and list_in_v_B == []: # Cas 1 : Emergence
                        self.droplet_dict[self.droplet_id] = Droplet(self.droplet_id, v, self.G.nodes[v]["size"], [Event('emergence', self.G.nodes[v]["in_flux"], None)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["local_density"])
                        self.G.nodes[v]["droplet"] = self.droplet_id
                        self.droplet_id += 1

                        self.G.edges[u, v]["event"] = 'emergence'

                elif u[1] > -.5: # Cas persistance standart
                        self.droplet_dict[self.G.nodes[u]["droplet"]].append(v, self.G.nodes[v]["size"], [Event('persistance', self.G.edges[v, u]["weight"], u )], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["local_density"])
                        self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('persistance', self.G.edges[v, u]["weight"], v))
                        self.G.nodes[v]["droplet"] = self.G.nodes[u]["droplet"]
                        self.G.edges[u, v]["event"] = 'persistance'

                        if self.G.nodes[u]["out_flux"] > 0:
                                self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('flux', self.G.nodes[u]["out_flux"], None))
                                self.G.edges[u, (v[0], -1)]["event"] = 'flux'

                        if self.G.nodes[v]["in_flux"] > 0:
                                self.droplet_dict[self.G.nodes[v]["droplet"]].append_in_event(v, Event('flux', self.G.nodes[v]["in_flux"], None))
                                self.G.edges[(v[0]-1, -1), v]["event"] = 'flux'

                        for in_v_B in list_in_v_B:
                                if in_v_B != u:
                                        in_v_B_id = self.G.nodes[in_v_B]["droplet"]
                                        self.droplet_dict[self.G.nodes[v]["droplet"]].append_in_event(v, Event('merge', self.G.edges[v, in_v_B]["weight"], in_v_B))
                                        self.droplet_dict[in_v_B_id].append_out_event(in_v_B, Event('merge', self.G.edges[v, in_v_B]["weight"], v))
                                        self.G.edges[v, in_v_B]["event"] = 'merge'

                                        if self.G.nodes[in_v_B]["out_flux"] > 0:
                                                self.droplet_dict[in_v_B_id].append_out_event(in_v_B, Event('flux', self.G.nodes[in_v_B]["out_flux"], None))
                                                self.G.edges[in_v_B, (v[0], -1)]["event"] = 'flux'
                    

        def Compute_edge_A(self, edge):

                u, v = edge
                
                list_in_v_B = self.node_neighbors_dict[v]["in_B"][:]
                list_in_v_A = self.node_neighbors_dict[v]["in_A"]

                for in_v_A in list_in_v_A:
                        if self.dict_max_give[in_v_A][1] < -.5:
                                if self.G.nodes[in_v_A]["droplet"] < -.5:
                                        list_in_v_B.append(in_v_A)
                                elif not(self.droplet_dict[self.G.nodes[in_v_A]["droplet"]].is_persistante_at_node(in_v_A)):
                                        list_in_v_B.append(in_v_A)

                if self.G.edges[edge]['event'] == "":

                        if u[1] < -.5 and list_in_v_B != []: # Cas 2 : persistance en cas de expand
                        
                                in_v_B_max = max({in_v_B : self.G.edges[(in_v_B, v)]["weight"] for in_v_B in list_in_v_B}.items(), key=lambda k: k[1])[0]
                                self.droplet_dict[self.G.nodes[in_v_B_max]["droplet"]].append(v, self.G.nodes[v]["size"], [Event('persistance', self.G.edges[v, in_v_B_max]["weight"], in_v_B_max)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["local_density"])
                                self.droplet_dict[self.G.nodes[in_v_B_max]["droplet"]].append_out_event(in_v_B_max, Event('persistance', self.G.edges[v, in_v_B_max]["weight"], v))
                                self.G.nodes[v]["droplet"] = self.G.nodes[in_v_B_max]["droplet"]
                                self.G.edges[in_v_B_max, v]["event"] = 'persistance'
                        
                                self.droplet_dict[self.G.nodes[in_v_B_max]["droplet"]].append_in_event(v, Event('expand', self.G.nodes[v]["in_flux"], None))
                                self.G.edges[u, v]["event"] = 'expand'

                                if self.G.nodes[in_v_B_max]["out_flux"] > 0:
                                        self.droplet_dict[self.G.nodes[in_v_B_max]["droplet"]].append_out_event(in_v_B_max, Event('flux', self.G.nodes[in_v_B_max]["out_flux"], None))
                                        self.G.edges[in_v_B_max, (v[0], -1)]["event"] = 'flux'

                                for in_v_B in list_in_v_B:
                                        if in_v_B != in_v_B_max and not(in_v_B in list_in_v_A):
                                                in_v_B_id = self.G.nodes[in_v_B]["droplet"]

                                                self.droplet_dict[self.G.nodes[v]["droplet"]].append_in_event(v, Event('amalgamate', self.G.edges[v, in_v_B]["weight"], in_v_B))
                                                self.droplet_dict[in_v_B_id].append_out_event(in_v_B, Event('amalgamate', self.G.edges[v, in_v_B]["weight"], v))
                                                self.G.edges[in_v_B, v]["event"] = 'amalgamate'

                                                if self.G.nodes[in_v_B]["out_flux"] > 0:
                                                        self.droplet_dict[in_v_B_id].append_out_event(in_v_B, Event('flux', self.G.nodes[in_v_B]["out_flux"], None))
                                                        self.G.edges[in_v_B, (v[0], -1)]["event"] = 'flux'

                        else :
                                if self.dict_max_give[u][1] < -.5: # g donne son maximum de particules au flux donc persistence ou splinter
                                       
                                        if self.droplet_dict[self.G.nodes[u]["droplet"]].is_persistante_at_node(u):
                                                out_node_B_max = None
                                        else:
                                                list_out_node_B = self.node_neighbors_dict[u]["out_B"]
                                                out_node_B_max = max({out_node_B : self.G.nodes[out_node_B]["size"] for out_node_B in list_out_node_B}.items(), key=lambda k: k[1])[0]
                                        if out_node_B_max == v:         # Cas de persistence

                                                self.droplet_dict[self.G.nodes[u]["droplet"]].append(v, self.G.nodes[v]["size"], [Event('persistance', self.G.edges[v, u]["weight"], u)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["local_density"])
                                                self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('persistance', self.G.edges[v, u]["weight"], v))
                                                self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('shrink', self.G.nodes[u]["out_flux"], None))
                                                self.G.nodes[v]["droplet"] = self.G.nodes[u]["droplet"]
                                                
                                                self.G.edges[u, v]["event"] = 'persistance'
                                                if self.G.nodes[v]["in_flux"] > 0:
                                                        self.droplet_dict[self.G.nodes[v]["droplet"]].append_in_event(v, Event('flux', self.G.nodes[v]["in_flux"], None))
                                                        self.G.edges[(v[0]-1, -1), v]["event"] = 'flux'
                                                self.G.edges[u, (v[0], -1)]["event"] = 'shrink'
                                                
                                        else:                           # Cas de splinter
                                                self.droplet_dict[self.droplet_id] = Droplet(self.droplet_id, v, self.G.nodes[v]["size"], [Event('splinter', self.G.edges[v, u]["weight"], u)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["local_density"])
                                                self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('splinter', self.G.edges[v, u]["weight"], v))
                                                self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('shrink', self.G.nodes[u]["out_flux"], None))
                                                self.G.nodes[v]["droplet"] = self.droplet_id
                                                self.droplet_id += 1

                                                self.G.edges[u, v]["event"] = 'splinter'
                                                if self.G.nodes[v]["in_flux"] > 0:
                                                        self.droplet_dict[self.G.nodes[v]["droplet"]].append_in_event(v, Event('flux', self.G.nodes[v]["in_flux"], None))
                                                        self.G.edges[(v[0]-1, -1), v]["event"] = 'flux'
                                                self.G.edges[u, (v[0], -1)]["event"] = 'shrink'

                                elif u == self.dict_max_receive[self.dict_max_give[u]]: # Cas de split
                                        self.droplet_dict[self.droplet_id] = Droplet(self.droplet_id, v, self.G.nodes[v]["size"], [Event('split', self.G.edges[v, u]["weight"], u)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["local_density"])
                                        self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('split', self.G.edges[v, u]["weight"], v))
                                        self.G.nodes[v]["droplet"] = self.droplet_id
                                        self.droplet_id += 1


                                        self.G.edges[u, v]["event"] = 'split'

                                        if self.G.nodes[v]["in_flux"] > 0:
                                                self.droplet_dict[self.G.nodes[v]["droplet"]].append_in_event(v, Event('flux', self.G.nodes[v]["in_flux"], None))
                                                self.G.edges[(v[0]-1, -1), v]["event"] = 'flux'
                                        if self.G.nodes[u]["out_flux"] > 0:
                                                self.droplet_dict[self.G.nodes[u]["droplet"]].append_in_event(u, Event('flux', self.G.nodes[u]["out_flux"], None))
                                                self.G.edges[u, (v[0], -1)]["event"] = 'flux'


                                else: # Cas de offload

                                        self.droplet_dict[self.droplet_id] = Droplet(self.droplet_id, v, self.G.nodes[v]["size"], [Event('offload', self.G.edges[v, u]["weight"], u)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["local_density"])
                                        self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('offload', self.G.edges[v, u]["weight"], v))
                                        self.G.nodes[v]["droplet"] = self.droplet_id
                                        
                                        self.droplet_id += 1

                                        self.G.edges[u, v]["event"] = 'offload'

                                        if self.G.nodes[v]["in_flux"] > 0:
                                                self.droplet_dict[self.G.nodes[v]["droplet"]].append_in_event(v, Event('flux', self.G.nodes[v]["in_flux"], None))
                                                self.G.edges[(v[0]-1, -1), v]["event"] = 'flux'
                                        if self.G.nodes[u]["out_flux"] > 0:
                                                self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('flux', self.G.nodes[u]["out_flux"], None))
                                                self.G.edges[u, (v[0], -1)]["event"] = 'flux'


                        list_in_v_B = self.node_neighbors_dict[v]["in_B"]
                                          

                for in_v_B in self.node_neighbors_dict[v]["in_B"]: # Cas de blend suite à un split quelconque
                        if in_v_B != u and self.G.edges[v, in_v_B]["event"] == "":
                                in_v_B_id = self.G.nodes[in_v_B]["droplet"]
                                self.droplet_dict[self.G.nodes[v]["droplet"]].append_in_event(v, Event('blend', self.G.edges[v, in_v_B]["weight"], in_v_B))
                                self.droplet_dict[in_v_B_id].append_out_event(in_v_B, Event('blend', self.G.edges[v, in_v_B]["weight"], v))
                                self.G.edges[v, in_v_B]["event"] = 'blend'

                                if self.G.nodes[in_v_B]["out_flux"] > 0:
                                        self.droplet_dict[in_v_B_id].append_out_event(in_v_B, Event('flux', self.G.nodes[in_v_B]["out_flux"], None))
                                        self.G.edges[in_v_B, (v[0], -1)]["event"] = 'flux'



        def Compute_edge_B(self, edge):  # On s'interesse maintenant à node comme à un g
                
                u, v = edge
                
                list_out_node_B = self.node_neighbors_dict[u]['out_B']

                if v[1] < -.5: # Cas 2 déjà en partie traité
                        if list_out_node_B == []: # Cas non traité d'évaporation

                                self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('evaporation', self.G.nodes[u]["out_flux"], None))
                                self.G.edges[u, v]["event"] = 'evaporation'

                        elif self.G.edges[u, v]["event"] == "":
                                print("\n -- \n")
                                print(u, v, self.G.nodes[u], self.G.nodes[v], self.G.edges[(u,v)])

                elif self.G.edges[u, v]["event"] == "":
                        print("\n -- -- \n")
                        print(u, v, self.G.nodes[u]["size"], self.G.nodes[u]["droplet"], self.G.nodes[v]["size"], self.G.nodes[v]["droplet"], self.dict_max_give[u], self.dict_max_receive[v], self.G.edges[(u,v)])


        def Compute_edge_C(self, edge):
                u, v = edge
                
                if self.G.edges[edge]["event"] == "":

                        u, v = edge

                        droplet_u = self.droplet_dict[self.G.nodes[u]["droplet"]]
                        droplet_v = self.droplet_dict[self.G.nodes[v]["droplet"]]
                        
                        bool_u = droplet_u.is_persistante_at_node(u)
                        bool_v = droplet_v.was_persistante_at_node(v)

                        if bool_u and bool_v:
                                droplet_u.append_out_event(u, Event("exchange", self.G.edges[edge]['weight'], v))
                                droplet_v.append_in_event(v, Event("exchange", self.G.edges[edge]['weight'], u))
                                self.G.edges[edge]["event"] = "exchange"
                        elif bool_u:
                                droplet_u.append_out_event(u, Event("coarsening", self.G.edges[edge]['weight'], v))
                                droplet_v.append_in_event(v, Event("coarsening", self.G.edges[edge]['weight'], u))
                                self.G.edges[edge]["event"] = "coarsening"
                        elif bool_v:
                                droplet_u.append_out_event(u, Event("ripening", self.G.edges[edge]['weight'], v))
                                droplet_v.append_in_event(v, Event("ripening", self.G.edges[edge]['weight'], u))
                                self.G.edges[edge]["event"] = "ripening"
                        else:
                                droplet_u.append_out_event(u, Event("flux", self.G.edges[edge]['weight'], v))
                                droplet_v.append_in_event(v, Event("flux", self.G.edges[edge]['weight'], u))
                                self.G.edges[edge]["event"] = "flux"
                
        
      

        def Compute_flux(self, cluster):

                self.p_flux_dict = {node : 0 for node in self.G.nodes}
                for v,w in self.G.edges:
                        if v[1] == cluster:
                                self.p_flux_dict[w] += (-1)**(v[0] > w[0]) * self.G.edges[v,w]["weight"]
                        elif w[1] == cluster:
                                self.p_flux_dict[v] += (-1)**(v[0] < w[0]) * self.G.edges[v,w]["weight"]

                self.flux_min = min([flux if cluster > -.5 and frame < self.N-1 and frame > 1.5 else 2*self.nLiq for (frame, cluster), flux in self.p_flux_dict.items()])
                self.flux_max = max([flux if cluster > -.5 and frame < self.N-1 and frame > 1.5 else -1 for (frame, cluster), flux in self.p_flux_dict.items()])

                self.node_flux = np.array([flux/max(abs(self.flux_min),abs(self.flux_max),1) if frame > 1.5 and frame < self.N-1 and cluster > -.5 else 0 for (frame, cluster), flux in self.p_flux_dict.items()])


                fig = plt.figure(figsize=(24,16))
                ax = fig.add_subplot()
                # network = nx.draw_networkx(self.G,pos = {(frame, cluster) : (frame, cluster) for frame, cluster in self.G.nodes()},ax = ax, node_size = 500, width = 5, node_color = self.node_flux, edge_color = self.edge_size, cmap = BrBG, edge_cmap = viridis, vmin = -1, vmax = 1, edge_vmin = -0.1, edge_vmax = 1.1, with_labels=False, linewidths=1, edgecolors="black")

                # cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=-1, vmax=1, clip=False), cmap=BrBG),  ax=ax, shrink = 0.7)

                sub_ax1 = inset_axes(
                        parent_axes=ax,
                        width="40%",
                        height="40%",
                        borderpad=1,
                        loc = "upper right")

                sub_ax1.scatter([flux if cluster > -.5 and frame < self.N - 1 and frame > 1.5 else 0 for (frame, cluster), flux in self.p_flux_dict.items()], [self.G.nodes[(frame, cluster)]["size"] if cluster > -.5 and frame > 1.5 and frame < self.N - 1 else 0 for (frame, cluster), _ in self.p_flux_dict.items()])
                sub_ax1.set_xlim((-max(abs(self.flux_min),abs(self.flux_max),1)-1, max(abs(self.flux_min),abs(self.flux_max),1)+1))
                sub_ax1.set_ylim((self.dropSizeMin-1, self.dropSizeMax+1))


                fig.savefig(f"/home/ppuel/data/flux_cloud_{self.filePath.split('data/')[1].split('/process.h5')[0].replace('/','_')}.png")



        def Print_film(self):
                self.reader = hdf5Reader(self.outputDir, "traj.h5", -1, readLiq=True, readPoly=False, backInBox=False)

                os.chdir(self.outputDir)
                exp = self.outputDir.split("data/")[1].split("/")[0]
                params = "_".join(self.outputDir.split("data/")[1].split("/N/")[0].split("/")[1:])
                self.outputPath = f"/home/ppuel/data/{exp}/traj_vtk"
                os.makedirs(self.outputPath, exist_ok=True)


                self.PrintBox()


                for i in range(self.reader.N):
                        data = next(self.reader)
                        self.PrintLiqFrame(data, i)

                        if (i+1) % 100 == 0:
                                print("Printed %d out of %d configurations" % (i+1, self.reader.N))



        def PrintLiqFrame(self, data, i):

                fileLiqName = 'liq{:05d}.vtp'.format(i-1)
                fileLiqPath = os.path.join(self.outputPath,fileLiqName)

                points = vtkPoints()
                liqDensity = vtkFloatArray()
                liqDisplacement = vtkFloatArray()
                liqDroplet = vtkFloatArray()
                liqLifeTime = vtkFloatArray()

                liqDensity.SetName("Density")
                liqDensity.SetNumberOfComponents(1)

                liqDisplacement.SetName("Displacement")
                liqDisplacement.SetNumberOfComponents(3)

                liqDroplet.SetName("Droplet")
                liqDroplet.SetNumberOfComponents(1)

                liqLifeTime.SetName("LifeTime")
                liqLifeTime.SetNumberOfComponents(1)

                for j in range(self.reader.nLiq):

                        aveDensity = data.liqDens[j]

                        x = (data.liqPos[j][0]-12)%24
                        y = (data.liqPos[j][1]-12)%24
                        z = data.liqPos[j][2]

                        dx = data.liqDisp[j][0]
                        dy = data.liqDisp[j][1]
                        dz = data.liqDisp[j][2]

                        points.InsertNextPoint(x, y, z)

                        liqDensity.InsertNextValue(aveDensity)
                        liqDisplacement.InsertNextTuple3(dx, dy, dz)

                        if self.liqInfo[i][j] > -.5 and i > .5:
                                # print('node : ', self.liqInfo[i][j],' | droplet : ', self.G.nodes[(i, int(self.liqInfo[i][j]))]["droplet"])
                                if self.G.nodes[(i, self.liqInfo[i][j])]["droplet"] != None and self.G.nodes[(i, self.liqInfo[i][j])]["droplet"] > -.5:
                                        # print('node : ', self.liqInfo[i][j], " ", i,' | droplet : ', self.G.nodes[(i, self.liqInfo[i][j])]["droplet"], ' | droplet : ', self.G.nodes[(i, self.liqInfo[i][j])]["droplet"])
                                        liqDroplet.InsertNextValue(self.G.nodes[(i, self.liqInfo[i][j])]["droplet"])
                                        liqLifeTime.InsertNextValue(self.droplet_dict[self.G.nodes[(i, self.liqInfo[i][j])]["droplet"]].tau)
                                else:
                                        liqDroplet.InsertNextValue(-2)
                                        liqLifeTime.InsertNextValue(-2)
                        else:
                                liqDroplet.InsertNextValue(-1)
                                liqLifeTime.InsertNextValue(-1)

                polyData = vtkPolyData()
                writer = vtkXMLPolyDataWriter()

                polyData.SetPoints(points)

                polyData.GetPointData().AddArray(liqDensity)
                polyData.GetPointData().AddArray(liqDisplacement)
                polyData.GetPointData().AddArray(liqDroplet)
                polyData.GetPointData().AddArray(liqLifeTime)

                writer.SetFileName(fileLiqPath)
                writer.SetInputData(polyData)

                writer.Write()

        def PrintBox(self):
                fileBoxName = 'box.vtp'
                fileBoxPath = os.path.join(self.outputPath,fileBoxName)

                cubeSource = vtkCubeSource()

                L = self.reader.boxDim[0]

                cubeSource.SetCenter((L-0.5)/2., (L-0.5)/2., (L-0.5)/2.)

                cubeSource.SetXLength(L+0.5)
                cubeSource.SetYLength(L+0.5)
                cubeSource.SetZLength(L+0.5)

                cubeSource.Update()

                writer = vtkXMLPolyDataWriter()

                writer.SetFileName(fileBoxPath)
                writer.SetInputConnection(cubeSource.GetOutputPort())

                writer.Write()




        def Print(self):

                font = font_manager.FontProperties(     weight='bold',
                                                        style='normal', size=60)
                fontlabel = {"labelsize" : 60}

                
                fig = plt.figure(figsize=(32,24))
                ax = fig.add_subplot()
                


                max_size = -1
                for droplet in self.droplet_dict.values():
                        max_size = max(max_size, np.max(droplet.sizes))

                rainbow = mpl.colormaps['plasma'].resampled(max_size)
                for droplet in self.droplet_dict.values():
                        if droplet.tau > 3 and max(droplet.sizes) > 4:
                        # volume = np.array(droplet.sizes)*20e-9**3/2**(1/2)

                        # ax.scatter([droplet.nodes[0][0]+ i for i in range(len(droplet.sizes))], droplet.sizes, color = rainbow(np.mean(droplet.local_density)), s = 100)
                                ax.plot([droplet.nodes[0][0] + i for i in range(len(droplet.sizes))], [dens*12 for dens in droplet.local_density], color = rainbow(np.max(droplet.sizes)), linewidth=8)

                exp = self.outputDir.split("data/")[1].split("/")[0]
                params = "_".join(self.outputDir.split("data/")[1].split("/N/")[0].split("/")[1:])

                ax.set_title("Film of the growth of the droplet", font = font)
                ax.set_ylabel("Size Droplet (Number of PRC1)", font = font)
                ax.set_xlabel("Time (kMCS)", font = font)
                
                ax.tick_params(axis="both", **fontlabel)
                

                # stats = (f'C_puncta mean = {mean_C_puncta:.2f}\nC_total mean = {mean_C_total:.2f}')
                # bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
                # ax.text(0.95, 0.95, stats, bbox=bbox,
                #         transform=ax.transAxes, horizontalalignment='right', verticalalignment='top', font = font)



                fig.savefig(f"/home/ppuel/data/{exp}/LifeTime_rho_loc_{params}.png")

        def Print_data(self):
                # for ids, droplet in self.droplet_dict.items():
                #         if ids == 0: 
                #                 print(ids, droplet.nodes[:])
                with open(os.path.join(self.outputDir, "liq_droplets.pickle"), "wb") as f:
                        pickle.dump(self.droplet_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
                        
                nx.to_graph6_bytes(self.G, os.path.join(self.outputDir, "liq_graph.graph6"))
                

if __name__ == "__main__":
        if len(sys.argv) != 2:
                print("\033[1;31mUsage is %s outputDir\033[0m" % sys.argv[0])
                sys.exit()

        outputDir = sys.argv[1]


        A = time.time()
        lifeTime = LifeTime(outputDir)
        B = time.time()
        print("init time : ",np.round((B-A),4))
        lifeTime.Compute_graph()
        C = time.time()
        print("graph time : ",np.round((C-B),4))
        lifeTime.Compute_droplet()
        D = time.time()
        print("droplet time : ",np.round((D-C),4))
        lifeTime.Print_data()
        E = time.time()
        print("data time : ",np.round((E-D),4))
        

# color_out_event = {'evaporation' : 'brown', "merge" : 'darkred', 'amalgamate' : 'red', 'blend' : 'pink'}
                # color_in_event = {'emergence' : 'purple', "split" : 'darkblue', "splintering" : 'blue', 'offload' : 'cyan' }
                

                        # in_color = 'black'
                        # out_color = 'black'


                        # for in_event in droplet.in_event[0]:
                        #         if in_event.name in color_in_event.keys():
                        #                 in_color = color_in_event[in_event.name]

                        # for out_event in droplet.out_event[-1]:
                        #         if out_event.name in color_out_event.keys():
                        #                 out_color = color_out_event[out_event.name]


                        # ax.plot([droplet.nodes[0][0]-0.2, droplet.nodes[0][0]-0.2], [droplet.sizes[0] - self.dropSizeMax*0.1, droplet.sizes[0] + self.dropSizeMax*0.1], color = in_color)
                        # ax.plot([droplet.nodes[-1][0]+0.2, droplet.nodes[-1][0]+0.2], [droplet.sizes[-1] - self.dropSizeMax*0.1, droplet.sizes[-1] + self.dropSizeMax*0.1], color = out_color)



        # for event_list in droplet.in_event:
                        #         for event in event_list:
                        #                 try :
                        #                         print(event.name, event.strength)
                        #                 except :
                        #                         print(event)
                        # for event_list in droplet.out_event:
                        #         for event in event_list:
                        #                 try :
                        #                         print(event.name, event.strength)
                        #                 except :
                        #                         print(event)

                # network = nx.draw_networkx(self.G,pos = {(frame, cluster) : (frame, cluster) for frame, cluster in self.G.nodes()},ax = ax, node_size = 500, width = 5, node_color = self.node_flux, edge_color = self.edge_size, cmap = BrBG, edge_cmap = viridis, vmin = -1, vmax = 1, edge_vmin = -0.1, edge_vmax = 1.1, with_labels=False, linewidths=1, edgecolors="black")

                # cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=-1, vmax=1, clip=False), cmap=BrBG),  ax=ax, shrink = 0.7)

                # sub_ax1 = inset_axes(
                #         parent_axes=ax,
                #         width="40%",
                #         height="40%",
                #         borderpad=1,
                #         loc = "upper right")

                # sub_ax1.scatter([flux if cluster > .5 and frame < self.N - 1 and frame > 1.5 else 0 for (frame, cluster), flux in self.p_flux_dict.items()], [self.G.nodes[(frame, cluster)]["size"] if cluster > -.5 and frame > 1.5 and frame < self.N - 1 else 0 for (frame, cluster), _ in self.p_flux_dict.items()])
                # sub_ax1.set_xlim((-max(abs(self.flux_min),abs(self.flux_max),1)-1, max(abs(self.flux_min),abs(self.flux_max),1)+1))
                # sub_ax1.set_ylim((self.dropSizeMin-1, self.dropSizeMax+1))
