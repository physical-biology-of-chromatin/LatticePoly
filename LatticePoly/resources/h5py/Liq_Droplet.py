##
##  LiqCluster.py
##  LatticePoly
##
##  Created by ppuel on 19/02/2025.
##  Copyright © 2020 ENS Lyon. All rights reserved.
##

import os
import sys
import h5py
import pickle

import numpy as np
import networkx as nx

from Reader import Reader


# from vtk import vtkXMLPolyDataWriter, vtkPoints, vtkFloatArray, vtkPolyData, vtkCellArray, vtkIntArray, vtkLine, vtkCubeSource

class SDroplet():
    def __init__(self, droplet_id, tau, nodes, sizes, center_of_mass, mean_number_of_neighbors):
        self.droplet_id = droplet_id
        self.tau = tau
        self.sizes = sizes
        self.center_of_mass = center_of_mass
        self.mean_number_of_neighbors = mean_number_of_neighbors
        self.frame_start = nodes[0]

class Droplet():
    def __init__(self, droplet_id, node, size, in_event, center_of_mass, mean_number_of_neighbors):
        self.droplet_id = droplet_id
        self.nodes = [node]
        self.tau = 1
        self.sizes = [size]
        self.in_event = [in_event]
        self.out_event = [[]]
        self.center_of_mass = [center_of_mass]
        self.mean_number_of_neighbors = [mean_number_of_neighbors]
        self.frame_start = node[0]

    def append(self, node, size, in_event, center_of_mass, mean_number_of_neighbors):
        self.nodes.append(node)
        self.tau += 1
        self.sizes.append(size)
        self.in_event.append(in_event)
        self.out_event.append([])
        self.center_of_mass.append(center_of_mass)
        self.mean_number_of_neighbors.append(mean_number_of_neighbors)
        
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
        return(str(self.name) + " " + str(self.strength) + " " + str(self.node))

    def __repr__(self):
        return(str(self.name) + " " + str(self.strength) + " " + str(self.node))

class LifeTime():
    def __init__(self, input_dir):
        print(f"LifeTime : Init {input_dir}")

        self.input_dir = input_dir     
        self.process_path = os.path.join(input_dir, "process.h5")        
        self.droplet_dict = {}
        self.droplet_id = 0
        self.node_neighbors_dict = {}
        self.simple_droplet_dict = {}
        self.G = nx.Graph()

        ### Get data
        
        # reader = Reader(input_dir)
        # self.box_dim = reader.box_dim
        # del reader
        
        with h5py.File(self.process_path, 'r') as processFile:
               
            self.liq_info = np.asarray(processFile["liq_info"], dtype = np.int32)
            self.liq_drop_info = np.array(processFile["liq_drop_info"])

        self.n_frame, self.n_liq = np.shape(self.liq_info)
        self.dict_edges_per_frame_and_type = {(frame, frame+1) : {0 : [], 'A' : [], 'B' : [], 'C' : []} for frame in range(self.n_frame)}
        
        # self.dropSizeMax = np.max(self.liq_drop_info['size'])
        # self.dropSizeMin = np.min(np.where(
        #     self.liq_drop_info['size'] > 0,
        #     self.liq_drop_info['size'],
        #     2*self.n_liq
        # ))
        
            
    def compute_graph(self):
        print("\n")

        self.G.add_node((0, int(-1)), size = self.n_liq)

        for frame in range(1, self.n_frame):
            if (frame) % 10 == 0:
                print("LifeTime[Compute_graph] : Process %d out of %d frames" % (frame, self.n_frame))
           
            n_liq_gaz = int((self.liq_info[frame] < -.5).sum())
            
            if n_liq_gaz > .5:
                self.G.add_node((frame, int(-1)), size = n_liq_gaz)
                
            cindex = 0
            
            while self.liq_drop_info[frame, cindex]['size'] > 0:
                self.G.add_node(
                    (frame, int(cindex)),
                    size = int(self.liq_drop_info[frame][cindex]['size']),
                    droplet = -1,
                    in_flux = 0,
                    out_flux = 0,
                    center_of_mass = tuple(self.liq_drop_info[frame][cindex]['center_of_mass']),
                    mean_number_of_neighbors = self.liq_drop_info[frame][cindex]['mean_number_of_neighbors'],
                    r_gyr = self.liq_drop_info[frame][cindex]['r_gyr'],
                    aniso = self.liq_drop_info[frame][cindex]['aniso'])

                cindex += 1
            
            if frame > 1:
                for particule in range(self.n_liq):
                    try:
                        self.G.edges[(frame-1, int(self.liq_info[frame-1][particule])),(frame, int(self.liq_info[frame][particule]))]["weight"] += 1
                    except Exception:
                        self.G.add_edge((frame-1, int(self.liq_info[frame-1][particule])),(frame, int(self.liq_info[frame][particule])), event = "")
                        self.G.edges[(frame-1, int(self.liq_info[frame-1][particule])),(frame, int(self.liq_info[frame][particule]))]["weight"] = 1
                    if int(self.liq_info[frame-1][particule]) < -.5 and int(self.liq_info[frame][particule]) > -.5:
                        self.G.nodes[(frame,int(self.liq_info[frame][particule]))]["in_flux"] += 1
                    elif int(self.liq_info[frame-1][particule]) > -.5 and int(self.liq_info[frame][particule]) < -.5:
                        self.G.nodes[(frame-1,int(self.liq_info[frame-1][particule]))]["out_flux"] += 1

            elif frame > 0:
                for particule in range(self.n_liq):
                    try:
                        self.G.edges[(0, -1),(frame, int(self.liq_info[frame][particule]))]["weight"] += 1
                    except Exception:
                        self.G.add_edge((0, -1),(frame, int(self.liq_info[frame][particule])))
                        self.G.edges[(0, -1),(frame, int(self.liq_info[frame][particule]))]["weight"] = 1
                    if int(self.liq_info[frame][particule]) > -.5:
                        self.G.nodes[(frame,int(self.liq_info[frame][particule]))]["in_flux"] += 1

        self.dict_max_give = {(frame, cluster) : max({(neigh_frame, neigh_cluster) : weight if neigh_frame > frame else -1 for _, (neigh_frame, neigh_cluster), weight in self.G.edges((frame, cluster), data = "weight")}.items(), key=lambda k: k[1])[0] if frame > -.5 else (frame, cluster) for frame, cluster in self.G.nodes()}
        self.dict_max_receive = {(frame, cluster) : max({(neigh_frame, neigh_cluster) : weight if neigh_frame < frame else -1 for _, (neigh_frame, neigh_cluster), weight in self.G.edges((frame, cluster), data = "weight")}.items(), key=lambda k: k[1])[0] if frame < self.n_frame else (frame, cluster) for frame, cluster in self.G.nodes()}

        for (u, v) in self.G.edges:
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

        print("LifeTime : Compute_graph completed")

        
    def compute_droplet(self):
        print("\n")

        for node in self.G.nodes():
            if node[1] > -.5:
                self.compute_node_neighbors(node)
        
        for frame in range(self.n_frame):
            if (frame) % 10 == 0:
                print("LifeTime[compute_droplet] : Process %d out of %d frames" % (frame, self.n_frame))
           
            for edges in self.dict_edges_per_frame_and_type[(frame, frame+1)][0]:
                self.compute_edge_0(edges)
                
            for edges in self.dict_edges_per_frame_and_type[(frame, frame+1)]['A']:
                self.compute_edge_A(edges)
                
            for edges in self.dict_edges_per_frame_and_type[(frame, frame+1)]['B']:
                self.compute_edge_B(edges)
                
            for edges in self.dict_edges_per_frame_and_type[(frame, frame+1)]['C']:
                self.compute_edge_C(edges)

        self.check_droplet()
        
        print("LifeTime : compute_droplet completed")


    def check_droplet(self):
        for node in self.G.nodes():
            if  node[1] > -.5 and node[0] < self.n_frame - 1 and self.G.nodes[node]["droplet"] == -1:
                raise NameError(f"node {node} as no droplet associated")
        for edge in self.G.edges():
            if (edge[0][1] > -0.5 or edge[1][1] > -0.5) and self.G.edges[edge]["event"] == "":
                raise NameError(f"edge {edge} as no event associated")
        for keys, droplet in self.droplet_dict.items():
            for i in range(droplet.tau):
                if droplet.nodes[i][0] != droplet.frame_start + i:
                    raise NameError(f"node {droplet.nodes[i]}'s frame in droplet {droplet.droplet_id} is not consistant")

            self.simple_droplet_dict[keys] = SDroplet(droplet.droplet_id, droplet.tau, droplet.nodes, droplet.sizes, droplet.center_of_mass, droplet.mean_number_of_neighbors)


    def compute_node_neighbors(self, node):
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

           
    def compute_edge_0(self, edge):

        u, v = edge
        list_in_v_B = self.node_neighbors_dict[v]["in_B"]

        if u[1] < -.5 and list_in_v_B == []: # Cas 1 : Emergence
            self.droplet_dict[self.droplet_id] = Droplet(self.droplet_id, v, self.G.nodes[v]["size"], [Event('emergence', self.G.nodes[v]["in_flux"], None)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["mean_number_of_neighbors"])
            self.G.nodes[v]["droplet"] = self.droplet_id
            self.droplet_id += 1

            self.G.edges[u, v]["event"] = 'emergence'

        elif u[1] > -.5: # Cas persistance standart
            self.droplet_dict[self.G.nodes[u]["droplet"]].append(v, self.G.nodes[v]["size"], [Event('persistance', self.G.edges[v, u]["weight"], u )], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["mean_number_of_neighbors"])
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
            

    def compute_edge_A(self, edge):

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
                self.droplet_dict[self.G.nodes[in_v_B_max]["droplet"]].append(v, self.G.nodes[v]["size"], [Event('persistance', self.G.edges[v, in_v_B_max]["weight"], in_v_B_max)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["mean_number_of_neighbors"])
                self.droplet_dict[self.G.nodes[in_v_B_max]["droplet"]].append_out_event(in_v_B_max, Event('persistance', self.G.edges[v, in_v_B_max]["weight"], v))
                self.G.nodes[v]["droplet"] = self.G.nodes[in_v_B_max]["droplet"]
                self.G.edges[in_v_B_max, v]["event"] = 'persistance'
            
                self.droplet_dict[self.G.nodes[in_v_B_max]["droplet"]].append_in_event(v, Event('expand', self.G.nodes[v]["in_flux"], None))
                self.G.edges[u, v]["event"] = 'expand'

                if self.G.nodes[in_v_B_max]["out_flux"] > 0:
                    self.droplet_dict[self.G.nodes[in_v_B_max]["droplet"]].append_out_event(in_v_B_max, Event('flux', self.G.nodes[in_v_B_max]["out_flux"], None))
                    self.G.edges[in_v_B_max, (v[0], -1)]["event"] = 'flux'

                for in_v_B in list_in_v_B:
                    if in_v_B != in_v_B_max and in_v_B not in list_in_v_A:
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
                    if out_node_B_max == v:     # Cas de persistence

                        self.droplet_dict[self.G.nodes[u]["droplet"]].append(v, self.G.nodes[v]["size"], [Event('persistance', self.G.edges[v, u]["weight"], u)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["mean_number_of_neighbors"])
                        self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('persistance', self.G.edges[v, u]["weight"], v))
                        self.droplet_dict[self.G.nodes[u]["droplet"]].append_out_event(u, Event('shrink', self.G.nodes[u]["out_flux"], None))
                        self.G.nodes[v]["droplet"] = self.G.nodes[u]["droplet"]
                        
                        self.G.edges[u, v]["event"] = 'persistance'
                        if self.G.nodes[v]["in_flux"] > 0:
                            self.droplet_dict[self.G.nodes[v]["droplet"]].append_in_event(v, Event('flux', self.G.nodes[v]["in_flux"], None))
                            self.G.edges[(v[0]-1, -1), v]["event"] = 'flux'
                        self.G.edges[u, (v[0], -1)]["event"] = 'shrink'
                        
                    else:               # Cas de splinter
                        self.droplet_dict[self.droplet_id] = Droplet(self.droplet_id, v, self.G.nodes[v]["size"], [Event('splinter', self.G.edges[v, u]["weight"], u)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["mean_number_of_neighbors"])
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
                    self.droplet_dict[self.droplet_id] = Droplet(self.droplet_id, v, self.G.nodes[v]["size"], [Event('split', self.G.edges[v, u]["weight"], u)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["mean_number_of_neighbors"])
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

                    self.droplet_dict[self.droplet_id] = Droplet(self.droplet_id, v, self.G.nodes[v]["size"], [Event('offload', self.G.edges[v, u]["weight"], u)], self.G.nodes[v]["center_of_mass"], self.G.nodes[v]["mean_number_of_neighbors"])
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


    def compute_edge_B(self, edge):  # On s'interesse maintenant à node comme à un g
        
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


    def compute_edge_C(self, edge):
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
        
    

    # def Print(self):

    #     font = font_manager.FontProperties(     weight='bold',
    #                         style='normal', size=60)
    #     fontlabel = {"labelsize" : 60}

        
    #     fig = plt.figure(figsize=(32,24))
    #     ax = fig.add_subplot()
        


    #     max_size = -1
    #     for droplet in self.droplet_dict.values():
    #         max_size = max(max_size, np.max(droplet.sizes))

    #     rainbow = mpl.colormaps['plasma'].resampled(max_size)
    #     for droplet in self.droplet_dict.values():
    #         if droplet.tau > 3 and max(droplet.sizes) > 4:
    #         # volume = np.array(droplet.sizes)*20e-9**3/2**(1/2)

    #         # ax.scatter([droplet.nodes[0][0]+ i for i in range(len(droplet.sizes))], droplet.sizes, color = rainbow(np.mean(droplet.mean_number_of_neighbors)), s = 100)
    #             ax.plot([droplet.nodes[0][0] + i for i in range(len(droplet.sizes))], [dens*12 for dens in droplet.mean_number_of_neighbors], color = rainbow(np.max(droplet.sizes)), linewidth=8)

    #     exp = self.input_dir.split("data/")[1].split("/")[0]
    #     params = "_".join(self.input_dir.split("data/")[1].split("/N/")[0].split("/")[1:])

    #     ax.set_title("Film of the growth of the droplet", font = font)
    #     ax.set_ylabel("Size Droplet (Number of PRC1)", font = font)
    #     ax.set_xlabel("Time (kMCS)", font = font)
        
    #     ax.tick_params(axis="both", **fontlabel)
        

    #     fig.savefig(f"/home/ppuel/data/{exp}/LifeTime_rho_loc_{params}.png")



    def print(self):
        print("\n")
        with open(os.path.join(self.input_dir, "liq_droplets.pickle"), "wb") as pfile:
            pickle.dump(self.droplet_dict, pfile, protocol=pickle.HIGHEST_PROTOCOL)
        
        with open(os.path.join(self.input_dir, "liq_simple_droplets.pickle"), "wb") as pfile:
            pickle.dump(self.simple_droplet_dict, pfile, protocol=pickle.HIGHEST_PROTOCOL)
        
        print("Pickle liq_droplets.pickle \n& liq_simple_droplets.pickle printed")
    
    # @staticmethod
    # @numba.njit
    # def centerOfMass(pts : np.ndarray[tuple[int, int], np.dtype[np.float32]], dims : np.ndarray[tuple[int], np.dtype[np.int32]]) -> tuple[float, float, float]: #, dims : list[np.int32]
    #     cord0 = pts[0]
        
    #     for i in range(1, len(pts)):
    #         cord = pts[i]
            
    #         cord0 = np.where(np.abs(cord-cord0) > dims/2,
    #                     (cord0 * i + cord - np.copysign(dims[0], cord - cord0))/(i+1),
    #                     (cord0 * i + cord)/(i+1))
        
    #     cm_x, cm_y, cm_z = np.where(np.logical_or(cord0 < 0, cord0 > dims[0]),
    #                     cord0 - np.copysign(dims, cord0),
    #                     cord0)
        
    #     return(cm_x, cm_y, cm_z)
        

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage is {sys.argv[0]} input_dir")
        sys.exit()

    input_dir = sys.argv[1]

    lifeTime = LifeTime(input_dir)
    
    lifeTime.compute_graph()
    lifeTime.compute_droplet()
    lifeTime.print()   