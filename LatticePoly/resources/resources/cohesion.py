import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib
import time


outputDir = sys.argv[1]

genome=np.arange(2*872)
intra1_contact=np.loadtxt(outputDir+"/cohesion_pattern_cis1.res")
intra2_contact=np.loadtxt(outputDir+"cohesion_pattern_cis2.res")
inter_contact=np.loadtxt(outputDir+"/cohesion_pattern_trans.res")
plt.figure(figsize=(30, 3))

G = nx.DiGraph()

G.add_nodes_from(genome)
pos2 = {}
chromatid1=int(len(genome)/2)
for i, node in enumerate(range(chromatid1)):
    pos2[node] = [i,1.01]
for i, node in enumerate(range(chromatid1)):
    pos2[node+chromatid1] = [i,1]
color_map=[]
for node in G:
    if node >= chromatid1:
        color_map.append('r')
    else: 
        color_map.append('b') 
    edge_set_2 = []
edge_set_inter = []
for i in range(len(inter_contact)):
    if(i%2==0):
        edge_set_inter.append((int(inter_contact[i]),int(inter_contact[i+1]+chromatid1)))
edge_set_intra1 = []
for i in range(len(intra1_contact)):
    if(i%2==0):
        if(int(intra1_contact[i])>int(intra1_contact[i+1])):
            edge_set_intra1.append((int(intra1_contact[i]),int(intra1_contact[i+1])))
edge_set_intra2 = []
for i in range(len(intra2_contact)):
    if(i%2==0):
        if(int(intra2_contact[i])<int(intra2_contact[i+1])):
            edge_set_intra2.append((int(chromatid1+intra2_contact[i]),int(chromatid1+intra2_contact[i+1])))
nx.draw(G, pos2,node_color=color_map, node_size=20)
nx.draw_networkx_edges(G, pos2, edge_set_inter, width=3,edge_color='gold',arrowstyle="-")
nx.draw_networkx_edges(G, pos2, edge_set_intra1,width=3, edge_color='gold',arrowstyle="-",connectionstyle="arc3,rad=0.4")
nx.draw_networkx_edges(G, pos2, edge_set_intra2, width=3,edge_color='gold',arrowstyle="-",connectionstyle="arc3,rad=0.4")
plt.savefig(outputDir+'/'+str(time.time())+'cohesion.png')
