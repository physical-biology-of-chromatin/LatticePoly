#!/bin/bash


##
##  slurm.sh
##  LatticePoly
## 
##  Created by mtortora on 20/06/2022
##  Copyright © 2022 ENS Lyon. All rights reserved.
##

#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

# Relative path to code root directory
ROOTDIR=${SCRIPTDIR}/../..

# Set working directory to root


 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_120_d0.05 N_120_d0.05
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_120_d0.075 N_120_d0.075
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_120_d0.1 N_120_d0.1
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_120_d0.2 N_120_d0.2
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_120_d0.3 N_120_d0.3
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_120_d0.4 N_120_d0.4
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_120_d0.5 N_120_d0.5
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_120_d0.6 N_120_d0.6
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_120_d0.7 N_120_d0.7

 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_110_d0.05 N_110_d0.05
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_110_d0.075 N_110_d0.075
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_110_d0.1 N_110_d0.1
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_110_d0.2 N_110_d0.2
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_110_d0.3 N_110_d0.3
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_110_d0.4 N_110_d0.4
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_110_d0.5 N_110_d0.5
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_110_d0.6 N_110_d0.6
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_110_d0.7 N_110_d0.7

 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_100_d0.05 N_100_d0.05
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_100_d0.075 N_100_d0.075
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_100_d0.1 N_100_d0.1
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_100_d0.2 N_100_d0.2
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_100_d0.3 N_100_d0.3
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_100_d0.4 N_100_d0.4
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_100_d0.5 N_100_d0.5
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_100_d0.6 N_100_d0.6
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_100_d0.7 N_100_d0.7


 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_90_d0.05 N_90_d0.05
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_90_d0.075 N_90_d0.075
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_90_d0.1 N_90_d0.1
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_90_d0.2 N_90_d0.2
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_90_d0.3 N_90_d0.3
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_90_d0.4 N_90_d0.4
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_90_d0.5 N_90_d0.5
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_90_d0.6 N_90_d0.6
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_90_d0.7 N_90_d0.7

 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_80_d0.05 N_80_d0.05
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_80_d0.075 N_80_d0.075
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_80_d0.1 N_80_d0.1
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_80_d0.2 N_80_d0.2
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_80_d0.3 N_80_d0.3
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_80_d0.4 N_80_d0.4
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_80_d0.5 N_80_d0.5
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_80_d0.6 N_80_d0.6
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_80_d0.7 N_80_d0.7
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.05 N_70_d0.05
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.075 N_70_d0.075
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.1 N_70_d0.1
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.2 N_70_d0.2
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.3 N_70_d0.3
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.4 N_70_d0.4
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.5 N_70_d0.5
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.6 N_70_d0.6
 python  /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.7 N_70_d0.7

