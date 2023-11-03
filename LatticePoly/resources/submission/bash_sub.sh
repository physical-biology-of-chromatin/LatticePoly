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
#rm -rf /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/*/*/r_1*
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_120_d0.05 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.05 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.075 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.2 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.3 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.4 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.5 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.6 5010 1
bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_all_frames_general.py /scratch/Bio/ddasaro/LatticeData/Tads_2L_0.0 500 2
 

