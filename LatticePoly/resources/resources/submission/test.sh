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

#python /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices.py /home/ddasaro/LatticePoly/LatticePoly/data/output/January23/segregation_25_long/null/ chrVII null
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/PolyGyration_repl_chr1.py  ${CURRENTDIR}  25000
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 10  
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 20  
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 40  
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 80  
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   2000 500 100  
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 200  
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 300  
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 400  
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 500
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 600  
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 700  
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 800  
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 900 
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${CURRENTDIR}   20000 500 999
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 50050
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 50150
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 50250
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 50350
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 50450
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 50550
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 50650
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 50750
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 50850
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 50950
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 51050
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 51150
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 51250
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 50000 2 51350
#python /home/ddasaro//LatticePoly/LatticePoly/resources//Poly_Rcmdiff_SCs.py ${CURRENTDIR} 20000
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${CURRENTDIR} 50000 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${CURRENTDIR} 25000 0 2 
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${CURRENTDIR} 25000 0 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${CURRENTDIR} 25000 50 2 
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${CURRENTDIR} 25000 50 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${CURRENTDIR} 25000 100 2 
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${CURRENTDIR} 25000 100 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${CURRENTDIR} 25000 200 2 
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${CURRENTDIR} 25000 200 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${CURRENTDIR} 25000 400 2 
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${CURRENTDIR} 25000 400 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${CURRENTDIR} 25000 800 2 
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${CURRENTDIR} 25000 800 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${CURRENTDIR} 25000 25 2 
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${CURRENTDIR} 25000 25 2

#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${CURRENTDIR} 20000 620 2 
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${CURRENTDIR} 20000 620 2

#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${CURRENTDIR} 20000 800 2 
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${CURRENTDIR} 20000 800 2

#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${CURRENTDIR} 20000 20 2 
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${CURRENTDIR} 20000 20 2

#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 51850 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 0 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 100 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 200 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 300 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 400 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 500 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 600 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 700 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 800 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 900 2
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 1000 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 0 4
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 100 4
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 200 4
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 300 4
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 400 4
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 500 4
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 600 4
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 700 4
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 800 4
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 900 4
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 1000 4

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 0 3
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 100 3
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 200 3
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 300 3
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 400 3
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 500 3
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 600 3
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 700 3
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 800 3
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 900 3
#python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_intervals.py ${CURRENTDIR} 1000 3

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_all_frames.py ${CURRENTDIR} 10000 3

python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_200_cool.py ${TMPDIR}   20000 2
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_200_cool.py ${TMPDIR}   20000 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${TMPDIR} 2000 0 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${TMPDIR} 2000 0 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${TMPDIR} 25000 25 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${TMPDIR} 25000 25 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${TMPDIR} 25000 50 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${TMPDIR} 25000 50 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${TMPDIR} 25000 100 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${TMPDIR} 25000 100 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${TMPDIR} 25000 200 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${TMPDIR} 25000 200 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${TMPDIR} 25000 400 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${TMPDIR} 25000 400 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${TMPDIR} 25000 800 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${TMPDIR} 25000 800 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${TMPDIR} 25000 1000 2
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${TMPDIR} 25000 1000 2
