##
##  LiqPolyCoincidence.py
##  LatticePoly
##Liq
##  Created by ppuel on ../../2023.
##  Copyright © 2023 ENS Lyon. All rights reserved.
##

import os, sys, re

import numpy as np

import h5py

class SimLiqPolyTime():

    def __init__(self, inputDir, refDir, outputDir, N, Nmeas, L, lDens, is_poly):

        self.inputDir = inputDir
        self.refDir = refDir
        self.outputDir = outputDir
        self.N = N
        self.Nmeas = Nmeas
        self.L = L
        self.lDens = lDens
        self.is_poly = is_poly

        if 'NINTER' in self.inputDir.split("/"):
            self.nInter = int(self.inputDir.split("/")[self.inputDir.split("/").index('NINTER')+1])
        else:
            try:
                with open(os.path.join("/".join(self.inputDir.split("/")[:6]), "input_slurm.cfg")) as cfile:
                    self.nInter = int(re.findall('Ninter = \d*', cfile.read())[0].split(' = ')[1])
            except:
                self.nInter = 100000


        if 'NMEAS' in self.inputDir.split("/"):
            self.Nmeas_ref = int(self.inputDir.split("/")[self.inputDir.split("/").index('NMEAS')+1])
        else:
            try:
                with open(os.path.join("/".join(self.inputDir.split("/")[:6]), "input_slurm.cfg")) as cfile:
                    self.Nmeas_ref = int(re.findall('Nmeas = \d*', cfile.read())[0].split(' = ')[1])
            except:
                self.Nmeas_ref = 100

        if 'NINTER' in self.refDir.split("/"):
            self.nInter_ref = int(self.refDir.split("/")[self.refDir.split("/").index('NINTER')+1])
        else:
            try:
                with open(os.path.join("/".join(self.refDir.split("/")[:6]), "input_slurm.cfg")) as cfile:
                    self.nInter_ref = int(re.findall('Ninter = \d*', cfile.read())[0].split(' = ')[1])
            except:
                self.nInter_ref = 100000

        try:
            with open(os.path.join("/".join(self.refDir.split("/")[:6]), "input_slurm.cfg")) as cfile:
                self.Nref = int(re.findall('Nstat = \d*', cfile.read())[0].split(' = ')[1])
        except:
            self.Nref = 1
        


        self.DPoly = 1e4 # nm²/s^(1/2)
        self.DLiq = 1e6 # nm²/s

    

    def Compute(self):

        # Get absolute reference 
        
        self.D_CoM_ref, self.D_liq_ref, self.D_poly_ref = self.ProcessPath_reference()

        # self.D_CoM_ref, self.D_liq_ref, self.D_poly_ref = 0.0046717,406.73,114.44

        print(self.D_CoM_ref, self.D_liq_ref, self.D_poly_ref)

        print("Relative references start")

        self.ProcessPath_experience()
            
        print("Relative references calculated")


    def ProcessPath_reference(self):
        print(self.refDir)

        liqMSD = np.zeros(self.Nmeas_ref//2-1)

        if self.is_poly:
            LiqPolyCoM = np.zeros(self.Nmeas_ref//2-1)
            polyMSD = np.zeros(self.Nmeas_ref//2-1)
        else:
            liqMSD_CoM = np.zeros(self.Nmeas_ref//2-1)
    
        for n in range(self.Nref):
            print(n)
            nFile = h5py.File(os.path.join(self.refDir,f"N/{n}/process.h5"), 'r')
            print(n, 'done')
            liqMSD += nFile["liqMSD"][1:self.Nmeas_ref//2]

            if self.is_poly:
                LiqPolyCoM += nFile["LiqPolyCoM"][1:self.Nmeas_ref//2]
                polyMSD += nFile["polyHomMSD"][1:self.Nmeas_ref//2] #approximation
            else: 
                liqMSD_CoM += nFile["liqMSD_CoM"][1:self.Nmeas_ref//2]
            
            nFile.close()

        D_liq_ref = np.mean(liqMSD/np.arange(1,self.Nmeas_ref//2))/self.nInter_ref*20**2*2/self.Nref # nm²/MCS
        
        if self.is_poly:
            D_CoM_ref = np.mean(LiqPolyCoM/np.arange(1,self.Nmeas_ref//2))/self.nInter_ref*20**2*2/self.Nref # nm²/MCS
            D_poly_ref = np.mean(polyMSD/(np.arange(1,self.Nmeas_ref//2)*self.nInter_ref)**(1/2))*20**2*2/self.Nref # nm²/MCS^(1/2)
        else: 
            D_CoM_ref = np.mean(liqMSD_CoM/np.arange(1,self.Nmeas_ref//2))/self.nInter_ref*20**2*2/self.Nref # nm²/MCS
            D_poly_ref = -1
        
        # Checking the slope of MSD, plus the ratio between polymer and particles diffusion ratio

        # p, _ = np.polyfit(np.log10(np.arange(1,self.Nmeas_ref+1)*self.nInter), np.log10(liqMSD), 1)
        
        # print(p[0], 1, np.abs(p[0]-1)/D_poly_ref)
        # print(10**(p[1]), D_liq_ref, np.abs(10**(p[1])-D_liq_ref)/D_liq_ref)
        
        # if self.is_poly:
            # p, _ = np.polyfit(np.log10(np.arange(1,self.Nmeas_ref+1)*self.nInter), np.log10(polyMSD), 1)
            
            # print(p[0], 0.5, np.abs(p[0]-0.5)/0.5)
            # print(10**(p[1]), D_liq_ref, np.abs(10**(p[1])-D_liq_ref)/D_liq_ref)

            # print(D_poly_ref/D_liq_ref, self.DPoly/self.DLiq, np.abs(D_poly_ref/D_liq_ref - self.DPoly/self.DLiq)/self.DPoly/self.DLiq)
       
        # MSD = D * MCS ** alpha

        # log(MSD) = log(D) + alpha * log(MCS)



        return(D_CoM_ref, D_liq_ref, D_poly_ref)

    def ProcessPath_experience(self):
        
        print("ProcessPath_experience")

        if self.is_poly:
            LiqPolyCoM = np.zeros(self.Nmeas//2-1)
        else:
            liqMSD_CoM = np.zeros(self.Nmeas//2-1)
    
        for n in range(N):
            print(n)
            nFile = h5py.File(os.path.join(self.inputDir,f"N/{n}/process.h5"), 'r')
            
            if self.is_poly:
                LiqPolyCoM += nFile["LiqPolyCoM"][1:self.Nmeas//2]
            else: 
                liqMSD_CoM += nFile["liqMSD_CoM"][1:self.Nmeas//2]
            
            nFile.close()
            print(n,' done')


        if self.is_poly:
            D_CoM_abs = np.mean(LiqPolyCoM/np.arange(1,self.Nmeas//2))/self.nInter*20**2*2/N # nm²/MCS
        else: 
            D_CoM_abs = np.mean(liqMSD_CoM/np.arange(1,self.Nmeas//2))/self.nInter*20**2*2/N # nm²/MCS


        print(D_CoM_abs, self.D_CoM_ref, self.lDens, self.D_liq_ref, self.DLiq, self.D_poly_ref, self.DPoly)

        ratio = D_CoM_abs/self.D_CoM_ref*min(self.D_liq_ref/self.DLiq, self.D_poly_ref/self.DPoly)

        SimTime = np.arange(self.Nmeas)*self.nInter*ratio

        outputFile = h5py.File(os.path.join(self.outputDir, 'post_process.h5'), 'a')
        
        if 'SimTime' in outputFile.keys():
            tmp = outputFile['SimTime']
            tmp[:] = SimTime
        else:
            outputFile.create_dataset("SimTime", data = SimTime)
    

if __name__ == "__main__":
    if len(sys.argv) != 9:
        print("\033[1;31mUsage is %s inputDir refDir outputDir N Nmeas L lDens is_poly\033[0m" % sys.argv[0])
        sys.exit()

    inputDir = sys.argv[1]
    refDir = sys.argv[2]
    outputDir = sys.argv[3]

    N = int(sys.argv[4])
    Nmeas = int(sys.argv[5])

    L = int(sys.argv[6])
    lDens = float(sys.argv[7])
    
    is_poly = sys.argv[8]=="1"

    print("Starting")

    SimTime = SimLiqPolyTime(inputDir, refDir, outputDir, N, Nmeas, L, lDens, is_poly)

    SimTime.Compute()
    
    print("\n")
    print("SimLiqPolyTime : Done\n\n")
