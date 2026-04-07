##
##  LiqPolyCoincidence.py
##  LatticePoly
##Liq
##  Created by ppuel on ../../2023.
##  Copyright © 2023 ENS Lyon. All rights reserved.
##

import os, subprocess, sys, h5py
import numpy as np

from utils import find_exp, exp_mapping, exp_pathList

from itertools import product, zip_longest

class LiqPolySimTime():

        def __init__(self, experience, reference):
                print(f"LiqPolySimTime : Init {experience} {reference}")
                                
                self.experience = experience
                self.reference = reference

                initDir = "/Xnfs/physbiochrom/ppuel/data/"

                exp_name = find_exp(self.experience, initDir)
                ref_name = find_exp(self.reference,  initDir)

                self.expDir = os.path.join(initDir, exp_name)
                self.refDir = os.path.join(initDir, ref_name)

                self.dict_parameters, self.isPoly, self.meta_parameter_N, self.meta_parameter_Nmeas = exp_mapping(self.expDir)
                self.ref_parameters , _          , self.ref_parameter_N , self.ref_parameter_Nmeas  = exp_mapping(self.refDir)

                self.pathList = exp_pathList(self.dict_parameters, -1, self.expDir)

                self.DPoly = 1e4 # nm²/s^(1/2)
                self.DLiq = 1e6 # nm²/s
                
        def coincidence_exp_ref(self):
                coincidence_ref_parameters = {}
                
                for keys, exp_values in self.dict_parameters.items():
                        ref_values = self.ref_parameters[keys]


                        if len(ref_values) == 1:
                                if len(exp_values) == 1 and ref_values[0] != exp_values[0]:
                                        print(f"missmatch with {keys} : {exp_values}, {ref_values}")
                        else:
                                if np.isin(exp_values, ref_values).all():
                                        coincidence_ref_parameters[keys] = exp_values
                                else:
                                        print(f"error with {keys} : {exp_values}, {ref_values}")
                                        sys.exit()

                return(coincidence_ref_parameters)                

        def Compute(self):
                print('\n')
                # Get absolute reference 
                
                self.coincidence_ref_parameters = self.coincidence_exp_ref()

                self.D_CoM_ref = {}
                self.D_liq_ref = {}
                self.D_poly_ref = {}

                ref_pathList = exp_pathList(self.coincidence_ref_parameters, -1, self.refDir)
                print(ref_pathList)
                for i, ref_path in enumerate(ref_pathList):
                        self.D_CoM_ref[ref_path], self.D_liq_ref[ref_path], self.D_poly_ref[ref_path] = 0, 0, 0 #self.ProcessPath_absolute(ref_path)
                        
                        if (i+1) % 10 == 0:
                                        print("Processed %d out of %d path" % (i+1, len(self.ref_pathList)))

                print("\nAbsolute references calculated\n")

                absoluteDiffusionList = ' '.join([f"{ref_path.replace("/", "_")} {self.D_CoM_ref[ref_path]:0.2e} {self.D_liq_ref[ref_path]:0.2e} {self.D_poly_ref[ref_path]:0.2e}" for ref_path in ref_pathList])
                refParamSlurmList = ' '.join([keys.upper() for keys in self.coincidence_ref_parameters.keys()])
                print(f".venv/bin/python3 resources/submission/submit_slurm_NewProcess.py {self.experience} resources/h5py/LiqPolyRelativeTime.py {refParamSlurmList} {absoluteDiffusionList}")
                #subprocess.run(f".venv/bin/python3 resources/submission/submit_slurm_NewProcess.py {self.experience} resources/h5py/LiqPolyRelativeTime.py {refParamSlurmList} {absoluteDiffusionList}", shell = True, executable = "\bin\bash")


        def ProcessPath_absolute(self, ref_path):
                
                Nmeas_MSD = self.ref_parameter_Nmeas//2

                if len(self.ref_parameters["Ninter"]) == 1:
                        nInter = self.ref_parameters["Ninter"][0]
                else:
                        nInter = int(ref_path.split("/")[ref_path.split("/").index('NINTER')+1])
                
                liqMSD = np.zeros(Nmeas_MSD)

                if self.isPoly:
                        LiqPolyCoM = np.zeros(Nmeas_MSD)
                        polyMSD = np.zeros(Nmeas_MSD)
                else:
                        liqMSD_CoM = np.zeros(Nmeas_MSD)
        
                for n in range(self.ref_parameter_N):
                        with h5py.File(os.path.join(ref_path,f"N/{n}/process.h5"), 'r') as nFile:
                        
                            liqMSD += nFile["liqMSD"][1:Nmeas_MSD+1]

                            if self.isPoly:
                                    LiqPolyCoM += nFile["LiqPolyCoM"][1:Nmeas_MSD+1]
                                    polyMSD += (nFile["polyHomMSD"][1:Nmeas_MSD+1] + nFile["polyHetMSD"][1:Nmeas_MSD+1])/2 #approximation
                            else: 
                                    liqMSD_CoM += nFile["liqMSD_CoM"][1:Nmeas_MSD+1]
                        
                D_liq_ref = np.mean(liqMSD/np.arange(1,Nmeas_MSD+1))/nInter*20**2*2/self.ref_parameter_N # nm²/MCS
                
                if self.isPoly:
                        D_CoM_ref = np.mean(LiqPolyCoM/np.arange(1,Nmeas_MSD+1))/nInter*20**2*2/self.ref_parameter_N # nm²/MCS
                        D_poly_ref = np.mean(polyMSD/(np.arange(1,Nmeas_MSD+1)*nInter)**(1/2))*20**2*2/self.ref_parameter_N # nm²/MCS^(1/2)
                else: 
                        D_CoM_ref = np.mean(liqMSD_CoM/np.arange(1,Nmeas_MSD+1))/nInter*20**2*2/self.ref_parameter_N # nm²/MCS
                        D_poly_ref = -1
                
                return(D_CoM_ref, D_liq_ref, D_poly_ref)


if __name__ == "__main__":
        if len(sys.argv)  != 3:
                print("\033[1;31mUsage is %s experience reference\033[0m" % sys.argv[0])
                sys.exit()

        experience = sys.argv[1]
        reference = sys.argv[2]

        SimTime = LiqPolySimTime(experience, reference)

        SimTime.Compute()
                
        print("\n")
        print("LiqPolySimTime : Done\n\n")