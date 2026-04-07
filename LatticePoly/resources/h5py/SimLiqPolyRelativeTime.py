##
##  LiqPolyCoincidence.py
##  LatticePoly
##Liq
##  Created by ppuel on ../../2023.
##  Copyright © 2023 ENS Lyon. All rights reserved.
##

import os
import re
import sys

import h5py
import numpy as np

class SimLiqPolyRelativeTime:
    def __init__(self, inDir, diffusionData):
        print(f"SimLiqPolyRelativeTime : Init {inDir} {diffusionData}")

        self.inDir = inDir
        self.diffusionData = diffusionData

        # processParam = {}

        # for param in re.findall('/\\D*/(\\.\\d)*'):
        #         processParam[param.split('/')[1]] = param.split('/')[2]

        # paramRefSlurm = {}

        # i = 0
        # while i < len(diffusionData):
        #         if re.fdiffusionData[i][-1] == '=':
        #                 paramRefSlurm[diffusionData[i][:-1]] = diffusionData[i+1]
        #                 i += 2
        #         else:
        #                 if all(f"{keys.upper()}_{values}" in diffusionData[i] for keys, values in paramRefSlurm.items()):
        #                         self.D_CoM_ref, self.D_liq_ref, self.D_poly_ref = float(diffusionData[i+1]), float(diffusionData[i+2]), float(diffusionData[i+3])
        #                         i += len(diffusionData)
        #                 else:
        #                         i += 4

        self.DPoly_experimental = 1e4  # nm²/s^(1/2)
        self.DLiq_experimental = 1e6  # nm²/s

    def Compute(self):

        nStat = 0

        isPoly, nMeas, nInter, lDens = "", "", "", ""

        for file in os.listdir(os.path.join(self.inDir, "N/")):
            if file != "6":
                if not nStat and file.isdigit():
                    with open(
                        os.path.join(self.inDir, f"N/{file}/input.cfg"), "r"
                    ) as cfile:
                        isPoly, nMeas, nInter, lDens = re.findall(
                            "domainPath = \\D* ; |Nmeas  = \\d*|Ninter = \\d*|Ldens = \\d*.\\d*",
                            cfile.read(),
                        )

                nStat += int(file.isdigit())

        isPoly = "toy_domain" not in isPoly
        nMeas = int(nMeas.split(" = ")[1])
        nInter = int(nInter.split(" = ")[1])
        lDens = float(lDens.split(" = ")[1])

        D_CoM_ref = float(self.diffusionData[0])
        D_Liq_ref = float(self.diffusionData[1])
        D_Poly_ref = float(self.diffusionData[2])  # / lDens

        if isPoly:
            D_CoM_ref_eff = D_CoM_ref / min(
                D_Liq_ref / self.DLiq_experimental, D_Poly_ref / self.DPoly_experimental
            )
        else:
            D_CoM_ref_eff = (
                D_CoM_ref
                / (float(self.diffusionData[1]) * lDens + float(self.diffusionData[2]))
                * self.DLiq_experimental
            )
        
        nMeasMSD = nMeas // 2

        MSDCoM = np.zeros(nMeasMSD)

        for n in range(nStat):
            if n != 6:
                with h5py.File(os.path.join(inDir, f"N/{n}/process.h5"), "r") as nFile:
                    if isPoly:
                        MSDCoM += nFile["LiqPolyCoM"][1 : nMeasMSD + 1]
                    else:
                        MSDCoM += nFile["liqMSD_CoM"][1 : nMeasMSD + 1]
            
        if isPoly:
            D_CoM_rel = (
                np.mean(MSDCoM / np.arange(1, nMeasMSD + 1))
                / nInter
                * 20**2
                * 2
                / nStat
            )  # nm²/MCS
        else:
            D_CoM_rel = (
                np.mean(MSDCoM / np.arange(1, nMeasMSD + 1))
                / nInter
                * 20**2
                * 2
                / nStat
            )  # nm²/MCS
        
        ratio = D_CoM_rel / D_CoM_ref_eff

        SimTime = np.arange(nMeas) * nInter * ratio
        
        with h5py.File(os.path.join(self.inDir, "post_process.h5"), "a") as hfile:
            self.PrintDataset(hfile, "SimTime", SimTime)

        with open(os.path.join(self.inDir, "time.txt"), "w") as timefile:
            timefile.write(str(nInter * ratio * nMeas))

    def PrintDataset(self, processFile, dataset_name, data):

        if dataset_name in processFile.keys():
            tmp = processFile[dataset_name]
            tmp[:] = data
        else:
            processFile.create_dataset(dataset_name, data=data)

        print(f"Dataset {dataset_name} printed")


if __name__ == "__main__":  # Usage is sys.argv[0] inDir outDir **diffusionData
    inDir = sys.argv[1]
    diffusionData = sys.argv[2:]

    RelativeTime = SimLiqPolyRelativeTime(inDir, diffusionData)

    RelativeTime.Compute()

    print("\n")
    print("SimLiqPolyRelativeTime : Done\n\n")

    RelativeTime.Compute()

    print("\n")
    print("SimLiqPolyRelativeTime : Done\n\n")
