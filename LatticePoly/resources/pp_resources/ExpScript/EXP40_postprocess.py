import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import timedelta
import h5py
import matplotlib as mpl
import scipy.optimize as op


plasma = mpl.colormaps["plasma"].resampled(5)

exp = "EXP40_3R_PRC1_premiere_calibration"

ldens_list = ["0.0008","0.0016","0.0032"]
Nmeas_list = ["10","100"]
Ninter_list = ["100","1000","10000","100000"]
polyType_list = ["MCHeteroPoly","MCLivingPoly"]

os.makedirs(f"/home/ppuel/data/{exp}/calibration",exist_ok=True)

def MSD_figure():

        fig = plt.figure()
        ax = fig.add_subplot()
        X = []
        for Ninter in Ninter_list:
                for Nmeas in Nmeas_list:
                        if not(Ninter == "100000" and Nmeas == "100"):
                                for i in range(1,int(Nmeas)+1):
                                        if not(int(Ninter)*i in X):
                                                X.append(int(Ninter)*i)
        Y = [0 for i in range(len(X))]
        cumul = [0 for i in range(len(X))]
        for polyType in polyType_list:
                for Nmeas in Nmeas_list:
                        for Ninter in Ninter_list:
                                if not(Ninter == "100000" and Nmeas == "100"):
                                        for ldens in ldens_list:
                                                file = h5py.File(f"/Xnfs/physbiochrom/ppuel/data/{exp}/LDENS/{ldens}/NMEAS/{Nmeas}/NINTER/{Ninter}/POLYTYPE/{polyType}/N/0/process.h5")
                                                MSDpoly = np.array(file["polyHetMSD"])
                                                # print(MSDpoly[0:5])
                                                # ax.plot(np.log10([i*int(Ninter) for i in range(1,int(Nmeas)+1)]),np.log10(MSDpoly[1:]*(20e-3)**2),color = plasma(int(np.log10(float(ldens)/0.0008))+1), label = f'{ldens}')
                                                for i in range(1,int(Nmeas)+1):
                                                        Y[X.index(int(Ninter)*i)] += MSDpoly[i]
                                                        cumul[X.index(int(Ninter)*i)] += 1
                                                # print(MSDpoly[0], polyType, Nmeas, Ninter)
                                                # Y[0] += MSDpoly[0]
                                                # cumul[0] += 1
        X = np.array(X)
        Y = np.array(Y)/(np.array(cumul))
        # alpha = np.zeros(len(X)-1)
        # for i in range(1,len(X)):
        #         # print(Y[i],Y[0],X[i])
        #         alpha[i-1] = np.mean(np.log10(Y[1:i])/np.log10(X[1:i]))

        # poly = np.polyfit(X,Y,1)
        # print(poly)
        logX = np.log10(X)
        logY = np.log10(Y)
        ax.scatter(logX,logY,color='red')
        ax1 = ax.twinx()
        ax1.scatter(logX[:-1],np.diff(logY)/np.diff(logX),color='blue')
        ax1.set_ylim(0,1)
        for i in range(4):
                poly = np.polyfit(logX[logX <= 3+i],logY[logX <= 3+i],1)
                ax.plot(logX[logX <= 3+i], np.polyval(poly, logX[logX <= 3+i]))
                # print(poly)


        # f = lambda x, a: x/2 + a
        # res = op.curve_fit(f, np.log10(X), np.log10(Y))
        # print(res[0][0])
        D_mcs = (np.mean((np.array(Y[np.logical_and(logX>=5, logX<=5.5)])*(20e-3*2**(1/2))**2)/(np.array(X[np.logical_and(logX>=5, logX<=5.5)]))**(1/2))/0.01)**(2)
        print(D_mcs)
        # ax.plot(np.log10(X),f(np.log10(X),res[0][0]))
        # print(T_mcs)


        fig.savefig(f"/home/ppuel/data/{exp}/calibration/Diffusion2.png")


def duration_figure():


        fig = plt.figure()
        ax = fig.add_subplot()


        for polyType in polyType_list:
                X = []
                Y = []
                for Nmeas in Nmeas_list:
                        for Ninter in Ninter_list:
                                for ldens in ldens_list:
                                        try :
                                                file = open(f"/Xnfs/physbiochrom/ppuel/data/{exp}/LDENS/{ldens}/NMEAS/{Nmeas}/NINTER/{Ninter}/POLYTYPE/{polyType}/N/0/log.out")

                                                for line in file.readlines():
                                                        if line[:5] == 'Total':
                                                                duration = float(line.split(" ")[2])
                                        except :
                                                duration = 1
                                                print("error")
                                        X.append(np.log10((int(Nmeas)+100)*int(Ninter)))
                                        Y.append(np.log10(duration))

                ax.scatter(X,Y,s=100)

                poly = np.polyfit(X,Y,1)
                print(poly)
                ax.plot(X,np.polyval(poly, X))

                print(timedelta(seconds = 60*10**float(np.polyval(poly, np.log10(172628170)))))


        fig.savefig(f"/home/ppuel/data/{exp}/calibration/duration.png")


MSD_figure()