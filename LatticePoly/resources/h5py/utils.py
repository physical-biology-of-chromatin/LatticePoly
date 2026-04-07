##
##  utils.py
##  LatticePoly
##
##  Modified by ppuel on 8/10/2025.
##  Copyright © 2025 ENS Lyon. All rights reserved.
##

# from ast import Dict
# from cProfile import label
# from logging import handlers
import os, sys, h5py
import numpy as np
from hdf5Reader import hdf5Reader
from itertools import product, zip_longest
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple, HandlerPatch

from matplotlib.colors import LogNorm, CenteredNorm
from scipy.spatial.distance import squareform

def AggregateHist(expPath, listDatasetName, pathList, metaParameterN, metaParameterNmeas):

        metaParameterN = int(metaParameterN)
        metaParameterNmeas = int(metaParameterNmeas)

        c_max = len(pathList)*metaParameterN
        c = 0
        c_problem = 0

        datasetnHist = {}

        with h5py.File(os.path.join(pathList[0], "N/0/process.h5"), 'r') as tmpfile:
                for datasetName in listDatasetName:
                        datasetnHist[datasetName] = np.shape(tmpfile[datasetName])[-1] if len(np.shape(tmpfile[datasetName]))>1 else 1


        for dataPath in pathList:
                Nb_sim = 0

                datasetArray = {}
                for datasetName in listDatasetName:
                        datasetArray[datasetName] = np.zeros((metaParameterNmeas, datasetnHist[datasetName], metaParameterN)) if datasetnHist[datasetName] > 1 else np.zeros((metaParameterNmeas, metaParameterN))

                for n in range(metaParameterN):
                        print(f"{c/c_max*100: 2.1f}%  {c_problem/c_max*100: 2.1f}%", end='\r')
                        try:
                                
                                with h5py.File(os.path.join(dataPath, f"N/{n}/process.h5"), 'r') as dataFile:
                                        for datasetName in listDatasetName:
                                                if datasetnHist[datasetName] > 1 :
                                                        datasetArray[datasetName][:,:,n] = dataFile[datasetName]
                                                else:
                                                        datasetArray[datasetName][:,n] = dataFile[datasetName]
                                                
                        except Exception as inst:
                                
                                with open(f"{expPath}/simu.txt", 'a') as tfile:
                                        tfile.write(dataPath+f" {n}\n"+str(inst)+'\n\n')
                                        c_problem += 1


                        c += 1

                with h5py.File(os.path.join(dataPath,"aggregated_process.h5"), 'a') as aggregatedFile:
                        for datasetName in listDatasetName:
                                PrintDataset(aggregatedFile, datasetName, datasetArray[datasetName])    


def FormatParameters(suffixe, listParameters, searchDictParam, dict_parameters):
        tmpDict = {}
        tmpNum = 1
        if len(listParameters)>0:
                for ids, aggParam in enumerate(listParameters):
                        if aggParam in dict_parameters.keys():
                                tmpDict[aggParam.upper()] = dict_parameters[aggParam]
                                searchDictParam[aggParam.upper()] = (suffixe,ids)
                                tmpNum *= len(dict_parameters[aggParam])
        else:
              tmpDict["None"] = [None]

        return(tmpDict, searchDictParam, tmpNum)


def str_to_slice(stringSlice):
        return(slice(*list(map(lambda x: int(x) if len(x)>0 else None, stringSlice.split(':')))))


def PrintAggregateHist(expPath, figDir, listDatasetName, dict_parameters, metaParameterN, metaParameterNmeas):
        
        metaParameterN = int(metaParameterN)
        metaParameterNmeas = int(metaParameterNmeas)

        tmp_dict_param = {}
        for keys, values in dict_parameters.items():
                if len(values) > 1:
                        tmp_dict_param[keys.upper()] = values
        
        dict_parameters["MCS"] = [str(i) for i in range(metaParameterNmeas)]

        # listAllParameters = ["MCS"] + list(tmp_dict_param.keys())

        figParameters = ["LDENS"]
        axsParameters = ["JLP", "JLL"]         #["colomn", "row"] or ["unique"]
        pltParameters = ["MCS"]                # size 1
        aggParameters = ["EV"]
        
        DictParametersStringSlice = {"MCS" : "::100", "JLL" : "::5", "JLP" : "::2"}
        DictParametersSlicing = {}

        for parameter, stringSlice in DictParametersStringSlice.items():
                DictParametersSlicing[parameter] = str_to_slice(stringSlice)

        for parameter, slicing in DictParametersSlicing.items():
                dict_parameters[parameter] = dict_parameters[parameter][slicing] 
        
        aggDataDir = os.path.join(expPath, "aggData")
        os.makedirs(aggDataDir, exist_ok=True)

        searchDictParam = {}

        figDictParam, searchDictParam, figNum = FormatParameters('fig', figParameters, searchDictParam, dict_parameters)
        axsDictParam, searchDictParam, axsNum = FormatParameters('axs', axsParameters, searchDictParam, dict_parameters)
        pltDictParam, searchDictParam, pltNum = FormatParameters('plt', pltParameters, searchDictParam, dict_parameters)
        aggDictParam, searchDictParam, aggNum = FormatParameters('agg', aggParameters, searchDictParam, dict_parameters)
        
        colormapNameList = ["Purples", "Greens"]
        
        colormapList = {}
        for ids, datasetName in enumerate(listDatasetName):
                colormapList[datasetName] = colormaps[colormapNameList[ids]].resampled(pltNum+2)


        for params in product(*tmp_dict_param.values()):
                tmp_path = "/".join(map(lambda x : "/".join(x), zip_longest(tmp_dict_param.keys(), params)))+"/"
                pathInit = os.path.join(expPath, tmp_path)
                break

        datasetnHist = {}

        with h5py.File(os.path.join(pathInit, "N/0/process.h5"), 'r') as tmpfile:
                for datasetName in listDatasetName:
                        datasetnHist[datasetName] = np.shape(tmpfile[datasetName])[1] if len(np.shape(tmpfile[datasetName])) > 1 else 1

        
        c_max = figNum*axsNum*pltNum*aggNum
        c = 0
        
        dictParametersValues = {"fig" : None, "axs" : None, "plt" : None, "agg" : None}

        for figParamValue in product(*figDictParam.values()):

                dictParametersValues["fig"] = figParamValue

                if len(axsParameters) == 1:
                        axsSeparation = int(np.sqrt(axsNum))
                        while axsNum%axsSeparation > 0:
                                axsSeparation -= 1
                                                              
                        fig = plt.figure(figsize=(4*axsNum//axsSeparation, 4*axsSeparation), layout='constrained')
                else:
                        fig = plt.figure(figsize=(4*len(axsDictParam[axsParameters[1]]), 4*len(axsDictParam[axsParameters[0]])), layout='constrained')
                        
                fig.suptitle(", ".join(map(lambda x : " = ".join(x), zip_longest(figDictParam.keys(), figParamValue))))
                
                for axsIds, axsParamValue in enumerate(product(*axsDictParam.values())):

                        dictParametersValues["axs"] = axsParamValue                
                        
                        listDatasetAxs = []
                        
                        for localAxsIds in range(len(listDatasetName)):
                                if localAxsIds == 0:
                                        if len(axsParameters) == 1:
                                                listDatasetAxs.append(fig.add_subplot(axsNum//axsSeparation, axsSeparation, axsIds+1))
                                        else:
                                                listDatasetAxs.append(fig.add_subplot(len(axsDictParam[axsParameters[0]]), len(axsDictParam[axsParameters[1]]), axsIds+1))
                                else:
                                        listDatasetAxs.append(listDatasetAxs[0].twinx())
                                
                                
                        listDatasetAxs[0].set_title(", ".join(map(lambda x : " = ".join(x), zip_longest(axsDictParam.keys(), axsParamValue))))
   
                        for pltIds, pltParamValue in enumerate(product(*pltDictParam.values())):

                                dictParametersValues["plt"] = pltParamValue   
                                
                                aggDataPath = aggDataDir
                                aggDataFileName = ""
                                for parameters in tmp_dict_param.keys():
                                        if parameters not in aggParameters:
                                                suffixe, ids = searchDictParam[parameters.upper()]
                                                aggDataFileName += f"_{parameters.upper()}_{dictParametersValues[suffixe][ids]}"

                                USE_SAVED_DATA = 1

                                datasetArray = {}
                                datasetArraySaved = [False for _ in range(len(listDatasetName))]
                                
                                for datasetIds, datasetName in enumerate(listDatasetName):
                                        if USE_SAVED_DATA and f"{datasetName}"+aggDataFileName+".npi" in os.listdir(aggDataPath):
                                                datasetArray[datasetName] = np.load(os.path.join(aggDataPath, f"{datasetName}"+aggDataFileName+".npi"))
                                                datasetArraySaved[datasetIds] = True
                                        else:
                                                datasetArray[datasetName] = np.zeros((pltNum, datasetnHist[datasetName], metaParameterN*aggNum)) if datasetnHist[datasetName] > 1 else np.zeros((metaParameterNmeas, metaParameterN*aggNum))
                                
                                if not all(datasetArraySaved):

                                        for aggIds, aggParamValue in enumerate(product(*aggDictParam.values())):

                                                dictParametersValues["agg"] = aggParamValue   

                                                dataPath = expPath
                                                for parameters in tmp_dict_param.keys():
                                                        suffixe, ids = searchDictParam[parameters.upper()]
                                                        dataPath = os.path.join(dataPath, f"{parameters.upper()}/{dictParametersValues[suffixe][ids]}")

                                                with h5py.File(os.path.join(dataPath, f"aggregated_process.h5"), 'r') as dataFile:
                                                        for datasetIds, datasetName in enumerate(listDatasetName):
                                                                if not datasetArraySaved[datasetIds]:
                                                                        if "MCS" in pltParameters:
                                                                                if datasetnHist[datasetName] > 1:
                                                                                        datasetArray[datasetName][:,:,metaParameterN*aggIds:metaParameterN*(aggIds+1)] = dataFile[datasetName][DictParametersSlicing["MCS"]]
                                                                                else:
                                                                                        datasetArray[datasetName][:,metaParameterN*aggIds:metaParameterN*(aggIds+1)] = dataFile[datasetName][DictParametersSlicing["MCS"]]
                                                                        # if "MCS" in aggParameters:
                                                                        #         if datasetnHist[datasetName] > 1:
                                                                        #                 datasetArray[datasetName][pltIds,:,metaParameterN*aggIds:metaParameterN*(aggIds+1)] = dataFile[datasetName]
                                                                        #         else:
                                                                        #                 datasetArray[datasetName][pltIds,metaParameterN*aggIds:metaParameterN*(aggIds+1)] = dataFile[datasetName]
                                                                        



                                                print(f"{c/c_max*100: 2.1f}%", end='\r')
                                                c += 1 
                                
                                else:
                                        print(f"{c/c_max*100: 2.1f}%", end='\r')
                                        c += aggNum 
                                
                                for datasetIds, datasetName in enumerate(listDatasetName):
                                        if not datasetArraySaved[datasetIds]:
                                                np.save(os.path.join(aggDataPath, f"{datasetName}"+aggDataFileName+".npi"), datasetArray[datasetName])
                                

                                for datasetIds, datasetName in enumerate(listDatasetName):
                                        hist, bin_edges = np.histogram(datasetArray[datasetName][pltIds].flatten(), bins='auto')
                                        listDatasetAxs[datasetIds].plot((bin_edges[1:]+bin_edges[:-1])/2, np.where(hist>0, np.log10(hist), np.nan), color = colormapList[datasetName](pltIds+1))
                                        
                pltLeg = fig.legend(handles = [tuple([Line2D([], [], color = colormapList[datasetName](pltIds+1)) for datasetName in listDatasetName]) for pltIds, pltParamValue in enumerate(product(*pltDictParam.values()))], labels = [', '.join(pltParamValue) for pltParamValue in product(*pltDictParam.values())], numpoints=1, handler_map={tuple: HandlerTuple(ndivide=None)}, loc="outside lower center", ncols = pltNum)
                pltLeg.set_title(", ".join(pltDictParam.keys()))

                fig.add_artist(pltLeg)

                datasetLeg = fig.legend(handles = [Patch(color = colormapList[datasetName](pltNum+2), label = datasetName) for datasetName in listDatasetName], loc="outside upper left")
                datasetLeg.set_title("Datasets")
                        
                # fig.subplots_adjust(wspace=0.5, hspace=0.5)
                
                fig.savefig(os.path.join(figDir, "hist_" + "_".join(listDatasetName) + "_" + "_".join(map(lambda x: "=".join(x), zip_longest(figDictParam.keys(), figParamValue)))+"_"+"_".join(map(lambda x: "%".join(x), DictParametersStringSlice.items()))+".png"))
                plt.close(fig=fig)

def AggregateData(expPath, listDatasetName, pathList, metaParameterN):

        metaParameterN = int(metaParameterN)

        c_max = len(pathList)*metaParameterN
        c = 0
        c_problem = 0

        datasetShape = {}

        with h5py.File(os.path.join(pathList[0], "N/0/process.h5"), 'r') as tmpfile:
                for datasetName in listDatasetName:
                        datasetShape[datasetName] = np.shape(tmpfile[datasetName])


        for dataPath in pathList:
                Nb_sim = 0

                datasetArray = {}
                for datasetName in listDatasetName:
                        datasetArray[datasetName] = np.zeros(datasetShape[datasetName])

                for n in range(metaParameterN):
                        print(f"{c/c_max*100: 2.1f}%  {c_problem/c_max*100: 2.1f}%", end='\r')
                        try:
                                with h5py.File(os.path.join(dataPath, f"N/{n}/process.h5"), 'r') as dataFile:
                                        for datasetName in listDatasetName:
                                                datasetArray[datasetName] += dataFile[datasetName]
                                        
                                        Nb_sim += 1
                        except:
                                with open(f"{expPath}/simu.txt", 'a') as tfile:
                                        tfile.write(dataPath+'\n')
                                        c_problem += 1
                                                
                        c += 1

                with h5py.File(os.path.join(dataPath,"aggregated_process.h5"), 'a') as aggregatedFile:
                        for datasetName in listDatasetName:
                                PrintDataset(aggregatedFile, datasetName, datasetArray[datasetName]/Nb_sim)


def PrintAggregateData(expPath, figDir, listDatasetName, dict_parameters, metaParameterN, metaParameterNmeas, plotType):
        
        if plotType == 'map' and len(listDatasetName) > 1:
                sys.exit()        

        metaParameterN = int(metaParameterN)
        metaParameterNmeas = int(metaParameterNmeas)

        tmp_dict_param = {}
        for keys, values in dict_parameters.items():
                if len(values) > 1:
                        tmp_dict_param[keys.upper()] = values
        
        dict_parameters["MCS"] = [str(i) for i in range(metaParameterNmeas)]

        figParameters = ["LDENS"]
        axsParameters = ["JLP", "JLL"]
        pltParameters = []
        aggParameters = ["EV"]

        if len(axsParameters) not in [1, 2]:
                sys.exit()

        if plotType == 'map' and len(pltParameters) > 0:
                sys.exit() 

        
        if plotType == 'hist' and len(pltParameters) > 1:
                sys.exit()

        DictParametersStringSlice = {"MCS" : "::20", "JLL" : "::4", "JLP" : "::2"}
        DictParametersSlicing = {}

        for parameter, stringSlice in DictParametersStringSlice.items():
                DictParametersSlicing[parameter] = str_to_slice(stringSlice)

        for parameter, slicing in DictParametersSlicing.items():
                dict_parameters[parameter] = dict_parameters[parameter][slicing] 
        
        aggDataDir = os.path.join(expPath, "aggData")
        os.makedirs(aggDataDir, exist_ok=True)

        searchDictParam = {}

        figDictParam, searchDictParam, figNum = FormatParameters('fig', figParameters, searchDictParam, dict_parameters)
        axsDictParam, searchDictParam, axsNum = FormatParameters('axs', axsParameters, searchDictParam, dict_parameters)
        pltDictParam, searchDictParam, pltNum = FormatParameters('plt', pltParameters, searchDictParam, dict_parameters)
        aggDictParam, searchDictParam, aggNum = FormatParameters('agg', aggParameters, searchDictParam, dict_parameters)
        
        if plotType == "map":
                colormapNameList = ["seismic"]
        elif plotType == "hist":
                colormapNameList = ["Purples", "Greens"]
                colorNameList = ["blue", "orange"]
        
                colormapList = {}
                for ids, datasetName in enumerate(listDatasetName):
                        colormapList[datasetName] = colormaps[colormapNameList[ids]].resampled(pltNum+2) if pltNum > 1 else colorNameList[ids]
        elif plotType == "plot":
                colormapNameList = ["Purples", "Greens"]
                colormapList = {}
                for ids, datasetName in enumerate(listDatasetName):
                        colormapList[datasetName] = colormaps[colormapNameList[ids]].resampled(pltNum+2) if pltNum > 1 else colorNameList[ids]
        

        for params in product(*tmp_dict_param.values()):
                tmp_path = "/".join(map(lambda x : "/".join(x), zip_longest(tmp_dict_param.keys(), params)))+"/"
                pathInit = os.path.join(expPath, tmp_path)
                break
        
        c_max = figNum*axsNum*pltNum*aggNum
        c = 0

        if plotType == "map":
                datasetShape = {}

                with h5py.File(os.path.join(pathInit, "aggregated_process.h5"), 'r') as tmpfile:
                        for datasetName in listDatasetName:
                                datasetShape[datasetName] = np.shape(tmpfile[datasetName])
                
                reader = hdf5Reader(os.path.join(pathInit, "N/0"), "traj.h5", -1)
                nTad = reader.nTad
                domains = reader.domains  
                reader.Close() 
        elif plotType == 'hist':
                datasetnHist = {}

                with h5py.File(os.path.join(pathInit, "N/0/process.h5"), 'r') as tmpfile:
                        for datasetName in listDatasetName:
                                datasetnHist[datasetName] = np.shape(tmpfile[datasetName])[1] if len(np.shape(tmpfile[datasetName])) > 1 else 1
        elif plotType == "plot":
                datasetShape = {}

                with h5py.File(os.path.join(pathInit, "aggregated_process.h5"), 'r') as tmpfile:
                        for datasetName in listDatasetName:
                                datasetShape[datasetName] = np.shape(tmpfile[datasetName])
                
        

        dictParametersValues = {"fig" : None, "axs" : None, "plt" : None, "agg" : None}

        for figParamValue in product(*figDictParam.values()):

                dictParametersValues["fig"] = figParamValue
                
                if len(axsParameters) == 1:
                        axsSeparation = int(np.sqrt(axsNum))
                        while axsNum%axsSeparation > 0:
                                axsSeparation -= 1
                                                              
                        fig = plt.figure(figsize=(4*axsNum//axsSeparation, 4*axsSeparation), layout='constrained')
                else:
                        fig = plt.figure(figsize=(4*len(axsDictParam[axsParameters[1]]), 4*len(axsDictParam[axsParameters[0]])), layout='constrained')
                        
                fig.suptitle(", ".join(map(lambda x : " = ".join(x), zip_longest(figDictParam.keys(), figParamValue)))+f"\n{plotType}")
                
                listFigureAxs = []

                for axsIds, axsParamValue in enumerate(product(*axsDictParam.values())):

                        dictParametersValues["axs"] = axsParamValue                
                        
                        listDatasetAxs = []
                        
                        for localAxsIds in range(len(listDatasetName)):
                                if localAxsIds == 0:
                                        if len(axsParameters) == 1:
                                                listDatasetAxs.append(fig.add_subplot(axsNum//axsSeparation, axsSeparation, axsIds+1))
                                        else:
                                                listDatasetAxs.append(fig.add_subplot(len(axsDictParam[axsParameters[0]]), len(axsDictParam[axsParameters[1]]), axsIds+1))
                                else:
                                        listDatasetAxs.append(listDatasetAxs[0].twinx())
                        
                        listFigureAxs.append(listDatasetAxs[0])
                                
                        listDatasetAxs[0].set_title(", ".join(map(lambda x : " = ".join(x), zip_longest(axsDictParam.keys(), axsParamValue))))
                        
                        if plotType == 'plot' and pltParameters[0] == "MCS":
                                listDatasetAxs[0].set_xlim([0, 109])
   
                        for pltIds, pltParamValue in enumerate(product(*pltDictParam.values())):
                                
                                # if pltParameters[0] not in ["MCS", "None"]:
                                #         dictParametersValues["plt"] = pltParamValue

                                aggDataPath = aggDataDir
                                aggDataFileName = ""
                                for parameters in tmp_dict_param.keys():
                                        if parameters not in aggParameters:
                                                suffixe, ids = searchDictParam[parameters.upper()]
                                                aggDataFileName += f"_{parameters.upper()}_{dictParametersValues[suffixe][ids]}"

                                USE_SAVED_DATA = 0

                                datasetArray = {}
                                datasetArraySaved = [False for _ in range(len(listDatasetName))]
                                
                                for datasetIds, datasetName in enumerate(listDatasetName):
                                        if USE_SAVED_DATA and f"{datasetName}"+aggDataFileName+".npi" in os.listdir(aggDataPath):
                                                datasetArray[datasetName] = np.load(os.path.join(aggDataPath, f"{datasetName}"+aggDataFileName+".npi"))
                                                datasetArraySaved[datasetIds] = True
                                        else:
                                                datasetArray[datasetName] = np.zeros(datasetShape[datasetName])
                                                
                                if not all(datasetArraySaved):

                                        for aggParamValue in product(*aggDictParam.values()):

                                                dictParametersValues["agg"] = aggParamValue   

                                                dataPath = expPath
                                                for parameters in tmp_dict_param.keys():
                                                        suffixe, ids = searchDictParam[parameters.upper()]
                                                        dataPath = os.path.join(dataPath, f"{parameters.upper()}/{dictParametersValues[suffixe][ids]}")

                                                with h5py.File(os.path.join(dataPath, f"aggregated_process.h5"), 'r') as dataFile:
                                                        for datasetIds, datasetName in enumerate(listDatasetName):
                                                                if not datasetArraySaved[datasetIds]:
                                                                        if plotType in ["map", "plot"]:
                                                                                datasetArray[datasetName] += dataFile[datasetName]

                                                print(f"{c/c_max*100: 2.1f}%", end='\r')
                                                c += 1 
                                else:
                                        print(f"{c/c_max*100: 2.1f}%", end='\r')
                                        c += aggNum
                                
                                for datasetIds, datasetName in enumerate(listDatasetName):
                                        if not datasetArraySaved[datasetIds]:
                                                np.save(os.path.join(aggDataPath, f"{datasetName}"+aggDataFileName+".npi"), datasetArray[datasetName])
                                
                                for datasetIds, datasetName in enumerate(listDatasetName):
                                        if plotType == "map":
                                                im = listDatasetAxs[datasetIds].imshow(squareform(datasetArray[datasetName][-1]/aggNum), extent=(0, nTad, 0, nTad), origin='lower', norm=LogNorm())
                                                                        
                                                for d in domains:
                                                        if len(d) > 0:
                                                                x = [d[0], d[-1], nTad]
                                                                
                                                                y1 = [d[0], d[-1], d[-1]]
                                                                y2 = [d[0], d[0], d[0]]
                                        
                                                                listDatasetAxs[datasetIds].fill_between(x=x, y1=y1, y2=y2,  color='red', alpha=0.5, lw=0)
                                                                listDatasetAxs[datasetIds].fill_between(x=x[:2], y1=y1[:2], color='red', alpha=0.5, lw=0)

                                                listDatasetAxs[datasetIds].set_xlim([0, nTad])
                                                listDatasetAxs[datasetIds].set_ylim([0, nTad])
                                        
                                        elif plotType == 'plot':
                                                if pltParameters[0] == "MCS":
                                                        tmpDataArray = datasetArray[datasetName][int(pltParamValue[0])]/aggNum
                                                        listDatasetAxs[datasetIds].plot(np.where(tmpDataArray > 0, np.log10(tmpDataArray), np.nan), color = colormapList[datasetName](pltIds+1))

                if plotType in ["hist", "plot"]:
                        pltLeg = fig.legend(handles = [tuple([Line2D([], [], color = colormapList[datasetName](pltIds+1)) for datasetName in listDatasetName]) for pltIds, pltParamValue in enumerate(product(*pltDictParam.values()))], labels = [', '.join(pltParamValue) for pltParamValue in product(*pltDictParam.values())], numpoints=1, handler_map={tuple: HandlerTuple(ndivide=None)}, loc="outside lower center", ncols = pltNum)
                        pltLeg.set_title(", ".join(pltDictParam.keys()))

                        fig.add_artist(pltLeg)

                        datasetLeg = fig.legend(handles = [Patch(color = colormapList[datasetName](pltNum+2), label = datasetName) for datasetName in listDatasetName], loc="outside upper left")
                        datasetLeg.set_title("Datasets")
                
                
                if plotType == "map":
                        fig.colorbar(im, ax = listFigureAxs, orientation='horizontal', label = listDatasetName[0])
                
                fig.savefig(os.path.join(figDir, f"{plotType}_" + "_".join(listDatasetName) + "_" + "_".join(map(lambda x: "=".join(x), zip_longest(figDictParam.keys(), figParamValue)))+"_"+"_".join(map(lambda x: "%".join(x), DictParametersStringSlice.items()))+".png"))
                plt.close(fig=fig)


def find_exp(experience, expDir):
        verif = 0

        for tmp_experience in os.listdir(expDir):
                if tmp_experience.split("_")[0] == f'EXP{experience}':
                        expName = tmp_experience
                        verif += 1

        if verif != 1:
                print(f"The experience {experience} is not found or found multiple time")
                sys.exit() 

        return(expName)


def exp_mapping(expDir):
        def string_to_list(input):
                if not ',' in input:
                        return([input.strip()])
                else:
                        return([i.strip() for i in input.split(",")])

        dict_parameters = {}

        is_poly = "1"

        with open(os.path.join(expDir,"input_slurm.cfg"),'r') as file:
                for line in file.readlines():
                        if line.split(' = ')[0] == "Nstat":
                                meta_parameter_N = line.split(' = ')[1].strip()
                        elif line.split(' = ')[1][0] == '*':
                                pass
                        else:
                                data = string_to_list(line.split(' = ')[1])
                                if "/" in data[0]:
                                        dict_parameters[line.split(' = ')[0].upper()] = [i.split("/")[-1] for i in data]
                                else:
                                        dict_parameters[line.split(' = ')[0].upper()] = data
                                if line.split(' = ')[1].strip() == 'data/toy_domain.in':
                                        is_poly = "0"
                                if line.split(' = ')[0] == "Nmeas":
                                        meta_parameter_Nmeas = str(int(line.split(' = ')[1].strip())+1)

        return(dict_parameters, is_poly, meta_parameter_N, meta_parameter_Nmeas)


def exp_pathList(dict_parameters, metaParameterN, expDir):

        metaParameterN = int(metaParameterN)

        tmp_dict_param = {}
        for keys, values in dict_parameters.items():
                if len(values) > 1:
                        tmp_dict_param[keys.upper()] = values
        
        if metaParameterN != -1:
                tmp_dict_param['N'] = [str(i) for i in range(metaParameterN)]
        
        pathList = []
                                
        for params in product(*tmp_dict_param.values()):
                tmp_path = "/".join(map(lambda x : "/".join(x), zip_longest(tmp_dict_param.keys(), params)))+"/"
                pathList.append(os.path.join(expDir, tmp_path))
        
        return(pathList)


def PrintDataset(processFile, dataset_name, data, groupFolder = None, verbosity = None):
                
               
                tmpFolder = processFile[groupFolder] if groupFolder else processFile
                if dataset_name in list(tmpFolder.keys()):
                        del tmpFolder[dataset_name]
                
                tmpFolder.create_dataset(dataset_name, data = data)
            
                if verbosity:
                        print(f"Dataset {dataset_name} printed")

        
def MeanDataset(processFile, dataset_name, data, groupFolder = None, verbosity = None):
                
               
                tmpFolder = processFile[groupFolder] if groupFolder else processFile
                if dataset_name in list(tmpFolder.keys()):
                        tmpdata = tmpFolder[dataset_name][:]
                        del tmpFolder[dataset_name]

                        tmpFolder.create_dataset(dataset_name, data = (tmpdata+data)/2)
                
                if verbosity:
                        print(f"Dataset {dataset_name} averaged")


def conversion_kBT_to_kJ_per_mol(J):
        J = np.array([float(i) for i in J])
        T = 300
        kB = 1.380e-23
        Na = 6.022e23
        res = J * kB * Na * T * 1e-3
        return([f"{i:0.1f}" for i in res])
        
        
# Adapted from Calandrini et al. (https://doi.org/10.1051/sfn/201112010)
def msdFFT(posHist):
        N = posHist.shape[0]
        
        sqDist = np.square(posHist).sum(axis=1)
        sqDist = np.append(sqDist, 0)
        
        S2 = sum([_autoCorrFFT(posHist[:,i]) for i in range(3)])
        Q = 2*sqDist.sum()
        
        S1 = np.zeros(N)
        
        for m in range(N):
                Q -= sqDist[m-1] + sqDist[N-m]
                S1[m] = Q / (N-m)
                
        return S1 - 2*S2


def _autoCorrFFT(x):
        N = x.shape[0]

        F = np.fft.fft(x, n=2*N)
        PSD = F * F.conjugate()

        autoCorr = np.fft.ifft(PSD)

        autoCorr = autoCorr[:N]
        autoCorr = autoCorr.real

        n = N*np.ones(N) - np.arange(N)

        return autoCorr / n
