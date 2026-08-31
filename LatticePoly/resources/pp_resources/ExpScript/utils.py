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
import os
import sys
from itertools import product, zip_longest
from locale import normalize

import h5py
import matplotlib.pyplot as plt
import numpy as np
from hdf5Reader import hdf5Reader
from matplotlib import colormaps
from matplotlib.colors import CenteredNorm, LogNorm
from matplotlib.legend_handler import HandlerLine2D, HandlerPatch, HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.optimize import curve_fit
from scipy.spatial.distance import squareform
from tqdm import tqdm


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


# for i in range(144):
#     if is_non_zero_file(
#         f"/Xnfs/physbiochrom/ppuel/data/EXP60_the_phase_diagram_with_log_nmeas/tmp/process_14510593_{i + 1}.err"
#     ):
#         with open(
#             f"/Xnfs/physbiochrom/ppuel/data/EXP60_the_phase_diagram_with_log_nmeas/tmp/process_14510593_{i + 1}.err",
#             'r'
#         ) as efile:
#             if any(list(map(lambda x : "IPMI failure detected" not in x, efile.readlines()))):
#                 print(i + 1)



def conversion_density_to_yM_PRC1(D, b=20e-9):
    V_site = b**3 / np.sqrt(2)
    print(V_site)
    Na = 6.022e23
    C = D / V_site / Na  # [mol/m**-3]
    yM = C * 1e3  # [ymol/L**-3]
    return yM

def format_parameters(
    suffixe, listParameters, searchDictParam, dict_parameters, MCSLocation
):
    tmpDict = {}
    tmpNum = 1
    tmpDim = 0
    for ids, param in enumerate(listParameters):
        if param in dict_parameters.keys():
            tmpDict[param.upper()] = dict_parameters[param]
            searchDictParam[param.upper()] = (suffixe, ids)
            tmpNum *= len(dict_parameters[param])
            tmpDim += 1

    if MCSLocation == suffixe:
        if suffixe in ["fig", "axs", "plt"]:
            tmpDict["MCS"] = dict_parameters["MCS"]
        tmpNum *= len(dict_parameters["MCS"])
        tmpDim += 1

    if tmpDict == {}:
        tmpDict["None"] = [None]

    return (tmpDict, searchDictParam, tmpNum, tmpDim)


def str_to_slice(stringSlice):
    return slice(
        *list(map(lambda x: int(x) if len(x) > 0 else None, stringSlice.split(":")))
    )


# def rmData(datasetPattern, pathList, metaParameterN):
#         metaParameterN = int(metaParameterN)

#         c_max = len(pathList)*metaParameterN
#         c = 0

#         for dataPath in pathList:

#                 for n in range(metaParameterN):
#                         print(f"{c/c_max*100: 2.1f}%", end='\r')

#                         with h5py.File(os.path.join(dataPath, f"N/{n}/process.h5"), 'w') as dataFile:
#                                 # for dataset_name in dataFile.keys():
#                                 print(dataFile.keys())
#                                         # if datasetPattern in dataset_name and dataset_name != "contProbOverGDistance":
#                         break
#                         c += 1


def AggregateData(list_dataset_name, pathList, metaParameterN, verify=False):
    from LiqDensity import LiqDensity

    if type(list_dataset_name) == str:
        list_dataset_name = [list_dataset_name]

    if verify:
        tmplist_dataset_name = []

        with h5py.File(
            os.path.join(pathList[0], "aggregated_process.h5"), "a"
        ) as aggregatedFile:
            for dataset_name in list_dataset_name:
                if dataset_name not in aggregatedFile.keys():
                    tmplist_dataset_name.append(dataset_name)
        if len(tmplist_dataset_name) == 0:
            return 0

    else:
        tmplist_dataset_name = list_dataset_name

    metaParameterN = 1  # int(metaParameterN)

    c_max = len(pathList) * metaParameterN
    c = 0

    datasetShape = {}

    density = LiqDensity(
        os.path.join(pathList[0], f"N/{0}/"), os.path.join(pathList[0], f"N/{0}/")
    )
    density.Compute()
    density.Print()

    with h5py.File(os.path.join(pathList[0], "N/0/process.h5"), "r") as tmpfile:
        for dataset_name in tmplist_dataset_name:
            datasetShape[dataset_name] = np.shape(tmpfile[dataset_name])

    for dataPath in pathList:
        datasetArray = {}
        for dataset_name in tmplist_dataset_name:
            datasetArray[dataset_name] = np.zeros(datasetShape[dataset_name])

        for n in range(metaParameterN):
            print(f"{c / c_max * 100: 2.1f}%", end="\r")

            density = LiqDensity(
                os.path.join(dataPath, f"N/{n}/"), os.path.join(dataPath, f"N/{n}/")
            )
            density.Compute()
            density.Print()

            with h5py.File(
                os.path.join(dataPath, f"N/{n}/process.h5"), "r"
            ) as dataFile:
                for dataset_name in tmplist_dataset_name:
                    datasetArray[dataset_name] += dataFile[dataset_name]

            c += 1

        with h5py.File(
            os.path.join(dataPath, "aggregated_process.h5"), "a"
        ) as aggregatedFile:
            for dataset_name in tmplist_dataset_name:
                PrintDataset(
                    aggregatedFile,
                    dataset_name,
                    datasetArray[dataset_name] / metaParameterN,
                )


def DataOverBaseline(
    expPath, list_dataset_name, baseline, pathList, metaParameterN, verify=True
):
    print("DataOverBaseline : start\n")

    for dataset_name in list_dataset_name:
        AggregateData(expPath, list_dataset_name, pathList, metaParameterN, verify=verify)

    AggregateData(expPath, baseline, pathList, metaParameterN, verify=verify)

    print("DataOverBaseline : aggregate done\n")

    metaParameterN = int(metaParameterN)

    c_max = len(pathList)
    c = 0

    datasetShape = {}

    with h5py.File(os.path.join(pathList[0], "aggregated_process.h5"), "r") as tmpfile:
        baselineShape = np.shape(tmpfile[baseline])

        for dataset_name in list_dataset_name:
            datasetShape[dataset_name] = np.shape(tmpfile[dataset_name])

    print("DataOverBaseline : shape done")

    for dataPath in pathList:
        print(f"DataOverBaseline : {c / c_max * 100: 2.1f}%", end="\r")

        datasetArray = {}
        for dataset_name in list_dataset_name:
            datasetArray[dataset_name] = np.zeros(datasetShape[dataset_name])

        with h5py.File(
            os.path.join(dataPath, f"aggregated_process.h5"), "r"
        ) as dataFile:
            for dataset_name in list_dataset_name:
                binNum = baselineShape[0] + 1

                if datasetShape[dataset_name][0] == (binNum) * (binNum - 1) // 2:
                    cnt = 0

                    for i in range(binNum - 1):
                        for j in range(i + 1, binNum):
                            datasetArray[dataset_name][cnt] = (
                                dataFile[dataset_name][cnt]
                                / dataFile[baseline][j - i - 1]
                                * (binNum - (j - i))
                            )
                            cnt += 1

                else:
                    print(
                        binNum,
                        datasetShape[dataset_name][0],
                        (binNum) * (binNum + 1) // 2,
                    )
                    raise Exception()

        c += 1

        with h5py.File(
            os.path.join(dataPath, "normalized_process.h5"), "a"
        ) as normalizedFile:
            for dataset_name in list_dataset_name:
                PrintDataset(
                    normalizedFile,
                    dataset_name,
                    datasetArray[dataset_name],
                    verbosity=True,
                )


def PrintAggregateData(
    fileName, list_dataset_name, plotType, expName, logPlot=False, SimTime=False
):

    dict_translate_dataset_names: dict[str, str] = {
        "liqMean": "Proteins local density [#prot/#neighbor]",
        "polyGyration": "Gyration radius of methylated chromatin [×20 nm]",
        "PolyVolume": "Density of methylated chromatin droplet",
        "SimTime": "Simulation duration [s]",
        "RatioLiqDensityInAndOut": "Ratio between protein concentration around and far from the chromatin",
        "PixelLiqHist": "Proteins local concentration",
        "PixelLiqHist_2de_axis": "Total Intensity",
        "PixelHetHist_2de_axis": "Heterochromatin local concentration",
        "PixelLiqHist_With_Het_2de_axis" : "Puncta Intensity",
        "PixelLiqHist_Without_Het_2de_axis" : "Nuclei Intensity",
        "PixelLiqHist_summed": "Total Intensity",
        "PixelLiqHist_With_Het_summed" : "Puncta Intensity",
    }
    dict_translate_parameter_names: dict[str, str] = {
        "JLL": "Self interaction [$k_BT$]",
        "JLP": "Cross interaction [$k_BT$]",
        "JLL_VALENCY": "HP1 valency",
        "JLP_VALENCY": "Cross valency",
        "JPL_VALENCY": "Availability",
        "LDENS": "Protein density [#prot/#lattice_sites]",
        "MCS": "Monte Carlo Step",
    }

    print("Initialization...", end="\r")
    
    plt.style.use("./resources/h5py/presentation.mplstyle")

    XnfsDir = "/Xnfs/physbiochrom/ppuel/data/"
    expPath = os.path.join(XnfsDir, expName)

    dict_parameters, _, metaParameterN, metaParameterNmeas = exp_mapping(expPath)

    if type(list_dataset_name) == str:
        list_dataset_name = [list_dataset_name]

    if len(list_dataset_name) > 1:
        figDir = f"/home/ppuel/figure/{expName}/composed"
        os.makedirs(figDir, exist_ok=True)
    elif plotType == "image":
        figDir = f"/home/ppuel/figure/{expName}/image/{list_dataset_name[0]}"
        os.makedirs(figDir, exist_ok=True)
    else:
        figDir = f"/home/ppuel/figure/{expName}/{list_dataset_name[0]}"
        os.makedirs(figDir, exist_ok=True)

    # datasetSegment = {"ContProbOverGDistance" : ["All contacts", "PCG/PCG contacts", "PCG/NonPCG contacts"]}

    metaParameterN = int(metaParameterN)
    metaParameterNmeas = int(metaParameterNmeas)

    tmp_dict_param = {}
    for keys, values in dict_parameters.items():
        if len(values) > 1:
            tmp_dict_param[keys.upper()] = values

    dict_parameters["MCS"] = [str(i) for i in range(metaParameterNmeas)]

    figParameters = ["JLL_VALENCY"]
    axsParameters = ["JPL_VALENCY", "JLP"]
    pltParameters = []
    axiParameters = []
    aggParameters = []

    MCSLocation = "axi"  # MCS c'est un paramètre très particulier à mettre à part

    if len(axsParameters) not in [1, 2]:
        raise Exception()

    elif plotType == "image" and (
        len(pltParameters) != 0
        or len(axiParameters) != 2
        or len(list_dataset_name) > 1
        or MCSLocation != "agg"
    ):
        raise Exception()

    elif plotType == "map" and (len(pltParameters) != 0 or len(list_dataset_name) > 1):
        raise Exception()

    elif plotType == "hist" and (
        len(axiParameters) != 0 or len(pltParameters) > 1 or MCSLocation != "agg"
    ):
        raise Exception()

    elif plotType == "plot" and len(axiParameters) == 0 and MCSLocation != "axi":
        raise Exception()

    DictParametersStringSlice = {"MCS": ":"}
    # {"JLL_VALENCY" : "::2"} #{"MCS" : "::20", }
    DictParametersSlicing = {}

    if "MCS" not in DictParametersStringSlice.keys():
        DictParametersStringSlice["MCS"] = ":"

    for parameter, stringSlice in DictParametersStringSlice.items():
        DictParametersSlicing[parameter] = str_to_slice(stringSlice)

    for parameter, slicing in DictParametersSlicing.items():
        dict_parameters[parameter] = dict_parameters[parameter][slicing]

    if len(aggParameters) > 0:
        aggDataDir = os.path.join(expPath, "aggData")
        os.makedirs(aggDataDir, exist_ok=True)

    searchDictParam = {}

    figDictParam, searchDictParam, figNum, figDim = format_parameters(
        "fig", figParameters, searchDictParam, dict_parameters, MCSLocation
    )
    axsDictParam, searchDictParam, axsNum, axsDim = format_parameters(
        "axs", axsParameters, searchDictParam, dict_parameters, MCSLocation
    )
    pltDictParam, searchDictParam, pltNum, pltDim = format_parameters(
        "plt", pltParameters, searchDictParam, dict_parameters, MCSLocation
    )
    aggDictParam, searchDictParam, aggNum, aggDim = format_parameters(
        "agg", aggParameters, searchDictParam, dict_parameters, MCSLocation
    )
    axiDictParam, searchDictParam, axiNum, axiDim = format_parameters(
        "axi", axiParameters, searchDictParam, dict_parameters, MCSLocation
    )

    for params in product(*tmp_dict_param.values()):
        tmp_path = (
            "/".join(
                map(lambda x: "/".join(x), zip_longest(tmp_dict_param.keys(), params))
            )
            + "/"
        )
        pathInit = os.path.join(expPath, tmp_path)
        if fileName not in os.listdir(pathInit):
            print(pathInit, os.listdir(pathInit))
            raise Exception()
        break

    #  Shape definition

    datasetShape = {}

    if plotType == "map":
        with h5py.File(os.path.join(pathInit, "normalized_process.h5"), "r") as tmpfile:
            for dataset_name in list_dataset_name:
                datasetShape[dataset_name] = np.shape(tmpfile[dataset_name])

        if list_dataset_name == ["contactHiC"]:
            reader = hdf5Reader(os.path.join(pathInit, "N/0"), "traj.h5", -1)
            nTad = reader.nTad
            domains = reader.domains
            reader.Close()

    elif plotType == "image":
        datasetMin = {dataset_name: np.inf for dataset_name in list_dataset_name}
        datasetMax = {dataset_name: -np.inf for dataset_name in list_dataset_name}

        for dataset_name in list_dataset_name:
            datasetShape[dataset_name] = tuple(
                len(axiValues) for axiValues in axiDictParam.values()
            )

        if SimTime:
            SimTimeMax = -np.inf

            for params in product(*tmp_dict_param.values()):
                tmp_path = os.path.join(
                    expPath,
                    "/".join(
                        map(
                            lambda x: "/".join(x),
                            zip_longest(tmp_dict_param.keys(), params),
                        )
                    )
                    + "/",
                )
                with h5py.File(
                    os.path.join(tmp_path, "post_process.h5"), "r"
                ) as tmpfile:
                    SimTimeMax = max(SimTimeMax, tmpfile["SimTime"][1])

            SimTimeMin = np.inf

            for params in product(*tmp_dict_param.values()):
                tmp_path = os.path.join(
                    expPath,
                    "/".join(
                        map(
                            lambda x: "/".join(x),
                            zip_longest(tmp_dict_param.keys(), params),
                        )
                    )
                    + "/",
                )
                with h5py.File(
                    os.path.join(tmp_path, "post_process.h5"), "r"
                ) as tmpfile:
                    if tmpfile["SimTime"][-1] >= SimTimeMax:
                        SimTimeMin = min(SimTimeMin, tmpfile["SimTime"][-1])

            pathMCS = {}

            for params in product(*tmp_dict_param.values()):
                tmp_path = os.path.join(
                    expPath,
                    "/".join(
                        map(
                            lambda x: "/".join(x),
                            zip_longest(tmp_dict_param.keys(), params),
                        )
                    )
                    + "/",
                )
                with h5py.File(
                    os.path.join(tmp_path, "post_process.h5"), "r"
                ) as tmpfile:
                    if tmpfile["SimTime"][-1] < SimTimeMax:
                        pathMCS[os.path.join(tmp_path, fileName)] = np.nan
                    else:
                        pathMCS[os.path.join(tmp_path, fileName)] = min(
                            max(
                                int(SimTimeMax // tmpfile["SimTime"][1]) + 1,
                                int(SimTimeMin // tmpfile["SimTime"][1]),
                            ),
                            int(metaParameterNmeas),
                        )

        for params in product(*tmp_dict_param.values()):
            tmp_path = os.path.join(
                expPath,
                "/".join(
                    map(
                        lambda x: "/".join(x),
                        zip_longest(tmp_dict_param.keys(), params),
                    )
                )
                + "/",
            )

            with h5py.File(os.path.join(tmp_path, fileName), "r") as tmpfile:
                if SimTime and not np.isnan(pathMCS[os.path.join(tmp_path, fileName)]):
                    for dataset_name in list_dataset_name:
                        datasetMin[dataset_name] = min(
                            datasetMin[dataset_name],
                            np.min(
                                tmpfile[dataset_name][
                                    pathMCS[os.path.join(tmp_path, fileName)]
                                ]
                            ),
                        )
                        datasetMax[dataset_name] = max(
                            datasetMax[dataset_name],
                            np.max(
                                tmpfile[dataset_name][
                                    pathMCS[os.path.join(tmp_path, fileName)]
                                ]
                            ),
                        )
                else:
                    for dataset_name in list_dataset_name:
                        if logPlot:
                            tmpdataset = tmpfile[dataset_name][
                                DictParametersSlicing["MCS"]
                            ][tmpfile[dataset_name][DictParametersSlicing["MCS"]] > 0]
                        else:
                            tmpdataset = tmpfile[dataset_name][
                                DictParametersSlicing["MCS"]
                            ]

                        datasetMin[dataset_name] = min(
                            datasetMin[dataset_name],
                            np.min(tmpdataset) if len(tmpdataset) > 0 else +np.inf,
                        )
                        datasetMax[dataset_name] = max(
                            datasetMax[dataset_name],
                            np.max(tmpdataset) if len(tmpdataset) > 0 else -np.inf,
                        )

    elif plotType == "hist":
        with h5py.File(os.path.join(pathInit, fileName), "r") as tmpfile:
            for dataset_name in list_dataset_name:
                datasetShape[dataset_name] = (
                    len(dict_parameters["MCS"])
                    * np.prod(np.shape(tmpfile[dataset_name])[1:])
                    if np.ndim(tmpfile[dataset_name]) > 1
                    else len(dict_parameters["MCS"]),
                )

    elif plotType == "plot":
        with h5py.File(os.path.join(pathInit, fileName), "r") as tmpfile:
            for dataset_name in list_dataset_name:
                datasetShape[dataset_name] = (
                    axiNum,
                    *np.shape(tmpfile[dataset_name])[1:],
                )  #

    # Colors definition

    colormapNameList = ["plasma_r", "Greens", "Reds"]
    colorNameList = [f"C{i}" for i in range(10)]

    if plotType == "map" or plotType == "image":
        colormapNameList = ["plasma"]

    elif plotType == "hist":
        colormapList = {}
        if len(pltParameters) == 0:
            for datasetIds, dataset_name in enumerate(list_dataset_name):
                colormapList[dataset_name] = colorNameList[datasetIds]
        else:
            for datasetIds, dataset_name in enumerate(list_dataset_name):
                colormapList[dataset_name] = (
                    colormaps[colormapNameList[datasetIds]].resampled(pltNum + 2)
                    if pltNum > 1
                    else colorNameList[datasetIds]
                )

    elif plotType == "plot":
        colormapList = {}
        if len(pltParameters) == 0:
            for datasetIds, dataset_name in enumerate(list_dataset_name):
                colormapList[dataset_name] = colorNameList[datasetIds]
        else:
            for datasetIds, dataset_name in enumerate(list_dataset_name):
                colormapList[dataset_name] = (
                    colormaps[colormapNameList[datasetIds]].resampled(pltNum + 2)
                    if pltNum > 1
                    else colorNameList[datasetIds]
                )

    c_max = figNum * axsNum * pltNum * aggNum * axiNum

    tqdmProduit = tqdm(total=c_max)

    dictParametersValues = {
        "fig": None,
        "axs": None,
        "plt": None,
        "agg": None,
        "axi": None,
    }

    # if MCSLocation == 'fig':
    #         figProduct = product(*figDictParam.values(), dict_parameters["MCS"])
    # else:
    #         figProduct = product(*figDictParam.values())

    # if MCSLocation == "fig":
    #     figParameters.append("MCS")
    # elif MCSLocation == "axs":
    #     axsParameters.append("MCS")
    # elif MCSLocation == "plt":
    #     pltParameters.append("MCS")
    # elif MCSLocation == "axi":
    #     axiParameters.append("MCS")
    # elif MCSLocation == "agg":
    #     aggParameters.append("MCS")

    for figParamValue in product(*figDictParam.values()):
        dictParametersValues["fig"] = (
            figParamValue if MCSLocation != "fig" else figParamValue[:-1]
        )
        if MCSLocation == "fig":
            MCSIds = dict_parameters["MCS"].index(figParamValue[-1])

        if axsDim == 1:
            if axsNum < 4:
                axsSeparation = axsNum
            else:
                axsSeparation = int(np.sqrt(axsNum)) + 1
                while axsNum % axsSeparation > 0:
                    axsSeparation += 1

            fig = plt.figure(
                figsize=(6 * axsSeparation, 6 * axsNum // axsSeparation),
                layout="constrained",
            )

        elif axsDim == 2:
            fig = plt.figure(
                figsize=(
                    4 * len(axsDictParam[axsParameters[1]]),
                    4 * len(axsDictParam[axsParameters[0]]),
                ),
                layout="constrained",
            )

        if len(figParameters) > 0:
            fig.suptitle(
                ", ".join(
                    map(
                        lambda x: " = ".join(x),
                        zip_longest(
                            map(
                                lambda x: dict_translate_parameter_names[x],
                                figDictParam.keys(),
                            ),
                            figParamValue,
                        ),
                    )
                )
                + f"\n{plotType}"
            )

        listFigureAxs = []

        for axsIds, axsParamValue in enumerate(product(*axsDictParam.values())):
            dictParametersValues["axs"] = (
                axsParamValue[:-1] if MCSLocation == "axs" else axsParamValue
            )
            if MCSLocation == "axs":
                MCSIds = dict_parameters["MCS"].index(axsParamValue[-1])

            listDatasetAxs = []

            for datasetIds, dataset_name in enumerate(list_dataset_name):
                if datasetIds == 0:
                    if axsDim == 1:
                        listDatasetAxs.append(
                            fig.add_subplot(
                                axsNum // axsSeparation, axsSeparation, axsIds + 1
                            )
                        )
                    else:
                        listDatasetAxs.append(
                            fig.add_subplot(
                                len(axsDictParam[axsParameters[0]]),
                                len(axsDictParam[axsParameters[1]]),
                                axsIds + 1,
                            )
                        )
                else:
                    listDatasetAxs.append(listDatasetAxs[0].twinx())

                if plotType == "image":  # datasetNum == 1
                    listDatasetAxs[datasetIds].set_xlabel(
                        dict_translate_parameter_names[axiParameters[1]], weight = 'bold', fontsize = 20
                    )
                    listDatasetAxs[datasetIds].set_ylabel(
                        dict_translate_parameter_names[axiParameters[0]], weight = 'bold', fontsize = 20
                    )
                    listDatasetAxs[datasetIds].set_xticks(
                        np.arange(len(axiDictParam[axiParameters[1]])),
                        axiDictParam[axiParameters[1]],
                        rotation=90, weight = 'bold', fontsize = 20
                    )
                    listDatasetAxs[datasetIds].set_yticks(
                        np.arange(len(axiDictParam[axiParameters[0]])),
                        axiDictParam[axiParameters[0]], weight = 'bold', fontsize = 20
                    )

                elif plotType == "plot":
                    listDatasetAxs[datasetIds].tick_params(
                        axis="y",
                        labelcolor=colormapList[dataset_name](pltNum + 1)
                        if pltNum > 1
                        else colormapList[dataset_name],
                    )
                    if datasetIds > 1:
                        listDatasetAxs[datasetIds].yaxis.set_major_formatter(
                            lambda x, _: " " * 2 * (datasetIds - 1) + f"{x:0.1f}"
                        )  # axsIds
                    else:
                        if MCSLocation == "axi":
                            if len(dict_parameters["MCS"]) < 10:
                                nb_ticks = len(
                                    dict_parameters["MCS"]
                                )  # len(dict_parameters["MCS"])
                                ticks_distance = 1
                            else:
                                nb_ticks = 10  # len(dict_parameters["MCS"]) // 10
                                ticks_distance = len(dict_parameters["MCS"]) // 10 + 1

                            listDatasetAxs[datasetIds].set_xticks(
                                np.arange(nb_ticks) * ticks_distance,
                                dict_parameters["MCS"][::ticks_distance][:nb_ticks],
                            )

                            listDatasetAxs[datasetIds].set_xlabel(
                                dict_translate_parameter_names["MCS"]
                            )

                        else:
                            listDatasetAxs[datasetIds].set_xticks(
                                np.arange(len(axiDictParam[axiParameters[0]])),
                                axiDictParam[axiParameters[0]],
                            )

                            listDatasetAxs[datasetIds].set_xlabel(
                                dict_translate_parameter_names[axiParameters[0]]
                            )

                elif plotType == "hist":
                    listDatasetAxs[datasetIds].tick_params(
                        axis="y",
                        labelcolor=colormapList[dataset_name](pltNum + 1)
                        if pltNum > 1
                        else colormapList[dataset_name],
                    )
                    if datasetIds > 1:
                        listDatasetAxs[datasetIds].yaxis.set_major_formatter(
                            lambda x, _: " " * 2 * (datasetIds - 1) + f"{x:0.1f}"
                        )  # axsIds
                    elif len(axiParameters) > 1:
                        if MCSLocation == "axi":
                            listDatasetAxs[datasetIds].set_xticks(
                                np.arange(len(dict_parameters["MCS"])),
                                dict_parameters["MCS"],
                            )

                            listDatasetAxs[datasetIds].set_xlabel(
                                dict_translate_parameter_names["MCS"]
                            )

                        else:
                            listDatasetAxs[datasetIds].set_xticks(
                                np.arange(len(axiDictParam[axiParameters[0]])),
                                axiDictParam[axiParameters[0]],
                            )

                            listDatasetAxs[datasetIds].set_xlabel(
                                dict_translate_parameter_names[axiParameters[0]]
                            )

            listFigureAxs.append(listDatasetAxs[0])

            if axsDim > 0:
                listDatasetAxs[0].set_title(
                    ",\n".join(
                        map(
                            lambda x: " = ".join(x),
                            zip_longest(
                                map(
                                    lambda x: dict_translate_parameter_names[x],
                                    axsDictParam.keys(),
                                ),
                                axsParamValue,
                            ),
                        )
                    ), weight = 'bold', fontsize = 20
                )

            for pltIds, pltParamValue in enumerate(product(*pltDictParam.values())):
                if len(pltParameters) > 0:
                    dictParametersValues["plt"] = (
                        pltParamValue[:-1] if MCSLocation == "plt" else pltParamValue
                    )
                if MCSLocation == "plt":
                    MCSIds = dict_parameters["MCS"].index(pltParamValue[-1])

                USE_SAVED_DATA = 0

                if USE_SAVED_DATA:
                    aggDataPath = aggDataDir
                    aggDataFileName = ""
                    for parameters in tmp_dict_param.keys():
                        if parameters not in aggParameters:
                            suffixe, ids = searchDictParam[parameters.upper()]
                            aggDataFileName += f"_{parameters.upper()}_{dictParametersValues[suffixe][ids]}"

                    datasetArraySaved = [False for _ in range(len(list_dataset_name))]

                    for datasetIds, dataset_name in enumerate(list_dataset_name):
                        if (
                            USE_SAVED_DATA
                            and f"{dataset_name}" + aggDataFileName + ".npi"
                            in os.listdir(aggDataPath)
                        ):
                            datasetArray[dataset_name] = np.load(
                                os.path.join(
                                    aggDataPath,
                                    f"{dataset_name}" + aggDataFileName + ".npi",
                                )
                            )
                            datasetArraySaved[datasetIds] = True
                        else:
                            datasetArray[dataset_name] = np.zeros(
                                datasetShape[dataset_name]
                            )
                else:
                    datasetArraySaved = [False for _ in range(len(list_dataset_name))]

                    datasetArray = {}
                    for datasetIds, dataset_name in enumerate(list_dataset_name):
                        datasetArray[dataset_name] = np.zeros(datasetShape[dataset_name])

                if not USE_SAVED_DATA or not all(datasetArraySaved):
                    if plotType == "image":
                        for axiParamValue in product(*axiDictParam.values()):
                            dictParametersValues["axi"] = axiParamValue

                            datasetMean = {
                                dataset_name: 0 for dataset_name in list_dataset_name
                            }

                            for aggParamValue in product(*aggDictParam.values()):
                                dictParametersValues["agg"] = aggParamValue

                                dataPath = expPath
                                for parameters in tmp_dict_param.keys():
                                    suffixe, ids = searchDictParam[parameters.upper()]
                                    dataPath = os.path.join(
                                        dataPath,
                                        f"{parameters.upper()}/{dictParametersValues[suffixe][ids]}",
                                    )

                                with h5py.File(
                                    os.path.join(dataPath, fileName), "r"
                                ) as dataFile:
                                    for datasetIds, dataset_name in enumerate(
                                        list_dataset_name
                                    ):
                                        if not datasetArraySaved[datasetIds]:
                                            if SimTime:
                                                if np.isnan(
                                                    pathMCS[
                                                        os.path.join(dataPath, fileName)
                                                    ]
                                                ):
                                                    datasetMean[dataset_name] += np.nan
                                                else:
                                                    datasetMean[dataset_name] += (
                                                        dataFile[dataset_name][
                                                            pathMCS[
                                                                os.path.join(
                                                                    dataPath, fileName
                                                                )
                                                            ]
                                                        ]
                                                    )
                                            else:
                                                datasetMean[dataset_name] += (
                                                    np.mean(
                                                        dataFile[dataset_name][
                                                            DictParametersSlicing["MCS"]
                                                        ]
                                                    )
                                                    if "MCS"
                                                    in DictParametersStringSlice
                                                    else np.mean(dataFile[dataset_name])
                                                )

                                tqdmProduit.update(1)

                            for datasetIds, dataset_name in enumerate(list_dataset_name):
                                if not datasetArraySaved[datasetIds]:
                                    datasetArray[dataset_name][
                                        axiDictParam[axiParameters[0]].index(
                                            axiParamValue[0]
                                        )
                                    ][
                                        axiDictParam[axiParameters[1]].index(
                                            axiParamValue[1]
                                        )
                                    ] = datasetMean[dataset_name]

                    else:
                        for axiIds, axiParamValue in enumerate(
                            product(*axiDictParam.values())
                        ):
                            dictParametersValues["axi"] = axiParamValue

                            for aggIds, aggParamValue in enumerate(
                                product(*aggDictParam.values())
                            ):
                                dictParametersValues["agg"] = aggParamValue

                                dataPath = expPath
                                for parameters in tmp_dict_param.keys():
                                    suffixe, ids = searchDictParam[parameters.upper()]
                                    dataPath = os.path.join(
                                        dataPath,
                                        f"{parameters.upper()}/{dictParametersValues[suffixe][ids]}",
                                    )

                                with h5py.File(
                                    os.path.join(dataPath, fileName), "r"
                                ) as dataFile:
                                    for datasetIds, dataset_name in enumerate(
                                        list_dataset_name
                                    ):
                                        if not datasetArraySaved[datasetIds]:
                                            if plotType == "map":
                                                datasetArray[dataset_name] += dataFile[
                                                    dataset_name
                                                ]
                                            elif plotType == "plot":
                                                if MCSLocation == "axi":
                                                    datasetArray[dataset_name] += (
                                                        dataFile[dataset_name][
                                                            DictParametersSlicing["MCS"]
                                                        ]
                                                    )
                                                else:
                                                    datasetArray[dataset_name][
                                                        axiIds
                                                    ] += (
                                                        np.sum(
                                                            dataFile[dataset_name][
                                                                DictParametersSlicing[
                                                                    "MCS"
                                                                ]
                                                            ]
                                                        )
                                                        if MCSLocation == "agg"
                                                        else dataFile[dataset_name][
                                                            MCSIds
                                                        ]
                                                    )
                                            elif plotType == "hist":
                                                datasetArray[dataset_name][
                                                    (0 + aggIds)
                                                    * datasetShape[dataset_name] : (
                                                        1 + aggIds
                                                    )
                                                    * datasetShape[dataset_name]
                                                ] = (
                                                    dataFile[dataset_name][
                                                        DictParametersSlicing["MCS"]
                                                    ].flatten()
                                                    if MCSLocation == "agg"
                                                    else dataFile[dataset_name]
                                                )

                                if aggNum != 1:
                                    tqdmProduit.update(aggNum)

                            if axiNum != 1:
                                tqdmProduit.update(axiNum)

                else:
                    tqdmProduit.update(aggNum)

                if USE_SAVED_DATA:
                    for datasetIds, dataset_name in enumerate(list_dataset_name):
                        if not datasetArraySaved[datasetIds]:
                            np.save(
                                os.path.join(
                                    aggDataPath,
                                    f"{dataset_name}" + aggDataFileName + ".npi",
                                ),
                                datasetArray[dataset_name],
                            )

                if plotType == "map":
                    if list_dataset_name == ["contactHiC"]:
                        im = listDatasetAxs[0].imshow(
                            squareform(datasetArray["contactHiC"]) / aggNum,
                            origin="lower",
                            extent=(0, nTad, 0, nTad),
                            norm=LogNorm(),
                        )

                        for d in domains:
                            if len(d) > 0:
                                x = [d[0], d[-1], nTad]

                                y1 = [d[0], d[-1], d[-1]]
                                y2 = [d[0], d[0], d[0]]

                                listDatasetAxs[0].fill_between(
                                    x=x, y1=y1, y2=y2, color="red", alpha=0.5, lw=0
                                )
                                listDatasetAxs[0].fill_between(
                                    x=x[:2], y1=y1[:2], color="red", alpha=0.5, lw=0
                                )

                        listDatasetAxs[datasetIds].set_xlim([0, nTad])
                        listDatasetAxs[datasetIds].set_ylim([0, nTad])

                    else:
                        for datasetIds, dataset_name in enumerate(list_dataset_name):
                            im = listDatasetAxs[datasetIds].imshow(
                                datasetArray[dataset_name] / aggNum, origin="lower"
                            )

                elif plotType == "image":
                    for datasetIds, dataset_name in enumerate(list_dataset_name):
                        im = listDatasetAxs[datasetIds].imshow(
                            np.where(
                                datasetArray[dataset_name] > 0,
                                np.log2(datasetArray[dataset_name]),
                                np.nan,
                            )
                            if logPlot
                            else 2
                            * np.logical_and(
                                datasetArray[dataset_name] > 3,
                                datasetArray[dataset_name] < 5,
                            )
                            + np.logical_and(
                                datasetArray[dataset_name] > 2,
                                datasetArray[dataset_name] < 6,
                            ),
                            origin="lower",
                            cmap=colormapNameList[0],
                            vmin=np.where(
                                datasetMin[dataset_name] > 0,
                                np.log2(datasetMin[dataset_name]),
                                -4,
                            )
                            if logPlot
                            else 0,  # datasetMin[dataset_name],
                            vmax=np.log2(datasetMax[dataset_name])
                            if logPlot
                            else 2,  # datasetMax[dataset_name]
                        )
                        # im = listDatasetAxs[datasetIds].imshow(np.where(datasetMax[dataset_name]/datasetArray[dataset_name] < 100, np.log10(datasetMax[dataset_name]/datasetArray[dataset_name]), np.nan), origin = 'lower', cmap = colormapNameList[0], vmin = np.log10(datasetMax[dataset_name]/datasetMax[dataset_name]), vmax = np.log10(100))#datasetMax[dataset_name]/datasetMin[dataset_name]))

                elif plotType == "plot":
                    if pltParameters == []:
                        for datasetIds, dataset_name in enumerate(list_dataset_name):
                            if len(datasetShape[dataset_name]) == 1:
                                tmpDataArray = datasetArray[dataset_name] / aggNum
                                listDatasetAxs[datasetIds].plot(
                                    np.where(
                                        tmpDataArray > 0, np.log10(tmpDataArray), np.nan
                                    )
                                    if logPlot
                                    else tmpDataArray,
                                    color=colormapList[dataset_name],
                                )
                            elif len(datasetShape[dataset_name]) == 2:
                                for i in range(datasetShape[dataset_name][1]):
                                    tmpDataArray = (
                                        datasetArray[dataset_name][:, i] / aggNum
                                    )
                                    listDatasetAxs[datasetIds].plot(
                                        np.where(
                                            tmpDataArray > 0,
                                            np.log10(tmpDataArray),
                                            np.nan,
                                        )
                                        if logPlot
                                        else tmpDataArray,
                                        color=colormapList[dataset_name](i),
                                    )

                    else:
                        for datasetIds, dataset_name in enumerate(list_dataset_name):
                            if len(datasetShape[dataset_name]) == 1:
                                tmpDataArray = datasetArray[dataset_name] / aggNum
                                listDatasetAxs[datasetIds].plot(
                                    np.where(
                                        tmpDataArray > 0, np.log10(tmpDataArray), np.nan
                                    )
                                    if logPlot
                                    else tmpDataArray,
                                    color=colormapList[dataset_name](pltIds + 1),
                                )

                            elif len(datasetShape[dataset_name]) == 2:
                                tmpDataArray = (
                                    np.mean(datasetArray[dataset_name][-10:, :], axis=0)
                                    / aggNum
                                )
                                listDatasetAxs[datasetIds].plot(
                                    np.where(
                                        tmpDataArray > 0, np.log10(tmpDataArray), np.nan
                                    )
                                    if logPlot
                                    else tmpDataArray,
                                    color=colormapList[dataset_name](pltIds + 1),
                                )

                            # Some analyse to be done systematically
                            #

                            # L = len(tmpDataArray)

                            # listDatasetAxs[datasetIds].plot(np.arange(L), np.log((max(tmpDataArray)-tmpDataArray)/(max(tmpDataArray)-min(tmpDataArray))), 'o')

                            # def f(x, a, b, c, d):
                            #     return(a / (1. + np.exp(-c * (x - d))) + b)
                            # try:
                            #     coef, _  = curve_fit(f, np.arange(len(tmpDataArray)), tmpDataArray, method="trf")
                            #     listDatasetAxs[datasetIds].plot(np.arange(len(tmpDataArray)), f(np.arange(len(tmpDataArray)), *coef), '--')
                            # except Exception:
                            #     pass
                elif plotType == "hist":
                    datasetHist = {}
                    datasetbinEdges = {}
                    for datasetIds, dataset_name in enumerate(list_dataset_name):
                        datasetHist[dataset_name], datasetbinEdges[dataset_name] = (
                            np.histogram(
                                datasetArray[dataset_name], density=True, bins=5
                            )
                        )

                    if pltParameters == []:
                        for datasetIds, dataset_name in enumerate(list_dataset_name):
                            listDatasetAxs[datasetIds].plot(
                                (
                                    datasetbinEdges[dataset_name][1:]
                                    + datasetbinEdges[dataset_name][:-1]
                                )
                                / 2,
                                np.where(
                                    datasetHist[dataset_name] > 0,
                                    np.log10(datasetHist[dataset_name]),
                                    np.nan,
                                )
                                if logPlot
                                else datasetHist[dataset_name],
                                color=colormapList[dataset_name],
                            )
                        for datasetIds, dataset_name in enumerate(list_dataset_name):
                            listDatasetAxs[datasetIds].plot(
                                (
                                    datasetbinEdges[dataset_name][1:]
                                    + datasetbinEdges[dataset_name][:-1]
                                )
                                / 2,
                                np.where(
                                    datasetHist[dataset_name] > 0,
                                    np.log10(datasetHist[dataset_name]),
                                    np.nan,
                                )
                                if logPlot
                                else datasetHist[dataset_name],
                                "o--",
                                color=colormapList[dataset_name](pltIds + 1),
                                markersize=12,
                            )

        if plotType in ["hist", "plot"]:
            if len(pltParameters) == 0:
                # pltLeg = fig.legend(handles = [tuple([Line2D([], [], color = colormapList[dataset_name]) for dataset_name in list_dataset_name]) for pltIds, pltParamValue in enumerate(product(*pltDictParam.values()))], labels = [', '.join(pltParamValue) for pltParamValue in product(*pltDictParam.values())], numpoints=1, handler_map={tuple: HandlerTuple(ndivide=None)}, loc="outside lower center", ncols = pltNum)
                # pltLeg.set_title(", ".join(map(lambda x: dict_translate_parameter_names[x], pltDictParam.keys())))

                # fig.add_artist(pltLeg)

                datasetLeg = fig.legend(
                    handles=[
                        Patch(
                            color=colormapList[dataset_name],
                            label=dict_translate_dataset_names[dataset_name],
                        ) for dataset_name in list_dataset_name
                    ],
                    loc="outside upper left",
                )
                datasetLeg.set_title("Datasets")
            else:
                pltLeg = fig.legend(
                    handles=[
                        tuple(
                            [
                                Line2D(
                                    [], [], color=colormapList[dataset_name](pltIds + 1)
                                )
                                for dataset_name in list_dataset_name
                            ]
                        )
                        for pltIds, pltParamValue in enumerate(
                            product(*pltDictParam.values())
                        )
                    ],
                    labels=[
                        ", ".join(pltParamValue)
                        for pltParamValue in product(*pltDictParam.values())
                    ],
                    numpoints=1,
                    handler_map={tuple: HandlerTuple(ndivide=None)},
                    loc="outside lower center",
                    ncols=pltNum,
                )
                pltLeg.set_title(
                    ", ".join(
                        map(
                            lambda x: dict_translate_parameter_names[x],
                            pltDictParam.keys(),
                        )
                    )
                )

                fig.add_artist(pltLeg)

                datasetLeg = fig.legend(
                    handles=[
                        Patch(
                            color=colormapList[dataset_name](pltNum + 2),
                            label=dict_translate_dataset_names[dataset_name],
                        )
                        for dataset_name in list_dataset_name
                    ],
                    loc="outside upper left",
                )
                datasetLeg.set_title("Datasets")

        if plotType == "map" or plotType == "image":
            cbar = fig.colorbar(
                im,
                ax=listFigureAxs,
                orientation="horizontal",
                label=dict_translate_dataset_names[list_dataset_name[0]]
                if list_dataset_name[0] in dict_translate_dataset_names.keys()
                else list_dataset_name[0],
            )
            

        tmpFigString = (
            "_"
            + "_".join(
                map(
                    lambda x: "=".join(x),
                    zip_longest(figDictParam.keys(), figParamValue),
                )
            )
            if len(figParameters) > 0
            else ""
        )
        tmpSliString = (
            "_"
            + "_".join(map(lambda x: "%".join(x), DictParametersStringSlice.items()))
            if len(DictParametersStringSlice) > 0
            else ""
        )
        tmpImaString = "_".join(axiParameters) + "_" if plotType == "image" else ""
        tmpTimeString = "_byTime" if SimTime else ""
        fig.savefig(
            os.path.join(
                figDir,
                f"{fileName}_{plotType}_"
                + tmpImaString
                + "_".join(list_dataset_name)
                + tmpFigString
                + tmpSliString
                + tmpTimeString
                + ".png",
            )
        )
        plt.close(fig=fig)

    tqdmProduit.close()


def find_exp(experience: int, expDir: str):
    expName = ""

    verif = 0

    for tmp_experience in os.listdir(expDir):
        if tmp_experience.split("_")[0] == f"EXP{experience}":
            expName = tmp_experience
            verif += 1

    if verif != 1:
        print(f"The experience {experience} is not found or found multiple time")
        sys.exit()

    return expName


def exp_mapping(expDir: str):
    def string_to_list(input: str):
        if "," not in input:
            return [input.strip()]
        else:
            return [i.strip() for i in input.split(",")]

    dict_parameters: dict[str, list[str]] = {}

    is_poly: str = "1"
    meta_parameter_N: str = "0"
    meta_parameter_Nmeas: str = "0"

    with open(os.path.join(expDir, "input_slurm.cfg"), "r") as file:
        for line in file.readlines():
            if line.split(" = ")[0] == "Nstat":
                meta_parameter_N = line.split(" = ")[1].strip()
            if line.split(" = ")[1][0] == "*":
                pass
            else:
                data = string_to_list(line.split(" = ")[1])
                if "/" in data[0]:
                    dict_parameters[line.split(" = ")[0].upper()] = [
                        i.split("/")[-1] for i in data
                    ]
                else:
                    dict_parameters[line.split(" = ")[0].upper()] = data
                if line.split(" = ")[1].strip() == "data/toy_domain.in":
                    is_poly = "0"
                if line.split(" = ")[0] == "Nmeas":
                    meta_parameter_Nmeas = str(int(line.split(" = ")[1].strip()) + 1)

    return (dict_parameters, is_poly, meta_parameter_N, meta_parameter_Nmeas)


def exp_path_list(
    dict_parameters: dict[str, list[str]], metaParameterN: str, expDir: str
):

    N = int(metaParameterN)

    tmp_dict_param: dict[str, list[str]] = {}
    for keys, values in dict_parameters.items():
        if len(values) > 1:
            tmp_dict_param[keys.upper()] = values

    if N != -1:
        tmp_dict_param["N"] = [str(i) for i in range(metaParameterN)]

    pathList: list[str] = []

    for params in product(*tmp_dict_param.values()):
        tmp_path = (
            "/".join(
                map(lambda x: "/".join(x), zip_longest(tmp_dict_param.keys(), params))
            )
            + "/"
        )
        pathList.append(os.path.join(expDir, tmp_path))

    return pathList


def PrintDataset(processFile, dataset_name, data, groupFolder=None, verbosity=None):

    tmpFolder = processFile[groupFolder] if groupFolder else processFile
    if dataset_name in list(tmpFolder.keys()):
        del tmpFolder[dataset_name]

    tmpFolder.create_dataset(dataset_name, data=data)

    if verbosity:
        print(f"Dataset {dataset_name} printed")


def conversion_kBT_to_kJ_per_mol(J):
    J = np.array([float(i) for i in J])
    T = 300
    kB = 1.380e-23
    Na = 6.022e23
    res = J * kB * Na * T * 1e-3
    return [f"{i:0.1f}" for i in res]


# Adapted from Calandrini et al. (https://doi.org/10.1051/sfn/201112010)
def msdFFT(posHist):
    N = posHist.shape[0]

    sqDist = np.square(posHist).sum(axis=1)
    sqDist = np.append(sqDist, 0)

    S2 = sum([_autoCorrFFT(posHist[:, i]) for i in range(3)])
    Q = 2 * sqDist.sum()

    S1 = np.zeros(N)

    for m in range(N):
        Q -= sqDist[m - 1] + sqDist[N - m]
        S1[m] = Q / (N - m)

    return S1 - 2 * S2


def _autoCorrFFT(x):
    N = x.shape[0]

    F = np.fft.fft(x, n=2 * N)
    PSD = F * F.conjugate()

    autoCorr = np.fft.ifft(PSD)

    autoCorr = autoCorr[:N]
    autoCorr = autoCorr.real

    n = N * np.ones(N) - np.arange(N)

    return autoCorr / n


# def AggregateHist(expPath, list_dataset_name, pathList, metaParameterN, metaParameterNmeas):

#         metaParameterN = int(metaParameterN)
#         metaParameterNmeas = int(metaParameterNmeas)

#         c_max = len(pathList)*metaParameterN
#         c = 0
#         c_problem = 0

#         datasetnHist = {}

#         with h5py.File(os.path.join(pathList[0], "N/0/process.h5"), 'r') as tmpfile:
#                 for dataset_name in list_dataset_name:
#                         datasetnHist[dataset_name] = np.shape(tmpfile[dataset_name])[-1] if len(np.shape(tmpfile[dataset_name]))>1 else 1


#         for dataPath in pathList:
#                 Nb_sim = 0

#                 datasetArray = {}
#                 for dataset_name in list_dataset_name:
#                         datasetArray[dataset_name] = np.zeros((metaParameterNmeas, datasetnHist[dataset_name], metaParameterN)) if datasetnHist[dataset_name] > 1 else np.zeros((metaParameterNmeas, metaParameterN))

#                 for n in range(metaParameterN):
#                         print(f"{c/c_max*100: 2.1f}%  {c_problem/c_max*100: 2.1f}%", end='\r')
#                         try:

#                                 with h5py.File(os.path.join(dataPath, f"N/{n}/process.h5"), 'r') as dataFile:
#                                         for dataset_name in list_dataset_name:
#                                                 if datasetnHist[dataset_name] > 1 :
#                                                         datasetArray[dataset_name][:,:,n] = dataFile[dataset_name]
#                                                 else:
#                                                         datasetArray[dataset_name][:,n] = dataFile[dataset_name]

#                         except Exception as inst:

#                                 with open(f"{expPath}/simu.txt", 'a') as tfile:
#                                         tfile.write(dataPath+f" {n}\n"+str(inst)+'\n\n')
#                                         c_problem += 1


#                         c += 1

#                 with h5py.File(os.path.join(dataPath,"aggregated_process.h5"), 'a') as aggregatedFile:
#                         for dataset_name in list_dataset_name:
#                                 PrintDataset(aggregatedFile, dataset_name, datasetArray[dataset_name])


# def PrintAggregateHist(expPath, figDir, list_dataset_name, dict_parameters, metaParameterN, metaParameterNmeas):

#         metaParameterN = int(metaParameterN)
#         metaParameterNmeas = int(metaParameterNmeas)

#         tmp_dict_param = {}
#         for keys, values in dict_parameters.items():
#                 if len(values) > 1:
#                         tmp_dict_param[keys.upper()] = values

#         dict_parameters["MCS"] = [str(i) for i in range(metaParameterNmeas)]

#         # listAllParameters = ["MCS"] + list(tmp_dict_param.keys())

#         figParameters = ["LDENS"]
#         axsParameters = ["JLP", "JLL"]         #["colomn", "row"] or ["unique"]
#         pltParameters = ["MCS"]                # size 1
#         aggParameters = ["EV"]

#         DictParametersStringSlice = {"MCS" : "::100", "JLL" : "::5", "JLP" : "::2"}
#         DictParametersSlicing = {}

#         for parameter, stringSlice in DictParametersStringSlice.items():
#                 DictParametersSlicing[parameter] = str_to_slice(stringSlice)

#         for parameter, slicing in DictParametersSlicing.items():
#                 dict_parameters[parameter] = dict_parameters[parameter][slicing]

#         aggDataDir = os.path.join(expPath, "aggData")
#         os.makedirs(aggDataDir, exist_ok=True)

#         searchDictParam = {}

#         figDictParam, searchDictParam, figNum = format_parameters('fig', figParameters, searchDictParam, dict_parameters)
#         axsDictParam, searchDictParam, axsNum = format_parameters('axs', axsParameters, searchDictParam, dict_parameters)
#         pltDictParam, searchDictParam, pltNum = format_parameters('plt', pltParameters, searchDictParam, dict_parameters)
#         aggDictParam, searchDictParam, aggNum = format_parameters('agg', aggParameters, searchDictParam, dict_parameters)

#         colormapNameList = ["Purples", "Greens"]

#         colormapList = {}
#         for ids, dataset_name in enumerate(list_dataset_name):
#                 colormapList[dataset_name] = colormaps[colormapNameList[ids]].resampled(pltNum+2)


#         for params in product(*tmp_dict_param.values()):
#                 tmp_path = "/".join(map(lambda x : "/".join(x), zip_longest(tmp_dict_param.keys(), params)))+"/"
#                 pathInit = os.path.join(expPath, tmp_path)
#                 break

#         datasetnHist = {}

#         with h5py.File(os.path.join(pathInit, "N/0/process.h5"), 'r') as tmpfile:
#                 for dataset_name in list_dataset_name:
#                         datasetnHist[dataset_name] = np.shape(tmpfile[dataset_name])[1] if len(np.shape(tmpfile[dataset_name])) > 1 else 1


#         c_max = figNum*axsNum*pltNum*aggNum
#         c = 0

#         dictParametersValues = {"fig" : None, "axs" : None, "plt" : None, "agg" : None}

#         for figParamValue in product(*figDictParam.values()):

#                 dictParametersValues["fig"] = figParamValue

#                 if len(axsParameters) == 1:
#                         axsSeparation = int(np.sqrt(axsNum))
#                         while axsNum%axsSeparation > 0:
#                                 axsSeparation -= 1
#                         axsSeparation = axsNum // axsSeparation

#                         fig = plt.figure(figsize=(4*axsNum//axsSeparation, 4*axsSeparation), layout='constrained')
#                 else:
#                         fig = plt.figure(figsize=(4*len(axsDictParam[axsParameters[1]]), 4*len(axsDictParam[axsParameters[0]])), layout='constrained')

#                 fig.suptitle(", ".join(map(lambda x : " = ".join(x), zip_longest(figDictParam.keys(), figParamValue))))

#                 for axsIds, axsParamValue in enumerate(product(*axsDictParam.values())):

#                         dictParametersValues["axs"] = axsParamValue

#                         listDatasetAxs = []

#                         for localAxsIds in range(len(list_dataset_name)):
#                                 if localAxsIds == 0:
#                                         if len(axsParameters) == 1:
#                                                 listDatasetAxs.append(fig.add_subplot(axsNum//axsSeparation, axsSeparation, axsIds+1))
#                                         else:
#                                                 listDatasetAxs.append(fig.add_subplot(len(axsDictParam[axsParameters[0]]), len(axsDictParam[axsParameters[1]]), axsIds+1))
#                                 else:
#                                         listDatasetAxs.append(listDatasetAxs[0].twinx())


#                         listDatasetAxs[0].set_title(", ".join(map(lambda x : " = ".join(x), zip_longest(axsDictParam.keys(), axsParamValue))))

#                         for pltIds, pltParamValue in enumerate(product(*pltDictParam.values())):

#                                 dictParametersValues["plt"] = pltParamValue

#                                 aggDataPath = aggDataDir
#                                 aggDataFileName = ""
#                                 for parameters in tmp_dict_param.keys():
#                                         if parameters not in aggParameters:
#                                                 suffixe, ids = searchDictParam[parameters.upper()]
#                                                 aggDataFileName += f"_{parameters.upper()}_{dictParametersValues[suffixe][ids]}"

#                                 USE_SAVED_DATA = 1

#                                 datasetArray = {}
#                                 datasetArraySaved = [False for _ in range(len(list_dataset_name))]

#                                 for datasetIds, dataset_name in enumerate(list_dataset_name):
#                                         if USE_SAVED_DATA and f"{dataset_name}"+aggDataFileName+".npi" in os.listdir(aggDataPath):
#                                                 datasetArray[dataset_name] = np.load(os.path.join(aggDataPath, f"{dataset_name}"+aggDataFileName+".npi"))
#                                                 datasetArraySaved[datasetIds] = True
#                                         else:
#                                                 datasetArray[dataset_name] = np.zeros((pltNum, datasetnHist[dataset_name], metaParameterN*aggNum)) if datasetnHist[dataset_name] > 1 else np.zeros((metaParameterNmeas, metaParameterN*aggNum))

#                                 if not all(datasetArraySaved):

#                                         for aggIds, aggParamValue in enumerate(product(*aggDictParam.values())):

#                                                 dictParametersValues["agg"] = aggParamValue

#                                                 dataPath = expPath
#                                                 for parameters in tmp_dict_param.keys():
#                                                         suffixe, ids = searchDictParam[parameters.upper()]
#                                                         dataPath = os.path.join(dataPath, f"{parameters.upper()}/{dictParametersValues[suffixe][ids]}")

#                                                 with h5py.File(os.path.join(dataPath, f"aggregated_process.h5"), 'r') as dataFile:
#                                                         for datasetIds, dataset_name in enumerate(list_dataset_name):
#                                                                 if not datasetArraySaved[datasetIds]:
#                                                                         if "MCS" in pltParameters:
#                                                                                 if datasetnHist[dataset_name] > 1:
#                                                                                         datasetArray[dataset_name][:,:,metaParameterN*aggIds:metaParameterN*(aggIds+1)] = dataFile[dataset_name][DictParametersSlicing["MCS"]]
#                                                                                 else:
#                                                                                         datasetArray[dataset_name][:,metaParameterN*aggIds:metaParameterN*(aggIds+1)] = dataFile[dataset_name][DictParametersSlicing["MCS"]]
#                                                                         # if "MCS" in aggParameters:
#                                                                         #         if datasetnHist[dataset_name] > 1:
#                                                                         #                 datasetArray[dataset_name][pltIds,:,metaParameterN*aggIds:metaParameterN*(aggIds+1)] = dataFile[dataset_name]
#                                                                         #         else:
#                                                                         #                 datasetArray[dataset_name][pltIds,metaParameterN*aggIds:metaParameterN*(aggIds+1)] = dataFile[dataset_name]


#                                                 print(f"{c/c_max*100: 2.1f}%", end='\r')
#                                                 c += 1

#                                 else:
#                                         print(f"{c/c_max*100: 2.1f}%", end='\r')
#                                         c += aggNum

#                                 for datasetIds, dataset_name in enumerate(list_dataset_name):
#                                         if not datasetArraySaved[datasetIds]:
#                                                 np.save(os.path.join(aggDataPath, f"{dataset_name}"+aggDataFileName+".npi"), datasetArray[dataset_name])


#                                 for datasetIds, dataset_name in enumerate(list_dataset_name):
#                                         hist, bin_edges = np.histogram(datasetArray[dataset_name][pltIds].flatten(), bins='auto')
#                                         listDatasetAxs[datasetIds].plot((bin_edges[1:]+bin_edges[:-1])/2, np.where(hist>0, np.log10(hist), np.nan), color = colormapList[dataset_name](pltIds+1))

#                 pltLeg = fig.legend(handles = [tuple([Line2D([], [], color = colormapList[dataset_name](pltIds+1)) for dataset_name in list_dataset_name]) for pltIds, pltParamValue in enumerate(product(*pltDictParam.values()))], labels = [', '.join(pltParamValue) for pltParamValue in product(*pltDictParam.values())], numpoints=1, handler_map={tuple: HandlerTuple(ndivide=None)}, loc="outside lower center", ncols = pltNum)
#                 pltLeg.set_title(", ".join(pltDictParam.keys()))

#                 fig.add_artist(pltLeg)

#                 datasetLeg = fig.legend(handles = [Patch(color = colormapList[dataset_name](pltNum+2), label = dataset_name) for dataset_name in list_dataset_name], loc="outside upper left")
#                 datasetLeg.set_title("Datasets")

#                 # fig.subplots_adjust(wspace=0.5, hspace=0.5)

#                 fig.savefig(os.path.join(figDir, "hist_" + "_".join(list_dataset_name) + "_" + "_".join(map(lambda x: "=".join(x), zip_longest(figDictParam.keys(), figParamValue)))+"_"+"_".join(map(lambda x: "%".join(x), DictParametersStringSlice.items()))+".png"))
#                 plt.close(fig=fig)
