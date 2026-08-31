import importlib
import itertools
import os
import pickle
import subprocess
import sys
from itertools import product, zip_longest

import h5py
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import seaborn as sns
from hdf5Reader import hdf5Reader
from matplotlib import axis, colormaps, markers
from matplotlib.colors import LogNorm
from mpl_toolkits import axisartist
from mpl_toolkits.axes_grid1 import host_subplot
from scipy.spatial.distance import squareform
from utils import *
plt.rcParams.update({"text.usetex": True, "font.family": "Cambria"})


plt.style.use("./resources/h5py/presentation.mplstyle")

expNum = 46
XnfsDir = "/Xnfs/physbiochrom/ppuel/data/"
expName = find_exp(expNum, XnfsDir)
expPath = os.path.join(XnfsDir, expName)


dict_parameters, _, metaParameterN, metaParameterNmeas = exp_mapping(expPath)

pathList = exp_pathList(dict_parameters, -1, expPath)

# with h5py.File(os.path.join(pathList[0], "N/0/process.h5"), 'r') as tmpFile:
#         print(tmpFile.keys())

listDatasetName = ["contactHiC"]
listDatasetRealName = ["Contact between PcG region and bulk", "Intra PcG contact"]

# baseline = "contactProb"
# DataOverBaseline(expPath, listDatasetName, baseline, pathList, metaParameterN)

# fileName = "normalized_process.h5"
fileName = "normalized_process.h5"

# for datasetName in listDatasetName:
#         PrintAggregateData(fileName, datasetName, "image", expName)


# with h5py.File(os.path.join(pathList[0], "normalized_process.h5"), 'r') as hfile:
#         for dataset in hfile.keys():
#                 print(hfile[dataset].shape, dataset)

# PrintAggregateData(fileName, listDatasetName, "map", expName, logPlot=True)


def supplementary_figure1():

    font_dirs = ["./resources/fonts/"]  # The path to the custom font file.
    font_files = fm.findSystemFonts(fontpaths=font_dirs)
    CambriaManager = fm.fontManager
    for font_file in font_files:
        print(font_file)
        CambriaManager.addfont(font_file)

    titleFont = fm.FontProperties(
        size=plt.rcParams["font.size"] * 2,
        fname="/home/ppuel/Simulation/LatticePoly/LatticePoly/resources/fonts/Cambria-Font-For-Linux-Bold.ttf",
        weight="bold",
    )
    print(titleFont.get_family())
    jll_list = [f"{i * 0.2 + 0.2:0.1f}" for i in range(12)]  #
    jlp_list = [f"{i * 0.2 + 0.8:0.1f}" for i in range(3)]
    ldens_list = [f"{i * 0.00135 + 0.00270:0.5f}" for i in range(3)]  #
    ev_list = [f"{i * 2 + 8:0d}" for i in range(3)]
    percentage_list = [r"0.5", r"0.75", r"1"]

    colormapNameList = ["Purples", "Reds", "Greens"]

    listDatasetName = ["maximeDropNum_PC", "maximeDropMean_PC", "maximeLiqFraction_PC"]
    listDatasetRealName = [
        "Number of foci",
        "Volume of foci [\\#PRC1]",
        "PRC1 free nucleoplasmic fraction",
    ]

    datasetMin = {datasetName: np.inf for datasetName in listDatasetName}
    datasetMax = {datasetName: -np.inf for datasetName in listDatasetName}

    c = 1
    cmax = 3 * 3 * 3 * 12

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(12):
                    with h5py.File(
                        os.path.join(
                            expPath,
                            f"JLL/{jll_list[l]}/JLP/{jlp_list[k]}/EV/{ev_list[j]}/LDENS/{ldens_list[i]}/aggregated_process.h5",
                        )
                    ) as hfile:
                        for datasetName in listDatasetName:
                            tmpValue = np.mean(hfile[datasetName][-10:])
                            if datasetName == "maximeLiqFraction_PC":
                                datasetMin[datasetName] = min(
                                    datasetMin[datasetName], 1 - tmpValue
                                )
                                datasetMax[datasetName] = max(
                                    datasetMax[datasetName], 1 - tmpValue
                                )
                            else:
                                datasetMin[datasetName] = min(
                                    datasetMin[datasetName], tmpValue
                                )
                                datasetMax[datasetName] = max(
                                    datasetMax[datasetName], tmpValue
                                )

    fileName = "aggregated_process.h5"

    fig = plt.figure(figsize=(12, 12), layout="constrained")
    subFigList = fig.subfigures(
        5, 4, width_ratios=[1, 4, 4, 4], height_ratios=[1, 4, 4, 4, 1]
    )

    for i in range(5):
        for j in range(4):
            localSubfig = subFigList[i][j]

            if i > 0 and j > 0 and i < 4:
                localSubplotList = localSubfig.subplots(3, 1, sharex=True)
                localSubfig.supylabel("$E_{P-H}$ [$k_BT$]")

                datasetArray = {}
                for datasetName in listDatasetName:
                    datasetArray[datasetName] = np.zeros((3, 12))

                for k in range(3):
                    for l in range(12):
                        with h5py.File(
                            os.path.join(
                                expPath,
                                f"JLL/{jll_list[l]}/JLP/{jlp_list[k]}/EV/{ev_list[j - 1]}/LDENS/{ldens_list[i - 1]}/aggregated_process.h5",
                            )
                        ) as hfile:
                            for datasetName in listDatasetName:
                                if datasetName == "maximeLiqFraction_PC":
                                    datasetArray[datasetName][k][l] = 1 - np.mean(
                                        hfile[datasetName][-10:]
                                    )
                                else:
                                    datasetArray[datasetName][k][l] = np.mean(
                                        hfile[datasetName][-10:]
                                    )
                        print(f"{c / cmax * 100:2.1f}%", end="\r")
                        c += 1
                datasetImShow = {}

                def tmp(x, _):
                    print(x)
                    return jll_list[int(x)]

                localSubplotList[2].xaxis.set_major_locator(
                    ticker.FixedLocator([1, 4, 7, 10])
                )
                localSubplotList[2].xaxis.set_minor_locator(
                    ticker.MultipleLocator(1, offset=0)
                )
                localSubplotList[2].xaxis.set_major_formatter(
                    ticker.FixedFormatter(np.array(jll_list)[1::3])
                )
                localSubplotList[2].xaxis.set_minor_formatter(ticker.NullFormatter())

                for datasetIds, datasetName in enumerate(listDatasetName):
                    datasetImShow[datasetIds] = localSubplotList[datasetIds].imshow(
                        datasetArray[datasetName],
                        origin="lower",
                        cmap=colormapNameList[datasetIds],
                        vmin=datasetMin[datasetName],
                        vmax=datasetMax[datasetName],
                    )
                    localSubplotList[datasetIds].set_yticks(
                        [i for i in range(3)], jlp_list
                    )
                    if datasetIds < 2:
                        localSubplotList[datasetIds].tick_params(
                            axis="x", which="both", bottom=False
                        )

                localSubplotList[2].set_xlabel("$E_{P-P}$ [$k_BT$]")
            elif i == 4 and j > 0:
                ax = localSubfig.add_subplot()  # axes(rect = (0.2,0.6,0.8,0.5))
                localSubfig.colorbar(
                    datasetImShow[j - 1],
                    cax=ax,
                    orientation="horizontal",
                    label=listDatasetRealName[j - 1],
                )
            elif i == 0 and j > 0:
                localSubfig.text(
                    0.5,
                    0.5,
                    "$E_{SH}$ = " + f"${ev_list[j - 1]}~k_BT$",
                    horizontalalignment="center",
                    verticalalignment="center",
                    font=titleFont,
                )
            elif i > 0 and i < 4 and j == 0:
                localSubfig.text(
                    0.5,
                    0.5,
                    "$R_{P/H}$ = " + percentage_list[i - 1],
                    rotation=+90,
                    horizontalalignment="center",
                    verticalalignment="center",
                    font=titleFont,
                )

    fig.savefig(f"/home/ppuel/figure/{expName}/supplementary_figure_paper_V1.svg")

    plt.close(fig=fig)


def figure_1():

    listDatasetName = ["maximeDropNum_PC", "maximeDropMean_PC", "maximeLiqFraction_PC"]
    listDatasetRealName = ["Number of droplets", "Size of droplets", "Gaz fraction"]

    fileName = "aggregated_process.h5"

    jll_list = ["1.0", "1.2", "1.4", "1.6", "1.8"]  #
    jlp_list = "0.8"
    ldens = "0.00270"
    ev = "8"

    ldens_star = "0.00540"

    yM_list = conversion_density_to_yM_PRC1(np.array([0.00270, 0.00540]))

    markerList = ["o", "*"]

    l = len(jll_list)

    # datasetLim = [(4,13.25), (125, 275), (0, 0.30)]
    datasetColor = ["blue", "red", "green"]

    fig = plt.figure(figsize=(8, 6), layout="constrained")

    host = host_subplot(111, axes_class=axisartist.Axes, figure=fig)

    par1 = host.twinx()
    par2 = host.twinx()

    par2.axis["right"] = par2.new_fixed_axis(loc="right", offset=(60, 0))

    par1.axis["right"].toggle(all=True)
    par2.axis["right"].toggle(all=True)

    datasetArray = {}

    for datasetName in listDatasetName:
        datasetArray[datasetName] = np.zeros(l)

    for i in range(l):
        with h5py.File(
            os.path.join(
                expPath,
                f"JLL/{jll_list[i]}/JLP/{jlp_list}/EV/{ev}/LDENS/{ldens}/aggregated_process.h5",
            )
        ) as hfile:
            for datasetName in listDatasetName:
                if datasetName == "maximeLiqFraction_PC":
                    datasetArray[datasetName][i] = 1 - np.mean(hfile[datasetName][-10])
                else:
                    datasetArray[datasetName][i] = np.mean(hfile[datasetName][-10])

    (p1,) = host.plot(datasetArray["maximeDropNum_PC"], "o-", color=datasetColor[0])
    (p2,) = par1.plot(
        16 * np.pi / 3 * datasetArray["maximeDropMean_PC"] ** 3,
        "o-",
        color=datasetColor[1],
    )
    (p3,) = par2.plot(datasetArray["maximeLiqFraction_PC"], "o-", color=datasetColor[2])

    with h5py.File(
        os.path.join(
            expPath,
            f"JLL/{jll_list[0]}/JLP/{jlp_list}/EV/{ev}/LDENS/{ldens_star}/aggregated_process.h5",
        )
    ) as hfile:
        par3 = host.scatter(
            [0],
            np.mean(hfile["maximeDropNum_PC"][-10]),
            marker=markerList[1],
            color=datasetColor[0],
        )
        par1.scatter(
            [0],
            16 * np.pi / 3 * np.mean(hfile["maximeDropMean_PC"][-10]) ** 3,
            marker=markerList[1],
            color=datasetColor[1],
        )
        par2.scatter(
            [0],
            1 - np.mean(hfile["maximeLiqFraction_PC"][-10]),
            marker=markerList[1],
            color=datasetColor[2],
        )

    host.set(ylim=(4, 13.25), xlabel="$E_{P-P}$", ylabel="Number of foci")
    par1.set(ylim=(125, 275), ylabel="Size of foci")
    par2.set(ylim=(0, 0.30), ylabel="Gaz fraction")

    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())
    par2.axis["right"].label.set_color(p3.get_color())

    host.set_xticks(np.arange(l), jll_list)
    datasetLeg = fig.legend(
        handles=[
            Patch(color=datasetColor[datasetIds], label=datasetName)
            for datasetIds, datasetName in enumerate(listDatasetRealName)
        ],
        loc="outside upper left",
    )
    datasetLeg.set_title("Datasets")

    fig.add_artist(datasetLeg)

    LdensLeg = fig.legend(
        handles=[par3, p1],
        labels=[f"{yM_list[0]:0.1f} $\\mu$M", f"{yM_list[1]:0.1f} $\\mu$M"],
        loc="outside upper right",
    )
    LdensLeg.set_title("PRC1 Concentration")

    fig.savefig(f"/home/ppuel/figure/{expName}/figure_paper3")

    plt.close(fig=fig)


def figure_1_V2():

    listDatasetName = [
        "maximeDropNum_PC",
        "maximeDropVolHist_PC",
        "maximeLiqFraction_PC",
    ]
    listDatasetRealName = [
        "Number of foci",
        "Size of foci",
        "PRC1 free\nnucleoplasmic fraction",
    ]

    jll_list = ["1.0", "1.2", "1.4", "1.6", "1.8"]
    jlp_list = "0.8"
    ldens = "0.00270"
    ev = "8"

    ldens_star = "0.00540"

    markerList = ["o", "*"]

    l_jll = len(jll_list)

    ldens_list = ["0.00270", "0.00405", "0.00540"]
    percentage_list = [r"0.5", r"0.75", r"1"]
    yM_list = conversion_density_to_yM_PRC1(
        np.fromiter(map(float, ldens_list), dtype=np.float32)
    )

    l_ldens = len(ldens_list)
    jll = "1.0"

    datasetColor = ["darkblue", "darkred", "darkgreen"]

    fig, axs_matrix = plt.subplots(
        nrows=3,
        ncols=2,
        sharex="col",
        figsize=(1.5 * 8, 1 * 7),
        # layout='constrained',
        width_ratios=[5, 3],
        height_ratios=[2, 3, 2],
    )

    ## ----------- First column

    datasetArray = {}

    for datasetName in listDatasetName:
        if datasetName == "maximeDropVolHist_PC":
            datasetArray[datasetName] = []
            sum_droplet = []
        else:
            datasetArray[datasetName] = np.zeros(l_jll)

    for i in range(l_jll):
        with h5py.File(
            os.path.join(
                expPath,
                f"JLL/{jll_list[i]}/JLP/{jlp_list}/EV/{ev}/LDENS/{ldens}/aggregated_process.h5",
            )
        ) as hfile:
            for datasetName in listDatasetName:
                if datasetName == "maximeLiqFraction_PC":
                    datasetArray[datasetName][i] = 1 - np.mean(hfile[datasetName][-10:])
                elif datasetName == "maximeDropVolHist_PC":
                    datasetArray[datasetName].append(
                        [
                            size + 1
                            for size, count in enumerate(
                                np.sum(hfile[datasetName][-10:], axis=0)
                            )
                            for _ in range(int(count * 40))
                        ]
                    )
                    sum_droplet.append(np.sum(hfile[datasetName][-10:]))
                else:
                    datasetArray[datasetName][i] = np.mean(hfile[datasetName][-10])

    widthList = [(sum_droplet[i]) / (max(sum_droplet)) * 0.75 for i in range(l_jll)]

    axs = axs_matrix[:, 0]

    (p1,) = axs[0].plot(datasetArray["maximeDropNum_PC"], "o-", color=datasetColor[0])

    # v1  = axs[1].violinplot(datasetArray["maximeDropVolHist_PC"], positions = [i for i in range(l_jll)], quantiles = [[0.25, 0.5, 0.75] for _ in range(l_jll)], showextrema = False, widths = widthList)
    sns.violinplot(
        datasetArray["maximeDropVolHist_PC"],
        ax=axs[1],
        color="darkred",
        density_norm="area",
        linecolor="red",
        width=1,
    )

    _ = axs[2].plot(datasetArray["maximeLiqFraction_PC"], "o-", color=datasetColor[2])

    # for violinBody in v1['bodies']:
    #         violinBody.set_facecolor("darkred")
    #         violinBody.set_edgecolor("red")
    #         violinBody.set_alpha(0.9)
    # violinQuantile = v1['cquantiles']
    # violinQuantile.set_color("black")
    # violinQuantile.set_alpha(0.9)

    lp = 0

    axs[0].set_title(
        "\\textbf{Mapping of $\\mathbf{E_{P-P}}$}\n PRC1/H3K27 ratio = 0.5"
    )
    axs[0].set(ylim=(4.5, 7.5))
    axs[0].set_ylabel(
        ylabel="Number of foci",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )
    axs[1].set(ylim=(-100, 1500))
    axs[1].set_ylabel(
        "Volume of foci\n[\\# PRC1]",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )
    axs[2].set(ylim=(0, 0.15), xlabel="$\\mathbf{E_{P-P}}$ [$k_BT$]")
    axs[2].set_ylabel(
        "PRC1\nFree nucleoplasmic\nfraction",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )

    axs[0].yaxis.label.set_color(datasetColor[0])
    axs[1].yaxis.label.set_color(datasetColor[1])
    axs[2].yaxis.label.set_color(datasetColor[2])

    axs[2].set_xticks([i for i in range(l_jll)], jll_list)

    ## ----------- second column

    axs = axs_matrix[:, 1]

    datasetArray = {}

    for datasetName in listDatasetName:
        if datasetName == "maximeDropVolHist_PC":
            datasetArray[datasetName] = []
            sum_droplet = []
        else:
            datasetArray[datasetName] = np.zeros(l_ldens)

    for i in range(l_ldens):
        with h5py.File(
            os.path.join(
                expPath,
                f"JLL/{jll}/JLP/{jlp_list}/EV/{ev}/LDENS/{ldens_list[i]}/aggregated_process.h5",
            )
        ) as hfile:
            for datasetName in listDatasetName:
                if datasetName == "maximeLiqFraction_PC":
                    datasetArray[datasetName][i] = 1 - np.mean(hfile[datasetName][-10:])
                elif datasetName == "maximeDropVolHist_PC":
                    datasetArray[datasetName].append(
                        [
                            size + 1
                            for size, count in enumerate(
                                np.sum(hfile[datasetName][-10:], axis=0)
                            )
                            for _ in range(int(count * 40))
                        ]
                    )
                    sum_droplet.append(np.sum(hfile[datasetName][-10:]))
                else:
                    datasetArray[datasetName][i] = np.mean(hfile[datasetName][-10:])

    widthList = [(sum_droplet[i]) / (max(sum_droplet)) * 0.75 for i in range(l_ldens)]

    (p1,) = axs[0].plot(datasetArray["maximeDropNum_PC"], "o-", color=datasetColor[0])

    sns.violinplot(
        datasetArray["maximeDropVolHist_PC"],
        ax=axs[1],
        color="darkred",
        density_norm="area",
        linecolor="red",
        width=1,
    )

    # v1  = axs[1].violinplot(datasetArray["maximeDropVolHist_PC"], positions = [i for i in range(l_ldens)], quantiles = [[0.25, 0.5, 0.75] for _ in range(l_ldens)], showextrema = False, widths = widthList)
    _ = axs[2].plot(datasetArray["maximeLiqFraction_PC"], "o-", color=datasetColor[2])

    # for violinBody in v1['bodies']:
    #         violinBody.set_facecolor("darkred")
    #         violinBody.set_edgecolor("red")
    #         violinBody.set_alpha(0.9)
    # violinQuantile = v1['cquantiles']
    # violinQuantile.set_color("black")
    # violinQuantile.set_alpha(0.9)

    axs[0].set_title("\\textbf{Mapping of PRC1/H3K27 ratio}\n$E_{P-P}$ = 1 $k_BT$")
    axs[0].set(ylim=(6.5, 13.5))
    axs[0].set_ylabel(
        "Number of foci",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )
    axs[1].set(ylim=(-100, 1500))
    axs[1].set_ylabel(
        "Volume of foci\n[\\# PRC1]",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )  # ylim=(125, 275),
    axs[2].set(ylim=(0.125, 0.30), xlabel="\\textbf{PRC1/H3K27 ratio}")
    axs[2].set_ylabel(
        "PRC1\nFree nucleoplasmic\nfraction",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )

    for i in range(3):
        axs[i].yaxis.label.set_color(datasetColor[i])
        axs[i].tick_params(
            axis="y", labelleft=False, left=False, labelright=True, right=True
        )
        axs[i].yaxis.set_label_position("right")

    axs[2].set_xticks([i for i in range(l_ldens)], percentage_list)

    # datasetLeg = fig.legend(handles = [Patch(color = datasetColor[datasetIds], label = datasetName) for datasetIds, datasetName in enumerate(listDatasetRealName)], loc="outside upper left")
    # fig.add_artist(datasetLeg)

    # LdensLeg = fig.legend(handles = [p1, s1], labels=[f"{yM_list[0]:0.1f} $\\mu$M", f"{yM_list[1]:0.1f} $\\mu$M"], loc="outside upper right")
    # LdensLeg.set_title("PRC1 Concentration")
    # fig.align_labels()
    for ax in fig.axes:
        ax.tick_params(direction="in")
    plt.subplots_adjust(wspace=0.5)
    fig.savefig(f"/home/ppuel/figure/{expName}/PRC1_paper_figure_1")

    plt.close(fig=fig)

    # with h5py.File(os.path.join(expPath, f"JLL/{jll_list[0]}/JLP/{jlp_list}/EV/{ev}/LDENS/{ldens_star}/aggregated_process.h5")) as hfile:

    #         maximeDropVolHist_PC = [size+1 for size, count in enumerate(np.sum(hfile["maximeDropVolHist_PC"][-10:], axis = 0)) for _ in range(int(count*40))]

    #         s1 = axs[0].scatter([0], np.mean(hfile["maximeDropNum_PC"][-10]),
    #                             marker=markerList[1], color = datasetColor[0])
    #         v2 = axs[1].violinplot(maximeDropVolHist_PC, positions = [0], showextrema = False, quantiles = [0.25, 0.5, 0.75])

    #         for violinBody in v2['bodies']:
    #                 violinBody.set_facecolor("orange")
    #                 violinBody.set_edgecolor("darkorange")
    #                 violinBody.set_alpha(0.5)
    #         violinQuantile = v2['cquantiles']
    #         violinQuantile.set_color("black")
    #         violinQuantile.set_alpha(0.5)

    #         _  = axs[2].scatter([0], 1-np.mean(hfile['maximeLiqFraction_PC'][-10]),
    #                             marker=markerList[1], color = datasetColor[2])


def figure_1_V3():

    listDatasetName = [
        "maximeDropNum_PC",
        "maximeDropVolHist_PC",
        "maximeLiqFraction_PC",
    ]
    listDatasetRealName = [
        "Number of foci",
        "Size of foci",
        "PRC1 free\nnucleoplasmic fraction",
    ]

    jll_list = ["1.0", "1.2", "1.4", "1.6", "1.8"]
    jlp = "0.8"
    ldens = "0.00270"
    ev = "8"

    ldens_star = "0.00540"

    markerList = ["o", "*"]

    l_jll = len(jll_list)

    jlp_list = ["1.2", "1.0", "0.8"]

    ldens_list = ["0.00540", "0.00405", "0.00270"]

    l_jlp = len(jlp_list)

    jll = "1.0"

    datasetColor = ["darkblue", "darkred", "darkgreen"]

    fig, axs_matrix = plt.subplots(
        nrows=3,
        ncols=2,
        sharex="col",
        figsize=(1.5 * 8, 1 * 7),
        # layout='constrained',
        width_ratios=[5, 3],
        height_ratios=[2, 3, 2],
    )

    ## ----------- First column

    datasetArray = {}

    for datasetName in listDatasetName:
        if datasetName == "maximeDropVolHist_PC":
            datasetArray[datasetName] = []
            sum_droplet = []
        else:
            datasetArray[datasetName] = np.zeros(l_jll)

    for i in range(l_jll):
        with h5py.File(
            os.path.join(
                expPath,
                f"JLL/{jll_list[i]}/JLP/{jlp}/EV/{ev}/LDENS/{ldens}/aggregated_process.h5",
            )
        ) as hfile:
            for datasetName in listDatasetName:
                if datasetName == "maximeLiqFraction_PC":
                    datasetArray[datasetName][i] = 1 - np.mean(hfile[datasetName][-10])
                elif datasetName == "maximeDropVolHist_PC":
                    datasetArray[datasetName].append(
                        [
                            size + 1
                            for size, count in enumerate(
                                np.sum(hfile[datasetName][-10:], axis=0)
                            )
                            for _ in range(int(count * 40))
                        ]
                    )
                    sum_droplet.append(np.sum(hfile[datasetName][-10:]))
                else:
                    datasetArray[datasetName][i] = np.mean(hfile[datasetName][-10])

    widthList = [(sum_droplet[i]) / (max(sum_droplet)) * 0.75 for i in range(l_jll)]

    axs = axs_matrix[:, 0]

    (p1,) = axs[0].plot(datasetArray["maximeDropNum_PC"], "o-", color=datasetColor[0])

    print(datasetArray["maximeDropNum_PC"]*400)
    
    # v1  = axs[1].violinplot(datasetArray["maximeDropVolHist_PC"], positions = [i for i in range(l_jll)], quantiles = [[0.25, 0.5, 0.75] for _ in range(l_jll)], showextrema = False, widths = widthList)
    sns.violinplot(
        datasetArray["maximeDropVolHist_PC"],
        ax=axs[1],
        color="darkred",
        density_norm="area",
        linecolor="red",
        width=1,
    )

    _ = axs[2].plot(datasetArray["maximeLiqFraction_PC"], "o-", color=datasetColor[2])

    # for violinBody in v1['bodies']:
    #         violinBody.set_facecolor("darkred")
    #         violinBody.set_edgecolor("red")
    #         violinBody.set_alpha(0.9)
    # violinQuantile = v1['cquantiles']
    # violinQuantile.set_color("black")
    # violinQuantile.set_alpha(0.9)

    lp = 0

    axs[0].set_title(
        "\\textbf{Mapping of $\\mathbf{E_{P-P}}$}\n PRC1/H3K27 ratio = 0.5"
    )
    axs[0].set(ylim=(4.5, 7.5))
    axs[0].set_ylabel(
        ylabel="Number of foci",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )
    axs[1].set(ylim=(-100, 1500))
    axs[1].set_ylabel(
        "Volume of foci\n[\\# PRC1]",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )
    axs[2].set(ylim=(0, 0.15), xlabel="$\\mathbf{E_{P-P}}$ [$k_BT$]")
    axs[2].set_ylabel(
        "PRC1\nFree nucleoplasmic\nfraction",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )

    axs[0].yaxis.label.set_color(datasetColor[0])
    axs[1].yaxis.label.set_color(datasetColor[1])
    axs[2].yaxis.label.set_color(datasetColor[2])

    axs[2].set_xticks([i for i in range(l_jll)], jll_list)

    ## ----------- second column

    axs = axs_matrix[:, 1]

    datasetArray = {}

    for datasetName in listDatasetName:
        if datasetName == "maximeDropVolHist_PC":
            datasetArray[datasetName] = []
            sum_droplet = []
        else:
            datasetArray[datasetName] = np.zeros(l_jlp)

    for i in range(len(ldens_list)):
        with h5py.File(
            os.path.join(
                expPath,
                f"JLL/{jll}/JLP/{jlp}/EV/{ev}/LDENS/{ldens_list[i]}/aggregated_process.h5",
            )
        ) as hfile:
            for datasetName in listDatasetName:
                if datasetName == "maximeLiqFraction_PC":
                    datasetArray[datasetName][i] = 1 - np.mean(hfile[datasetName][-10])
                elif datasetName == "maximeDropVolHist_PC":
                    datasetArray[datasetName].append(
                        [
                            size + 1
                            for size, count in enumerate(
                                np.sum(hfile[datasetName][-10:], axis=0)
                            )
                            for _ in range(int(count * 40))
                        ]
                    )
                    sum_droplet.append(np.sum(hfile[datasetName][-10:]))
                else:
                    datasetArray[datasetName][i] = np.mean(hfile[datasetName][-10])

    widthList = [(sum_droplet[i]) / (max(sum_droplet)) * 0.75 for i in range(l_jlp)]

    (p1,) = axs[0].plot(datasetArray["maximeDropNum_PC"], "o-", color=datasetColor[0])

    print(datasetArray["maximeDropNum_PC"]*400)
    
    sns.violinplot(
        datasetArray["maximeDropVolHist_PC"],
        ax=axs[1],
        color="darkred",
        density_norm="area",
        linecolor="red",
        width=1,
    )

    # v1  = axs[1].violinplot(datasetArray["maximeDropVolHist_PC"], positions = [i for i in range(l_ldens)], quantiles = [[0.25, 0.5, 0.75] for _ in range(l_ldens)], showextrema = False, widths = widthList)
    _ = axs[2].plot(datasetArray["maximeLiqFraction_PC"], "o-", color=datasetColor[2])

    # for violinBody in v1['bodies']:
    #         violinBody.set_facecolor("darkred")
    #         violinBody.set_edgecolor("red")
    #         violinBody.set_alpha(0.9)
    # violinQuantile = v1['cquantiles']
    # violinQuantile.set_color("black")
    # violinQuantile.set_alpha(0.9)

    axs[0].set_title("\\textbf{Mapping of PRC1/H3K27 ratio}\n$E_{P-P}$ = 1 $k_BT$")
    # axs[0].set(ylim=(3,22))
    axs[0].set_ylabel(
        "Number of foci",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )
    # axs[1].set(ylim=(-100, 1100))
    axs[1].set_ylabel(
        "Volume of foci\n[\\# PRC1]",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )  # ylim=(125, 275),
    # axs[2].set(ylim=(0.005, 0.15), xlabel="\\textbf{PRC1/H3K27 ratio}")
    axs[2].set_ylabel(
        "PRC1\nFree nucleoplasmic\nfraction",
        horizontalalignment="center",
        verticalalignment="center",
        labelpad=lp,
    )

    for i in range(3):
        axs[i].yaxis.label.set_color(datasetColor[i])
        axs[i].tick_params(
            axis="y", labelleft=False, left=False, labelright=True, right=True
        )
        axs[i].yaxis.set_label_position("right")

    axs[2].set_xticks([i for i in range(l_jlp)], jlp_list, rotation=-90)

    # datasetLeg = fig.legend(handles = [Patch(color = datasetColor[datasetIds], label = datasetName) for datasetIds, datasetName in enumerate(listDatasetRealName)], loc="outside upper left")
    # fig.add_artist(datasetLeg)

    # LdensLeg = fig.legend(handles = [p1, s1], labels=[f"{yM_list[0]:0.1f} $\\mu$M", f"{yM_list[1]:0.1f} $\\mu$M"], loc="outside upper right")
    # LdensLeg.set_title("PRC1 Concentration")
    # fig.align_labels()
    for ax in fig.axes:
        ax.tick_params(direction="in")
    plt.subplots_adjust(wspace=0.5)
    ## fig.savefig(f"/home/ppuel/figure/{expName}/PRC1_paper_figure_1_V3")

    plt.close(fig=fig)

    # with h5py.File(os.path.join(expPath, f"JLL/{jll_list[0]}/JLP/{jlp_list}/EV/{ev}/LDENS/{ldens_star}/aggregated_process.h5")) as hfile:

    #         maximeDropVolHist_PC = [size+1 for size, count in enumerate(np.sum(hfile["maximeDropVolHist_PC"][-10:], axis = 0)) for _ in range(int(count*40))]

    #         s1 = axs[0].scatter([0], np.mean(hfile["maximeDropNum_PC"][-10]),
    #                             marker=markerList[1], color = datasetColor[0])
    #         v2 = axs[1].violinplot(maximeDropVolHist_PC, positions = [0], showextrema = False, quantiles = [0.25, 0.5, 0.75])

    #         for violinBody in v2['bodies']:
    #                 violinBody.set_facecolor("orange")
    #                 violinBody.set_edgecolor("darkorange")
    #                 violinBody.set_alpha(0.5)
    #         violinQuantile = v2['cquantiles']
    #         violinQuantile.set_color("black")
    #         violinQuantile.set_alpha(0.5)

    #         _  = axs[2].scatter([0], 1-np.mean(hfile['maximeLiqFraction_PC'][-10]),
    #                             marker=markerList[1], color = datasetColor[2])

figure_1_V3()
# c_max = 12*3*3
# c = 0


# percentile_list = np.array([99, 98, 96, 94, 92, 90, 85, 80, 75, 70, 65, 60, 50])


# for jll in jll_list:
#         i = jll_list.index(jll)

#         for ldens in ldens_list:
#                 l = ldens_list.index(ldens)
#                 fig = plt.figure(figsize=(12,12))
#                 ax = fig.add_subplot()

#                 for jlp in jlp_list:
#                         j = jlp_list.index(jlp)
#                         with h5py.File(f"{expPath}/aggregated_data.h5", 'r') as dfile:
#                                 liqHist = dfile[f"liqHist_JLL_{jll}_JLP_{jlp}_LDENS_{ldens}"]
#                                 weights = np.mean(liqHist[-10:,:], axis=0)

#                         print(f"{c/c_max*100: 2.1f}%", end='\r')
#                         c+=1

#                         ax.plot(percentile_list, np.percentile(np.arange(Nsite_per_cut)-1, percentile_list, method='inverted_cdf', weights = weights), color = plasma(i+1), label = f"{jll}")
#                 fig.legend()
#                 fig.suptitle(f"JLL : {jll}, LDENS : {ldens}")
#                 fig.savefig(f"{figDir}/LiqPercentile_JLL_{jll}_LDENS_{ldens}.png")
#                 plt.close(fig)


# lenght, width = 24,12

# fig1 = plt.figure(figsize=(lenght, width))
# fig2 = plt.figure(figsize=(lenght, width))
# fig3 = plt.figure(figsize=(lenght, width))
# fig4 = plt.figure(figsize=(lenght, width))

# liqFractionEquilibrium_matrix  = np.load(f"/home/ppuel/data/{exp}/liqFractionEquilibrium_matrix.npy" )
# liqDropNumEquilibrium_matrix   = np.load(f"/home/ppuel/data/{exp}/liqDropNumEquilibrium_matrix.npy"  )
# liqDropSizeEquilibrium_matrix  = np.load(f"/home/ppuel/data/{exp}/liqDropSizeEquilibrium_matrix.npy" )
# polyGyrationEquilibrium_matrix = np.load(f"/home/ppuel/data/{exp}/polyGyrationEquilibrium_matrix.npy")


# liqFraction_matrix = np.load(f"/home/ppuel/data/{exp}/liqFraction_matrix.npy")
# liqDropNum_matrix = np.load(f"/home/ppuel/data/{exp}/liqDropNum_matrix.npy")
# liqDropSize_matrix = np.load(f"/home/ppuel/data/{exp}/liqDropSize_matrix.npy")
# polyGyration_matrix = np.load(f"/home/ppuel/data/{exp}/polyGyration_matrix.npy")

# fig1.suptitle(f'{np.min(liqFractionEquilibrium_matrix[:,:,:,:,3]):0.2f} - {np.max(liqFractionEquilibrium_matrix[:,:,:,:,3]):0.2f}')
# fig2.suptitle(f'{np.min(liqDropNumEquilibrium_matrix[:,:,:,:,3]):0.2f} - {np.max(liqDropNumEquilibrium_matrix[:,:,:,:,3]):0.2f}')
# fig3.suptitle(f'{np.min(liqDropSizeEquilibrium_matrix[:,:,:,:,3]):0.2f} - {np.max(liqDropSizeEquilibrium_matrix[:,:,:,:,3]):0.2f}')
# fig4.suptitle(f'{np.min(polyGyrationEquilibrium_matrix[:,:,:,:,3]):0.2f} - {np.max(polyGyrationEquilibrium_matrix[:,:,:,:,3]):0.2f}')

# LiqPolySimTime_matrix = np.load(f"/home/ppuel/data/{exp}/LiqPolySimTime_matrix.npy")


# np.save(f"/home/ppuel/data/{exp}/liqFractionEquilibrium_matrix.npy" , liqFractionEquilibrium_matrix )
# np.save(f"/home/ppuel/data/{exp}/liqDropNumEquilibrium_matrix.npy"  , liqDropNumEquilibrium_matrix  )
# np.save(f"/home/ppuel/data/{exp}/liqDropSizeEquilibrium_matrix.npy" , liqDropSizeEquilibrium_matrix )
# np.save(f"/home/ppuel/data/{exp}/polyGyrationEquilibrium_matrix.npy", polyGyrationEquilibrium_matrix)


# LiqPolySimTime_matrix[i,j,k,l] = np.array(pfile["SimTime"])[-1]

# LiqPolySimTime_matrix_reshape = np.zeros((12*3*3,3))

# for i in range(3):
#         LiqPolySimTime_matrix_reshape[:,i] = LiqPolySimTime_matrix[:,i,:,:].flatten()

# for i in range(12):
#         ax.boxplot(LiqPolySimTime_matrix_reshape//60)
#         ax.set_xticks([1,2,3],jlp_list_conv)

# ax.set_xlabel('P <-> H')
# ax.set_ylabel('Simulated time')
# fig.suptitle('Simulated time\nper P-H energy of interaction')
# fig.savefig(f"/home/ppuel/data/{exp}/LiqPolySimTime_matrix_P-H.png")
# plt.close('all')


# ax1 = fig1.add_subplot(3,3,k+3*l+1)
# ax2 = fig2.add_subplot(3,3,k+3*l+1)
# ax3 = fig3.add_subplot(3,3,k+3*l+1)
# ax4 = fig4.add_subplot(3,3,k+3*l+1)

# im1 = ax1.imshow(np.transpose(liqFractionEquilibrium_matrix[:, :, k, l, 3]), origin = 'lower', cmap = 'plasma', vmin = np.min(liqFractionEquilibrium_matrix[:,:,:,:,3]), vmax = np.max(liqFractionEquilibrium_matrix[:,:,:,:,3]))
# im2 = ax2.imshow(np.transpose(liqDropNumEquilibrium_matrix[:, :, k, l, 3]), origin = 'lower', cmap = 'plasma', vmin = np.min(liqDropNumEquilibrium_matrix[:,:,:,:,3]), vmax = np.max(liqDropNumEquilibrium_matrix[:,:,:,:,3]))
# im3 = ax3.imshow(np.transpose(liqDropSizeEquilibrium_matrix[:, :, k, l, 3]), origin = 'lower', cmap = 'plasma', vmin = np.min(liqDropSizeEquilibrium_matrix[:,:,:,:,3]), vmax = np.max(liqDropSizeEquilibrium_matrix[:,:,:,:,3]))
# im4 = ax4.imshow(np.transpose(polyGyrationEquilibrium_matrix[:, :, k, l, 3]), origin = 'lower', cmap = 'plasma', vmin = np.min(polyGyrationEquilibrium_matrix[:,:,:,:,3]), vmax = np.max(polyGyrationEquilibrium_matrix[:,:,:,:,3]))

# ax1.set_xlabel("P <-> P")
# ax1.set_ylabel("P <-> H")
# ax1.set_xticks([i for i in range(12)], np.array(jll_list_conv))
# ax1.set_yticks([i for i in range(3)], np.array(jlp_list_conv))
# ax1.set_title(f"P% = {percentage_list[l]}, EV = {EV_list_conv[k]}")

# ax2.set_xlabel("P <-> P")
# ax2.set_ylabel("P <-> H")
# ax2.set_xticks([i for i in range(12)], np.array(jll_list_conv))
# ax2.set_yticks([i for i in range(3)], np.array(jlp_list_conv))
# ax2.set_title(f"P% = {percentage_list[l]}, EV = {EV_list_conv[k]}")

# ax3.set_xlabel("P <-> P")
# ax3.set_ylabel("P <-> H")
# ax3.set_xticks([i for i in range(12)], np.array(jll_list_conv))
# ax3.set_yticks([i for i in range(3)], np.array(jlp_list_conv))
# ax3.set_title(f"P% = {percentage_list[l]}, EV = {EV_list_conv[k]}")

# ax4.set_xlabel("P <-> P")
# ax4.set_ylabel("P <-> H")
# ax4.set_xticks([i for i in range(12)], np.array(jll_list_conv))
# ax4.set_yticks([i for i in range(3)], np.array(jlp_list_conv))
# ax4.set_title(f"P% = {percentage_list[l]}, EV = {EV_list_conv[k]}")


#         # ax_divider = make_axes_locatable(axs[i,2])
#         # Add an Axes to the right of the main Axes.
#         # cax1 = ax_divider.append_axes("right", size="7%", pad="2%")
#         # cb1 = fig.colorbar(im, cax=cax1)


# fig1.savefig(f"/home/ppuel/data/{exp}/liqFractionEquilibriumBeta_matrix.png")
# fig2.savefig(f"/home/ppuel/data/{exp}/liqDropNumEquilibriumBeta_matrix.png")
# fig3.savefig(f"/home/ppuel/data/{exp}/liqDropSizeEquilibriumBeta_matrix.png")
# fig4.savefig(f"/home/ppuel/data/{exp}/polyGyrationEquilibriumBeta_matrix.png")

plt.close("all")
plt.close('all')
