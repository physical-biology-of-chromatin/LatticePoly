import importlib
import itertools
import os
import pickle
import re
import subprocess
import sys
from typing import cast

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas
from hdf5Reader import hdf5Reader
from matplotlib import axis, colormaps, layout_engine
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from networkx import rescale_layout
from scipy.optimize import Bounds, curve_fit
from scipy.spatial.distance import squareform
from seaborn import reset_defaults
from tqdm import tqdm
from utils import PrintAggregateData, exp_mapping, exp_pathList, find_exp

plt.style.use("./resources/h5py/presentation.mplstyle")
plasma = colormaps["plasma"].resampled(256)
viridis = colormaps["viridis"].resampled(256)

expNum = 59
XnfsDir = "/Xnfs/physbiochrom/ppuel/data/"
expName = find_exp(expNum, XnfsDir)
expPath = os.path.join(XnfsDir, expName)
fig_dir = f"/home/ppuel/figure/{expName}/"


os.makedirs(fig_dir, exist_ok=True)

dict_parameters, _, meta_parameter_n, meta_parameter_n_meas = exp_mapping(expPath)

pathList = exp_pathList(dict_parameters, -1, expPath)

listDatasetName = [
    "PixelLiqHist_summed",
    "PixelLiqHist_With_Het_summed"
]

plotType = "plot"

fileName = "aggregated_process.h5"
print(fileName)
with h5py.File(os.path.join(pathList[0], fileName), "r") as hfile:
    for dataset in hfile.keys():
        print(hfile[dataset].shape, dataset)

print("\n\nN/0/process.h5")

with h5py.File(os.path.join(pathList[0], "N/0/process.h5"), "r") as hfile:
    for dataset in hfile.keys():
        print(hfile[dataset].shape, dataset)

PrintAggregateData(fileName, listDatasetName, plotType, expName, logPlot=False)


def figure_compared_WT():

    WT_list = []
    liqMean_list = []
    polyGyration_list = []

    borne_inf = 5
    borne_sup = 6.5

    for path in pathList:
        with h5py.File(os.path.join(path + fileName)) as h5file:
            tmp_dataset = cast(list[float], h5file["RatioLiqDensityInAndOut"])[-1]
            if tmp_dataset > borne_inf and tmp_dataset < borne_sup:
                WT_list.append(path)
                liqMean_list.append(h5file["liqMean"][-1])
                polyGyration_list.append(h5file["polyGyration"][-1])

    max_liqMean_list = max(liqMean_list)
    min_liqMean_list = min(liqMean_list)
    max_polyGyration_list = max(polyGyration_list)
    min_polyGyration_list = min(polyGyration_list)

    l_fig = 4
    fig = plt.figure(figsize=(8 * l_fig, 3 * l_fig), layout="constrained")
    axlist = fig.subplots(3, 6)
    re_path = re.compile(r"/[\D_]+/(\d+\.?\d*)")

    WT_param = [re.findall(re_path, path) for path in WT_list]
    for id_path, [jlp, jll_valency, jpl_valency] in enumerate(WT_param):
        axlist[0, int(jll_valency) - 1].scatter(
            [float(jlp)],
            [int(jpl_valency)],
            c=plasma(
                (liqMean_list[id_path] - min_liqMean_list)
                / (max_liqMean_list - min_liqMean_list)
            ),
        )
        axlist[1, int(jll_valency) - 1].scatter(
            [float(jlp)],
            [int(jpl_valency)],
            c=viridis(
                (polyGyration_list[id_path] - min_polyGyration_list)
                / (max_polyGyration_list - min_polyGyration_list)
            ),
        )
        for id_ax, ax in enumerate(axlist[0]):
            ax.set_xlim([0 - 0.1, 1 + 0.1])
            ax.set_ylim([1 - 0.5, 6 + 0.5])
            ax.set_title(f"Self interaction valency = {id_ax + 1}")
            ax.set_xticks(
                list(map(float, dict_parameters["JLP"])),
                dict_parameters["JLP"],
                rotation=-90,
            )
            ax.set_yticks(
                list(map(int, dict_parameters["JPL_VALENCY"])),
                dict_parameters["JPL_VALENCY"],
            )
            ax.set_xlabel("Cross interaction energy")
            ax.set_ylabel("Cross interaction valency")

        for id_ax, ax in enumerate(axlist[1]):
            ax.set_xlim([0, 1 + 0.1])
            ax.set_ylim([1 - 0.5, 6 + 0.5])
            ax.set_title(f"Self interaction valency = {id_ax + 1}")
            ax.set_xticks(
                list(map(float, dict_parameters["JLP"])),
                dict_parameters["JLP"],
                rotation=-90,
            )
            ax.set_yticks(
                list(map(int, dict_parameters["JPL_VALENCY"])),
                dict_parameters["JPL_VALENCY"],
            )
            ax.set_xlabel("Cross interaction energy")
            ax.set_ylabel("Cross interaction valency")

    norm = Normalize(vmin=min_liqMean_list, vmax=max_liqMean_list)

    fig.colorbar(
        ScalarMappable(norm=norm, cmap=plasma),
        cax=axlist[2, 1],
        orientation="horizontal",
        label="protein local density",
    )

    norm = Normalize(vmin=min_polyGyration_list, vmax=max_polyGyration_list)

    fig.colorbar(
        ScalarMappable(norm=norm, cmap=viridis),
        cax=axlist[2, 4],
        orientation="horizontal",
        label="chromatin gyration radius",
    )

    fig.savefig(os.path.join(fig_dir, "compared_WT"))


# figure_compared_WT()


def exp_diffusion(params, meta_parameter_n, meta_parameter_n_meas):

    txt_params = "_" + "_".join(map(str, params))

    meta_parameter_n = int(meta_parameter_n) - 1
    n_meas_MSD = int(meta_parameter_n_meas) - 1
    n_inter = 100000

    d_liq = 1e6
    d_poly = 1e4

    fig = plt.figure(figsize=(18, 6), layout="constrained")
    ax_list = []
    ax_list = fig.subplots(1, 3)

    print(ax_list)

    d_liq_mean = 0
    d_liq_var = 0
    d_com_mean = 0
    d_com_var = 0

    list_file = os.listdir(fig_dir)

    if f"liq_MSD_matrix_Karpen{txt_params}.npy" in list_file:
        liq_MSD_matrix_Karpen = np.load(
            os.path.join(fig_dir, f"liq_MSD_matrix_Karpen{txt_params}.npy")
        )
    else:
        liq_MSD_matrix_Karpen = np.zeros((meta_parameter_n, n_meas_MSD))

    if f"poly_MSD_matrix_Karpen{txt_params}.npy" in list_file:
        poly_MSD_matrix_Karpen = np.load(
            os.path.join(fig_dir, f"poly_MSD_matrix_Karpen{txt_params}.npy")
        )
    else:
        poly_MSD_matrix_Karpen = np.zeros((meta_parameter_n, n_meas_MSD))

    if f"com_MSD_matrix_Karpen{txt_params}.npy" in list_file:
        com_MSD_matrix_Karpen = np.load(
            os.path.join(fig_dir, f"com_MSD_matrix_Karpen{txt_params}.npy")
        )
    else:
        com_MSD_matrix_Karpen = np.zeros((meta_parameter_n, n_meas_MSD))

    if (
        f"liq_MSD_matrix_Karpen{txt_params}.npy" not in list_file
        or f"poly_MSD_matrix_Karpen{txt_params}.npy" not in list_file
        or f"com_MSD_matrix_Karpen{txt_params}.npy" not in list_file
    ):
        for n in range(meta_parameter_n):
            if n != 6:
                with h5py.File(
                    os.path.join(
                        expPath,
                        f"JLP/{params[0]}/JLL_VALENCY/{params[1]}/JPL_VALENCY/{params[2]}/N/{n}/process.h5",
                    ),
                    "r",
                ) as nFile:
                    if f"liq_MSD_matrix_Karpen{txt_params}.npy" not in list_file:
                        liq_MSD_matrix_Karpen[n] = nFile.require_dataset(
                            "liqMSD", n_meas_MSD + 1, np.float64
                        )[1 : n_meas_MSD + 1]
                    if f"poly_MSD_matrix_Karpen{txt_params}.npy" not in list_file:
                        poly_MSD_matrix_Karpen[n] = nFile.require_dataset(
                            "polyHetMSD", n_meas_MSD + 1, np.float64
                        )[1 : n_meas_MSD + 1]
                    if f"com_MSD_matrix_Karpen{txt_params}.npy" not in list_file:
                        com_MSD_matrix_Karpen[n] = nFile.require_dataset(
                            "LiqPolyCoM", n_meas_MSD + 1, np.float64
                        )[1 : n_meas_MSD + 1]

    if f"liq_MSD_matrix_Karpen{txt_params}.npy" not in list_file:
        np.save(
            os.path.join(fig_dir, f"liq_MSD_matrix_Karpen{txt_params}.npy"),
            liq_MSD_matrix_Karpen,
        )
    if f"poly_MSD_matrix_Karpen{txt_params}.npy" not in list_file:
        np.save(
            os.path.join(fig_dir, f"poly_MSD_matrix_Karpen{txt_params}.npy"),
            poly_MSD_matrix_Karpen,
        )
    if f"com_MSD_matrix_Karpen{txt_params}.npy" not in list_file:
        np.save(
            os.path.join(fig_dir, f"com_MSD_matrix_Karpen{txt_params}.npy"),
            com_MSD_matrix_Karpen,
        )

    liq_MSD_matrix_Karpen_Mean = np.mean(liq_MSD_matrix_Karpen, axis=0)
    liq_MSD_matrix_Karpen_Var = np.var(liq_MSD_matrix_Karpen, axis=0)

    poly_MSD_matrix_Karpen_Mean = np.mean(poly_MSD_matrix_Karpen, axis=0)
    poly_MSD_matrix_Karpen_Var = np.var(poly_MSD_matrix_Karpen, axis=0)

    com_MSD_matrix_Karpen_Mean = np.mean(com_MSD_matrix_Karpen, axis=0)
    com_MSD_matrix_Karpen_Var = np.var(com_MSD_matrix_Karpen, axis=0)

    ax_list[0].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        liq_MSD_matrix_Karpen_Mean * 20**2 * 2,
        "darkred",
    )
    ax_list[0].tick_params(axis="y", labelcolor="darkred")

    ax_twinx = ax_list[0].twinx()
    ax_twinx.loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        np.sqrt(liq_MSD_matrix_Karpen_Var * 20**2 * 2 / meta_parameter_n),
        "darkblue",
    )
    ax_twinx.tick_params(axis="y", labelcolor="darkblue")

    ax_list[0].set_xlabel("MCS")
    ax_list[0].set_ylabel("MSD of particles", color="darkred")
    ax_twinx.set_ylabel("Var of MSD", color="darkblue")

    ax_list[1].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        poly_MSD_matrix_Karpen_Mean * 20**2 * 2,
        "darkred",
    )
    ax_list[1].tick_params(axis="y", labelcolor="darkred")

    ax_twinx = ax_list[1].twinx()
    ax_twinx.loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        np.sqrt(poly_MSD_matrix_Karpen_Var * 20**2 * 2 / meta_parameter_n),
        "darkblue",
    )
    ax_twinx.tick_params(axis="y", labelcolor="darkblue")

    ax_list[1].set_xlabel("MCS")
    ax_list[1].set_ylabel("MSD of polymer", color="darkred")
    ax_twinx.set_ylabel("Var of MSD", color="darkblue")

    ax_list[2].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        com_MSD_matrix_Karpen_Mean * 20**2 * 2,
        "darkred",
    )
    ax_list[2].tick_params(axis="y", labelcolor="darkred")

    ax_twinx = ax_list[2].twinx()
    ax_twinx.loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        np.sqrt(com_MSD_matrix_Karpen_Var * 20**2 * 2 / meta_parameter_n),
        "darkblue",
    )
    ax_twinx.tick_params(axis="y", labelcolor="darkblue")

    ax_list[2].set_xlabel("MCS")
    ax_list[2].set_ylabel("MSD of CoM", color="darkred")
    ax_twinx.set_ylabel("Var of MSD", color="darkblue")

    ax_list[0].set_title(f"Particles density : {0.019}", color="darkgreen")
    ax_list[1].set_title(f"Particles density : {0.019}", color="darkgreen")
    ax_list[2].set_title(f"Particles density : {0.019}", color="darkgreen")

    liq_coef = np.polyfit(
        np.log10(np.arange(1, n_meas_MSD + 1) * n_inter),
        np.log10(liq_MSD_matrix_Karpen_Mean * 20**2 * 2),
        deg=1,
        w=1 / np.sqrt(liq_MSD_matrix_Karpen_Var),
        cov="unscaled",
    )
    poly_coef = np.polyfit(
        np.log10(np.arange(1, n_meas_MSD + 1) * n_inter),
        np.log10(poly_MSD_matrix_Karpen_Mean * 20**2 * 2),
        deg=1,
        w=1 / np.sqrt(poly_MSD_matrix_Karpen_Var),
        cov="unscaled",
    )
    com_coef = np.polyfit(
        np.log10(np.arange(1, n_meas_MSD + 1) * n_inter),
        np.log10(com_MSD_matrix_Karpen_Mean * 20**2 * 2),
        deg=1,
        w=1 / np.sqrt(com_MSD_matrix_Karpen_Var),
        cov="unscaled",
    )

    alpha_liq_mean = liq_coef[0][0]
    d_liq_mean = 10 ** liq_coef[0][1]
    alpha_liq_var = np.sqrt(liq_coef[1][0][0]) / meta_parameter_n
    d_liq_var = np.sqrt(liq_coef[1][1][1]) / meta_parameter_n
    alpha_poly_mean = poly_coef[0][0]
    d_poly_mean = 10 ** poly_coef[0][1]
    alpha_poly_var = np.sqrt(poly_coef[1][0][0]) / meta_parameter_n
    d_poly_var = np.sqrt(poly_coef[1][1][1]) / meta_parameter_n
    alpha_CoM_mean = com_coef[0][0]
    d_com_mean = 10 ** com_coef[0][1]
    alpha_CoM_var = np.sqrt(com_coef[1][0][0]) / meta_parameter_n
    d_com_var = np.sqrt(com_coef[1][1][1]) / meta_parameter_n

    ax_list[0].text(
        0.9,
        0.075,
        f"$\\alpha = ${alpha_liq_mean:0.2f}\n$D = ${d_liq_mean:0.2e}",
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax_list[0].transAxes,
    )

    ax_list[1].text(
        0.9,
        0.075,
        f"$\\alpha = ${alpha_poly_mean:0.2f}\n$D = ${d_poly_mean:0.2e}",
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax_list[1].transAxes,
    )

    ax_list[2].text(
        0.9,
        0.075,
        f"$\\alpha = ${alpha_CoM_mean:0.2f}\n$D = ${d_com_mean:0.2e}",
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax_list[2].transAxes,
    )

    print(f"{alpha_liq_mean = }")
    print(f"{alpha_liq_var = }")
    print(f"{d_liq_mean = }")
    print(f"{d_liq_var = }")
    print(f"{alpha_poly_mean = }")
    print(f"{alpha_poly_var = }")
    print(f"{d_poly_mean = }")
    print(f"{d_poly_var = }")
    print(f"{alpha_CoM_mean = }")
    print(f"{alpha_CoM_var = }")
    print(f"{d_com_mean = }")
    print(f"{d_com_var = }")

    ax_list[0].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        10 ** liq_coef[0][1]
        * (np.arange(1, n_meas_MSD + 1) * n_inter) ** liq_coef[0][0],
        "--",
        color="black",
    )
    ax_list[1].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        10 ** poly_coef[0][1]
        * (np.arange(1, n_meas_MSD + 1) * n_inter) ** poly_coef[0][0],
        "--",
        color="black",
    )
    ax_list[2].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        10 ** com_coef[0][1]
        * (np.arange(1, n_meas_MSD + 1) * n_inter) ** com_coef[0][0],
        "--",
        color="black",
    )

    fig.savefig(
        os.path.join(fig_dir, f"DiffusionCoefficientRef_fullMSD{txt_params}.png")
    )


# exp_diffusion(["1.00", "6", "6"], meta_parameter_n, meta_parameter_n_meas)


def ref_diffusion_49(dict_parameters, meta_parameter_n, meta_parameter_n_meas):

    meta_parameter_n = int(meta_parameter_n)
    n_meas_MSD = int(meta_parameter_n_meas) - 1
    n_inter = 100000

    d_liq = 1e6
    d_poly = 1e4

    fig = plt.figure(figsize=(18, 6), layout="constrained")
    ax_list = []
    ax_list = fig.subplots(1, 3)

    print(ax_list)

    d_liq_mean = 0
    d_liq_var = 0
    d_com_mean = 0
    d_com_var = 0

    list_file = os.listdir(fig_dir)

    if "liq_MSD_matrix_Karpen.npy" in list_file:
        liq_MSD_matrix_Karpen = np.load(
            os.path.join(fig_dir, "liq_MSD_matrix_Karpen.npy")
        )
    else:
        liq_MSD_matrix_Karpen = np.zeros((meta_parameter_n, n_meas_MSD))

    if "poly_MSD_matrix_Karpen.npy" in list_file:
        poly_MSD_matrix_Karpen = np.load(
            os.path.join(fig_dir, "poly_MSD_matrix_Karpen.npy")
        )
    else:
        poly_MSD_matrix_Karpen = np.zeros((meta_parameter_n, n_meas_MSD))

    if "com_MSD_matrix_Karpen.npy" in list_file:
        com_MSD_matrix_Karpen = np.load(
            os.path.join(fig_dir, "com_MSD_matrix_Karpen.npy")
        )
    else:
        com_MSD_matrix_Karpen = np.zeros((meta_parameter_n, n_meas_MSD))

    if (
        "liq_MSD_matrix_Karpen.npy" not in list_file
        or "poly_MSD_matrix_Karpen.npy" not in list_file
        or "com_MSD_matrix_Karpen.npy" not in list_file
    ):
        for n in range(meta_parameter_n):
            with h5py.File(
                os.path.join(
                    expPath,
                    f"LDENS/0.019/n_inter/100000/N/{n}/process.h5",
                ),
                "r",
            ) as nFile:
                if "liq_MSD_matrix_Karpen.npy" not in list_file:
                    liq_MSD_matrix_Karpen[n] = nFile.require_dataset(
                        "liqMSD", n_meas_MSD + 1, np.float64
                    )[1 : n_meas_MSD + 1]
                if "poly_MSD_matrix_Karpen.npy" not in list_file:
                    poly_MSD_matrix_Karpen[n] = nFile.require_dataset(
                        "polyHetMSD", n_meas_MSD + 1, np.float64
                    )[1 : n_meas_MSD + 1]
                if "com_MSD_matrix_Karpen.npy" not in list_file:
                    com_MSD_matrix_Karpen[n] = nFile.require_dataset(
                        "LiqPolyCoM", n_meas_MSD + 1, np.float64
                    )[1 : n_meas_MSD + 1]

    if "liq_MSD_matrix_Karpen.npy" not in list_file:
        np.save(
            os.path.join(fig_dir, "liq_MSD_matrix_Karpen.npy"), liq_MSD_matrix_Karpen
        )
    if "poly_MSD_matrix_Karpen.npy" not in list_file:
        np.save(
            os.path.join(fig_dir, "poly_MSD_matrix_Karpen.npy"), poly_MSD_matrix_Karpen
        )
    if "com_MSD_matrix_Karpen.npy" not in list_file:
        np.save(
            os.path.join(fig_dir, "com_MSD_matrix_Karpen.npy"),
            com_MSD_matrix_Karpen,
        )

    liq_MSD_matrix_Karpen_Mean = np.mean(liq_MSD_matrix_Karpen, axis=0)
    liq_MSD_matrix_Karpen_Var = np.var(liq_MSD_matrix_Karpen, axis=0)

    poly_MSD_matrix_Karpen_Mean = np.mean(poly_MSD_matrix_Karpen, axis=0)
    poly_MSD_matrix_Karpen_Var = np.var(poly_MSD_matrix_Karpen, axis=0)

    com_MSD_matrix_Karpen_Mean = np.mean(com_MSD_matrix_Karpen, axis=0)
    com_MSD_matrix_Karpen_Var = np.var(com_MSD_matrix_Karpen, axis=0)

    ax_list[0].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        liq_MSD_matrix_Karpen_Mean * 20**2 * 2,
        "darkred",
    )
    ax_list[0].tick_params(axis="y", labelcolor="darkred")

    ax_twinx = ax_list[0].twinx()
    ax_twinx.loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        np.sqrt(liq_MSD_matrix_Karpen_Var * 20**2 * 2 / meta_parameter_n),
        "darkblue",
    )
    ax_twinx.tick_params(axis="y", labelcolor="darkblue")

    ax_list[0].set_xlabel("MCS")
    ax_list[0].set_ylabel("MSD of particles", color="darkred")
    ax_twinx.set_ylabel("Var of MSD", color="darkblue")

    ax_list[1].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        poly_MSD_matrix_Karpen_Mean * 20**2 * 2,
        "darkred",
    )
    ax_list[1].tick_params(axis="y", labelcolor="darkred")

    ax_twinx = ax_list[1].twinx()
    ax_twinx.loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        np.sqrt(poly_MSD_matrix_Karpen_Var * 20**2 * 2 / meta_parameter_n),
        "darkblue",
    )
    ax_twinx.tick_params(axis="y", labelcolor="darkblue")

    ax_list[1].set_xlabel("MCS")
    ax_list[1].set_ylabel("MSD of polymer", color="darkred")
    ax_twinx.set_ylabel("Var of MSD", color="darkblue")

    ax_list[2].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        com_MSD_matrix_Karpen_Mean * 20**2 * 2,
        "darkred",
    )
    ax_list[2].tick_params(axis="y", labelcolor="darkred")

    ax_twinx = ax_list[2].twinx()
    ax_twinx.loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        np.sqrt(com_MSD_matrix_Karpen_Var * 20**2 * 2 / meta_parameter_n),
        "darkblue",
    )
    ax_twinx.tick_params(axis="y", labelcolor="darkblue")

    ax_list[2].set_xlabel("MCS")
    ax_list[2].set_ylabel("MSD of CoM", color="darkred")
    ax_twinx.set_ylabel("Var of MSD", color="darkblue")

    ax_list[0].set_title(f"Particles density : {0.019}", color="darkgreen")
    ax_list[1].set_title(f"Particles density : {0.019}", color="darkgreen")
    ax_list[2].set_title(f"Particles density : {0.019}", color="darkgreen")

    liq_coef = np.polyfit(
        np.log10(np.arange(1, n_meas_MSD + 1) * n_inter),
        np.log10(liq_MSD_matrix_Karpen_Mean * 20**2 * 2),
        deg=1,
        w=1 / np.sqrt(liq_MSD_matrix_Karpen_Var),
        cov="unscaled",
    )
    poly_coef = np.polyfit(
        np.log10(np.arange(1, n_meas_MSD + 1) * n_inter),
        np.log10(poly_MSD_matrix_Karpen_Mean * 20**2 * 2),
        deg=1,
        w=1 / np.sqrt(poly_MSD_matrix_Karpen_Var),
        cov="unscaled",
    )
    com_coef = np.polyfit(
        np.log10(np.arange(1, n_meas_MSD + 1) * n_inter),
        np.log10(com_MSD_matrix_Karpen_Mean * 20**2 * 2),
        deg=1,
        w=1 / np.sqrt(com_MSD_matrix_Karpen_Var),
        cov="unscaled",
    )

    alpha_liq_mean = liq_coef[0][0]
    d_liq_mean = 10 ** liq_coef[0][1]
    alpha_liq_var = np.sqrt(liq_coef[1][0][0]) / meta_parameter_n
    d_liq_var = np.sqrt(liq_coef[1][1][1]) / meta_parameter_n
    alpha_poly_mean = poly_coef[0][0]
    d_poly_mean = 10 ** poly_coef[0][1]
    alpha_poly_var = np.sqrt(poly_coef[1][0][0]) / meta_parameter_n
    d_poly_var = np.sqrt(poly_coef[1][1][1]) / meta_parameter_n
    alpha_CoM_mean = com_coef[0][0]
    d_com_mean = 10 ** com_coef[0][1]
    alpha_CoM_var = np.sqrt(com_coef[1][0][0]) / meta_parameter_n
    d_com_var = np.sqrt(com_coef[1][1][1]) / meta_parameter_n

    ax_list[0].text(
        0.9,
        0.075,
        f"$\\alpha = ${alpha_liq_mean:0.2f}\n$D = ${d_liq_mean:0.2e}",
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax_list[0].transAxes,
    )

    ax_list[1].text(
        0.9,
        0.075,
        f"$\\alpha = ${alpha_poly_mean:0.2f}\n$D = ${d_poly_mean:0.2e}",
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax_list[1].transAxes,
    )

    ax_list[2].text(
        0.9,
        0.075,
        f"$\\alpha = ${alpha_CoM_mean:0.2f}\n$D = ${d_com_mean:0.2e}",
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax_list[2].transAxes,
    )

    print(f"{alpha_liq_mean = }")
    print(f"{alpha_liq_var = }")
    print(f"{d_liq_mean = }")
    print(f"{d_liq_var = }")
    print(f"{alpha_poly_mean = }")
    print(f"{alpha_poly_var = }")
    print(f"{d_poly_mean = }")
    print(f"{d_poly_var = }")
    print(f"{alpha_CoM_mean = }")
    print(f"{alpha_CoM_var = }")
    print(f"{d_com_mean = }")
    print(f"{d_com_var = }")

    ax_list[0].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        10 ** liq_coef[0][1]
        * (np.arange(1, n_meas_MSD + 1) * n_inter) ** liq_coef[0][0],
        "--",
        color="black",
    )
    ax_list[1].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        10 ** poly_coef[0][1]
        * (np.arange(1, n_meas_MSD + 1) * n_inter) ** poly_coef[0][0],
        "--",
        color="black",
    )
    ax_list[2].loglog(
        np.arange(1, n_meas_MSD + 1) * n_inter,
        10 ** com_coef[0][1]
        * (np.arange(1, n_meas_MSD + 1) * n_inter) ** com_coef[0][0],
        "--",
        color="black",
    )

    fig.savefig(os.path.join(fig_dir, "DiffusionCoefficientRef_fullMSD.png"))


# ref_diffusion_49(dict_parameters, meta_parameter_n, meta_parameter_n_meas)
