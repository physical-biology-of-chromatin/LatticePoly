##
##  hdf5Reader.py
##  LatticePoly
##
##  Created by ppuel on 18/10/2024.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
from typing import cast

import h5py
import numba
import numpy as np


class hdf5Reader:
    def __init__(
        self,
        output_directory,
        file_name,
        init_frame=-1,
        read_liq=False,
        read_poly=True,
        back_in_box=False,
    ):
        self.output_directory = os.path.realpath(output_directory)
        self.file_name = file_name
        self.file_path = os.path.join(self.output_directory, self.file_name)

        if os.path.exists(self.output_directory):
            self._h5_file = os.path.join(self.output_directory, "test.h5")
            # print("\033[1;34mParsing directory '%s'\033[0m" % self.output_directory)

        else:
            raise IOError(f"Directory '{self.output_directory}' does not exist")

        # self.liq_pos = None
        # self.poly_pos = None

        # self.liq_disp = None
        # self.liq_dens = None

        # self.poly_type = None
        # self.poly_painter = None

        # self.box_dim = []

        self.n_frame = 0
        self.frame = self.init_frame = init_frame
        self.n_liq = self.n_tad = self.n_euc = self.n_het = 0

        self._read_liq = read_liq
        self._read_poly = read_poly
        self._back_in_box = back_in_box

        self.file = h5py.File(self.file_path, "r")

        self.init_reader(init_frame)

    def __iter__(self):
        return self

    def __next__(self):
        if self.frame < self.n_frame + self.init_frame:
            if self._read_liq:
                self._read_liq_frame()

            if self._read_poly:
                self._read_poly_frame()

            self.frame += 1

            return self

        else:
            raise StopIteration

    def __len__(self):
        return self.n_frame

    def close(self):
        self.file.close()

    def init_reader(self, init_frame):
        try:
            self._check_range(init_frame)

            self.box_dim = np.array(
                [
                    cast(
                        np.ndarray[tuple[int, ...], np.dtype[np.int32]],
                        self.file["Pol"].attrs["L"],
                    )[0]
                    for _ in range(3)
                ]
            )  # /!\

            # print("Box linear dimensions: (%.0f,%.0f,%.0f)" % tuple(self.box_dim))

            if self._read_liq:
                self._read_liq_frame()

                self.n_liq = self.liq_dens.size

                # print("Initial liquid state: %d occupied sites" % self.n_liq)

            if self._read_poly:
                self._read_poly_frame()

                self.n_tad = self.poly_type.size

                # print("Initial chromatin state: %d TADs inc. %d heterochromatic loci" % (self.n_tad, self.n_het))

        except IOError:
            raise

    def _read_liq_frame(self):
        grp_liq = cast(h5py.Group, self.file["Liq"])

        self.liq_pos = np.array(
            cast(h5py.Dataset, grp_liq["{:05d}_position_dataset".format(self.frame)]),
            dtype=np.float32,
        )

        self.liq_dens = np.array(
            cast(h5py.Dataset, grp_liq["{:05d}_density_dataset".format(self.frame)]),
            dtype=np.float32,
        )
        self.liq_disp = np.array(
            cast(
                h5py.Dataset, grp_liq["{:05d}_displacement_dataset".format(self.frame)]
            ),
            dtype=np.float32,
        )

    def _read_poly_frame(self):

        grp_poly = cast(h5py.Group, self.file["Pol"])

        self.poly_pos = np.array(
            cast(h5py.Dataset, grp_poly["{:05d}_position_dataset".format(self.frame)]),
            dtype=np.float32,
        )

        self.poly_type = np.array(
            cast(h5py.Dataset, grp_poly["{:05d}_type_dataset".format(self.frame)]),
            dtype=np.float32,
        )
        self.poly_painter = np.array(
            cast(h5py.Dataset, grp_poly["{:05d}_painter_dataset".format(self.frame)]),
            dtype=np.float32,
        )

        self.n_euc = np.count_nonzero(self.poly_type == 0)
        self.n_het = np.count_nonzero(self.poly_type == 1)

        het_domains = np.nonzero(self.poly_type == 1)[0]
        self.domains = np.split(het_domains, np.where(np.diff(het_domains) != 1)[0] + 1)

        self.n_dom = len(self.domains)

        if self._back_in_box:
            self._fix_PBCs(self.box_dim, self.poly_pos)

    def _check_range(self, init_frame):
        self._parse_file_seqs()

        if self._read_poly & self._read_liq:
            min_frame = max(self._min_frame_liq, self._min_frame_poly)
            max_frame = min(self._max_frame_liq, self._max_frame_poly)

            self.frame = self.init_frame = min_frame if init_frame == -1 else init_frame

            if not min_frame <= self.init_frame <= max_frame:
                raise IOError(
                        f"Frame not in range ({min_frame}, {max_frame})"
                )

            self.n_frame = max_frame - self.init_frame + 1

        elif self._read_poly:
            self.frame = self.init_frame = (
                self._min_frame_poly if init_frame == -1 else init_frame
            )

            if not self._min_frame_poly <= self.init_frame <= self._max_frame_poly:
                raise IOError(
                    f"Frame not in range ({self._min_frame_poly}, {self._max_frame_poly})"
                )

            self.n_frame = self._max_frame_poly - self.init_frame + 1

        elif self._read_liq:
            self.frame = self.init_frame = (
                self._min_frame_liq if init_frame == -1 else init_frame
            )

            if not self._min_frame_liq <= self.init_frame <= self._max_frame_liq:
                raise IOError(
                    f"Frame not in range ({self._min_frame_liq}, {self._max_frame_liq})"
                )

            self.n_frame = self._max_frame_liq - self.init_frame + 1

    def _parse_file_seqs(self):
        if self._read_poly:
            try:
                grp_pol = cast(h5py.Group, self.file["Pol"])

                poly_seq = list(grp_pol.keys())
                poly_seq.sort()

                self._min_frame_poly = int(poly_seq[0].split("_")[0])
                self._max_frame_poly = int(poly_seq[-1].split("_")[0])

            except Exception:
                raise IOError(
                    f"Could not locate any polymer configuration dataset in '{self.file_path}'"
                )
        if self._read_liq:
            try:
                grp_liq = cast(h5py.Group, self.file["Pol"])

                liq_seq = list(grp_liq.keys())
                liq_seq.sort()

                self._min_frame_liq = int(liq_seq[0].split("_")[0])
                self._max_frame_liq = int(liq_seq[-1].split("_")[0])

            except Exception:
                raise IOError(
                    f"Could not locate any liquid configuration dataset in '{self.file_path}'"
                )

    @staticmethod
    @numba.jit("void(i4[:], f4[:,:])", nopython=True)
    def _fix_PBCs(
        dims: np.ndarray[tuple[int], np.dtype[np.int32]],
        pts: np.ndarray[tuple[int, int], np.dtype[np.float32]],
    ):
        n_points = pts.shape[0]

        for i in range(n_points):
            for j in range(3):
                while pts[i, j] < 0:
                    pts[i, j] += dims[j]

                while pts[i, j] >= dims[j]:
                    pts[i, j] -= dims[j]
