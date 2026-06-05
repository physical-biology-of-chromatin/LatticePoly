##
##  Reader.py
##  LatticePoly
##
##  Created by ppuel on 18/10/2024.
##  Copyright © 2024 ENS Lyon. All rights reserved.
##

import os

import h5py

import numpy as np


class Reader: 
    """
    Reader is a context managers, iterator object on a traj file in h5 format.
    We put the emphasis on when the h5 file is open. When the reader is created,
    _init_reader open once the h5 file to gather the initial configuration and
    the metadata. Then this reader object can be "open" with a :
        "with reader as iterator"
    to open the h5 file. Then it behave as an interator that return data for each
    frame of the simulation. Iterator can be used in a for loop with the syntax:
        "for findex, fdata in iterator"
    The use of a context manager allows the h5 file to be properly closed if an
    error occurs during the loop.  
    
    Attributes
    ----------
    path_traj_dir : str
        path to the directory where the trajectory file is located.
        
    name_traj_file : str, optional
        name of the .h5 trajectory file (default is traj.h5).
        
    init_frame : int, optional
        first frame of the iterator 
        (default is -1 to start from the begining).
        
    read_liq : bool, optional
        iterate on liquid frame (default is False).
        
    read_poly : bool, optional
        iterate on polymer frame (default is True).
        
    back_in_box : bool, optional
        compute polymer positions to account for 
        periodic boundaries conditions (default is False).
    
    Methods
    -------
    _init_reader()
        check frame range and initialize the first frame.
        
    _read_liq_frame()
        fetch the liquid position, density and displacement
        of the current frame.
        
    _read_poly_frame()
        fetch the polymer position, type and painter
        of the current frame and fix PBCs if necessary.
        
    _check_range()
        check initial frame with fetched frame 
        boundaries for liquid and/or poly.
        
    _parse_file_datasets()
        fetch frame boundaries for liquid and/or poly
        from trajectory file.
        
    _fix_PBCs()
        update polymer positions to account
        for periodic boundaries conditions.
    """

    def __init__(
        self,
        path_traj_dir,
        name_traj_file='traj.h5',
        init_frame=-1,
        read_liq=False,
        read_poly=True,
        back_in_box=False,
    ):
        """
            Parameters
            ----------
            path_traj_dir : str
                path to the directory where the traj file is located.
                
            name_traj_file : str, optional
                name of the .h5 traj file (default is traj.h5).
                
            init_frame : int, optional
                first frame of the iterator 
                (default is -1 to start from the begining).
                
            read_liq : bool, optional
                iterate on liquid frame (default is False).
                
            read_poly : bool, optional
                iterate on polymer frame (default is True).
                
            back_in_box : bool, optional
                compute polymer positions to account for periodic 
                boundaries conditions (default is False).
        """

        self.path_traj_dir = os.path.realpath(path_traj_dir)

        if not os.path.exists(self.path_traj_dir):
            raise IOError(f"Directory '{self.path_traj_dir}' does not exist")
        
        self.name_traj_file = name_traj_file
        self.path_traj_file = os.path.join(
            self.path_traj_dir, 
            self.name_traj_file
        )

        self.n_frame = 0
        self.frame = self.init_frame = init_frame
        self.n_liq = self.n_tad = self.n_euc = self.n_het = 0

        self._read_liq = read_liq
        self._read_poly = read_poly
        self._back_in_box = back_in_box

        self._init_reader()
        
    
    def __enter__(self):
        self.file = h5py.File(self.path_traj_file, "r")
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()
    
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

    def _init_reader(self):
        """Check frame range and initialize the first frame.

        Raises
        ------
        IOError
            Catch IOError(s) found during _check_range.
        """
        with h5py.File(self.path_traj_file, "r") as self.file:
    
            try:
                self._check_range()
    
                self.box_dim = self.file["Pol"].attrs["L"]
    
                print(f"Box linear dimensions: {tuple(self.box_dim)}.")
    
                if self._read_liq:
                    self._read_liq_frame()
    
                    self.n_liq = self.liq_dens.size
    
                    print(f"Initial liquid state: {self.n_liq:0d} occupied sites.")
    
                if self._read_poly:
                    self._read_poly_frame()
    
                    self.n_tad = self.poly_type.size
    
                    print(f"Initial chromatin state: {self.n_tad} monomer")
                    print(f"inc. {self.n_het} heterochromatic loci.")
    
            except IOError:
                raise

    def _read_liq_frame(self):
        """Fetch the liquid position, density and displacement
        of the current frame.
        """
        
        grp_liq = self.file["Liq"]

        self.liq_pos = np.array(
            grp_liq[f"{self.frame:05d}_position_dataset"],
            dtype=np.float32,
        )

        self.liq_dens = np.array(
            grp_liq[f"{self.frame:05d}_density_dataset"],
            dtype=np.float32,
        )
        self.liq_disp = np.array(
            grp_liq[f"{self.frame:05d}_displacement_dataset"],
            dtype=np.float32,
        )

    def _read_poly_frame(self):
        """Fetch the polymer position, type and painter
        of the current frame and fix PBCs if necessary.
        """
        
        grp_poly = self.file["Pol"]

        self.poly_pos = np.array(
            grp_poly[f"{self.frame:05d}_position_dataset"],
            dtype=np.float32,
        )

        self.poly_type = np.array(
            grp_poly[f"{self.frame:05d}_type_dataset"],
            dtype=np.float32,
        )
        
        self.poly_painter = np.array(
            grp_poly[f"{self.frame:05d}_painter_dataset"],
            dtype=np.float32,
        )

        self.n_euc = np.count_nonzero(self.poly_type == 0)
        self.n_het = np.count_nonzero(self.poly_type == 1)

        het_domains = np.nonzero(self.poly_type == 1)[0]
        self.domains = np.split(het_domains, np.where(np.diff(het_domains) != 1)[0] + 1)

        self.n_dom = len(self.domains)

        if self._back_in_box:
            self._fix_PBCs()

    def _check_range(self):
        """Check initial frame with fetched frame 
        boundaries for liquid and/or polymer.

        Raises
        ------
        IOError
            If initial frame required is not in range.
        """
        
        self._parse_file_datasets()

        if self._read_poly & self._read_liq:
            min_frame = max(self._min_frame_liq, self._min_frame_poly)
            max_frame = min(self._max_frame_liq, self._max_frame_poly)

            self.frame = self.init_frame = min_frame if self.init_frame == -1 else self.init_frame

            if not min_frame <= self.init_frame <= max_frame:
                raise IOError(
                        f"Frame not in range ({min_frame}, {max_frame})"
                )

            self.n_frame = max_frame - self.init_frame + 1

        elif self._read_poly:
            self.frame = self.init_frame = (
                self._min_frame_poly if self.init_frame == -1 else self.init_frame
            )

            if not self._min_frame_poly <= self.init_frame <= self._max_frame_poly:
                raise IOError(
                    f"Frame not in range ({self._min_frame_poly}, {self._max_frame_poly})"
                )

            self.n_frame = self._max_frame_poly - self.init_frame + 1

        elif self._read_liq:
            self.frame = self.init_frame = (
                self._min_frame_liq if self.init_frame == -1 else self.init_frame
            )

            if not self._min_frame_liq <= self.init_frame <= self._max_frame_liq:
                raise IOError(
                    f"Frame not in range ({self._min_frame_liq}, {self._max_frame_liq})"
                )

            self.n_frame = self._max_frame_liq - self.init_frame + 1

    def _parse_file_datasets(self):
        """Fetch frame boundaries for liquid and/or poly.

        Raises
        ------
        IOError
            If there was an output error during trajectory computation.
        """
        
        if self._read_poly:
            try:
                grp_pol = self.file["Pol"]

                poly_seq = list(grp_pol.keys())
                poly_seq.sort()

                self._min_frame_poly = int(poly_seq[0].split("_")[0])
                self._max_frame_poly = int(poly_seq[-1].split("_")[0])

            except Exception:
                raise IOError(
                    f"Could not locate any polymer configuration dataset in '{self.path_traj_file}'/Pol."
                )
        if self._read_liq:
            try:
                grp_liq = self.file["Pol"]

                liq_seq = list(grp_liq.keys())
                liq_seq.sort()

                self._min_frame_liq = int(liq_seq[0].split("_")[0])
                self._max_frame_liq = int(liq_seq[-1].split("_")[0])

            except Exception:
                raise IOError(
                    f"Could not locate any liquid configuration dataset in '{self.path_traj_file}'/Liq."
                )

    def _fix_PBCs(self):
        """Update polymer positions to account
        for periodic boundaries conditions.
        """
        np.mod(self.poly_pos, self.box_dim, out=self.poly_pos)

    # @staticmethod
    # @numba.jit("void(i4[:], f4[:,:])", nopython=True)
    # def _fix_PBCs(
    #     dims: np.ndarray[tuple[int], np.dtype[np.int32]],
    #     pts: np.ndarray[tuple[int, int], np.dtype[np.float32]],
    # ):
    #     n_points = pts.shape[0]

    #     for i in range(n_points):
    #         for j in range(3):
    #             while pts[i, j] < 0:
    #                 pts[i, j] += dims[j]

    #             while pts[i, j] >= dims[j]:
    #                 pts[i, j] -= dims[j]
