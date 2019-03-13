"""
Tauren-MD trajectory objects.
"""
# Copyright © 2018-2019 Tauren-MD Project
#
# Tauren-MD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Tauren-MD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Tauren-MD. If not, see <http://www.gnu.org/licenses/>.
#
# Contributors to this file:
# - João M.C. Teixeira (https://github.com/joaomcteixeira)
import sys
import string
from collections import namedtuple
import numpy as np

from abc import ABC, abstractmethod

import mdtraj
import MDAnalysis as mda
from MDAnalysis.analysis import align as mdaalign
from MDAnalysis.analysis.rms import RMSD as mdaRMSD

from tauren import logger

log = logger.get_log(__name__)

StorageKey = namedtuple(
    "StorageKey",
    [
        "datatype",
        "identifier",
        "filenaming",
        ],
    )

StorageData = namedtuple(
    "StorageData",
    [
        "columns",
        "data",
        ],
    )


class TaurenTraj(ABC):
    """
    Base class for Tauren-MD sub Traj classes.
    """
    
    _err_frame_index = \
        "    frame '{}' does NOT exist in trajectory, ignoring..."
    
    def __init__(self):
        
        self._set_full_frames_list()
        self._update_traj_slicer(
            start=0,
            end=len(self.full_frames_list),
            step=1
            )
        
        self.observables = None
        self._rmsds_counter = 0
        
        return
    
    # @property
    # @abstractmethod
    # def trajectory(self):
        # """
        # Returns the trajectory object.
        # """
        # pass
    
    @property
    @abstractmethod
    def trajectory(self):
        pass
    
    @trajectory.setter
    def trajectory(self):
        log.info("* CANT manually set a sliced traj. Ignoring...")
        return None
    
    @property
    def observables(self):
        """A dictionary containing all data obtained from the system."""
        return self._observables
    
    @observables.setter
    def observables(self, obs=None):
        self._observables = obs or TrajObservables()
    
    @abstractmethod
    def _set_full_frames_list(self, num_frames):
        """
        Defines the list of frames in the input trajectory.
        """
        self._full_frames_list = list(range(1, num_frames + 1))
    
    @property
    def full_frames_list(self):
        """
        The list of frames in the input trajectory.
        """
        return self._full_frames_list
    
    @property
    def sliced_frames_list(self):
        """
        The list of frames in the current frame slicing.
        """
        return self.full_frames_list[self._fslicer]
    
    @property
    def slice_tuple(self):
        """
        A 3 element tuple identifing the current slicing
        (start, end, step).
        """
        return self._slice_tuple
    
    @slice_tuple.setter
    def slice_tuple(self, tuple_):
        """
        Tuple (start, end, step)
        """
        start = tuple_[0]
        end = tuple_[1]
        step = tuple_[2]
        
        if step == 0:
            raise ValueError("step can NOT be zero.")
        
        elif step > 0 and not(start < end) \
                or step < 0 and not(start > end):
            raise ValueError(
                f"This tuple combination '{tuple_}'"
                " will render an empty selection."
                " Use start > end for step > 0 or"
                " start < end for steps < 0."
                )
        
        log.info(f"slice_tuple: {tuple_}")
        
        self._slice_tuple = tuple_
        
        return
    
    @property
    def n_frames(self):
        """
        The number of frames in the trajectory considering the
        current slicing as defined by the :attr:`slice_tuple`.
        """
        return len(self.full_frames_list[self._fslicer])
    
    @property
    @abstractmethod
    def totaltime(self):
        """
        The total time in the current slicing.
        """
        pass
    
    @property
    @abstractmethod
    def timestep(self):
        """
        The time step in the current slicing.
        """
        pass
    
    @property
    def atom_selection(self):
        try:
            return self._atom_selection
        
        except AttributeError:
            return "all"
    
    @atom_selection.setter
    def atom_selection(self, selector):
        """
        String that defines the current atom selection.
        
        Parameters
        ----------
        selector : str
            The string the defines the selector.
        """
        self._atom_selection = selector
    
    @property
    @abstractmethod
    def n_residues(self):
        """
        The number of residues in the current slicing.
        """
        pass
    
    @property
    @abstractmethod
    def n_atoms(self):
        """
        The number of atoms in the current slicing.
        """
        pass
    
    def _update_traj_slicer(
            self,
            start,
            end,
            step,
            ):
        """
        Updates the current frame slicing object.
        """
        
        self._check_correct_slice(end)
        
        self._fslicer = slice(start, end, step)
        
        self.slice_tuple = (start, end, step)
        
        log.debug(f"<_fslicer> updated: {self._fslicer}")
        
        return
    
    def _check_correct_slice(self, end):
        """
        Checks if slicing end is less or equal than traj length.
        """
        
        if not isinstance(end, int):
            raise TypeError(f"<end> must be int, {type(end)} given.")
        
        _err = (
            "* WARNING! *"
            "\n"
            f"Your slicing range '{end}' "
            "goes beyond the trajectory length. "
            "The maximum length of your trajectory "
            f"is {len(self.full_frames_list)} frames."
            )
        
        if end > len(self.full_frames_list):
            log.warning(_err)
            raise ValueError(_err)
        
        return
    
    def report(self):
        """Prints general information about trajectory."""
        
        info = (
            "* Trajectory details:\n"
            "\n"
            f"    num of frames: {self.n_frames}\n"
            f"    num of residues: {self.n_residues}\n"
            f"    num of atoms: {self.n_atoms}\n"
            "\n"
            f"    total time: {self.totaltime}\n"
            "\n"
            f"    time_step: {self.timestep} (usually ps)\n"
            f"        or {self.timestep / 1000} (usually ns)\n"
            "*\n"
            )
        
        log.info(info)
        
        return
    
    def remove_solvent(self, **kwargs):
        """
        Removes solvent from trajectory.
        """
        
        self._remove_solvent(**kwargs)
    
    @abstractmethod
    def _remove_solvent(self):
        pass
    
    def align_traj(
            self,
            *,
            weights="mass",
            file_name="aligned_traj.dcd",
            inplace=True,
            ):
        """
        Aligns trajectory to the topology structure (reference).
        """
        
        log.info("* Aligning trajectory... ")
        
        if not(isinstance(weights, str)):
            raise TypeError("weights is NOT str type.")
        
        if not(isinstance(file_name, str)):
            raise TypeError("file_name is NOT str type.")
        
        if not(isinstance(inplace, bool)):
            raise TypeError("inplace is NOT bool.")
        
        self._align_traj(
            weights,
            file_name,
            inplace,
            )
        
        log.info("    done")
        
        return
    
    @abstractmethod
    def _align_traj(self):
        return
    
    @abstractmethod
    def image_molecules(self):
        """
        Images molecules.
        """
        # Do whatever needs to be done to start using the trajectory
        # without solvent. The implementation may differ drastically
        # from library used.
        pass
    
    def frame_slice(
            self,
            start=1,
            end=None,
            step=1,
            ):
        """
        Slices trajectory in frames for the subsequent operations.
        
        Parameters
        ----------
        start : int
            The starting frame.
            Frame index starts at 1.
            Defaults to None, slices from the start of trajectory.
        
        end : int
            The end frame for the new slicing (INCLUSIVE).
            END should be lower or equal than the traj length.
            Defaults to None, slices until the end of the trajectory.
        
        step : int
            Integer value which determines the increment between
            each index for slicing.
            Defaults to 1.
        
        Exceptions
        ----------
        TypeError
            If any parameter is not integer type.
        
        ValueError
            Terminates if any parameter equals 0.
        """
        
        # None values can be received from the config.json file
        start = start or 1
        end = end or len(self.full_frames_list)
        step = step or 1
        
        log.info(f"* Slicing trajectory [{start}, {end}, {step}]...")
        
        _err = "parameter should be integer"
        
        if not isinstance(start, int):
            raise TypeError(f"start {_err}: '{type(start)}'")
        
        if not isinstance(end, int):
            raise TypeError(f"end {_err}: '{type(end)}'")
        
        if not isinstance(step, int):
            raise TypeError(f"step {_err}: '{type(step)}'")
        
        if start == 0 or end == 0 or step == 0:
            _err = (
                "* ERROR * "
                "Traj frames indexes start from 1. "
                "You can NOT slice from 0. "
                "Remember 'end' is INCLUSIVE. "
                "Review your configuration file and rerun."
                )
            log.info(_err)
            sys.exit("* EXIT *")
        
        self._update_traj_slicer(start - 1, end, step)
        
        log.info("    done.")
        
        return
    
    def set_atom_selection(self, selector, **kwargs):
        """
        Sets the current atom selection.
        
        Atom selection will be used in subsequent operations.
        
        Parameters
        ----------
        selector : str
            The selection string. This may deppend on the
            Trajectory library type chosen for the calculation.
        
        Raises
        ------
        TypeError
            If selection is not string.
        """
        
        log.debug(f"<selection>: {selector}")
        
        if isinstance(selector, str):
            self.atom_selection = selector
        
        elif selector is None:
            self.atom_selection = "all"
        
        else:
            raise TypeError(
                "<selection> parameter must be STRING or None types."
                f"'{type(selector)}' given."
                )
        
        return
    
    def frames2file(
            self,
            frames="all",
            prefix="_",
            ext="pdb",
            ):
        """
        Extracts trajectory frames to PDB files using prefix name.
        
        Parameters
        ----------
        frames : str, optional
            Frame or range of frames to extract.
            Defaults to "all": extract all frames from current slicing.
            
            A frame range can be defined as follows:
            End values are INCLUSIVE.
                - "1"           -> the first frame
                - "1,5,100"     -> frames 1, 5 and 100
                - "10:"         -> from frame 10 to the end
                - ":500"        -> from first frame to frame 500 (inclusive)
                - "10:50"       -> from frame 10 to 50 (inclusive)
                - "10:100:2"    -> from frame 10 to 100 in steps of 2
        
        prefix : str, optional
            The prefix name for extracted PDBs.
            Defaults to "_".
        
        ext : str, optional ["pdb"]
            The file extention to which the frames are saved.
            Deppending on the trajectory type used the allowed
            file types may differ. Reffer to the documentation of
            the MD analysis library you are using.
        
        Exceptions
        ----------
        TypeError
            If frames is not of string type.
            
        ValueError
            If frames string is not consistent with parameter
            description.
        """
        
        log.info("* Extracting frames...")
        
        log.debug(f"<frames>: {frames}")
        log.debug(f"<prefix>: {prefix}")
        log.debug(f"<ext>: {ext}")
        
        # frames_to_extract is a list of the frames number
        # starting at 1.
        if frames == "all":
            frames_to_extract = self.sliced_frames_list
        
        elif isinstance(frames, str):
            
            if frames.isdigit():
                frames_to_extract = [int(frames)]
            
            elif "," in frames and frames.replace(",", "").isdigit():
                frames_to_extract = frames.split(",")
                
            elif ":" in frames and frames.replace(":", "").isdigit():
                frames_to_extract = \
                    self.full_frames_list[
                        self._gen_frame_slicer_from_string(frames)
                        ]
            
            else:
                raise ValueError(
                    "<frames> not of valid format see: "
                    f"{self.frames2file.__doc__}"
                    )
        
        else:
            raise TypeError(
                f"<frames> should be string type: '{type(frames)}' given."
                )
                
        pdb_name_fmt = \
            prefix \
            + self._gen_pdb_name_format(len(self.full_frames_list), ext)
        
        assert isinstance(frames_to_extract, list), (
            "<frames_to_extract> should be LIST! "
            f"{type(frames_to_extract)} given"
            )
        
        assert isinstance(pdb_name_fmt, str), (
            "<pdb_name_format> should be string type. "
            f"{type(pdb_name_fmt)} given"
            )
        
        self._frames2file(frames_to_extract, pdb_name_fmt)
        
        log.info("    frames extracted successfully!")
        
        return
    
    def _gen_frame_slicer_from_string(self, s):
        """
        Returns a slicer object.
        
        Does not check for s integrity.
        """
        
        if s.isdigit():
            start = int(s) - 1
            end = int(s)
            step = 1
        
        elif s.endswith(":") and s.count(":") == 1:
            ss = s.split(":")
            start = int(ss[0]) - 1
            end = len(self.full_frames_list)
            step = 1
        
        elif s.startswith(":") and s.count(":") == 1:
            ss = s.split(":")
            start = 0
            end = int(ss[0])
            step = 1
        
        elif s.count(":") == 1:
            ss = s.split(":")
            start = int(ss[0]) - 1
            end = int(ss[1])
            step = 1
        
        elif s.count(":") == 2:
            ss = s.split(":")
            start = int(ss[0] - 1)
            end = int(ss[1])
            step = int(ss[2])
        
        else:
            raise ValueError("slice string not valid")
        
        slicer = slice(start, end, step)
        
        assert isinstance(slicer, slice), "NOT A SLICE OBJECT"
        return slicer
    
    @staticmethod
    def _gen_pdb_name_format(num_of_frames, ext):
        """
        Creates a formatting string to apply leading zeros.
        
        Parameters
        ----------
        num_of_frames : int or str
            The number of frames.
        
        ext : str
            The file extension.
        """
        
        leading_zeros = str(len(str(num_of_frames)))
        pdb_name_fmt = "{:0>" + leading_zeros + "}." + f"{ext}"
        
        return pdb_name_fmt
    
    @abstractmethod
    def _frames2file(self, frames_list, pdb_name_fmt):
        """
        frames_list is a list of integers with the frames to extract.
        frames_list should be indexed at 1 (human way not python way)
        
        pdb_name_fmt is a .format() prepared string where the number
        of the extracted frame will fit in.
        """
        return
    
    def save_traj(
            self,
            file_name="traj_output.dcd",
            **kwargs
            ):
        """
        Saves trajectory to file.
        
        Overwrites existing files.
    
        Parameters
        ----------
        file_name : str
            Name of the output trajectory file.
            File extention is taken from file_name.
        """
        log.info(f"* Exporting trajectory to: {file_name}")
    
        self._save_traj(file_name)
    
        log.info("    ... saved")
        
        return
    
    @abstractmethod
    def _save_traj(self, file_name):
        """The subclass algorithm to save a trajectory."""
        pass
    
    def calc_rmsds_combined_chains(
            self,
            *,
            chains="all",
            ref_frame=0,
            storage_key="rmsds_combined_chains",
            **kwargs
            ):
        """
        Calculates combined RMSDs for a set of chains.
        
        Calculated RMSDs are stored in the form of np.ndarray
        in trajectory's observables attribute.
        Storage key is the tuple of strings: (storage_key, chains).
        
        Parameters
        ----------
        chains : str or list of identifiers, optional
            Defaults to "all", all chains are used.
            BUT, previous selection is considered, therefore, chains
            will subselect over a previous
            :meth:`~TaurenTraj.set_atom_selection`.
            
            With str, use: "1" or comma separated identifers, "1,2,4".
            
            With list, use a list of identifiers: [1,2,4] or ["A", "D"].
            
            Remember that:
            # when using **MDAnalysis**, identifiers are the segid
            characters.
            # when using **MDTraj**, identifiers are digits that
            represent chain order.
        
        ref_frame : int, optional
            The reference frame for the RMSD calculation.
        
        storage_key : str, optional
            The first element of the key tuple with which the
            calculated RMSD data will be stored in the trajectory's
            observables' dictionary. Defaults to "rmsds_combined_chains".
        
        Returns
        -------
        key : tuple
        The key with which data was stored in observables
        attribute dictionary.
        
        Exceptions
        ----------
        TypeError
            If chains is not list or string.
        
        ValueError
            If chains str or list is not of valid format (read above).
        
        IndexError
            If trajectory can not be sliced according to chains.
        """
        
        log.info("* Calculating RMSDs for combined chains...")
        
        storage_key = storage_key or "rmsds_combined_chains"
        
        self._check_chains_argument(chains)
                
        chain_list = self._gen_chain_list(chains)  # abstractmethod
        
        combined_rmsds = self._calc_rmsds_combined_chains(
            chain_list,
            ref_frame,
            )
        
        assert combined_rmsds.ndim == 1, (
            "<combined_rmsds> array should have only one dimension. "
            f"Detected array with {combined_rmsds.ndim}."
            )
        
        frames_array = np.array(self.sliced_frames_list)
        
        assert combined_rmsds.size == frames_array.size, (
            "combined_rmsds and frames_array size does not match. "
            f"{combined_rmsds.size} vs. {frames_array.size}"
            )
        
        data = np.vstack((frames_array, combined_rmsds)).T
        
        chain_name_export = chains.replace(',', '-')
        key = StorageKey(
            datatype=storage_key,
            identifier=(
                f"{self.atom_selection} "
                f"for '{chain_name_export}' chains"
                ),
            filenaming=(
                f"{storage_key}"
                f"_{self.atom_selection.replace(' ', '-')}"
                f"_{chain_name_export.replace(' ', '-')}"
                ),
            )
        
        datatuple = StorageData(
            columns=["frames", key.identifier],
            data=data,
            )
        
        self.observables.store(key, datatuple)
        
        return key
    
    @staticmethod
    def _check_chains_argument(chains):
        """
        Checks validity of <chains> input argument.
        chains should be string of comma separated chars or digits
        or list of chars or digits. Raise TypeError otherwise.
        """
        
        def valid_chain_id(c):
            
            return (isinstance(c, int)
                    or isinstance(c, str) and (c.isdigit() or c.isalpha()))
        
        if isinstance(chains, str):
            commafree = chains.replace(",", "")
            if commafree.isdigit() or commafree.isalpha():
                return
        
            else:
                _err = (
                    "chains identifiers should be letters or digits, "
                    "separated by comma ','. "
                    f"Wrong input: {chains}"
                    )
                log.debug(_err)
                raise ValueError(_err)
        
        elif isinstance(chains, list):
            if all(valid_chain_id(c) for c in chains):
                return
            
            else:
                _err = (
                    "Chains identifiers in list should be letters or digits. "
                    f"Wrong input found: {chains}"
                    )
                log.debug(_err)
                raise TypeError(_err)
            
        else:
            _err = (
                "chains arguments should be of type str or list. "
                f"'{type(chains)}' given."
                )
            log.debug(_err)
            raise TypeError(_err)
        
    @abstractmethod
    def _calc_rmsds_combined_chains(self):
        """
        MD analysis library specific implementation.
        
        Parameters
        ----------
        chain_list : list of strs
            List of strings containing the chains to operate with
        
        ref_frame : int
            The reference frame to which calculate te rmsds
        
        Returns
        -------
        numpy.array of shape=(X,)
            The returned array should be sliced according to the
            current frame slicer (self._fslicer).
        """
        
        pass
        
    @abstractmethod
    def _gen_chain_list(self, chains):
        """
        Genereates a list of chains based on a <chains> string value.
        """
        pass
    
    def calc_rmsds_separated_chains(
            self,
            *,
            chains="all",
            ref_frame=0,
            storage_key="rmsds_separated_chains",
            **kwargs
            ):
        """
        Calculates RMSDs for each chain separately.
        
        Calculated RMSDs are stored in the form of np.ndarray
        in trajectory's observables attribute.
        Storage key is the tuple of strings:
        (storage_key, chains_headers), where chains_headers is a "-"
        separated list of the chainids parsed (A-B-C).
        
        Parameters
        ----------
        chains : str or list of idetifiers, optional
            Defaults to "all", all chains are used.
            
            With str, use: "1" or comma separated identifers, "1,2,4".
            
            With list, use a list of identifiers: [1,2,4] or ["A", "D"].
            
            Remember that:
            # when using **MDAnalysis**, identifiers are the segid
            characters.
            # when using **MDTraj**, identifiers are digits that
            represent chain order.
        
        ref_frame : int, optional
            The reference frame for the RMSD calculation.
        
        storage_key : str, optional
            The first element of the key tuple with which the
            calculated RMSD data will be stored in the trajectory's
            observables' dictionary. Defaults to "rmsds_combined_chains".
        
        Returns
        -------
        key : tuple
            The key with which data was stored in observables
            attribute dictionary.
        
        Exceptions
        ----------
        TypeError
            If chains is not list or string.
        
        ValueError
            If chains str or list is not of valid format (read above).
        
        IndexError
            If trajectory can not be sliced according to chains.
        """
        
        log.info("* Calculating RMSDs for each chain separately")
        
        self._check_chains_argument(chains)
        
        chain_list = self._gen_chain_list(chains)  # abstractmethod
        
        rmsds, chains_headers = self._calc_rmsds_separated_chains(
            chain_list,
            ref_frame=ref_frame,
            )
        
        frames_array = np.array(self.sliced_frames_list)
        
        assert rmsds.shape[0] == frames_array.size, (
            "RMSDs array does not match frames_array size. "
            f"{rmsds.shape[0]} vs. {frames_array.size}."
            )
        
        data = np.concatenate(
            (
                frames_array.reshape(frames_array.size, 1),
                rmsds
                ),
            axis=1,
            )
        
        chains_columns = list(map(
            lambda x: f"{self.atom_selection}_{x}",
            chains_headers
            ))
        
        key = StorageKey(
            datatype=storage_key,
            identifier=",".join(chains_columns),
            filenaming=(
                f"{self.atom_selection.replace(' ','-')}"
                f"_{'-'.join(chains_headers)}"
                ),
            )
        
        datatuple = StorageData(
            columns=["frames", *chains_columns],
            data=data,
            )
        
        self.observables.store(key, datatuple)
        return key
    
    @abstractmethod
    def _calc_rmsds_separated_chains(self):
        """
        MD analysis library specific implementation.
        
        Parameters
        ----------
        chain_list : list of strs
            List of strings containing the chains to operate with
        
        ref_frame : int
            The reference frame to which calculate te rmsds
        
        Returns
        -------
        numpy.array of shape=(Y,X)
            Where Y is the number of frames and X number of chains.
            The returned array should be sliced according to the
            current frame slicer (self._fslicer).
        
        list of strs
            List of chains IDs evaluated.
        """
        pass
    
    def _gen_selector(
            self,
            identifiers,
            selection="segid",
            boolean="or",
            ):
        """
        Generates a selector for chains, atoms, residues parsing.
        """
        
        log.debug(f"identifier: {identifiers}")
        
        lid = list(map(lambda x: f"{selection} {x}", identifiers))
        
        selector = f" {boolean} ".join(lid)
        
        log.debug(f"selector: {selector}")
        
        assert isinstance(selector, str), "selector should be string"
        return selector
    
    def export_data(
            self,
            key,
            file_name="table.csv",
            sep=",",
            header="",
            ):
        """
        Exports data arrays to file.
        
        Parameters
        ----------
        key : tuple
            The key with which the data is stored in the
            observables attribute dictionary.
        
        file_name : str
            The name of the file.
            Defaults to "table.csv".
            
        sep : str
            The column separator.
            Defaults to comma ",".
        
        header : str
            Any text you wish to add as comment as file header.
            Headers are identified by "#".
            Defaults to nothing.
        """
        
        log.info(f"* Exporting {file_name} data")
        
        header = f"{key}\n{header}\n{','.join(self.observables[key].columns)}"
        
        np.savetxt(
            file_name,
            self.observables[key].data,
            delimiter=sep,
            header=header,
            )
        
        log.info(f"    saved {file_name}")
        
        return


class TaurenMDAnalysis(TaurenTraj):
    
    def __init__(self, trajectory, topology):
        
        self.universe = mda.Universe(topology, trajectory)
        self.topology = mda.Universe(topology)
        self.original_traj = self.universe.trajectory
        
        super().__init__()
        
        return
    
    def _set_full_frames_list(self):
        super()._set_full_frames_list(self.original_traj.n_frames)
    
    @property
    def original_traj(self):
        return self._original_traj
    
    @original_traj.setter
    def original_traj(self, traj):
        self._original_traj = traj
    
    @TaurenTraj.trajectory.getter
    def trajectory(self):
        return self.universe.select_atoms(self.atom_selection)
    
    @TaurenTraj.totaltime.getter
    def totaltime(self):
        return self.original_traj[self._fslicer][-1].time
    
    @TaurenTraj.timestep.getter
    def timestep(self):
        return (
            self.original_traj[self._fslicer][1].time
            - self.original_traj[self._fslicer][0].time
            )
    
    @TaurenTraj.n_residues.getter
    def n_residues(self):
        return self.universe.atoms.n_residues
    
    @TaurenTraj.n_atoms.getter
    def n_atoms(self):
        return len(self.universe.atoms)
    
    def _remove_solvent(self, exclude=None, inplace=True, **kwargs):
        """
        Removes solvent
        
        NOT IMPLEMENTED
        """
        
        log.info("remove solvent is not implemented for MDAnalysis routines")
        
        return
    
    def image_molecules(self):
        log.info("image_molecules method not implemented for MDAnalaysis")
        return
    
    def _align_traj(
            self,
            weights,
            file_name,
            inplace,
            ):
        
        # https://www.mdanalysis.org/docs/documentation_pages/analysis/align.html#MDAnalysis.analysis.align.AlignTraj
        alignment = mdaalign.AlignTraj(
            self.universe,
            self.topology,
            filename=file_name,
            in_memory=inplace,
            verbose=True,
            start=self.slice_tuple[0],
            end=self.slice_tuple[1],
            step=self.slice_tuple[2],
            )
        
        alignment.run()
        
        return
    
    def _frames2file(
            self,
            frames_to_extract,
            pdb_name_fmt,
            ):
        """
        frames_to_extract, list of frames index (int)
        pdb_name_fmt, "prefix_{FORMATTING CONDITION}.extension"
        """
        
        # frames are treated separatelly to allow exception capture and log
        for frame in frames_to_extract:
            
            try:
                self.trajectory.write(
                    filename=pdb_name_fmt.format(frame),
                    frames=[frame],
                    file_format="PDB",
                    bonds=None,
                    )
            
            except IndexError as e:
                log.info(self._err_frame_index.format(frame))
                log.debug(e)
                continue
            
            log.info(f"    extracted {pdb_name_fmt.format(frame)}")
        
        return
    
    def _save_traj(
            self,
            file_name,
            ):
        
        # https://www.mdanalysis.org/MDAnalysisTutorial/writing.html#trajectories
        selection = self.universe.select_atoms(self.atom_selection)
        with mda.Writer(file_name, selection.n_atoms) as W:
            for ts in self.original_traj[self._fslicer]:
                W.write(selection)
                log.info(f"    exported {ts}")
        
        return
    
    def _calc_rmsds_combined_chains(
            self,
            chain_list,
            ref_frame,
            ):
        
        absolute_selector = self._gen_selector(chain_list)
        
        chain_selectors = absolute_selector.split(" or ")
        
        filtered_selectors = self._filter_existent_selectors(chain_selectors)
        
        log.debug(f"chain_selector: {filtered_selectors}")
        
        if len(filtered_selectors) == 0:
            _err = (
                "* ERROR *"
                " The chain list does not match any selection:"
                f" {chain_list}."
                )
            log.info(_err)
            sys.exit(1)
        
        final_selection = (
            f"{self.atom_selection}"
            f" and ({' or '.join(filtered_selectors)})"
            )
        
        log.debug(f"<final_selection>: {final_selection}")
        
        # checks for empty selection
        if not(self._filter_existent_selectors([final_selection])):
            log.info(
                "   * EMPTY SELECTION ERROR *"
                f" The atom selection provided '{final_selection}'"
                " gives an empty selection.\n"
                "* Aborting calculation..."
                )
            sys.exit(1)
        
        # https://www.mdanalysis.org/docs/documentation_pages/analysis/align.html#rms-fitting-tutorial
        # https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html#MDAnalysis.analysis.rms.RMSD
        R = mdaRMSD(
            self.universe,
            self.topology,
            select=final_selection,
            groupselection=None,
            ref_frame=ref_frame,
            )
        
        R.run()
        
        return R.rmsd[:, 2][self._fslicer]  # numpy array
    
    def _calc_rmsds_separated_chains(
            self,
            chain_list,
            ref_frame,
            ):
        
        absolute_selector = self._gen_selector(chain_list)
        
        chain_selectors = absolute_selector.split(" or ")
        
        filtered_selectors = self._filter_existent_selectors(chain_selectors)
        
        rmsds = np.empty((self.n_frames, len(filtered_selectors)))
        
        subplot_has_data = []
        
        for ii, chain_selector in enumerate(filtered_selectors):
            
            final_selection = (
                f"{self.atom_selection}"
                f" and ({chain_selector})"
                )
            
            atoms = self.universe.select_atoms(final_selection)
            
            if len(atoms) == 0:
                log.debug("len of atoms is 0. Continuing...")
                subplot_has_data.append(False)
                continue
            
            atoms_top = self.topology.select_atoms(final_selection)
            
            R = mdaRMSD(
                atoms,
                atoms_top,
                groupselection=None,
                ref_frame=ref_frame,
                verbose=False,
                )
            
            R.run(verbose=False)
            
            rmsds[:, ii] = R.rmsd[:, 2][self._fslicer]
            
            subplot_has_data.append(True)
        
        column_headers = list(map(
            lambda x: x.replace("segid ", ""),
            filtered_selectors
            ))
        
        assert isinstance(column_headers, list), "c_selectors NOT list!"
        
        return (
            rmsds[:, subplot_has_data],
            np.array(column_headers)[subplot_has_data],
            )
    
    def _gen_chain_list(
            self,
            chains,
            ):
        
        log.debug(f"input chains: {chains}")
        
        if chains == "all":
            chain_list = list(string.ascii_letters + string.digits)
        
        elif isinstance(chains, str) and chains.count(",") == 0:
            chain_list = [chains]
        
        elif isinstance(chains, str) and chains.count(",") > 0:
            chain_list = chains.split(",")
        
        else:
            raise ValueError("could NOT read chains.")
        
        log.debug(f"return chain_list: {chain_list}")
        assert isinstance(chain_list, list), "Not a list!"
        
        return chain_list
    
    def _filter_existent_selectors(self, selectors_list):
        
        # https://www.mdanalysis.org/docs/documentation_pages/selections.html#simple-selections
        
        selectors = []
        for selector in selectors_list:
        
            atoms = self.topology.select_atoms(selector)
            if len(atoms) == 0:
                log.debug(f"chain {selector} does not exist")
                continue
            
            else:
                selectors.append(selector)
        
        return selectors
    

class TaurenMDTraj(TaurenTraj):
    
    def __init__(self, trajectory, topology):
        
        traj_ = mdtraj.load(trajectory, top=topology)
        self.original_traj = traj_
        self.topology = traj_.topology
        
        super().__init__()
        
        return
    
    def _set_full_frames_list(self):
        super()._set_full_frames_list(self.original_traj.n_frames)
    
    @property
    def original_traj(self):
        return self._trajectory
    
    @original_traj.setter
    def original_traj(self, traj):
        self._trajectory = traj
    
    @TaurenTraj.trajectory.getter
    def trajectory(self):
        
        slicer = self.original_traj.topology.select(self.atom_selection)
        
        try:
            sliced_traj = self.original_traj.atom_slice(slicer, inplace=False)
        
        except IndexError:
            log.exception("Could not slice traj")
            sys.exit(1)
        
        log.debug(
            f"returning sliced traj for atoms '{self.atom_selection}'"
            f" in frames '{self._fslicer}'"
            )

        return sliced_traj[self._fslicer]
    
    @TaurenTraj.totaltime.getter
    def totaltime(self):
        return self.trajectory[self._fslicer].time[-1]
    
    @TaurenTraj.timestep.getter
    def timestep(self):
        return self.trajectory[self._fslicer].timestep
    
    @TaurenTraj.n_residues.getter
    def n_residues(self):
        return self.trajectory.n_residues
    
    @TaurenTraj.n_atoms.getter
    def n_atoms(self):
        return self.trajectory.n_atoms
    
    def _remove_solvent(
            self,
            *,
            exclude=None,
            inplace=True,
            **kwargs
            ):
        """
        Removes solvent from Trajectory.
        
        Performs: MDTraj.Trajectory.remove_solvent()
        
        Parameters
        ----------
        exclude : :obj:`list`
            List of solvent residue names to retain in the new
            trajectory.
            Defaults to None.
        
        inplace : :obj:`bool`
            Whether trajectory is modified in place or a copy
            is created (returned).
            If ``True``, the :att:`~original_traj` is replaced.
            Defaults to True.
        
        Return
        ------
        Trajectory.
            **If** ``inplace`` is ``True``. Otherwise, returns ``None``.
        """
        log.info("* Removing solvent...")
        log.info(f"    received trajectory: {self.trajectory}")
        
        new_traj = self.trajectory.remove_solvent(
            inplace=False,
            exclude=exclude
            )
        
        log.info(f"    solventless trajectory: {self.trajectory}")
        
        if inplace:
            self.original_traj = new_traj
            return None
        
        else:
            return new_traj
    
    def image_molecules(
            self,
            *,
            anchor_molecules=None,
            other_molecules=None,
            sorted_bonds=None,
            make_whole=True,
            inplace=True,
            ):
        """
        Performs MDTraj.Trajectory.image_molecules, accepts same arguments.
        
        Other Parameters
        ----------------
        inplace : bool
            Whether trajectory is modified in place or a copy
            is created (returned).
            Defaults to True
        """
        log.info("* Trying imaging molecules... this can take a while...")
        
        new_traj = self.trajectory.image_molecules(
            inplace=False,
            anchor_molecules=anchor_molecules,
            other_molecules=other_molecules,
            sorted_bonds=sorted_bonds,
            make_whole=make_whole,
            )
        
        log.info("    completed.")
        
        if inplace:
            self.trajectory = new_traj
            return None
        
        else:
            return new_traj
    
    def _align_traj(self, **kwargs):
        
        log.info(
            "* IMPORTANT *"
            "align_traj method is NOT implemented in Tauren-MD"
            " for MDTraj library."
            " You should use other library option."
            " IGNORING..."
            )
        
        return
    
    def _frames2file(
            self,
            frames_to_extract,
            pdb_name_fmt,
            ):
        """
        frames_to_extract, list of frames index (int)
        pdb_name_fmt, "prefix_{FORMATTING CONDITION}.extension"
        """
        
        for frame in frames_to_extract:
        
            try:
                slice_ = self.trajectory.slice(frame, copy=True)
            
            except IndexError as e:
                log.info(self._err_frame_index.format(frame))
                log.debug(e)
                continue
        
            pdb_name = pdb_name_fmt.format(frame)
            slice_.save_pdb(pdb_name)
            log.info(f"    extracted {pdb_name}")
        
        return
    
    def _save_traj(
            self,
            file_name,
            ):
        
        self.trajectory.save(file_name, force_overwrite=True)
        
        return
    
    def _gen_chain_list(
            self,
            chains,
            ):
        
        log.debug(chains)
        
        if chains == "all":
            
            chain_list = list(range(self.trajectory.n_chains))
        
        elif isinstance(chains, str) \
                and chains.isalpha() or chains.isdigit():
            
            chain_list = [chains]
        
        elif isinstance(chains, str) and chains.count(",") > 0:
            
            chain_list = chains.split(",")
        
        elif isinstance(chains, list):
            
            try:
                chain_list = list(int(i) for i in chains.split(","))
            
            except TypeError as e:
                _ = f"Chainid values must be integers: {chains}"
                log.debug(e)
                log.info(_)
                raise TypeError(_)
        
        try:
            log.debug(chain_list)
        
        except UnboundLocalError as e:
            log.debug(e)
            log.info("* ERROR * <chain_list> not defined")
            sys.exit("* Aborting *")
        
        assert isinstance(chain_list, list), "Should be list type!"
        
        return chain_list
    
    def _calc_rmsds_combined_chains(
            self,
            chain_list,
            ref_frame,
            ):
        
        if not(all(str(s).isdigit() for s in chain_list)):
            raise ValueError(
                "MDTraj requires chainid as integer values: "
                f"given: {chain_list}")
        
        log.debug(f"chain_list: {chain_list}")
        
        chain_selector = self._gen_selector(
            chain_list,
            selection="chainid",
            boolean="or",
            )
        
        sliced_traj = self._atom_slice_traj(chain_selector)
        
        log.debug(f"len sliced_traj: {len(sliced_traj)}")
        
        combined_rmsds = self._calc_rmsds(
            sliced_traj,
            ref_frame=ref_frame,
            )
        
        log.debug(f"combined_rmsds: {combined_rmsds.shape}")
        
        assert combined_rmsds.size == self.n_frames, (
            f"combined_rmsds size '{combined_rmsds.size}' NOT matching"
            f" n_frames '{self.n_frames}'."
            )
        return combined_rmsds
    
    def _atom_slice_traj(self, selector):
        """
        Slices trajectory according to selector.
        Returns a sliced_traj.
        """
        
        slicer = self.original_traj.topology.select(selector)
        
        try:
            sliced_traj = self.trajectory.atom_slice(slicer, inplace=False)
        
        except IndexError as e:
            log.debug(e)
            log.info(
                f"* ERROR * Could not slice traj using slicer '{selector}'."
                f" Most likely the slicer does NOT share selection"
                f" with the general atom selection '{self.atom_selection}'.\n"
                "* Aborting calculation *"
                )
            sys.exit(1)
        
        return sliced_traj
    
    @staticmethod
    def _calc_rmsds(
            trajectory,
            *,
            ref_frame=0,
            **kwargs
            ):
        """
        Calculates trajectory RMSDs.
        
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The trajectory upon which operate.
        
        ref_frame : int, optional
            Defaults to 0.
            The reference frame for RMSDs calculation.
        
        Return
        ------
        np.ndarray : float
            The calculated RMSDs.
        """
        rmsds = mdtraj.rmsd(
            trajectory,
            trajectory,
            frame=ref_frame,
            parallel=True,
            precentered=False,
            )
    
        log.debug(
            f"<rmsds>: max {rmsds.max()},"
            f" min {rmsds.min()},"
            f" average {rmsds.mean()}",
            )
        
        return rmsds
    
    def _calc_rmsds_separated_chains(
            self,
            chain_list,
            ref_frame,
            ):
        
        rmsds = np.empty((self.n_frames, len(chain_list)))
        
        for index, chain in enumerate(chain_list):
            
            sliced_traj_single_chain = \
                self._atom_slice_traj(f"chainid {chain}")
            
            rmsds[:, index] = self._calc_rmsds(
                sliced_traj_single_chain,
                ref_frame=ref_frame,
                )
        
        chain_list = list(map(lambda x: str(x), chain_list))
        
        assert isinstance(rmsds, np.ndarray), "rmsds it NOT ndarray"
        assert isinstance(chain_list, list), "chain_list is NOT list"
        assert all(isinstance(c, str) for c in chain_list), (
            "all items in chains_list should be string"
            )
        return rmsds, chain_list


class TrajObservables(dict):
    """
    Stores observables obtained from traj analysis.
    """
    
    def store(self, key, data):
        """
        Stored data with key.
        
        Warns user and ignores execution if key exists.
        """
        
        if not(isinstance(key, StorageKey)):
            raise TypeError(f"key should be namedtuple, '{type(key)}' given.")
        
        if not(isinstance(data, StorageData)):
            raise TypeError(f"data sould be SorageData, '{type(data)}' given.")
        
        if key in self:
            
            log.warning(
                "* WARNING *"
                f" The storage keyword {key} already exists in the"
                " trajectory observables data base."
                )
            
            key.datatype = f"{key.datatype}_"
            
            log.warning(
                f" CHANGING KEY TO... {key}"
                )
            
        self.setdefault(key, data)
        
        return
