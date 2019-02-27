"""
TAUREN OBJECTS

Copyright © 2018-2019 Tauren-MD Project

Contributors to this file:
- João M.C. Teixeira (https://github.com/joaomcteixeira)

Tauren-MD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Tauren-MD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tauren-MD. If not, see <http://www.gnu.org/licenses/>.
"""
import sys
import string

import numpy as np
from abc import ABC, abstractmethod

import mdtraj
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD as mdarmsd

from tauren import logger

log = logger.get_log(__name__)


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
    def n_frames(self):
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
    
    @abstractmethod
    def image_molecules(self):
        """
        Images molecules.
        """
        # Do whatever needs to be done to start using the trajectory
        # without solvent. The implementation may differ drastically
        # from library used.
        pass
    
    def slice(
            self,
            start=1,
            end=None,
            step=1,
            ):
        """
        Slices trajectory for the subsequent operations.
        
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
            
            With str, use: "1" or comma separated identifers, "1,2,4".
            
            With list, use a list of identifiers: [1,2,4] or ["A", "D"].
            
            Remember that:
            # when using **MDAnalysis**, identifiers are the segid
            characters.
            # when using **MDTraj**, identifiers are digits that
            represent chain order.
        
        ref_frame : int
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
        
        key = (storage_key, chains)
        columns = f"frames,combined_{chains}"
        self.observables.store(key, (columns, data))
        
        assert isinstance(key, tuple), "key should be tuple!"
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
        
        ref_frame : int
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
        
        key = (storage_key, ",".join(chains_headers))
        columns = "frames," + ",".join(chains_headers)
        self.observables.store(key, (columns, data))
        
        assert isinstance(key, tuple), "key should be tuple!"
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
        
        log.info(f"* Exporting {key} data")
        
        header = f"{key}\n{header}\n{self.observables[key][0]}"
        
        np.savetxt(
            file_name,
            self.observables[key][1],
            delimiter=sep,
            header=header,
            )
        
        log.info(f"    saved {file_name}")
        
        return


class TaurenMDAnalysis(TaurenTraj):
    
    def __init__(self, trajectory, topology):
        
        self.universe = mda.Universe(topology, trajectory)
        self.trajectory = self.universe.trajectory
        self.topology = mda.Universe(topology)
        
        super().__init__()
        
        return
    
    def _set_full_frames_list(self):
        super()._set_full_frames_list(self.trajectory.n_frames)
    
    @TaurenTraj.totaltime.getter
    def totaltime(self):
        return self.trajectory[self._fslicer][-1].time
    
    @TaurenTraj.timestep.getter
    def timestep(self):
        return (
            self.trajectory[self._fslicer][1].time
            - self.trajectory[self._fslicer][0].time
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
                self.universe.atoms.write(
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
        with mda.Writer(file_name, self.universe.atoms.n_atoms) as W:
            for ts in self.trajectory[self._fslicer]:
                W.write(self.universe.atoms)
                log.info(f"    exported {ts}")
        
        return
    
    def _calc_rmsds_combined_chains(
            self,
            chain_list,
            ref_frame,
            ):
        
        chain_selector = self._gen_selector(chain_list)
        
        log.debug(f"chain_selector: {chain_selector}")
        
        # https://www.mdanalysis.org/docs/documentation_pages/analysis/align.html#rms-fitting-tutorial
        # https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html#MDAnalysis.analysis.rms.RMSD
        R = mdarmsd(
            self.universe,
            self.topology,
            select=chain_selector,
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
        
        for ii, chain_selector in enumerate(filtered_selectors):
            
            atoms = self.universe.select_atoms(chain_selector)
            atoms_top = self.topology.select_atoms(chain_selector)
            
            R = mdarmsd(
                atoms,
                atoms_top,
                groupselection=None,
                ref_frame=ref_frame,
                verbose=False,
                )
            
            R.run(verbose=False)
            
            rmsds[:, ii] = R.rmsd[:, 2][self._fslicer]
        
        column_headers = list(map(
            lambda x: x.replace("segid ", ""),
            filtered_selectors
            ))
        
        assert isinstance(column_headers, list), "c_selectors NOT list!"
        
        return rmsds, column_headers
    
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
        
        self.trajectory = mdtraj.load(trajectory, top=topology)
        self.topology = self.trajectory.topology
        
        super().__init__()
        
        return
    
    def _set_full_frames_list(self):
        super()._set_full_frames_list(self.trajectory.n_frames)
    
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
        exclude : list
            List of solvent residue names to retain in the new
            trajectory.
            Defaults to None.
        
        inplace : bool
            Whether trajectory is modified in place or a copy
            is created (returned).
            Defaults to True
        """
        log.info("* Removing solvent...")
        log.info(f"    received trajectory: {self.trajectory}")
        
        new_traj = self.trajectory.remove_solvent(
            inplace=False,
            exclude=exclude
            )
        
        log.info(f"    solventless trajectory: {self.trajectory}")
        
        if inplace:
            self.trajectory = new_traj
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
        
        self.trajectory[self._fslicer].save(file_name, force_overwrite=True)
        
        return
    
    def _gen_chain_list(
            self,
            chains,
            ):
        
        if chains == "all":
            
            chain_list = list(range(self.trajectory.n_chains))
        
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
        
        combined_rmsds = self._calc_rmsds(
            sliced_traj[self._fslicer],
            ref_frame=ref_frame,
            )
        
        assert combined_rmsds.size == self.n_frames, (
            "combined_rmsds array size NOT valid"
            )
        return combined_rmsds
    
    def _atom_slice_traj(self, selector):
        """
        Slices trajectory according to selector.
        Returns a sliced_traj.
        """
        
        slicer = self.trajectory.topology.select(selector)
        
        try:
            sliced_traj = self.trajectory.atom_slice(slicer, inplace=False)
        
        except IndexError:
            log.exception("Could not slice traj")
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
                sliced_traj_single_chain[self._fslicer],
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
        
        if not(isinstance(key, tuple)):
            raise TypeError(f"key should be tuple, '{type(key)}' given.")
        
        if key in self:
            log.warning(
                f"The storage keyword {key} already exists in the "
                "trajectory observables data base. "
                "Ignoring..."
                )
                
        else:
            self.setdefault(key, data)
        
        return
