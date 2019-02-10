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

import numpy as np
import mdtraj

from tauren import logger

log = logger.get_log(__name__)

_wrong_traj_type = "Trajectory argument should be a {} trajectory type."
_msg_key_exists = (
    "The storage keyword {} already exists in the "
    "trajectory observables data base. "
    "Ignoring..."
    )
_msg_wrong_chains_list = (
    "The list of chains provided '{}' "
    "contains values which are not of integer type"
    )
_msg_wrong_param_type = (
    "The parameter '{}' is of wrong type '{}' ",
    "should be '{}'."
    )


class TaurenTraj:
    """
    Stores and handles trajectory data.
    
    Parameters
    ----------
    trajectory : Molecular Dynamics trajectory
    
    topology : topology used to initiate the trajectory
        
    Exceptions:
        raise TypeError if trajectory is not of trajtype.type.
    """
    
    def __init__(self):
        
        self.observables = TrajObservables()
        self._rmsds_counter = 0
        
        return
    
    @staticmethod
    def _get_pdb_name_format(frames):
        """
        Creates a formatting string to apply leading zeros.
        
        Parameters
        ----------
        frames : int or str
            The number of frames.
        """
        
        leading_zeros = str(len(str(frames)))
        pdb_name_fmt = "{:0>" + leading_zeros + "}.pdb"
        
        return pdb_name_fmt
    
    @staticmethod
    def _int_slicer(slicer, top, *args, **kwargs):
        """
        Returns a list of (integer) frames from a slicer string.
        
        Slicer can be of the following formats (INCLUSIVE):
            - "0"
            - "0,1,5,100"
            - "0:"
            - ":500"
            - "10:50"
            - "10:100:2"
        """
        
        if not isinstance(slicer, str):
            raise ValueError(f"<slicer> '{slicer}' should be of string type.")
        
        log.debug(f"<slicer>: {slicer}")
        
        list_of_slices = 0
        
        possible_inputs = (
            slicer.isdigit(),
            (slicer.find(":") < 0 and slicer.count(",") > 0),
            slicer.endswith(":"),
            slicer.startswith(":"),
            len(slicer.split(":")) == 2,
            slicer.count(":") == 2,
            )
        
        dict_of_actions = {
            0: lambda x: [int(x)],
            1: lambda x: [int(i) for i in x.split(",")],
            2: lambda x: list(range(int(x[:-1]), top + 1)),
            3: lambda x: list(range(0, int(x[1:]) + 1)),
            4: lambda x: list(range(
                int(x.split(":")[0]),
                int(x.split(":")[-1]) + 1)
                ),
            5: lambda x: list(range(
                int(x.split(":")[0]),
                int(x.split(":")[1]) + 1,
                int(x.split(":")[-1]),
                ))
            }
        
        log.debug(possible_inputs)
        
        slicer_not_valid = f"<slicer> input not valid: '{slicer}'"
        
        try:
            which_action = possible_inputs.index(True)
        
        except ValueError:
            log.exception(slicer_not_valid)
            sys.exit(1)
        
        log.debug(f"<index>: {which_action}")
        
        try:
            list_of_slices = dict_of_actions[which_action](slicer)
        
        except ValueError:
            log.exception(slicer_not_valid)
            sys.exit(1)
        
        except Exception:
            log.exception(slicer_not_valid)
            sys.exit(1)
        
        assert all([isinstance(i, int) for i in list_of_slices]), \
            "Not all items in slicer are integers"
        
        return list_of_slices
    
    @staticmethod
    def _get_items_list(slicer, max_num_of_items):
        """
        Generates a list of indexes for the object_
        based on slicer.
        
        Parameters
        ----------
        slicer : str or list of integers
        
        max_num_of_items : int
            The maximum number of items
        """
        
        if isinstance(slicer, list) \
                and all([isinstance(i, int) for i in slicer]):
            
            items_list = slicer
        
        elif slicer == "all":
            
            items_list = list(range(max_num_of_items))
        
        else:
            
            items_list = TaurenTraj._int_slicer(slicer, max_num_of_items)
        
        log.debug(f"<intems_list>: {items_list}")
        
        return items_list


class TaurenMDTraj(TaurenTraj):
    
    def __init__(self, trajectory, topology):
        
        super().__init__()
        
        self.trajectory = mdtraj.load(trajectory, top=topology)
        
        self.topology = self.trajectory.topology
        
        self._update_frames_list(self.trajectory.n_frames, 1)
        
        return
    
    def _update_frames_list(
            self,
            end,
            step,
            ):
        """
        Updates frames range.
        
        Useful to represent accurate plot axis.
        """
        
        index_update = list(range(0, end, step))
        
        self.frames_list = list(map(
            lambda x: x + 1,
            index_update,
            ))
        
        log.debug(f"frame_range updated: {self.frames_list}")
        
        return
    
    def slice_in_chains(
            self,
            chain_list,
            operator="or",
            ):
        """
        Slices the trajectory according to a list of chains.
        
        Parameters
        ----------
        chain_list : list of ints or int convertable strings
            Example: [1, 4, 6, 7]
        
        operator : str, optional ["or", "and"], defaults "or"
            The logical operator to perform slicing.
            According to MDTraj.
        
        Returns
        -------
        The sliced trajectory.
        
        Exceptions
        ----------
        IndexError
            If trajectory can not be sliced according to chains.
        """
        
        selector_traj_str = [f"chainid {chain}" for chain in chain_list]
        selector = f" {operator} ".join(selector_traj_str)
        slicer = self.trajectory.topology.select(selector)
        
        try:
            sliced_traj = self.trajectory.atom_slice(slicer, inplace=False)
        
        except IndexError:
            log.expection("Could not slice traj")
            sys.exit(1)
        
        return sliced_traj
    
    def slice(
            self,
            start,
            stop,
            step=1,
            inplace=True,
            ):
        """
        Slices trajectory.
        
        Parameters
        ----------
        start : int
            The starting frame
        
        stop : int
            The ending frame (EXCLUDED)
        
        step : int, optional
            Slicing step.
            Defaults to 1.
        
        inplace : bool
            Whether trajectory is modified in place or a copy
            is created (returned).
            Defaults to True
        """
        
        log.info("* Slicing trajectory")
        log.info(f"    current traj: {self.trajectory}")
        log.info(f"    start: {start}, stop: {stop}, step: {step}")
        
        new_traj = self.trajectory[start:stop:step]
        
        if inplace:
            self._update_frames_list(start, stop, step)
            self.trajectory = new_traj
            log.info(f"    new traj: {self.trajectory}")
            return None
        
        else:
            log.info(f"    returned traj: {new_traj}")
            return new_traj
    
    def reduce_equidistant(
            self,
            step=1,
            inplace=True,
            ):
        """
        Reduces trajectory in equidistant frames separated by
        step frames.
        
        Parameters
        ----------
        step : int
            The number of steps between equidistant frames.
            Defaults to 1.
        
        inplace : bool
            Whether trajectory is modified in place or a copy
            is created (returned).
            Defaults to True
        """
        log.info("* Reducing Trajectory to equidistant frames...")
            
        if not isinstance(step, int):
            raise ValueError(f"<step> must be integer type: '{step}'")
        
        log.info(f"    current trajectory status: {self.trajectory}")
        
        new_traj = self.trajectory[::step]
        
        log.info(f"    reduced trajectory status: {new_traj}")
        
        if inplace:
            self._update_frames_list(self.trajectory.n_frames, step)
            self.trajectory = new_traj
            return None
        
        else:
            return new_traj
    
    def remove_solvent(
            self,
            *,
            exclude=None,
            inplace=True,
            ):
        """
        Removes solvent from Trajectory.
        
        Performs: MDTraj.Trajectory.remove_solvent()
        
        Parameters
        ----------
        exclude : list
            List of solvent residue names to retain
                in the new trajectory.
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
    
    @staticmethod
    def calc_rmsds(
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
    
    def calc_rmsds_separated_chains(
            self,
            *,
            chains="all",
            storage_key="rmsds_separated_chains",
            **kwargs
            ):
        """
        Calculates RMSDs for each chain separately.
        
        Parameters
        ----------
        chains : str or list of int, optional
            Defaults to "all", all chains are used.
            
            If str type
            -----------
            The chains slicer identifying the chains upon which
            perform the combined RMSD calculation.
            
            Slicer the standard Python slicing formats
            BUT these are INCLUSIVE:
            - "0"
            - "0,1,5,100"
            - "0:"
            - ":500"
            - "10:50"
            - "10:100:2"
            
            If list if int
            --------------
            The exact chain numbers.
            
        
        ref_frame : int
            The reference frame for the RMSD calculation.
        
        storage_key : str, optional
            The first element of the key tuple with which the
            calculated RMSD data will be stored in the trajectory's
            observables' dictionary. Defaults to "rmsds_separated_chains".
        
        Exceptions
        ----------
        TypeError
            If chains is list of ints and list contains
            other then int type elements.
        
        IndexError
            If trajectory can not be sliced according to chains.
        """
        
        self._check_chains(chains)
        
        chain_list = self._get_items_list(
            chains,
            self.get_num_of_chains(),
            )
        
        rmsds = np.empty((len(chain_list), self.trajectory.n_frames))
        
        for index, chain in enumerate(chain_list):
            
            sliced_traj_single_chain = self.slice_in_chains([chain])
            
            rmsds[index, :] = self.calc_rmsds(sliced_traj_single_chain)
        
        key = (storage_key, ",".join(map(str, chain_list)))
        
        self.observables.store(key, rmsds)
        
        return key
    
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
        chains : str or list of int, optional
            Defaults to "all", all chains are used.
            
            If str type
            -----------
            The chains slicer identifying the chains upon which
            perform the combined RMSD calculation.
            
            Slicer the standard Python slicing formats
            BUT these are INCLUSIVE:
            - "0"
            - "0,1,5,100"
            - "0:"
            - ":500"
            - "10:50"
            - "10:100:2"
            
            If list of integers
            --------------
            The exact chain numbers.
            
        
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
            If chains is list of ints and list contains
            other then int type elements.
        
        IndexError
            If trajectory can not be sliced according to chains.
        """
        
        self._check_chains(chains)
        
        chain_list = self._get_items_list(
            chains,
            self.get_num_of_chains(),
            )
        
        sliced_traj = self.slice_in_chains(chain_list)
        
        combined_rmsds = self.calc_rmsds(
            sliced_traj,
            )
        
        key = (storage_key, ",".join(map(str, chain_list)))
        
        self.observables.store(key, combined_rmsds)
        
        return key
    
    def frames2PDB(
            self,
            frames="all",
            prefix="_",
            **kwargs
            ):
        """
        Extracts trajectory frames to PDB files using prefix name.
        
        Parameters
        ----------
        frames : str, optional
            Frame or range of frames to extract.
            Defaults to "all", all chains are used.
            
            If str type
            -----------
            The chains slicer identifying the chains upon which
            perform the combined RMSD calculation.
            
            Slicer the standard Python slicing formats
            BUT these are INCLUSIVE:
            - "0"
            - "0,1,5,100"
            - "0:"
            - ":500"
            - "10:50"
            - "10:100:2"
            
            If list of integers
            --------------
            The exact chain numbers.
            
        prefix  : str, optional
            The prefix name for extracted PDBs.
            Defaults to "_".
            
        """
        
        log.info("* Extracting frames...")
        
        log.debug(f"<frames>: {frames}")
        log.debug(f"<prefix>: {prefix}")
        
        frames_to_extract = self._get_items_list(
            frames,
            self.trajectory.n_frames,
            )
        
        pdb_name_fmt = \
            prefix \
            + self._get_pdb_name_format(self.trajectory.n_frames)
        
        for frame in frames_to_extract:
        
            try:
                slice_ = self.trajectory.slice(frame, copy=True)
            
            except IndexError:
                log.exception(
                    f"    frame '{frame}' does NOT exist in trajectory, ",
                    "ignoring...",
                    )
                continue
        
            pdb_name = pdb_name_fmt.format(frame)
            slice_.save_pdb(pdb_name)
            log.info(f"    extracted {pdb_name}")
        
        return
    
    def export_data(
            self,
            key,
            file_name="table.csv",
            sep=",",
            header="",
            ):
        """
        Exports data arrays to file.
        """
        
        header = f"{header}\nframes,{key[1]}"
        
        data = np.vstack((
            np.array(self.frames_list),
            self.observables[key]
            ))
        
        np.savetxt(
            file_name,
            data.T,
            fmt=["%i"] + ["%.18e"] * (data.shape[0] - 1),
            delimiter=sep,
            header=header,
            )
        
        log.info(f"saved {file_name}")
        
        return
    
    def save_traj(
            self,
            file_name="traj_output.dcd",
            overwrite=True,
            ):
        """
        Saves trajectory to file.
        
        Trajectory format is given by extension name,
        must be compatible with mdtraj.Trajectory.save().
    
        Parameters
        ----------
        file_name : str
            Name of the output trajectory file.
        
        overwrite : bool
            If file_name already exists, overwrites it.
        """
        log.info(f"* Exporting trajectory to: {file_name}")
    
        self.trajectory.save(file_name, force_overwrite=overwrite)
    
        log.info("    ... saved")
        
        return
    
    def report(self):
        
        info = f"""
* Trajectory details:

n_frames: {self.trajectory.n_frames}
n_residues: {self.trajectory.n_residues}
n_atmos: {self.trajectory.n_atoms}

time_step: {self.trajectory.timestep} ps
    or {self.trajectory.timestep / 1000} ns

total_time: {self.trajectory.time[-1]} ps
    or {self.trajectory.time[-1] / 1000} ns
"""
        
        log.info(info)
        
        return
    
    def get_num_of_chains(self):
        """
        Returns the number (integer) of chains in topology.
        """
        return sum(1 for _ in self.topology.chains)
    
    @staticmethod
    def _check_chains(chains):
        
        if isinstance(chains, str):
            return
        
        elif isinstance(chains, list):
            
            if all([isinstance(i, int) for i in chains]):
                return
            
            else:
                raise TypeError(_msg_wrong_chains_list.format(chains))
            
        else:
            raise TypeError(
                _msg_wrong_param_type.format(
                    "chains",
                    type(chains),
                    "string or list",
                    )
                )
        
        return


class TrajObservables(dict):
    """
    Stores observables obtained from traj analysis.
    """
    
    def store(self, key, data):
        """
        Stored data with key.
        
        Warns user and ignores execution if key exists.
        """
        
        if key in self:
            log.warning(_msg_key_exists)
                
        else:
            self.setdefault(key, data)
        
        return
