"""
Functions that perform combined operations that are related.
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
from tauren import logger
from tauren import plot

log = logger.get_log(__name__)


def rmsds_separated_chains(
        taurentraj,
        calc_rmsds_separated_chains,
        *,
        export_data=False,
        plot_rmsd_chain_per_subplot=False,
        plot_rmsd_individual_chains_one_subplot=False,
        **kwargs
        ):
    """
    Execute routines related to RMSDs of separated chains.
    """
        
    key = taurentraj.calc_rmsds_separated_chains(**calc_rmsds_separated_chains)
    
    if export_data:
        
        _update_export_data(
            export_data,
            key,
            )
        
        taurentraj.export_data(key, **export_data)
    
    if plot_rmsd_chain_per_subplot:
        
        _update_multiple_plot_config(
            plot_rmsd_chain_per_subplot,
            key,
            "plot_rmsd_chain_per_subplot",
            taurentraj.observables[key],
            )

        plot.rmsd_chain_per_subplot(
            taurentraj.observables[key][1][:, 0],
            taurentraj.observables[key][1][:, 1:].T,
            **plot_rmsd_chain_per_subplot,
            )
    
    if plot_rmsd_individual_chains_one_subplot:
        
        _update_multiple_plot_config(
            plot_rmsd_individual_chains_one_subplot,
            key,
            "plot_rmsd_individual_chains_one_subplot",
            taurentraj.observables[key],
            )
        
        plot.rmsd_individual_chains_one_subplot(
            taurentraj.observables[key][1][:, 0],
            taurentraj.observables[key][1][:, 1:].T,
            **plot_rmsd_individual_chains_one_subplot
            )
        
    return


def rmsds_combined_chains(
        taurentraj,
        calc_rmsds_combined_chains,
        *,
        export_data=False,
        plot_rmsd_combined_chains=False,
        **kwargs
        ):
    """
    Execute routines related to RMSDs of combined chains.
    """
        
    key = taurentraj.calc_rmsds_combined_chains(**calc_rmsds_combined_chains)
    
    if export_data:
        
        _update_export_data(
            export_data,
            key,
            )
        
        taurentraj.export_data(key, **export_data)
    
    if plot_rmsd_combined_chains:
        
        _update_single_plot_config(
            plot_rmsd_combined_chains,
            key,
            "plot",
            taurentraj.observables[key],
            )
        
        plot.rmsd_combined_chains(
            taurentraj.observables[key][1][:, 0],
            taurentraj.observables[key][1][:, 1],
            **plot_rmsd_combined_chains,
            )
        
    return


def _get_key_list(key):
    """
    .. deprecated:: 0.6.0
    """
    kl = key.split(",")
    log.debug(kl)
    return kl


def _update_export_data(
        kwargs,
        key
        ):
        
    if kwargs["file_name"] is None:
        kwargs["file_name"] = f"{key.datatype}_{key.filenaming}.csv"


def _update_single_plot_config(
        kwargs,
        key,
        name,
        data,
        ):
        
    if kwargs["label"] is None:
        kwargs["label"] = data.columns[1]
    
    if kwargs["fig_name"] is None:
        kwargs["fig_name"] = f"{name}_{key.filenaming}.pdf"


def _update_multiple_plot_config(
        kwargs,
        key,
        name,
        data,
        ):
    
    if kwargs["labels"] is None:
        kwargs["labels"] = data.columns[1:]  # key_list
    
    if kwargs["fig_name"] is None:
        kwargs["fig_name"] = f"{name}_{key.filenaming}.pdf"
    
    if kwargs["colors"] is None:
        kwargs.pop("colors")
