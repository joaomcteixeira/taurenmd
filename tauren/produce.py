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
        
    key = taurentraj.calc_rmsds_separated_chains(**calc_rmsds_separated_chains)
    
    if export_data:
        
        _update_export_data(
            export_data,
            key[1],
            "calc_rmsd_separated_chains",
            )
        
        taurentraj.export_data(key, **export_data)
    
    if plot_rmsd_chain_per_subplot:
        
        _update_multiple_plot_config(
            plot_rmsd_chain_per_subplot,
            key[1],
            "plot_rmsd_chain_per_subplot",
            )

        plot.rmsd_chain_per_subplot(
            taurentraj.frames_list,
            taurentraj.observables[key],
            **plot_rmsd_chain_per_subplot,
            )
    
    if plot_rmsd_individual_chains_one_subplot:
        
        _update_multiple_plot_config(
            plot_rmsd_individual_chains_one_subplot,
            key[1],
            "plot_rmsd_individual_chains_one_subplot",
            )
        
        plot.rmsd_individual_chains_one_subplot(
            taurentraj.frames_list,
            taurentraj.observables[key],
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
        
    key = taurentraj.calc_rmsds_combined_chains(**calc_rmsds_combined_chains)
    
    if export_data:
        
        _update_export_data(
            export_data,
            key[1],
            "calc_rmsd_combined_chains",
            )
        
        taurentraj.export_data(key, **export_data)
    
    if plot_rmsd_combined_chains:
        
        _update_single_plot_config(
            plot_rmsd_combined_chains,
            key[1],
            "plot_rmsd_combined_chains",
            )
        
        plot.rmsd_combined_chains(
            taurentraj.frames_list,
            taurentraj.observables[key],
            **plot_rmsd_combined_chains,
            )
        
    return


def _get_key_list(key):
    return key.split(",")


def _update_export_data(
        kwargs,
        key,
        name,
        ):
    
    key_list = _get_key_list(key)
    
    if kwargs["file_name"] is None:
        kwargs["file_name"] = f"{name}_{'-'.join(key_list)}.csv"


def _update_single_plot_config(
        kwargs,
        key,
        name,
        ):
    
    key_list = _get_key_list(key)
    
    if kwargs["label"] is None:
        kwargs["label"] = key_list
    
    if kwargs["fig_name"] is None:
        kwargs["fig_name"] = f"{name}_{'-'.join(key_list)}.pdf"


def _update_multiple_plot_config(
        kwargs,
        key,
        name,
        ):
    
    key_list = _get_key_list(key)
    
    if kwargs["labels"] is None:
        kwargs["labels"] = key_list
    
    if kwargs["fig_name"] is None:
        kwargs["fig_name"] = f"{name}_{'-'.join(key_list)}.pdf"
    
    if kwargs["colors"] is None:
        kwargs.pop("colors")
