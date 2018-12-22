"""
Defines executables

Copyright © 2018-2019 Tauren-MD

THIS FILE WAS ADAPTED FROM TREE-OF-LIFE PROJECT (version 1.0.0 - LGPLv3)
AND MODIFIED ACCORDINGLY TO THE NEEDS OF THE TAUREN-MD PROJECT.

Visit the original Tree-of-Life project at:

https://github.com/joaomcteixeira/Tree-of-Life

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

from install import system

# interesting readings:
# https://stackoverflow.com/questions/6943208/activate-a-virtualenv-with-a-python-script
# https://halotis.com/running-python-code-in-windows-batch-file-trick/
# https://docs.python.org/3.3/using/windows.html
# finally, shebangs can be used on Windows10
# allow double click execution

# changed from Tree-of-Life v1.0.0 to fit Tauren-MD
#shebang
shebang = "#! {}\n"

# define your executable scripts
exec1_code = r"""'''
MANAGES TAUREN WORKFLOW

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
'''
import sys
import os
import argparse
from pathlib import Path

software_folder = os.path.abspath(
    os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.pardir
        )
    )

sys.path.append(software_folder)

from tauren import system, logger
from tauren.core import openlib

log_path = Path(logger.log_file_name)
if log_path.exists():
    log_path.unlink()

log = logger.get_log(__name__)
log.info("* Tauren-MD initiated!")

ap = argparse.ArgumentParser(description=__doc__)

path_to_default_config = os.path.abspath(
    os.path.join(
        software_folder,
        'tauren_config.json'
        )
    )

# Mandatory
ap.add_argument(
    '-c',
    '--config',
    help=(
        "Tauren-MD configuration JSON file. "
        "If not provided runs the default config: "
        "simply loads trajectory and shows information"
        ),
    default=path_to_default_config
    )

# Options
ap.add_argument(
    '-traj',
    '--trajectory',
    default=None,
    help="Trajectory file ({})".format(system.trajectory_types)
    )

ap.add_argument(
    '-top',
    '--topology',
    default=None,
    help="Topology file ({})".format(system.topology_types)
    )

cmd = ap.parse_args()

# set path configuration
conf = openlib.load_json_config(cmd.config)

if cmd.trajectory:
    trajectory_path = cmd.trajectory

elif conf.input_data["trajectory"]:
    trajectory_path = conf.input_data["trajectory"]

else:
    log.info("* ERROR * No trajectory file provided")
    sys.exit(1)

if cmd.topology:
    topology_path = cmd.topology

elif conf.input_data["topology"]:
    topology_path = conf.input_data["topology"]

else:
    log.info("* ERROR * No topology file provided")
    sys.exit(1)

traj = openlib.load_traj(trajectory_path, topology_path)

for action, arguments in conf.actions.items():
    if arguments[0]:
        action_name = action.rstrip("_")
        log_msg = "Performing '{}' with args: '{}'"
        log.info(log_msg.format(action_name, arguments[1]))
        
        traj = system.actions_dict[action_name](traj, **arguments[1])

log.info("* Tauren-MD completed!")

"""

update_script_code = r"""
import sys
import os
import importlib
import pathlib

software_folder = os.path.abspath(
    os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.pardir
        )
    )

sys.path.append(software_folder)

if sys.version_info[0] != 3:
    sys.stderr.write("Python 3 is required to run Updater")
    sys.exit(1)

from install import logger
from install import messages
from install import system
from install import executables

try:
    import installation_vars
except ModuleNotFoundError as e:
    print(e)
    print("* ERROR * installation_vars.py file not found")
    print("* ERROR * this file has created during installation")
    print("* ERROR * and is required for UPDATING")
    print("* ERROR * Reinstall the SOFTWARE to repair updating functions")
    print(messages.additional_help)
    print(messages.abort)
    input(messages.terminate)
    sys.exit(1)

try:
    install_dir = installation_vars.install_dir
    python_exec = installation_vars.python_exec
    install_option = installation_vars.install_option
    conda_exec = installation_vars.conda_exec
    installed_env_file = installation_vars.installed_env_file
    installed_env_name = installation_vars.installed_env_name
    installed_env_version = installation_vars.installed_env_version
    miniconda_folder = installation_vars.miniconda_folder

except AttributeError as e:
    print(messages.update_var_missing)
    print()
    print(messages.consider_reinstall)
    print(messages.additional_help)
    print(e)
    print(messages.abort)
    input(messages.terminate)
    sys.exit(1)

list_of_paths = [
    install_dir,
    python_exec,
    conda_exec,
    miniconda_folder,
    miniconda_folder
    ]

for _path in list_of_paths:
    if isinstance(_path, pathlib.Path) and not _path.exists():
        print(messages.update_var_missing)
        print(os.fspath(_path) + " path does NOT exists")
        print()
        print()
        print(messages.consider_reinstall)
        print(messages.additional_help)
        print(messages.abort)
        input(messages.terminate)
        sys.exit(1)

update_log = install_dir.joinpath(system.update_log_name)

if update_log.exists():
    update_log.unlink()

log = logger.InstallLogger(__name__, log_file_name=update_log).gen_logger()

log.debug("start updating")

from install import updater
from install import commons
from install import condamanager

upf = updater.Updater(install_dir)
upf.run()

# reloads the libs updated version
importlib.reload(system)
importlib.reload(executables)

log.info("* Checking Conda environment...")

if install_option == 1:

    if system.latest_env_version > installed_env_version:

        log.info("* A NEW Python environment version is available")
        log.info("* Software's dependencies must be updated")
    
        if os.path.exists(conda_exec):
        
            log.info("* Miniconda installation found")
            log.info("   ... starting env update")
            
            upc = condamanager.CondaManager(cwd=install_dir)
            upc.set_conda_exec(conda_exec)
            upc.set_env_name(installed_env_name)
            upc.remove_env()
            upc.set_env_file(system.latest_env_file)
            upc.install_env()
            upc.logs_env_information()
            log.info("... Conda env UPDATED")
            
            # registers installation variables
            install_option = 1
            conda_exec = upc.get_conda_exec()
            python_exec = upc.get_env_python_exec()
            installed_env_file = upc.get_env_file()
            installed_env_name = upc.get_env_name()
            installed_env_version = upc.get_env_version()
            miniconda_folder = upc.get_miniconda_install_folder()
        
        else:
            log.info("* ERROR * Could not find the CONDA executable")
            log.info(messages.something_wrong)
            log.info(messages.additional_help)
            log.info(messages.update_continues)
            log.info(messages.consider_reinstall)
    else:
        log.info("   ...Conda env already in latest version")
        log.info("")

elif install_option == 2:
    log.info(
        "* You have previously configured Python libraries manually.\n"
        "* Please check if it's necessary to update the software's \n"
        "* Python dependencies, consult .yml file in 'install' folder."
        )
    
else:
    log.info("* ERROR* We couldn't access install information")
    log.info(messages.something_wrong)
    log.info(messages.additional_help)
    log.info(messages.consider_reinstall)
    log.info(messages.abort)
    sys.exit(1)

log.info("* Updating executable files...")

commons.create_executables(install_dir, python_exec)

commons.register_install_vars(
    install_dir=install_dir,
    python_exec=python_exec,
    install_option=install_option,
    conda_exec=conda_exec,
    env_file=installed_env_file,
    env_name=installed_env_name,
    env_version=installed_env_version,
    miniconda_folder=miniconda_folder
    )

log.info(messages.update_completed)
commons.sys_exit()
"""

# executable scripts file names and extensions
exec1 = "taurenmd{}".format(system.exec_file_extension)
updatescript = "update{}".format(system.exec_file_extension)

# dictionary listing the executable scripts
# keys are file names, values the string with code
executable_files = {
    exec1: exec1_code,
    updatescript: update_script_code
    }
