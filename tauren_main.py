#! /home/joao/Programming/Tauren-MD/miniconda/envs/taurenmd/bin/python
"""
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
"""
import sys
import os
import argparse
from pathlib import Path

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
        os.path.dirname(os.path.realpath(__file__)),
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
