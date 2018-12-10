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
import argparse

from tauren import system, logger
from tauren.core import openlib

log = logger.get_log(__name__)
ap = argparse.ArgumentParser(description=__doc__)

# Mandatory
ap.add_argument(
    'config',
    help="Tauren-MD configuration JSON file"
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
        log_msg = "Performing '{}' with args: '{}'"
        log.info(log_msg.format(action, arguments[1]))
        traj = system.actions_dict[action](traj, **arguments[1])
