"""
Specific variables for the host project.

THIS FILE WAS ADAPTED FROM TREE-OF-LIFE PROJECT (version 1.1.1 - LGPLv3)
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
# CONFIGURE ACCORDING TO THE HOST PROJECT

# host project name, this will be displayed in messages and throughout
software_name = "Tauren-MD"

# host project version
software_version = (1, 0, 1)

# provide a link and e-mail with further documentation on the install process
install_wiki = "https://github.com/joaomcteixeira/Tauren-MD/wiki"
mailist = "https://github.com/joaomcteixeira/Tauren-MD/issues"

# min GB required to install the host project
# usually depends on the Miniconda ENV
min_space_allowed = 5

# name of the LOG files
installation_log_name = "install.log"
update_log_name = "update.log"

# the Miniconda ENV file specific of the host project
env_file = "taurenmd.yml"

# Miniconda installation folder
miniconda_folder = "miniconda"

# the path where the host project is hosted
# serves updating purposes
new_version_url = \
    "https://github.com/joaomcteixeira/Tauren-MD/archive/master.zip"
    
# temporary ZIP file for the new version during update
new_version_zip = "master.zip"

# folders to remove during the update
folders_to_remove = ["install", "tauren", ".github"]
