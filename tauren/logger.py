"""
Provides a logger for Tauren-MD. Uses Python logging system.
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
import logging
import logging.config

log_file_name = "tauren-md.log"

tauren_log_config = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "debug_format": {
            "format": (
                "%(asctime)s - "
                "%(levelname)s - "
                "%(filename)s:%(name)s:%(funcName)s:%(lineno)d - "
                "%(message)s"
                )
            },
        "info_format": {}
        },
    
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "level": "INFO",
            "formatter": "info_format",
            "stream": "ext://sys.stdout"
            },
        "debug_file_handler": {
            "class": "logging.handlers.RotatingFileHandler",
            "level": "DEBUG",
            "formatter": "debug_format",
            "filename": log_file_name,
            "maxBytes": 10485760,
            "backupCount": 20,
            "encoding": "utf8"
            }
        },
    
    "loggers": {},
    
    "root": {
        "level": "DEBUG",
        "handlers": ["console", "debug_file_handler"]
        }
    }
"""
The Tauren-MD logger configuration.
"""


def get_log(name):
    """
    Returns a configured logger according to :const:`~tauren_log_config`.
    """
    
    logging.config.dictConfig(tauren_log_config)
    return logging.getLogger(name)
