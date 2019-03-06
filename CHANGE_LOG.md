# Index

- [v0.5.3](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v053)
- [v0.5.2](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v052)
- [v0.5.1](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v051)
- [v0.5.0](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v050)
- [v0.4.1](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v041)
- [v0.4.0](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v040)
- [v0.3.2](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v032)
- [v0.3.1](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v031)
- [v0.3.0](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v030)
- [v0.2.0](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v020)
- [v0.1.1](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v011)
- [v0.0.0](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v000)

# v0.5.3
[Back to Index][1]

- improved README and Documentation files.
- added "mdanalysis" in `taurenmd` executable options.

# v0.5.2
[Back to Index][1]

- corrected bug in `bin/update` script
  - now updates `docs` folder.

# v0.5.1
[Back to Index][1]

- implemented `align_traj` method
- improved `slicing` methods and distribution

# v0.5.0
[Back to Index][1]

- Improve code architecture in `tauren.tauren`:
  - `TaurenTraj` is not ABC class.
- Implemented initial interface for MDAnalysis
- changed method `frames2PDB` to `frames2file`.
- `overwrite` in `.save_traj` is deprecated.
- `reduce_equidistant` deprecated in favor of `slice`.

# v0.4.1
[Back to Index][1]

- Actions in config JSON file can be commented with `#` for deactivation.
- Updated `core` module with docstrings.

# v0.4.0
[Back to Index][1]

- Refactored Tauren module's architecture:
  - prepared polymorphism-based interface,
  - Improved plotting configuration.
- Improved configuration file.
- Version numbering shifted to 0:
  - the current version is considered beta version until publication, therefore it was restarted to 0 (zero) and does not follow the [_major_[_minor_[_bug_]] structure, instead, [_0_[_major_[_minor/bug_]].

# v0.3.2
[Back to Index][1]

- Added data export function `communicate.export.save_data_array_to_file()`:
  - now each plotting routine exports used data in .csv files.

# v0.3.1
[Back to Index][1]

- improved names of RMSD plotting functions:
- Added RMSD plotting functions to configuration file.


# v0.3.0
[Back to Index][1]

-Added three plotting routines to plot RMSDs:
    - combined RMSDs,
    - RMSDs per chain in different subplots,
    - RMSDs per chain in single subplots.

# v0.2.0
[Back to Index][1]

Improves and establishes general architecture and design [#5](https://github.com/joaomcteixeira/Tauren-MD/pull/5):

- function argument input is validated via decorators.
- in order to allow compatibility with future implementations:
    - all functions receive *argsand **kwargsparameters,
    - functions return (traj, ) tuple.
- corrected implementation errors, for example, checking argument validity via assert.
- improved library organization and names.
- adopts usage of f-strings.
- improves implementation of some functions.
- improves logging configuration and messages systemwide.
- adds VISION: problem identification, solution and implementation.
- adds CHANGE_LOG.

# v0.1.1
[Back to Index][1]

- Updates installer [#2](https://github.com/joaomcteixeira/Tauren-MD/pull/2)

# v0.0.0
[Back to Index][1]

v0 is a simple yet production ready version, though open to refactor and redesign.

**Hello World**


[1]: https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#Index
