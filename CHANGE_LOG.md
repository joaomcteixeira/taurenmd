- [v0.4.0](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v0.4.0)
- [v0.3.2](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v0.3.2)
- [v0.3.1](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v0.3.1)
- [v0.3.0](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v0.3.0)
- [v0.2.0](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v0.2.0)
- [v0.1.1](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v0.1.1)
- [v0.0.0](https://github.com/joaomcteixeira/Tauren-MD/blob/master/CHANGE_LOG.md#v0.0.0)

# v0.4.0

- Refactore of the Tauren module architecture:
  - prepared polymorphism-based interface,
  - Improved plotting configuration.
- Improved configuration file.
- Version numbering shifted to 0:
  - the current version is considered beta version until publication, therefore it was restarted to 0 (zero) and does not follow the [_major_[_minor_[_bug_]] structure, instead, [_0_[_major_[_minor/bug_]].

# v0.3.2

- Added data export function `communicate.export.save_data_array_to_file()`:
  - now each plotting routine exports used data in .csv files.

# v0.3.1

- improved names of RMSD plotting functions:
- Added RMSD plotting functions to configuration file.


# v0.3.0

-Added three plotting routines to plot RMSDs:
    - combined RMSDs,
    - RMSDs per chain in different subplots,
    - RMSDs per chain in single subplots.

# v0.2.0

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

- Updates installer [#2](https://github.com/joaomcteixeira/Tauren-MD/pull/2)

# v0.0.0

v0 is a simple yet production ready version, though open to refactor and redesign.

**Hello World**