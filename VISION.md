# Problem definition

Powerful libraries for the analysis of Molecular Dynamics' trajectories and related data are being actively developed, with a [Python](https://www.python.org/) interface the most known are [MDTraj](https://github.com/mdtraj/mdtraj), [MDAnalysis](https://www.mdanalysis.org/) and [OpenMM](https://github.com/pandegroup/openmm).

Many researchers with strong Biochemical and Biophysics background use Molecular Dynamics (MD) to generate _ab initio_ models or modulate their _in vitro_ or _in vivo_ acquired data. These collective of researchers use MD as a tool and not necessary have (or need to have) programming skills.

*Problem 1:* There is a growing gap between the libraries being developed and those users that must use them and lack programming skills.

*Problem 2:* Data analysis rapidly reach routine. Providing pre-defined routines for data mining, parameter calculation and plotting, considerably enhances analysis speed of multiple datasets.

# Proposed Solutions

The gap described in *problem 1* can be fulfilled either by teaching programming skills to this research collective, or asking them to be self-taught, or by providing an interface for them to use the MD analysis libraries without the need to learn how to program.

*Solution 1:* Tauren-MD aims to provide such interface for users without programming skills.

Also for users with programming skills, it is time consuming, and even challenging, to assemble routines that automate analysis steps.

*Solution 2:* Tauren-MD aims to provide predefined set routines to stream-line analysis procedures and, therefore, boost deliverables' rate.

# Implementation

## User Requirements

- Users should be able to access all Tauren-MD routines, hereafter named `actions`, through a configurable text file, hereafter named `config`.
- to configure the `config` file NO programming skills should be required. The minimum requirement being the syntax of the `config` file itself.
- the `config` file should be a list-like menu of the `actions` that the user wants to perform.
- `actions` entry should contain two fields: 1) its name and 2) its arguments (options).
- `actions` can be flagged (`true` or `false`) to (de)activate its execution.
- `actions` can be deactivated by removing them from the list.
- `actions` can be repeated.
- `actions` can be reordered.
- `config` files can optionally provide a `path` to the trajectory file.
- `config` files can optionally provide a `path` to the topology file.
- Tauren-MD execution is provided by a `bin` executable file.
- Tauren-MD executable takes three optional arguments:
    - a `path` to a config file
    - a `path` to a trajectory file
    - a `path` to a topology file
- for optional arguments not provided, a default value is used.

## Architecture

It should be possible to use Tauren-MD as an user interface or as an independent and organized, well documented, library. With that in mind:

### Installation and dependency management

Regardless of whatever additional methods are implemented, Tauren-MD installation should be possible via the [Tree-of-Life project](https://github.com/joaomcteixeira/Tree-of-Life).

### Tauren-MD executable file

Tauren-MD executable file should:

- Use `paths` from `--trajectory` and `--topology` arguments when provided;
    - otherwise use `paths` provided in the `config` file,
    - if none provided, `raise error`.
- use `config` file provided in arguments;
    - if none provided use `config` distributed by default with package,
    - default `config` should just load a `trajectory`.
- read over the `config`'s list of `actions` and execute them in order with the respective `action`'s arguments.
    - communication between `config` file `action`'s name and the `action`'s function should be provided by a dictionary in `tauren.core.system`.

### Tauren-MD actions

Tauren-MD `actions`:

- are functions defined in Tauren-MD libs.
- should always take as first positional argument the trajectory object.
- can be configured with `kwargs`;
    - these `kwargs` can be used in the `config` file `action` entries.
- should `return` the trajectory enclosed in a tuple `(traj, )`;
    - other values can be passed inside the tuple provided `traj` is indexed at `0`,
    - returned values should be `assert`ed.
- to allow compatibility with future implementations:
    - `action`'s parameter list should contain `*args` and `**kwargs`.

### Lib organization

General project organization:

- `tauren/`: _main_ lib folder
    - `core/`: libs that are used system wide, commons, decorators, helpers...
    - `communicate/`: deals with input (_read_) and output (_export_) of data
    - `transform/`: transforms trajectory or data
    - `plot/`: plotting templates
