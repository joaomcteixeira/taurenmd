# Tauren-MD

An interface that streamlines analisis routines for Molecular Dynamics.

Designed mostly for users without programming skills, but that serves also those proficient Python programmers.

Tauren-MD wraps around high performance MD analysis libraries, currently: [MDTraj](https://github.com/mdtraj/mdtraj), [MDAnalysis](https://www.mdanalysis.org/) and [OpenMM](https://github.com/pandegroup/openmm); and other Python libraries, such as [matplotlib](https://matplotlib.org/), to allow simple yet dynamic analysis workflows.

Tauren-MD contains sets of predefined functions for _loading_, _transforming_ and _exporting_ trajectories, for calculating restraints from trajectories and also predifined plotting routines to maximize the output quality. Tauren-MD is an _À la carte_ menu where the user can easily chose what she/he wants from hers/his trajectories.

Tauren-MD does not require any programming skills to be used, you can easily configure the `tauren_config.json` file with the routines you wish to execute and run it via `tauren_main.py` (created after installation).

```
python tauren_main.py -c your_config.json
```

You can access further functionalities via:

```
python tauren_main.py -h
```

For programmers, Tauren-MD libraries can also be imported and used as such, `import tauren`.

## Installation

Tauren-MD installation is supported by the [Tree-of-Life project](https://github.com/joaomcteixeira/Tree-of-Life), which was designed with a very specific [vision](https://github.com/joaomcteixeira/Tree-of-Life/blob/master/VISION.md) aimed at users without programming skills. To install a stand-alone Taren-MD, simply:

```
python install_tauren-md.py
```

## Documentation

- Read the Tauren-MD documentation [here](https://github.com/joaomcteixeira/Tauren-MD/wiki), where all the functionalities are explained, see how easy it is to use it! `:-)`

## License

The entire Tauren-MD code comes with no liability and is licensed under the [GPL-3.0](https://github.com/joaomcteixeira/Tauren-MD/blob/master/LICENSE).

<a href="https://www.gnu.org/licenses/gpl-3.0.en.html"><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/9/93/GPLv3_Logo.svg/1200px-GPLv3_Logo.svg.png" width="75" height="37"></a>

