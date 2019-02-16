<!--
   ease4lmp documentation master file, created by
   sphinx-quickstart on Tue Feb 12 09:26:07 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
-->

# Welcome to ease4lmp's documentation!

```eval_rst
.. toctree::
   :maxdepth: 1
   :caption: Contents:
   :hidden:

   apis/modules
```

`ease4lmp` aims at preparing [Lammps](https://github.com/lammps/lammps)' simulations using [ASE (Atomic Simulation Environment)](https://gitlab.com/ase/ase).
ASE has a lot of functionalities to manipulate atoms, and can be interfaced with many molecular dynamics (MD) and quantum mechanics (QM) simulators (including Lammps).
However, it lacks two abilities essential for conducting classical MD simulations: (1) considering bonds connecting atoms and angles/dihedrals/impropers defined by these bonds, and (2) setting different forcefield types to atoms with the same atomic number.
To make up the lack of abilities, `ease4lmp` has `BondedAtoms` class extending `ase.Atoms` class with 'bonds' data and 'types' data.
In addition, after creating and customizing atomic data (positions, velocities, etc.) using a `BondedAtoms` instance, a Lammps' data (or molecule) file required for simulation can be created from the instance by using a instance of `ease4lmp`'s `LammpsWriter` class.

## Prerequisites

* [ase/ase](https://gitlab.com/ase/ase)
* [numpy/numpy](https://github.com/numpy/numpy)

## Install

Clone and install.

```bash
git clone https://github.com/irisTa56/ease4lmp.git
cd ease4lmp
python setup.py install
```

Download and install from this repository using pip.

```bash
pip install git+https://github.com/irisTa56/ease4lmp.git
```

## Example Notebooks

* [Making Diamond](https://nbviewer.jupyter.org/github/irisTa56/ease4lmp/blob/master/examples/making_diamond.ipynb)
* [Simple Demo for `ease4lmp.LammpsWriter`](https://nbviewer.jupyter.org/github/irisTa56/ease4lmp/blob/master/examples/lammps_writer.ipynb)

## Acknowledgement

This project would not be possible without the following great open-source projects.

* [ase/ase](https://gitlab.com/ase/ase)
* [numpy/numpy](https://github.com/numpy/numpy)

## Indices and tables

* [Index](genindex)
* [Module Index](modindex)
* [Search](search)
