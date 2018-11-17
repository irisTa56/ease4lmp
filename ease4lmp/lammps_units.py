"""
@file ease4lmp/lammps_units.py
@brief This file contains unit data for each *Lammps unit*
(currently only `real` and `si` are supported).
@author Takayuki Kobayashi
@date 2018/05/19

Pease see also the followings for details:
* lammps/src/force.h
* lammps/src/update.cpp
* http://lammps.sandia.gov/doc/units.html
* https://physics.nist.gov/cuu/Constants/Table/allascii.txt
  * (Avogadro constant) = 6.022140857 * 10^23

From 2019/05/20:
* Avogadro constant is *exactly* 6.02214076*10^23 [mol^-1]
* Boltzmann constant is *exactly* 1.380649*10^−23 [J/K]
* Planck constant is *exactly* 6.626070150*10^-34 [J*s]
* Elementary charge is *exactly* 1.602176634*10^−19 [C]
"""

## Values of physical constants in each *unit* used in Lammps*.
lmp_units = {
  "real": {
    "boltz": 0.0019872067,                                              # Boltzmann constant (eng/degree-K)
    "hplanck": 95.306976368,                                            # Planck's constant (energy*time)
    "mvv2e": 48.88821291 * 48.88821291,                                 # conversion of mv^2 to energy
    "ftm2v": 1.0 / 48.88821291 / 48.88821291,                           # conversion of ft/m to velocity
    "mv2d": 1.0 / 0.602214129,                                          # conversion of mass/volume to density
    "nktv2p": 68568.415,                                                # conversion of NkT/V to pressure
    "qqr2e": 332.06371,                                                 # conversion of q^2/r to energy
    "qe2f": 23.060549,                                                  # conversion of qE to force
    "vxmu2f": 1.4393264316e4,                                           # conversion of vx dynamic-visc to force
    "xxt2kmu": 0.1,                                                     # conversion of xx/t to kinematic-visc
    "e_mass": 1.0 / 1836.1527556560675,                                 # electron mass
    "hhmrr2e": 0.0957018663603261,                                      # conversion of (hbar)^2/(mr^2) to energy
    "mvh2r": 1.5339009481951,                                           # conversion of mv/hbar to distance
    "angstrom": 1.0,
    "femtosecond": 1.0,
    "qelectron": 1.0,
    "_2si": {                                                           # conversion to SI units
      "mass": 1 / 6.022140857 * 1e-26,
      "distance": 1.0e-10,
      "time": 1.0e-15,
      "energy": 4.184 / 6.022140857 * 1e-20,
      "velocity": 1.0e5,
      "force": 4.184 / 6.022140857 * 1e-10,
      "torque": 4.184 / 6.022140857 * 1e-20,
      "temperature": 1.0,
      "pressure": 101325,
      "dynamic viscosity": 0.1,
      "charge": 1.6021765e-19,
      "dipole": 1.6021765e-29,
      "electric field": 1e10,
      "density": 1e3,                                                   # dim = 3
    }
  },
  "si": {
    "boltz": 1.3806504e-23,
    "hplanck": 6.62606896e-34,
    "mvv2e": 1.0,
    "ftm2v": 1.0,
    "mv2d": 1.0,
    "nktv2p": 1.0,
    "qqr2e": 8.9876e9,
    "qe2f": 1.0,
    "vxmu2f": 1.0,
    "xxt2kmu": 1.0,
    "e_mass": 0.0,                                                      # not yet set
    "hhmrr2e": 0.0,
    "mvh2r": 0.0,
    "angstrom": 1.0e-10,
    "femtosecond": 1.0e-15,
    "qelectron": 1.6021765e-19,
    "_2si": {
      "mass": 1.0,
      "distance": 1.0,
      "time": 1.0,
      "energy": 1.0,
      "velocity": 1.0,
      "force": 1.0,
      "torque": 1.0,
      "temperature": 1.0,
      "pressure": 1.0,
      "dynamic viscosity": 1.0,
      "charge": 1.0,
      "dipole": 1.0,
      "electric field": 1.0,
      "density": 1.0,
    }
  }
}
