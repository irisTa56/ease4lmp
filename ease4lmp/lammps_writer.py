"""
This file is for LammpsWriter class.

create: 2018/05/30 by Takayuki Kobayashi
"""

from .bonded_atoms import BondedAtoms
from .lammps_data import LammpsAtoms, LammpsSequences, LammpsDataLines
from .lammps_units import lmp_units

from ase.calculators.lammpslib import convert_cell
from ase.calculators.lammpsrun import Prism

import ase.units as au
import datetime as dt
import decimal as dec
import numpy as np

seq_keys = ["bond", "angle", "dihedral", "improper"]

class LammpsWriter:
  """
  This class ...
  """

  def __init__(
    self, atoms, atom_style,
    lammps_unit="real", tag_is_type=True, **kwargs
  ):
    """
    This constructor ...
    [Arguments]
    * atoms: <BondedAtoms/ase.Atoms>
    * atom_style: <str>
    * lammps_unit: <str>
    * tag_is_type: <bool>
    * kwargs: contains ...
      1. 'class2' for improper
    """

    lmp_unit = lmp_units[lammps_unit]

    pos_unit = 1 / au.m / lmp_unit["_2si"]["distance"]

    is_skewed, transform = self._parse_box(pos_unit*atoms.get_cell())

    to_str = np.vectorize(lambda x: format(x, "15.8e"))

    positions = \
      to_str(np.dot(transform, pos_unit*atoms.get_positions().T).T) \
      if is_skewed else to_str(pos_unit*atoms.get_positions())

    velocities = atoms.get_velocities()
    if velocities is not None:
      vel_unit = pos_unit * au.s * lmp_unit["_2si"]["time"]
      velocities = \
        to_str(np.dot(transform, vel_unit*velocities.T).T) \
        if is_skewed else to_str(vel_unit*velocities)

    self._atom_types = list(atoms.get_tags()) if tag_is_type else None

    self._atoms = LammpsAtoms(
      atom_style, positions, self._atom_types, velocities)

    if isinstance(atoms, BondedAtoms):
      bonds = [
        [i +b[0] for b in bs if b[0] != 0]
        for i, bs in enumerate(atoms.get_bonds())]
      self._seqs = {
        k: LammpsSequences.create(k, bonds, **kwargs) for k in seq_keys}

    if self._atom_types is not None:
      mass_unit = 1 / au.kg / lmp_unit["_2si"]["mass"]
      masses = atoms.get_masses() * mass_unit
      mass_dict = {t: m for m, t in zip(masses, self._atom_types)}
      self._mass_lines = LammpsDataLines("Masses", "{type} {mass}")
      self._mass_lines.set_data(
        type=mass_dict.keys(), mass=mass_dict.values())

  def get_required_atom_datanames(self):
    """
    This method ...
    """
    return self._atoms.get_required_datanames()

  def set_atom_data(self, **kwargs):
    """
    This method ...
    [Arguments]
    * kwargs:
      * (key) name of atom-value
      * (val) <list/tuple/numpy.ndarray>
    """
    if "type" in kwargs.keys():
      self._atom_types = list(kwargs["type"])
    self._atoms.set_data(**kwargs)

  def set_atom_types(self, types):
    self.set_atom_data(type=types)

  def set_bond_types(self, seq_to_type):
    self.set_sequence_types(bond=seq_to_type)

  def set_angle_types(self, seq_to_type):
    self.set_sequence_types(angle=seq_to_type)

  def set_dihedral_types(self, seq_to_type):
    self.set_sequence_types(dihedral=seq_to_type)

  def set_improper_types(self, seq_to_type):
    self.set_sequence_types(improper=seq_to_type)

  def set_sequence_types(self, **kwargs):
    """
    This method ...
    [Arguments]
    * kwargs:
      * (key) concrete name of sequence. 'Bonds', 'Bond', 'bonds',
        'bond', ... are OK.
      * (val) <dict> (<tuple> => <int>): from sequence to type
    """

    if self._atom_types is None:
      RuntimeError("Please set atom's type first")

    for k, v in kwargs.items():
      k_lower = k.lower()
      name = k_lower[:-1] if k_lower[-1] == "s" else k_lower
      self._seqs[name].set_types(v, self._atom_types)

  def set_masses(self, mass_dict):
    """
    This method ...
    [Arguments]
    * mass_dict: <dict>; <int> => <float> (in lammps' unit)
    """
    self._mass_lines = LammpsDataLines("Masses", "{type} {mass}")
    self._mass_lines.set_data(
      type=mass_dict.keys(), mass=mass_dict.values())

  def write_lammps_data(self, path, mass=False, **kwargs):
    """
    This method ...
    [Arguments]
    * path: <str>
    * mass: <bool>; to on/off Masses section
    * kwargs: contains ...
      1. 'velocity' to on/off Velocities section
      2. 'num_*_types'
    """

    num_atoms = self._atoms.get_num()
    num_atom_types = kwargs["num_atom_types"] \
      if "num_atom_types" in kwargs.keys() else self._atoms.get_num_type()

    num_seqs = [self._seqs[k].get_num() for k in seq_keys]
    num_seq_types = [self._seqs[k].get_num_type() for k in seq_keys]
    seq_is_set = [
      True if isinstance(n, int) else False for n in num_seq_types]

    for i, k in enumerate(seq_keys):
      n = "num_{}_types".format(k)
      if n in kwargs.keys():
        num_seq_types[i] = kwargs[n]

    with open(path, "w") as f:
      f.write("# written by ease4lmp.LammpsWriter at {}\n".format(
        dt.datetime.now().strftime("%Y/%m/%d %H:%M:%S")))

      f.write("\n{} atoms\n".format(num_atoms))
      f.write("\n".join([
        "{} {}s".format(num_seqs[i], k)
        for i, k in enumerate(seq_keys) if seq_is_set[i]]))

      f.write("\n{} atom types\n".format(num_atom_types))
      f.write("\n".join([
        "{} {} types".format(num_seq_types[i], k)
        for i, k in enumerate(seq_keys) if seq_is_set[i]]))

      f.write("\n{}\n".format("\n".join(self._box_data)))

    if hasattr(self, "_mass_lines") and mass:
      self._mass_lines.write(path)

    self._atoms.write_lines(path, **kwargs)

    for i, k in enumerate(seq_keys):
      if seq_is_set[i]:
        self._seqs[k].write_lines(path)

  def write_lammps_molecule(self, path):
    pass

  def _parse_box(self, cell):
    """
    This method ...
    [Arguments]
    * cell: <numpy.ndarray>; unit-converted ase's cell
    """

    is_skewed = Prism(cell).is_skewed()
    box, coord_transform = convert_cell(cell)

    self._box_data = [
      "{0:14.8e} {1:14.8e} {2}lo {2}hi".format(0.0, box[i,i], x)
      for i, x in enumerate(["x","y","z"])]

    if is_skewed:
      self._box_data.append(
        "\n{:14.8e} {:14.8e} {:14.8e} xy xz yz"
        .format(box[0,1], box[0,2], box[1,2]))

    return is_skewed, coord_transform
