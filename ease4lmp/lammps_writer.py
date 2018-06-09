"""
This file is for LammpsWriter class.

create: 2018/05/30 by Takayuki Kobayashi
"""

from .bonded_atoms import BondedAtoms
from .lammps_data \
  import LammpsAtoms, LammpsDataLines, LammpsTopology, LammpsSpecialBonds
from .lammps_units import lmp_units

from ase.calculators.lammpsrun import Prism

import ase.units as au
import datetime as dt
import numpy as np

topo_keys = ("bond", "angle", "dihedral", "improper")

class LammpsWriter:
  """
  This class ...
  """

  def __init__(
    self, atoms, atom_style, lammps_unit="real", special_bonds=False,
    tag_is_type=True, **kwargs):
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

    self._prism = ExtendedPrism(pos_unit*atoms.get_cell())

    positions = self._prism.transform_to_lammps(
      pos_unit*atoms.get_positions())

    velocities = atoms.get_velocities()
    if velocities is not None:
      vel_unit = pos_unit * au.s * lmp_unit["_2si"]["time"]
      velocities = self._prism.transform_to_lammps(vel_unit*velocities)

    self._atom_types = list(atoms.get_tags()) if tag_is_type else None

    self._lmp_atoms = LammpsAtoms(
      atom_style, positions, velocities, self._atom_types)

    topo_data = dict(zip(topo_keys, (
      atoms.get_bonded_bonds(), atoms.get_bonded_angles(),
      atoms.get_bonded_dihedrals(), atoms.get_bonded_impropers()))) \
      if isinstance(atoms, BondedAtoms) else {}

    self._topo = {
      k: LammpsTopology.create(k, v, **kwargs)
      for k, v in topo_data.items()}

    self._special_bonds = LammpsSpecialBonds(atoms.get_bonds_per_atom()) \
      if isinstance(atoms, BondedAtoms) else None

    if self._atom_types is not None:
      mass_unit = 1 / au.kg / lmp_unit["_2si"]["mass"]
      mass_dict = dict(sorted(dict(zip(
        self._atom_types, mass_unit*atoms.get_masses())).items()))
      self._mass_lines = LammpsDataLines(
        "Masses", "{type:4d} {mass:10.6f}").set_data({
        "type": mass_dict.keys(), "mass": mass_dict.values()})

  def get_required_datanames(self):
    """
    This method ...
    """
    return self._lmp_atoms.get_required_datanames()

  def get_required_datanames_for_molecule(self):
    """
    This method ...
    """
    return self._lmp_atoms.get_required_datanames(molecule=True)

  def get_sequence_patterns(self, key):
    """
    This method ...
    [Arguments]
    * key: <str>
    """
    if self._atom_types is None:
      RuntimeError("Please set atom's type in advance")
    return self._topo[key].get_sequence_patterns(self._atom_types)

  def get_bond_patterns(self):
    return self.get_sequence_patterns("bond")

  def get_angle_patterns(self):
    return self.get_sequence_patterns("angle")

  def get_dihedral_patterns(self):
    return self.get_sequence_patterns("dihedral")

  def get_improper_patterns(self):
    return self.get_sequence_patterns("improper")

  def get_maximum_per_atom(self, key):
    """
    This method ...
    [Arguments]
    * key: <str>
    """
    print("You might need to set 'extra/{}/per/atom' to: {}".format(
      key, self._topo[key].get_maximum_per_atom()))

  def get_max_bonds_per_atom(self):
    self.get_maximum_per_atom("bond")

  def get_max_angles_per_atom(self):
    self.get_maximum_per_atom("angle")

  def get_max_dihedrals_per_atom(self):
    self.get_maximum_per_atom("dihedral")

  def get_max_impropers_per_atom(self):
    self.get_maximum_per_atom("improper")

  def get_max_specials_per_atom(self):
    print("You might need to set 'extra/special/per/atom' to: {}".format(
      self._special_bonds.get_maximum_per_atom()))

  def set_atom_data(self, **kwargs):
    """
    This method ...
    [Arguments]
    * kwargs:
      * (key) name of atom-value
      * (val) <list/tuple/numpy.ndarray>
    """
    if "type" in kwargs:
      self._atom_types = list(kwargs["type"])
    self._lmp_atoms.set_data(**kwargs)

  def set_atom_types(self, types):
    self.set_atom_data(type=types)

  def set_masses(self, mass_dict):
    """
    This method ...
    [Arguments]
    * mass_dict: <dict>; <int> => <float> (in lammps' unit)
    """
    mass_dict = dict(sorted(mass_dict.items()))
    self._mass_lines = LammpsDataLines(
      "Masses", "{type:4d} {mass:10.6f}").set_data({
      "type": mass_dict.keys(), "mass": mass_dict.values()})

  def set_topology_types(self, **kwargs):
    """
    This method ...
    [Arguments]
    * kwargs:
      * (key) concrete name (key) of sequence.
      * (val) <dict> (<tuple> => <int>): from sequence to type
    """
    if self._atom_types is None:
      RuntimeError("Please set atom's type in advance")
    for k, v in kwargs.items():
      self._topo[k].set_types(v, self._atom_types)

  def set_bond_types(self, seq_to_type):
    self.set_topology_types(bond=seq_to_type)

  def set_angle_types(self, seq_to_type):
    self.set_topology_types(angle=seq_to_type)

  def set_dihedral_types(self, seq_to_type):
    self.set_topology_types(dihedral=seq_to_type)

  def set_improper_types(self, seq_to_type):
    self.set_topology_types(improper=seq_to_type)

  def write_lammps_data(self, path, mass=False, **kwargs):
    """
    This method ...
    [Arguments]
    * path: <str>
    * mass: <bool>; on/off Masses section
    * kwargs: contains ...
      1. 'velocity' on/off Velocities section
      2. 'num_*_type'
    """

    num_atom = self._lmp_atoms.get_num()
    num_atom_type = kwargs["num_atom_type"] \
      if "num_atom_type" in kwargs else self._lmp_atoms.get_num_type()

    num_topo = {
      k: v.get_num() if v.get_num_type() else 0
      for k, v in self._topo.items()} \
      if self._topo else dict.fromkeys(topo_keys, 0)
    num_topo_type = {k: v.get_num_type() for k, v in self._topo.items()} \
      if self._topo else dict.fromkeys(topo_keys, 0)

    for k in topo_keys:
      name = "num_{}_type".format(k)
      if name in kwargs:
        num_topo_type[k] = kwargs[name]

    with open(path, "w") as f:

      f.write("# written by ease4lmp.LammpsWriter at {}\n".format(
        dt.datetime.now().strftime("%Y/%m/%d %H:%M:%S")))

      f.write("\n{} atoms\n".format(num_atom))
      f.write("".join([
        "{} {}s\n".format(num_topo[k], k) for k in topo_keys]))

      f.write("\n{} atom types\n".format(num_atom_type))
      f.write("".join([
        "{} {} types\n".format(num_topo_type[k], k) for k in topo_keys]))

      xhi, yhi, zhi, xy, xz, yz = self._prism.get_lammps_prism()

      f.write("\n{}".format("".join([
        "{0:.8e} {1:.8e}  {2}lo {2}hi\n".format(0, hi, x)
        for hi, x in zip([xhi, yhi, zhi], ["x", "y", "z"])])))

      if self._prism.is_skewed():
        f.write("\n{:.8e} {:.8e} {:.8e}  xy xz yz\n".format(
          *map(float, (xy, xz, yz))))

    if mass:
      self._mass_lines.write(path)

    self._lmp_atoms.write_lines(path, **kwargs)

    for v in self._topo.values():
      if 0 < v.get_num() and 0 < v.get_num_type():
        v.write_lines(path)

  def write_lammps_molecule(self, path, special_bonds=True):
    """
    This method ...
    [Arguments]
    * path: <str>
    """

    num_atom = self._lmp_atoms.get_num()

    num_topo = {
      k: v.get_num() if v.get_num_type() else 0
      for k, v in self._topo.items()} \
      if self._topo else dict.fromkeys(topo_keys, 0)

    with open(path, "w") as f:

      f.write("# written by ease4lmp.LammpsWriter at {}\n".format(
        dt.datetime.now().strftime("%Y/%m/%d %H:%M:%S")))

      f.write("\n{} atoms\n".format(num_atom))
      f.write("".join([
        "{} {}s\n".format(num_topo[k], k) for k in topo_keys]))

    self._lmp_atoms.write_lines_for_molecule(path)

    for v in self._topo.values():
      if 0 < v.get_num() and 0 < v.get_num_type():
        v.write_lines(path)

    if special_bonds and hasattr(self, "_special_bonds"):
      self._special_bonds.write_lines(path)

class ExtendedPrism(Prism):

  def transform_to_lammps(self, vectors):
    """
    This method ...
    [Arguments]
    * path: <str>
    """
    if self.is_skewed():
      return np.dot(vectors, self.R.round(int(-self.car_prec.log10())))
    else:
      return vectors