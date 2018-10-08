"""
@file ease4lmp/lammps_writer.py
@brief This file is for `LammpsWriter` class, which is an interface to
write whole *Lammps data file* (or *molecule file*).
@author Takayuki Kobayashi
@date 2018/05/30
"""

from .bonded_atoms import BondedAtoms
from .lammps_data \
  import LammpsAtoms, LammpsDataLines, LammpsTopology, LammpsSpecialBonds
from .lammps_units import lmp_units

from ase.calculators.lammpsrun import Prism

import ase.units as au
import datetime as dt
import numpy as np


## Names of topology component.
topo_keys = ("bond", "angle", "dihedral", "improper")


class LammpsWriter:
  """
  This class is an interface to write whole  *Lammps data file*
  (or *molecule file*).
  """

  def __init__(
    self, atoms, atom_style, lammps_unit="real", tag_is_type=True,
    **kwargs):
    """
    @param atoms
      An instance of @BondedAtoms or `ase.Atoms`.
    @param atom_style
      A string denoting an *atom style* used in Lammps.
    @param lammps_unit
      A string denoting a *unit* used in Lammps.
    @param tag_is_type
      Whether to use tags of `atoms` as atom types (*tags* is one of
      properties of `ase.Atoms`).
    @param **kwargs
      A variable number of named arguments.
      One can use `class2` keyword.
    """

    lmp_unit = lmp_units[lammps_unit]
    pos_unit = 1 / au.m / lmp_unit["_2si"]["distance"]

    ## An ExtendedPrism instance.
    self._prism = ExtendedPrism(pos_unit*atoms.get_cell())

    positions = self._prism.transform_to_lammps(
      pos_unit*atoms.get_positions())

    velocities = atoms.get_velocities()
    if velocities is not None:
      vel_unit = pos_unit * au.s * lmp_unit["_2si"]["time"]
      velocities = self._prism.transform_to_lammps(vel_unit*velocities)

    ## List containing type for each atom.
    self._atom_types = list(atoms.get_tags()) if tag_is_type else None

    ## A @LammpsAtoms instance.
    self._lmp_atoms = LammpsAtoms(
      atom_style, positions, velocities, self._atom_types)

    topo_data = dict(zip(topo_keys, [
        atoms.get_bonded_bonds(), atoms.get_bonded_angles(),
        atoms.get_bonded_dihedrals(), atoms.get_bonded_impropers()
      ])) if isinstance(atoms, BondedAtoms) else {}

    ## Dictionary containing @LammpsTopology instances
    # for each topology components.
    self._topo = {
      k: LammpsTopology.create(k, v, **kwargs)
      for k, v in topo_data.items()
    }

    ## A @LammpsSpecialBonds instance.
    self._special_bonds = (
      LammpsSpecialBonds(atoms.get_bonds_per_atom())
      if isinstance(atoms, BondedAtoms) else None)

    if self._atom_types is not None:

      mass_unit = 1 / au.kg / lmp_unit["_2si"]["mass"]
      mass_dict = dict(sorted(dict(zip(
        self._atom_types, mass_unit*atoms.get_masses())).items()))

      ## A @LammpsDataLines instance for *Masses* section.
      self._mass_lines = LammpsDataLines(
        "Masses", "{type:4d} {mass:10.6f}").set_data({
        "type": mass_dict.keys(), "mass": mass_dict.values()})


  def get_required_datanames(self):
    """
    @return
      A set of required data names:
      *Lammps data file* needs those data,
      but #_lmp_atoms does not have them yet.
    """

    return self._lmp_atoms.get_required_datanames()


  def get_required_datanames_for_molecule(self):
    """
    @return
      A set of required data names:
      *Lammps molecule file* needs those data,
      but #_lmp_atoms does not have them yet.
    """

    return self._lmp_atoms.get_required_datanames(molecule=True)


  def get_sequence_patterns(self, key):
    """
    @param key
      Name of topology component:
      `bond`, `angle`, `dihedral` or `improper`.

    @return
      A set of unique sequences of atom types
      appearing in given topology components.
      Note that the number of elements in the set is not necessarily
      equal to the number of types of topology component:
      for example, sequence `(1, 2, 3)` and `(3, 2, 1)` are
      reduced to one angle type.
    """

    if self._atom_types is None:
      RuntimeError("Please set atom's type in advance")

    return self._topo[key].get_sequence_patterns(self._atom_types)


  def get_bond_patterns(self):
    """
    @return
      A set of unique sequences of atom types
      appearing in all the sequences of atoms forming bonds.
    """

    return self.get_sequence_patterns("bond")


  def get_angle_patterns(self):
    """
    @return
      A set of unique sequences of atom types
      appearing in all the sequences of atoms forming angles.
    """

    return self.get_sequence_patterns("angle")


  def get_dihedral_patterns(self):
    """
    @return
      A set of unique sequences of atom types
      appearing in all the sequences of atoms forming dihedrals.
    """

    return self.get_sequence_patterns("dihedral")


  def get_improper_patterns(self):
    """
    @return
      A set of unique sequences of atom types
      appearing in all the sequences of atoms forming impropers.
    """

    return self.get_sequence_patterns("improper")


  def print_maximum_per_atom(self, key):
    """
    @param key
      Name of topology component:
      `bond`, `angle`, `dihedral` or `improper`.

    @return
      The maximum number of given topology components per atom.
    """

    print(
      "You might need to set 'extra/{}/per/atom' to: {}"
      .format(key, self._topo[key].get_maximum_per_atom()))


  def print_max_bonds_per_atom(self):
    """
    @return
      The maximum number of bonds per atom.
    """

    self.print_maximum_per_atom("bond")


  def print_max_angles_per_atom(self):
    """
    @return
      The maximum number of angles per atom.
    """

    self.print_maximum_per_atom("angle")


  def print_max_dihedrals_per_atom(self):
    """
    @return
      The maximum number of dihedrals per atom.
    """

    self.print_maximum_per_atom("dihedral")


  def print_max_impropers_per_atom(self):
    """
    @return
      The maximum number of impropers per atom.
    """

    self.print_maximum_per_atom("improper")


  def print_max_specials_per_atom(self):
    """
    @return
      The maximum number of special bonds per atom.
    """

    print(
      "You might need to set 'extra/special/per/atom' to: {}"
      .format(self._special_bonds.get_maximum_per_atom()))


  def set_atom_data(self, **kwargs):
    """
    @param **kwargs
      A variable number of named arguments.
      Keywords are property names, and arguments are lists
      (or other array-like objects) containing property values.

    This method sets data to be written in *Atoms* (and *Velocities*)
    section to #_lmp_atoms. If keyword `type` exists, this method
    assigns the argument to #_atom_types.
    """

    if "type" in kwargs:
      self._atom_types = list(kwargs["type"])

    self._lmp_atoms.set_data(**kwargs)


  def set_masses(self, mass_dict):
    """
    @param mass_dict
      A dictionary from atom types to atom masses in Lammps unit.

    This method creates a @LammpsDataLines instance
    for *Masses* section and assigns it to #_mass_lines.
    """

    mass_dict = dict(sorted(mass_dict.items()))

    self._mass_lines = LammpsDataLines(
      "Masses", "{type:4d} {mass:10.6f}").set_data({
        "type": mass_dict.keys(), "mass": mass_dict.values()
      })


  def set_topology_types(self, **kwargs):
    """
    @param **kwargs
      A variable number of named arguments.
      Keywords are names of topology component (`bond`, `angle`,
      `dihedral` and `improper`), and arguments are dictionary
      mapping sequences of atom types to corresponding types of
      the topology component.
    """

    if self._atom_types is None:
      RuntimeError("Please set atom's type in advance")
    for k, v in kwargs.items():
      self._topo[k].set_types(v, self._atom_types)


  def set_bond_types(self, seq_to_type):
    """
    @param seq_to_type
      Dictionary mapping two-element tuples of atom types
      to corresponding types of bonds.

    *Example*

    ```py
    writer.set_bond_types({
      (1, 1): 1,
      (1, 2): 2,
      (2, 2): 3
    })
    ```

    Bond type for `(2, 1)` is always assumed to be the same as `(1, 2)`.
    """

    self.set_topology_types(bond=seq_to_type)


  def set_angle_types(self, seq_to_type):
    """
    @param seq_to_type
      Dictionary mapping three-element tuples of atom types
      to corresponding types of angles.
      Note that second atom of each tuple is assumed to be
      at the center of angle.

    *Example*

    ```py
    writer.set_angle_types({
      (1, 1, 1): 1,
      (1, 1, 2): 1,
      (1, 2, 1): 2,
      (1, 2, 2): 2,
      (2, 1, 2): 1,
      (2, 2, 2): 2,
    })
    ```

    Angle types for `(2, 1, 1)` and `(2, 2, 1) are always assumed
    to be the same as `(1, 1, 2)` and `(1, 2, 2), respectively.
    """

    self.set_topology_types(angle=seq_to_type)


  def set_dihedral_types(self, seq_to_type):
    """
    @param seq_to_type
      Dictionary mapping four-element tuples of atom types
      to corresponding types of dihedrals.
      Note that the four atoms are assumed to be connected linearly
      by three bonds in that order.
    """

    self.set_topology_types(dihedral=seq_to_type)


  def set_improper_types(self, seq_to_type):
    """
    @param seq_to_type
      Dictionary mapping four-element tuples of atom types
      to corresponding types of impropers.
      Note that, if *class2* forcefield is used, the second atom of
      each tuple is assumed to be at the center of improper;
      the first atom is assumed to be at the center
      for other forcefields.
    """

    self.set_topology_types(improper=seq_to_type)


  def write_lammps_data(
    self, path, mass=False, centering=False, **kwargs):
    """
    @param path
      File path to *Lammps data file*.
    @param mass
      Whether to write *Masses* section.
    @param centering
      Whether to shift a simulation box
      so that its center position becomes `(0, 0, 0)`.
    @param **kwargs
      A variable number of named arguments.
      * `velocity`: Whether to write *Velocities* section or not.
      * `num_atom_type`: Number of atom types overwriting the number of
      atom types in #_lmp_atoms. It is useful when new atoms will be
      added to the simulation later.
      * `num_{name}_type`: Number of types of topology component
      (`name` = `bond`, `angle`, `dihedral` or `improper`) overwriting
      the number of types of the component in #_topo. It is useful
      when new atoms will be added to the simulation later.

    This method writes *Lammps data file*.
    """

    num_atom = self._lmp_atoms.get_num()

    num_atom_type = (
      kwargs["num_atom_type"]
      if "num_atom_type" in kwargs else self._lmp_atoms.get_num_type())

    num_topo = ({
        k: v.get_num() if v.get_num_type() else 0
        for k, v in self._topo.items()
      } if self._topo else dict.fromkeys(topo_keys, 0))

    num_topo_type = ({
        k: v.get_num_type() for k, v in self._topo.items()
      } if self._topo else dict.fromkeys(topo_keys, 0))

    for k in topo_keys:
      name = "num_{}_type".format(k)
      if name in kwargs:
        num_topo_type[k] = kwargs[name]

    with open(path, "w") as f:

      f.write("# written by ease4lmp.LammpsWriter at {}\n".format(
        dt.datetime.now().strftime("%Y/%m/%d %H:%M:%S")))

      f.write("\n{} atoms\n".format(num_atom))
      f.write("".join([
          "{} {}s\n".format(num_topo[k], k)
          for k in topo_keys
        ]))

      f.write("\n{} atom types\n".format(num_atom_type))
      f.write("".join([
          "{} {} types\n".format(num_topo_type[k], k)
          for k in topo_keys
        ]))

      xlo, ylo, zlo, = 0.0, 0.0, 0.0
      xhi, yhi, zhi, xy, xz, yz = self._prism.get_lammps_prism()

      if centering:

        xhi *= 0.5
        yhi *= 0.5
        zhi *= 0.5
        xlo = -xhi
        ylo = -yhi
        zlo = -zhi

        self._lmp_atoms.shift_positions((xlo, ylo, zlo))

      f.write("\n{}".format("".join([
          "{0:.8e} {1:.8e}  {2}lo {2}hi\n".format(lo, hi, x)
          for lo, hi, x in zip(
            [xlo, ylo, zlo], [xhi, yhi, zhi], ["x", "y", "z"])
        ])))

      if self._prism.is_skewed():
        f.write(
          "\n{:.8e} {:.8e} {:.8e}  xy xz yz\n"
          .format(*map(float, (xy, xz, yz))))

    if mass:
      self._mass_lines.write(path)

    self._lmp_atoms.write_lines(path, **kwargs)

    for v in self._topo.values():
      if 0 < v.get_num() and 0 < v.get_num_type():
        v.write_lines(path)


  def write_lammps_molecule(self, path, special_bonds=True):
    """
    @param path
      File path to *Lammps molecule file*.
    @param special_bonds
      Whether to write *Special Bond Counts*
      and *Special Bonds** section.

    This method writes *Lammps molecule file*.
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
  """
  This class inherits from `ase.calculators.lammpsrun.Prism` class.
  """

  def transform_to_lammps(self, vectors):
    """
    @param vectors
      A `numpy.ndarray` for positions or velocities of atoms.

    @returns
      Transposed `vectors`. Note that transposing is done
      only if the simulation box is skewed.

    This method is required
    because Lammps might use a skewed coordinate system
    whereas Ase uses a cartesian coordinate system.
    """

    if self.is_skewed():
      return np.dot(vectors, self.R.round(int(-self.car_prec.log10())))
    else:
      return vectors