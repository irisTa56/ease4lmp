"""Submodule for an interface to write Lammps' data file."""

from .bonded_atoms import BondedAtoms
from .lammps_sectionwriter import (
  create_topology, LammpsAtoms, LammpsSpecialBonds)
from .lammps_units import lmp_units

import ase
from ase.calculators.lammpsrun import Prism

import ase.units as au
import datetime as dt
import numpy as np

#=======================================================================

class LammpsWriter:
  """An interface to write Lammps' data file (or molecule file)."""

  def __init__(self, atoms, atom_style, lammps_unit="real", **kwargs):
    """
    Parameters:

    atoms: instance of BondedAtoms or ase.Atoms
      Atoms to be written to Lammps' data file.

    atom_style: str
      Specifies an *atom style* used in Lammps.

    lammps_unit: str
      Specifies a *unit* used in Lammps.

    kwargs:
      A variable number of named arguments.

      * ``class2`` (bool) : If Lammps' data file written by this instance is
        used with CLASS2 forcefield, ``class2=True`` must be set.

    """
    if type(atoms) == ase.Atoms:
      atoms = BondedAtoms(atoms)

    lmp_unit = lmp_units[lammps_unit]

    dunit = 1 / au.m / lmp_unit["_2si"]["distance"]

    try:
      self._prism = ExtendedPrism(dunit*atoms.get_cell())
    except OverflowError:
      raise RuntimeError("Unit cell size must be set")

    pos = self._prism.transform_to_lammps(dunit*atoms.get_positions())

    vel = atoms.get_velocities()
    if vel is not None:
      vunit = dunit * au.s * lmp_unit["_2si"]["time"]
      vel = self._prism.transform_to_lammps(vunit*vel)

    mas = atoms.get_masses()
    if 0 < mas.size:
      mas *= 1 / au.kg / lmp_unit["_2si"]["mass"]
    else:
      mas = None

    self._lmpatoms = LammpsAtoms(
      atom_style, atoms.get_types(), pos, velocities=vel, masses=mas)

    # dictionary containing LammpsTopology instances
    # for each topology components.
    self._lmptopo_dict = {
      k: create_topology(k, v, **kwargs)
      for k, v in {
        "bond": atoms.get_bonded_bonds(),
        "angle": atoms.get_bonded_angles(),
        "dihedral": atoms.get_bonded_dihedrals(),
        "improper": atoms.get_bonded_impropers(),
      }.items()
    }

    self._lmpsbonds = LammpsSpecialBonds(atoms.get_bonds_per_atom())

  def get_required_datanames(self):
    """Return a set of required data names (keys).

    Names of data are returned if the data is required to
    write Lammps' data file and has not been set yet.

    """
    return self._lmpatoms.get_required_datanames()

  def get_required_datanames_for_molecule(self):
    """Return a set of required data names (keys).

    Names of data are returned if the data is required to
    write Lammps' molecule file and has not been set yet.

    """
    return self._lmpatoms.get_required_datanames(molecule=True)

  def get_sequence_patterns(self, key):
    """Return a set of unique sequences of atom types
    appearing in the given topology component.

    Parameters:

    key: str
      Name of topology component:
      'bond', 'angle', 'dihedral' or 'improper'.

    """
    atom_types = self._lmpatoms.get_atom_types()
    return self._lmptopo_dict[key].get_sequence_patterns(atom_types)

  def get_bond_patterns(self):
    """Return a set of unique sequences of atom types
    appearing in all the sequences describing bonds."""
    return self.get_sequence_patterns("bond")

  def get_angle_patterns(self):
    """Return a set of unique sequences of atom types
    appearing in all the sequences describing angles."""
    return self.get_sequence_patterns("angle")

  def get_dihedral_patterns(self):
    """Return a set of unique sequences of atom types
    appearing in all the sequences describing dihedrals."""
    return self.get_sequence_patterns("dihedral")

  def get_improper_patterns(self):
    """Return a set of unique sequences of atom types
    appearing in all the sequences describing impropers."""
    return self.get_sequence_patterns("improper")

  def print_maximum_per_atom(self, key):
    """Return the maximum number of given topology components per atom.

    Parameters:

    key: str
      Name of topology component:
      'bond', 'angle', 'dihedral' or 'improper'.

    """
    print(
      "You might need to set 'extra/{}/per/atom' to: {}"
      .format(key, self._lmptopo_dict[key].get_maximum_per_atom()))

  def print_max_bonds_per_atom(self):
    """Return the maximum number of bonds per atom."""
    self.print_maximum_per_atom("bond")

  def print_max_angles_per_atom(self):
    """Return the maximum number of angles per atom."""
    self.print_maximum_per_atom("angle")

  def print_max_dihedrals_per_atom(self):
    """Return the maximum number of dihedrals per atom."""
    self.print_maximum_per_atom("dihedral")

  def print_max_impropers_per_atom(self):
    """Return the maximum number of impropers per atom."""
    self.print_maximum_per_atom("improper")

  def print_max_specials_per_atom(self):
    """Return the maximum number of special bonds per atom."""
    print(
      "You might need to set 'extra/special/per/atom' to: {}"
      .format(self._lmpsbonds.get_maximum_per_atom()))

  def set_atom_data(self, **kwargs):
    """Set data for *Atoms* (and *Velocities*) section.

    Parameters:

    kwargs:
      A variable number of named arguments.
      Keywords are data names, and arguments are data values
      (list-like objects).

    """
    self._lmpatoms.set_data(**kwargs)

  def set_masses(self, type2mass):
    """Set data for *Masses* section.

    Parameters:

    type2mass: dict
      Dictionary mapping atom's type (int) to its mass (float).

    """
    self._lmpatoms.set_masses(type2mass)

  def set_topology_types(self, **kwargs):
    """Set types of the given topology components.

    Parameters:

    kwargs:
      A variable number of named arguments.
      Keywords are names of topology component ('bond', 'angle',
      'dihedral', 'improper'), and arguments are dictionary
      mapping sequences of atom types to corresponding types of
      the topology component.

    """
    for k, v in kwargs.items():
      self._lmptopo_dict[k].set_types(v, self._lmpatoms.get_atom_types())

  def set_bond_types(self, seq_to_type):
    """Set types of bonds.

    Parameters:

    seq_to_type: dict form tuple to int
      Mapping two-element tuples of atom types
      to corresponding types of bond.

    Examples:

    >>> from ease4lmp import BondedAtoms, LammpsWriter
    >>> from ase.build import molecule
    >>> methanol = BondedAtoms(molecule("CH3OH"))
    >>> methanol.get_atomic_numbers()
    array([6, 8, 1, 1, 1, 1])
    >>> methanol.get_distance(1, 3)  # O-H bond
    0.9700009076665858
    >>> methanol.set_types([1, 2, 3, 4, 3, 3])
    >>> bond_list = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 3)]
    >>> for t in bond_list:
    ...   methanol.add_bond(*t)
    >>> methanol.set_cell([[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]])
    >>> methanol.center()
    >>> writer = LammpsWriter(methanol, atom_style="molecular", special_bonds=True)
    LammpsAtoms: 'id' have been set
    LammpsAtoms: 'type' have been set
    LammpsAtoms: 'x' have been set
    LammpsAtoms: 'y' have been set
    LammpsAtoms: 'z' have been set
    LammpsAtoms: 'mass' have been set
    LammpsBonds: 'id' have been set
    LammpsBonds: 'atom1' have been set
    LammpsBonds: 'atom2' have been set
    LammpsAngles: 'id' have been set
    LammpsAngles: 'atom1' have been set
    LammpsAngles: 'atom2' have been set
    LammpsAngles: 'atom3' have been set
    LammpsDihedrals: 'id' have been set
    LammpsDihedrals: 'atom1' have been set
    LammpsDihedrals: 'atom2' have been set
    LammpsDihedrals: 'atom3' have been set
    LammpsDihedrals: 'atom4' have been set
    LammpsImpropers: 'id' have been set
    LammpsImpropers: 'atom1' have been set
    LammpsImpropers: 'atom2' have been set
    LammpsImpropers: 'atom3' have been set
    LammpsImpropers: 'atom4' have been set
    >>> writer.get_bond_patterns()
    {(1, 2), (1, 3), (2, 4)}
    >>> writer.set_bond_types({
    ...   (1, 2): 1,
    ...   (1, 3): 2,
    ...   (2, 4): 3,
    ... })
    LammpsBonds: 'type' have been set

    """
    self.set_topology_types(bond=seq_to_type)

  def set_angle_types(self, seq_to_type):
    """Set types of angles.

    Parameters:

    seq_to_type: dict form tuple to int
      Mapping three-element tuples of atom types
      to corresponding types of angle.
      Note that second atom of each tuple should be
      at the center of angle.

    Examples:

    >>> from ease4lmp import BondedAtoms, LammpsWriter
    >>> from ase.build import molecule
    >>> methanol = BondedAtoms(molecule("CH3OH"))
    >>> methanol.get_atomic_numbers()
    array([6, 8, 1, 1, 1, 1])
    >>> methanol.get_distance(1, 3)  # O-H bond
    0.9700009076665858
    >>> methanol.set_types([1, 2, 3, 4, 3, 3])
    >>> bond_list = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 3)]
    >>> for t in bond_list:
    ...   methanol.add_bond(*t)
    >>> methanol.set_cell([[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]])
    >>> methanol.center()
    >>> writer = LammpsWriter(methanol, atom_style="molecular", special_bonds=True)
    LammpsAtoms: 'id' have been set
    LammpsAtoms: 'type' have been set
    LammpsAtoms: 'x' have been set
    LammpsAtoms: 'y' have been set
    LammpsAtoms: 'z' have been set
    LammpsAtoms: 'mass' have been set
    LammpsBonds: 'id' have been set
    LammpsBonds: 'atom1' have been set
    LammpsBonds: 'atom2' have been set
    LammpsAngles: 'id' have been set
    LammpsAngles: 'atom1' have been set
    LammpsAngles: 'atom2' have been set
    LammpsAngles: 'atom3' have been set
    LammpsDihedrals: 'id' have been set
    LammpsDihedrals: 'atom1' have been set
    LammpsDihedrals: 'atom2' have been set
    LammpsDihedrals: 'atom3' have been set
    LammpsDihedrals: 'atom4' have been set
    LammpsImpropers: 'id' have been set
    LammpsImpropers: 'atom1' have been set
    LammpsImpropers: 'atom2' have been set
    LammpsImpropers: 'atom3' have been set
    LammpsImpropers: 'atom4' have been set
    >>> writer.get_angle_patterns()
    {(1, 2, 4), (2, 1, 3), (3, 1, 3)}
    >>> writer.set_angle_types({
    ...   (1, 2, 4): 1,
    ...   (2, 1, 3): 2,
    ...   (3, 1, 3): 3,
    ... })
    LammpsAngles: 'type' have been set

    """
    self.set_topology_types(angle=seq_to_type)

  def set_dihedral_types(self, seq_to_type):
    """Set types of dihedrals.

    Parameters:

    seq_to_type: dict form tuple to int
      Mapping four-element tuples of atom types
      to corresponding types of dihedral.
      Note that the four atoms should be connected linearly
      by three bonds in that order.

    Examples:

    >>> from ease4lmp import BondedAtoms, LammpsWriter
    >>> from ase.build import molecule
    >>> methanol = BondedAtoms(molecule("CH3OH"))
    >>> methanol.get_atomic_numbers()
    array([6, 8, 1, 1, 1, 1])
    >>> methanol.get_distance(1, 3)  # O-H bond
    0.9700009076665858
    >>> methanol.set_types([1, 2, 3, 4, 3, 3])
    >>> bond_list = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 3)]
    >>> for t in bond_list:
    ...   methanol.add_bond(*t)
    >>> methanol.set_cell([[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]])
    >>> methanol.center()
    >>> writer = LammpsWriter(methanol, atom_style="molecular", special_bonds=True)
    LammpsAtoms: 'id' have been set
    LammpsAtoms: 'type' have been set
    LammpsAtoms: 'x' have been set
    LammpsAtoms: 'y' have been set
    LammpsAtoms: 'z' have been set
    LammpsAtoms: 'mass' have been set
    LammpsBonds: 'id' have been set
    LammpsBonds: 'atom1' have been set
    LammpsBonds: 'atom2' have been set
    LammpsAngles: 'id' have been set
    LammpsAngles: 'atom1' have been set
    LammpsAngles: 'atom2' have been set
    LammpsAngles: 'atom3' have been set
    LammpsDihedrals: 'id' have been set
    LammpsDihedrals: 'atom1' have been set
    LammpsDihedrals: 'atom2' have been set
    LammpsDihedrals: 'atom3' have been set
    LammpsDihedrals: 'atom4' have been set
    LammpsImpropers: 'id' have been set
    LammpsImpropers: 'atom1' have been set
    LammpsImpropers: 'atom2' have been set
    LammpsImpropers: 'atom3' have been set
    LammpsImpropers: 'atom4' have been set
    >>> writer.get_dihedral_patterns()
    {(3, 1, 2, 4)}
    >>> writer.set_dihedral_types({
    ...   (3, 1, 2, 4): 1,
    ... })
    LammpsDihedrals: 'type' have been set

    """
    self.set_topology_types(dihedral=seq_to_type)

  def set_improper_types(self, seq_to_type):
    """Set types of impropers.

    Parameters:

    seq_to_type: dict form tuple to int
      Mapping four-element tuples of atom types
      to corresponding types of improper.
      Note that, if CLASS2 forcefield is used, the second atom of
      each tuple should be at the center of improper;
      the first atom should be at the center for other forcefields.

    Examples:

    >>> from ease4lmp import BondedAtoms, LammpsWriter
    >>> from ase.build import molecule
    >>> methanol = BondedAtoms(molecule("CH3OH"))
    >>> methanol.get_atomic_numbers()
    array([6, 8, 1, 1, 1, 1])
    >>> methanol.get_distance(1, 3)  # O-H bond
    0.9700009076665858
    >>> methanol.set_types([1, 2, 3, 4, 3, 3])
    >>> bond_list = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 3)]
    >>> for t in bond_list:
    ...   methanol.add_bond(*t)
    >>> methanol.set_cell([[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]])
    >>> methanol.center()
    >>> writer = LammpsWriter(methanol, atom_style="molecular", special_bonds=True)
    LammpsAtoms: 'id' have been set
    LammpsAtoms: 'type' have been set
    LammpsAtoms: 'x' have been set
    LammpsAtoms: 'y' have been set
    LammpsAtoms: 'z' have been set
    LammpsAtoms: 'mass' have been set
    LammpsBonds: 'id' have been set
    LammpsBonds: 'atom1' have been set
    LammpsBonds: 'atom2' have been set
    LammpsAngles: 'id' have been set
    LammpsAngles: 'atom1' have been set
    LammpsAngles: 'atom2' have been set
    LammpsAngles: 'atom3' have been set
    LammpsDihedrals: 'id' have been set
    LammpsDihedrals: 'atom1' have been set
    LammpsDihedrals: 'atom2' have been set
    LammpsDihedrals: 'atom3' have been set
    LammpsDihedrals: 'atom4' have been set
    LammpsImpropers: 'id' have been set
    LammpsImpropers: 'atom1' have been set
    LammpsImpropers: 'atom2' have been set
    LammpsImpropers: 'atom3' have been set
    LammpsImpropers: 'atom4' have been set
    >>> writer.get_improper_patterns()
    {(1, 2, 3, 3), (1, 3, 3, 3)}
    >>> writer.set_improper_types({
    ...   (1, 2, 3, 3): 1,
    ...   (1, 3, 3, 3): 2,
    ... })
    LammpsImpropers: 'type' have been set

    """
    self.set_topology_types(improper=seq_to_type)

  def write_lammps_data(
    self, path, centering=False, **kwargs):
    """Write Lammps' data file.

    Parameters:

    path: str
      File path to Lammps' data file.

    centering: bool
      Whether to shift a simulation box
      so that its center position becomes ``(0, 0, 0)``.

    kwargs:
      A variable number of named arguments.

      * ``velocity`` (bool) : Whether to write *Velocities* section or not.
      * ``mass`` (bool) : Whether to write *Masses* section or not.
      * ``num_atom_type`` (bool) : Number of atom types; this overwrites
        the number of atom types stored in this instance. It is useful
        when extra atoms will be added to the simulation box later.
      * ``num_{name}_type`` (bool) : Number of types of topology component
        (``name`` is one of *bond*, *angle*, *dihedral* and *improper*);
        this overwrites the number of the component types
        stored in this instance. It is useful when extra atoms
        will be added to the simulation box later.

    """
    num_atom = self._lmpatoms.get_num()

    num_atom_type = kwargs.get(
      "num_atom_type", self._lmpatoms.get_num_type())

    num_topo = {
      k: v.get_num() for k, v in self._lmptopo_dict.items()
    }

    num_topo_type = {
      k: kwargs.get(kw, v.get_num_type())
      for k, v, kw in (
        (k, v, "num_{}_type".format(k))
        for k, v in self._lmptopo_dict.items())
    }

    with open(path, "w") as f:
      f.write("# written by ease4lmp.LammpsWriter at {}\n".format(
        dt.datetime.now().strftime("%Y/%m/%d %H:%M:%S")))

      f.write("\n{} atoms\n".format(num_atom))
      f.write("".join(
        "{} {}s\n".format(num_topo[k], k)
        for k in self._lmptopo_dict.keys()))

      f.write("\n{} atom types\n".format(num_atom_type))
      f.write("".join(
        "{} {} types\n".format(num_topo_type[k], k)
        for k in self._lmptopo_dict.keys()))

      xlo, ylo, zlo, = 0.0, 0.0, 0.0
      xhi, yhi, zhi, xy, xz, yz = self._prism.get_lammps_prism()

      if centering:
        xhi *= 0.5
        yhi *= 0.5
        zhi *= 0.5
        xlo = -xhi
        ylo = -yhi
        zlo = -zhi
        self._lmpatoms.shift_positions((xlo, ylo, zlo))

      f.write("\n{}".format("".join(
          "{0:.8e} {1:.8e}  {2}lo {2}hi\n".format(lo, hi, x)
          for lo, hi, x in zip(
            [xlo, ylo, zlo], [xhi, yhi, zhi], ["x", "y", "z"]))))

      if self._prism.is_skewed():
        f.write(
          "\n{:.8e} {:.8e} {:.8e}  xy xz yz\n"
          .format(*map(float, (xy, xz, yz))))

    self._lmpatoms.write_lines(path, **kwargs)

    for v in self._lmptopo_dict.values():
      if 0 < v.get_num() and 0 < v.get_num_type():
        v.write_lines(path)

  def write_lammps_molecule(self, path, special_bonds=True, **kwargs):
    """Write Lammps' molecule file.

    Parameters:

    path: str
      File path to Lammps' molecule file.

    special_bonds: bool
      Whether to write *Special Bond Counts*
      and *Special Bonds** section.

    kwargs:
      A variable number of named arguments.

      * ``mass`` (bool) : Whether to write *Masses* section or not.

    """
    num_atom = self._lmpatoms.get_num()

    num_topo = {
      k: v.get_num() for k, v in self._lmptopo_dict.items()
    }

    with open(path, "w") as f:
      f.write("# written by ease4lmp.LammpsWriter at {}\n".format(
        dt.datetime.now().strftime("%Y/%m/%d %H:%M:%S")))

      f.write("\n{} atoms\n".format(num_atom))
      f.write("".join(
        "{} {}s\n".format(num_topo[k], k)
          for k in self._lmptopo_dict.keys()))

    self._lmpatoms.write_lines_for_molecule(path, **kwargs)

    for v in self._lmptopo_dict.values():
      if 0 < v.get_num() and 0 < v.get_num_type():
        v.write_lines(path)

    if special_bonds:
      self._lmpsbonds.write_lines(path)

#=======================================================================

class ExtendedPrism(Prism):
  """Extended ``ase.calculators.lammpsrun.Prism`` class."""

  def transform_to_lammps(self, vectors):
    """Return transposed *vectors*.

    This method is required to convert vectors
    from ASE's cartesian coordinate system
    to Lammps' skewed coordinate system.

    Note that transposing occurs only if the simulation box is skewed.

    Parameters:

    vectors: numpy.ndarray
      Positions or velocities of atoms.

    """
    if self.is_skewed():
      return np.dot(vectors, self.R.round(int(-self.car_prec.log10())))
    else:
      return vectors