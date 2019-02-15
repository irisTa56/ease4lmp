"""Submodule for classes writing each section in Lammps' data file."""

from .lammps_dataformats import (
  lmp_dataformats_atoms, lmp_dataformats_velocities)

import itertools
import numpy as np


def create_topology(name, sequences, **kwargs):
  """A factory function for LammpsTopology's subclass.

  Parameters:

  name: str
    Name of topology component:
    'bond', 'angle', 'dihedral' or 'improper'.

  sequences: numpy.ndarray
    Two-dimensional array describing of topology components.
    Shapes of the array for 'bond', 'angle', 'dihedral'
    and 'improper' are (*N*, 2), (*N*, 3), (*N*, 4) and (*N*, 4),
    respectively. Here *N* is the number of topology components.
    Elements of the array are zero-based indices for atoms.

  kwargs:
    A variable number of named arguments.

    * ``class2`` (bool) : Set to True when using CLASS2 forcefield.

  """
  if name == "bond":
    return LammpsBonds(
      sequences, ("id", "type", "atom1", "atom2"))
  elif name == "angle":
    return LammpsAngles(
      sequences, ("id", "type", "atom1", "atom2", "atom3"))
  elif name == "dihedral":
    return LammpsDihedrals(
      sequences, ("id", "type", "atom1", "atom2", "atom3", "atom4"))
  elif name == "improper":
    return LammpsImpropers(
      sequences, ("id", "type", "atom1", "atom2", "atom3", "atom4"), **kwargs)
  else:
    raise RuntimeError(
      "Invalid name for LammpsTopology '{}'".format(name))

def _write_section_lines(path, data, section_header, line_format):
    """Writes lines of a specified section to a specified file.

    Parameters:

    path: str
      File path to Lammps' data file (or molecule file).

    data: dict
      Dictionary from strings for data name to lists
      containing data values, which describes data to be written.

    section_header: str
      Header of a section to be written.

    line_format: str
      A string representing a format of each line in the section.
      Data values are assigned by Python's ``str.format()`` method.

    """
    line_format += "\n"

    keys, values = data.keys(), data.values()
    dicts = [dict(zip(keys, t)) for t in zip(*values)]

    with open(path, 'a') as f:
      f.write("\n{}\n\n".format(section_header))
      f.write("".join(line_format.format(**dic) for dic in dicts))

    print("'{}' section was written".format(section_header))

#=======================================================================

class LammpsAtoms:
  """An interface to write *Atoms* (and *Velocities*) section
  in Lammps' data file (or molecule file)."""

  def __init__(
    self, atom_style, types, positions, velocities=None, masses=None):
    """
    Parameters:

    atom_style: str
      Specify *atom style* of Lammps.

    types: numpy.ndarray
      Returned value of ``BondedAtoms.get_types()``.

    positions: numpy.ndarray
      Returned value of ``ase.Atoms.get_positions()``.

    velocities: None or numpy.ndarray
      Returned value of ``ase.Atoms.get_velocities()``.
      If None, the *Velocities* section will not be written.
      Note that the *i* th atom in ``velocities`` must correspond to
      the *i* th atom in ``types``.

    masses: None or list or tuple or numpy.ndarray
      A list whose length is equal to the number of atoms.
      Each element is a *mass* of each atom used in Lammps.
      Note that the *i* th atom in ``masses`` must correspond to
      the *i* th atom in ``types``.

    """
    # number of atoms.
    self._num = len(types)

    self._dataformats = lmp_dataformats_atoms[atom_style]

    # dictionary mapping data names to data values
    # (list-like objects) for *Atoms* (and *Masses*) section.
    self._data = {
      k: None
      for k in (s.split(":")[0] for s in self._dataformats)
    }

    self.set_data(id=range(1, self._num+1))  # `id` starts from 1 in Lammps
    self.set_data(type=types)
    self.set_data(**dict(zip(["x", "y", "z"], positions.T)))

    if velocities is not None:
      if len(velocities) != self._num:
        raise RuntimeError("Inconsistent length of velocities")

      try:
        self._dataformats_vel = lmp_dataformats_velocities[atom_style]
      except KeyError:
        self._dataformats_vel = lmp_dataformats_velocities["*"]

      # dictionary mapping data names to data values
      # (list-like objects) for *Velocities* section.
      self._data_vel = {
        k: self._data["id"] if k == "id" else None
        for k in (s.split(":")[0] for s in self._dataformats_vel)
      }

      self.set_data(**dict(zip(["vx", "vy", "vz"], velocities.T)))

    if masses is not None:
      self.set_data(mass=masses)

  def get_num(self):
    """Returns the number of atoms."""
    return self._num

  def get_num_type(self):
    """Returns the number of atom types."""
    return len(set(self._data["type"]))

  def get_atom_types(self):
    return self._data["type"]

  def get_required_datanames(self, molecule=False):
    """Returns a set of required data names (keys).

    Names of data are returned if the data is required to
    write Lammps' data file and has not been set yet.

    Parameters:

    molecule: bool
      Whether to return data names required for Lammps' molecule file.
      This method returns data names for Lammps' data file by default.

    """
    names = set(k for k, v in self._data.items() if v is None)

    if molecule:
      names &= {"id", "type", "q", "x", "y", "z"}  # "diameter" is currently not supported
    elif hasattr(self, "_data_vel"):
      names |= set(k for k, v in self._data_vel.items() if v is None)

    return names

  def set_data(self, **kwargs):
    """Stores given data in ``self._data``.

    Parameters:

    kwargs:
      A variable number of named arguments.
      Keywords are data names, and arguments are data values.

    """
    for k, v in kwargs.items():
      data = list(v)  # ensure data is `list`

      if len(data) != self._num:
        raise RuntimeError("Inconsistent length of data")

      if k in self._data:
        end = "" if self._data[k] is None else " again"
        self._data[k] = data
        print("LammpsAtoms: '{}' have been set{}".format(k, end))
      elif hasattr(self, "_data_vel") and k in self._data_vel:
        end = "" if self._data[k] is None else " again"
        self._data_vel[k] = data
        print("LammpsAtoms: '{}' have been set{}".format(k, end))
      elif k == "mass":
        end = "" if k not in self._data else " again"
        self._data[k] = data
        print("LammpsAtoms: '{}' have been set{}".format(k, end))

  def set_masses(self, type2mass):
    """Stores masses in ``self._data``.

    Though the parameter ``type2mass`` is a dictionary mapping
    atom's type to its mass, masses will be stored as a per-atom array.

    Parameters:

    type2mass: dict
      Dictionary mapping atom's type (int) to its mass (float).

    """
    if self._data["type"] is None:
      raise RuntimeError("You need to set types")

    self.set_data(mass=[type2mass[t] for t in self._data["type"]])

  def shift_positions(self, shift=(0.0, 0.0, 0.0)):
    """Shifts atom positions stored in ``self._data``.

    Parameters:

    shift: tuple or list
      A tuple or list defining a Cartesian vector by which
      atom positions shift.

    """
    for dim, d in zip(["x", "y", "z"], shift):
      self._data[dim] = list(np.array(self._data[dim]) + d)

  def write_lines(self, path, velocity=False, mass=False, **kwargs):
    """Writes lines of *Atoms* and *Velocities* section.

    Parameters:

    path: str
      File path to Lammps' data file.

    velocity: bool
      Whether to write *Velocities* section or not.

    mass: bool
      Whether to write *Masses* section or not.

    kwargs:
      A variable number of named arguments (currently not in use).

    """
    nones = [k for k, v in self._data.items() if v is None]

    if velocity:
      if hasattr(self, "_data_vel"):
        nones += [k for k, v in self._data_vel.items() if v is None]
      else:
        raise RuntimeError("You need to set velocities")

    if nones:
      raise RuntimeError(
        "You have not set '{}'".format("', '".join(nones)))
    else:
      if mass:
        if "mass" in self._data:
          massdct = dict(sorted(dict(zip(
            self._data["type"], self._data["mass"])).items()))
          _write_section_lines(
            path, {
            "type": list(massdct.keys()), "mass": list(massdct.values())
            }, "Masses", "{type:4d} {mass:9.6f}")
        else:
          raise RuntimeError("You need to set masses")

      _write_section_lines(
        path, self._data, "Atoms",
        "{{{}}}".format("} {".join(self._dataformats)))

      if velocity:
        _write_section_lines(
          path, self._data_vel, "Velocities",
          "{{{}}}".format("} {".join(self._dataformats_vel)))

  def write_lines_for_molecule(self, path, mass=False, **kwargs):
    """Writes lines of *Coords*, *Types*, *Charges* and *Masses* section
    to Lammps' molecule file.

    This method creates a new LammpsDataLines instance to write
    *Coords* and *Types* sections (and *Charges* section)
    in Lammps' molecule file.

    Parameters:

    path: str
      File path to Lammps' molecule file.

    mass: bool
      Whether to write *Masses* section or not.

    kwargs:
      A variable number of named arguments (currently not in use).

    """
    section_formats = {
      "Coords": ("id:4d", "x:15.8e", "y:15.8e", "z:15.8e"),
      "Types": ("id:4d", "type:4d")
    }

    if "q" in self._data:
      section_formats["Charges"] = ("id:4d", "q:9.6f")

    if mass:
      if "mass" in self._data:
        section_formats["Masses"] = ("id:4d", "q:9.6f")
      else:
        raise RuntimeError("You need to set masses")

    for h, dfs in section_formats.items():
      _write_section_lines(
        path, {
          k: v for k, v in self._data.items()
          if k in (s.split(":")[0] for s in dfs)
        }, h, "{{{}}}".format("} {".join(dfs)))

#=======================================================================

class LammpsTopology:
  """An (abstract) interface to write a section of topology component
  (*Bonds*, *Angles*, *Dihedrals*, *Impropers*)
  in Lammps' data file (or molecule file)."""

  def __init__(self, sequences, datanames):
    """
    Parameters:

    sequences: numpy.ndarray
      Two-dimensional array describing topology components.

    datanames: tuple
      Tuple of strings specifying data names written in a line.

    """
    # `numpy.ndarray` for topology components
    # whose elements are `atom-id`s used in Lammps
    self._sequences = sequences + 1

    # number of topology components.
    self._num = len(self._sequences)

    self._datanames = datanames

    # dictionary mapping data names to data values
    # (list-like objects) for topology components.
    self._data = {k: None for k in datanames}

    if 0 < self._num:
      self._set_data(id=range(1, self._num+1))  # `id` starts from 1 in Lammps
      self._set_data(**dict(zip(self._datanames[2:], self._sequences.T)))  # set `atom-id`s

  def get_num(self):
    """Returns the number of topology components.

    The number of topology components is considered as 0
    until types of the topology components are set.

    """
    return self._num if hasattr(self, "_num_type") else 0

  def get_num_type(self):
    """Returns the number of topology types."""
    try:
      return self._num_type
    except AttributeError:
      return 0

  def get_sequence_patterns(self, atom_types):
    """Returns a set of unique sequences of atom types
    appearing in the topology components.

    Parameters:

    atom_types: list or tuple or numpy.ndarray
      A one dimensional array-like object whose *i* th element is
      a type of the *i* th atom.

    """
    order_matters_set = set(
      tuple(atom_types[i-1] for i in seq) for seq in self._sequences)

    tmp = []

    for seq in order_matters_set:
      if seq not in tmp and seq[::-1] not in tmp:
        tmp.append(seq)

    return set(tmp)

  def get_maximum_per_atom(self):
    """Returns the maximum number of topology components per atom.

    For more details see also Lammps' source code of
    ``Atom::data_angles()``, ``Atom::data_dihedrals()``
    and ``Atom::data_impropers()``.

    """
    _, counts = np.unique(self._sequences.T[1], return_counts=True)

    return max(counts)

  def set_types(self, seq_to_type, atom_types):
    """Sets type of the topology components.

    Parameters:

    seq_to_type: dict
      Mapping from sequences of atom types (tuple)
      to types of topology components consisting of those atoms (int).
      Note that the sequences can be arbitrary order; for example,
      two sequences of angle ``(1, 2, 3)`` and ``(3, 2, 1)``
      have the same meaning.

    """
    # number of types of the topology components.
    self._num_type = len(set(seq_to_type.values()))

    full_dict = self._make_full_seq_to_type(seq_to_type)

    type_seqs = [
      tuple(atom_types[i-1] for i in seq) for seq in self._sequences
    ]

    self._set_data(type=[full_dict[seq] for seq in type_seqs])

  def write_lines(self, path):
    """Writes lines of a topology section
    in Lammps' data file (or molecule file).

    Parameters:

    path: str
      File path to Lammps' data file (or molecule file).

    """
    _write_section_lines(
      path, self._data, self.__class__.__name__.replace("Lammps", ""),
      "{{{}}}".format("} {".join(self._datanames)))

  def _make_full_seq_to_type(self, seq_to_type):
    """Returns *full* mapping from sequences of atom types
    to types of topology component.

    Returned dictionary (mapping) has both original
    and reverse ordered sequence of atom types as its keys.

    Parameters:

    seq_to_type: dict
      Mapping from sequences of atom types to types of topology
      components consisting of those atoms.

    """
    return {t: v for k, v in seq_to_type.items() for t in [k, k[::-1]]}

  def _set_data(self, **kwargs):
    """Stores given data in ``self._data``.

    Parameters:

    kwargs:
      A variable number of named arguments.
      Keywords are data names, and arguments are data values.

    """
    for k, v in kwargs.items():
      data = list(v)  # ensure data is <list>
      if k in self._data:
        end = "" if self._data[k] is None else " again"
        self._data[k] = data
        print(
          "{}: '{}' have been set{}"
          .format(self.__class__.__name__, k, end))

#=======================================================================

class LammpsBonds(LammpsTopology):
  """An interface to write *Bonds* section
  in Lammps' data file (or molecule file)."""

  def get_maximum_per_atom(self):
    """Returns the maximum number of topology components per atom.

    This method overrides ``LammpsTopology.get_maximum_per_atom()``
    because data for a bond is linked to the first atom in the bond
    whereas data for the other topology components is linked to
    the second atom in the component.
    See also Lammps' source code of ``Atom::data_bonds()``.

    """
    _, counts = np.unique(self._sequences.T[0], return_counts=True)

    return max(counts)

#=======================================================================

class LammpsAngles(LammpsTopology):
  """An interface to write *Angles* section
  in Lammps' data file (or molecule file)."""

#=======================================================================

class LammpsDihedrals(LammpsTopology):
  """An interface to write *Dihedrals* section
  in Lammps' data file (or molecule file)."""

#=======================================================================

class LammpsImpropers(LammpsTopology):
  """An interface to write *Impropers* section
  in Lammps' data file (or molecule file)."""

  def __init__(self, sequences, datanames, class2=False, **kwargs):
    """
    Parameters:

    sequences: numpy.ndarray
      Two-dimensional array describing topology components.
      Note that the first atom of each sequence is a center atom
      of improper (all the three bonds connect to the atom).

    datanames: tuple
      Tuple of strings specifying data names written in a line.

    class2: bool
      Whether forcefiled is CLASS2 or not.
      If forcefield is CLASS2, a center atom of each improper
      must be the second atom in a sequence.
      For other forcefields, the first atom is the center one.

    kwargs:
      A variable number of named arguments
      (for receiving excessive arguments).

    """
    if class2:
      sequences[:, 0:2] = sequences[:, 0:2][:, ::-1]

    super().__init__(sequences, datanames)

    self._class2 = class2

  def get_sequence_patterns(self, atom_types):
    """Returns a set of unique sequences of atom types
    appearing in the topology components.

    This method overrides ``LammpsTopology.get_sequence_patterns()``
    where original and reverse ordered sequence are considered as
    equivalent; it is not the case with improper.

    Parameters:

    atom_types: list or tuple or numpy.ndarray
      A one dimensional array-like object whose *i* th element is
      a type of the *i* th atom.

    """
    order_matters_set = set(
      tuple(atom_types[i-1] for i in seq) for seq in self._sequences)

    c = 1 if self._class2 else 0

    tmp = []

    for seq in order_matters_set:

      equivalents = [
        p[:c] + (seq[c],) + p[c:]
        for p in itertools.permutations(seq[:c]+seq[c+1:])
      ]

      if all(e not in tmp for e in equivalents):
        tmp.append(seq)

    return set(tmp)

  def _make_full_seq_to_type(self, seq_to_type):
    """Returns *full* mapping from sequences of atom types
    to types of topology component.

    Returned dictionary (mapping) has all sequences equivalent
    with the original sequence of atom types as its keys.

    This method overrides ``LammpsTopology._make_full_seq_to_type()``
    to be applicable to improper.
    Different from other topology component where atoms are
    connected *linearly*, improper has a center atom
    to which all the other three atoms connect. For that reason,
    sequence of atom types ``(1, 2, 3, 4)``, ``(1, 2, 4, 3)``,
    ``(1, 3, 2, 4)``, ``(1, 3, 4, 2)``, ``(1, 4, 2, 3)``
    and ``(1, 4, 3, 2)`` are considered as equivalent.
    In the above example, the first atom in the sequences is
    the center atom of improper.

    Parameters:

    seq_to_type: dict
      Mapping from sequences of atom types to types of topology
      components consisting of those atoms.

    """
    c = 1 if self._class2 else 0

    return {
      t[:c] + (k[c],) + t[c:]: v
      for k, v in seq_to_type.items()
      for t in itertools.permutations(k[:c]+k[c+1:])
    }

#=======================================================================

class LammpsSpecialBonds:
  """An interface to write *Special Bond Counts*
  and *Special Bonds* section in Lammps' molecule file."""

  def __init__(self, bonds_per_atom):
    """
    Parameters:

    bonds_per_atom: list
      Two-dimensional list containing bond data for each atom:
      each element of the first list corresponds to each atom,
      and the second list consists of *absolute* indices for atoms
      bonded with the each atom
      (Returned value of ``BondedAtoms.get_bonds_per_atom()``).

    """
    # number of atoms.
    self._num = len(bonds_per_atom)

    # dictionary mapping data names to data values
    # (list-like objects) for special bonds.
    self._data = [
      {"1-2": [], "1-3": [], "1-4": []} for i in range(self._num)]

    for i, bs in enumerate(bonds_per_atom):
      for j in bs:
        self._data[i]["1-2"].append(j+1)
        for k in bonds_per_atom[j]:
          if k == i:
            continue
          self._data[i]["1-3"].append(k+1)
          for l in bonds_per_atom[k]:
            if l == j or l == i:
              continue
            self._data[i]["1-4"].append(l+1)


  def get_maximum_per_atom(self):
    """Returns the maximum number of special bonds per atom
    (sum over *1-2*, *1-3*, and *1-4* special bonds)."""
    return max(sum(len(v) for v in d.values()) for d in self._data)

  def write_lines(self, path):
    """Writes lines of *Special Bond Counts*
    and *Special Bonds* section in Lammps' molecule file.

    Parameters:

    path: str
      File path to Lammps' molecule file.

    """
    with open(path, 'a') as f:

      f.write("\nSpecial Bond Counts\n\n")

      for i in range(self._num):
        data = self._data[i]
        f.write(
          "{} {} {} {}\n".format(
            i+1,
            len(data["1-2"]), len(data["1-3"]), len(data["1-4"])))

      f.write("\nSpecial Bonds\n\n")

      for i in range(self._num):
        data = self._data[i]
        f.write(
          "{} {}\n".format(
            i+1,
            " ".join(map(str, data["1-2"]+data["1-3"]+data["1-4"]))))

    print(
      "'Special Bond Counts' and 'Special Bonds' section was written")
