"""
This file is for LammpsData class and its subclasses.

create: 2018/05/30 by Takayuki Kobayashi
"""

from .lammps_dataformats \
  import lmp_dataformats_atoms, lmp_dataformats_velocities

import itertools as it
import numpy as np

datanames_molecule = {
  "id", "type", "q", "x", "y", "z", "diameter", "mass"}

class LammpsAtoms:
  """
  This class ...
  """

  def __init__(self, atom_style, positions, types=None, velocities=None):
    """
    This constructor ...
    [Arguments]
    * atom_style: <str>
    * positions: <numpy.ndarray>
    * types: <list>
    * velocities: <numpy.ndarray>
    """

    self._num = len(positions)

    #-------------------------------#
    # preparation for Atoms section #
    #-------------------------------#

    dataformats = lmp_dataformats_atoms[atom_style]
    datanames = [s.split(":")[0] for s in dataformats]
    self._lines = LammpsDataLines(
      "Atoms", "{{{}}}".format("} {".join(dataformats)))

    self._data = {k: None for k in datanames}

    self.set_data(id=range(1, self._num+1))
    self.set_data(**dict(zip(["x", "y", "z"], positions.T)))

    if types is not None:
      self._num_type = len(set(types))
      self.set_data(type=types)

    #------------------------------------#
    # preparation for Velocities section #
    #------------------------------------#

    if velocities is not None:

      try:
        dataformats_vel = lmp_dataformats_velocities[atom_style]
      except KeyError:
        dataformats_vel = lmp_dataformats_velocities["*"]

      datanames_vel = [s.split(":")[0] for s in dataformats_vel]
      self._lines_vel = LammpsDataLines(
        "Velocities", "{{{}}}".format("} {".join(datanames_vel)))

      self._data_vel = {k: None for k in datanames_vel}
      self._data_vel["id"] = self._data["id"]

      self.set_data(**dict(zip(["vx", "vy", "vz"], velocities.T)))

  def get_num(self):
    return self._num

  def get_num_type(self):
    return self._num_type

  def get_required_datanames(self, molecule=False):

    names = set([k for k, v in self._data.items() if v is None])

    if molecule:
      names &= datanames_molecule
    elif hasattr(self, "_data_vel"):
      names |= set([k for k, v in self._data_vel.items() if v is None])

    return names

  def set_data(self, **kwargs):
    """
    This method ...
    """

    if "type" in kwargs:
      self._num_type = len(set(kwargs["type"]))

    for k, v in kwargs.items():
      data = list(v)  # ensure data is <list>
      if k in self._data:
        if self._data[k] is None:
          self._data[k] = data
          print("LammpsAtoms: '{}' have been set".format(k))
        else:
          self._data[k] = data
          print("LammpsAtoms: '{}' have been set again".format(k))
      elif hasattr(self, "_data_vel") and k in self._data_vel:
        if self._data_vel[k] is None:
          self._data_vel[k] = data
          print("LammpsAtoms: '{}' have been set".format(k))
        else:
          self._data_vel[k] = data
          print("LammpsAtoms: '{}' have been set again".format(k))

  def write_lines(self, path, velocity=False, **kwargs):
    """
    This method ...
    [Arguments]
    * path: <str>
    * velocity: <bool>
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
      self._lines.write(path, self._data)
      if velocity:
        self._lines_vel.write(path, self._data_vel)

class LammpsDataLines:
  """
  This class ...
  """

  def __init__(self, section_header, line_format):
    """
    This constructor ...
    [Arguments]
    * section_header: <str>
    * line_format: <str>
    """
    self._section_header = section_header
    self._line_format = line_format + "\n"

  def write(self, path, data):
    """
    This method ...
    [Arguments]
    * path: <str>
    * data: <dict>
    """

    keys = data.keys()
    dicts = [dict(zip(keys, t)) for t in zip(*data.values())]

    with open(path, 'a') as f:
      f.write("\n{}\n\n".format(self._section_header))
      f.write("".join([
        self._line_format.format(**dic) for dic in dicts]))

    print("'{}' section was written".format(self._section_header))

class LammpsTopology:
  """
  This class ...
  """

  @staticmethod
  def create(name, bonds, **kwargs):
    """
    This method ...
    [Arguments]
    * name: <str>
    * bonds: <list>; not containing image flags
    * kwargs: contains ...
      1. 'class2' for improper
    [Return]
    <LammpsBonds/LammpsAngles/LammpsDihedrals/LammpsImpropers>
    """
    if name == "bond":
      return LammpsBonds(
        bonds, ("id", "type", "atom1", "atom2"), **kwargs)
    elif name == "angle":
      return LammpsAngles(
        bonds, ("id", "type", "atom1", "atom2", "atom3"), **kwargs)
    elif name == "dihedral":
      return LammpsDihedrals(
        bonds, ("id", "type", "atom1", "atom2", "atom3", "atom4"), **kwargs)
    elif name == "improper":
      return LammpsImpropers(
        bonds, ("id", "type", "atom1", "atom2", "atom3", "atom4"), **kwargs)
    else:
      raise RuntimeError(
        "Invalid name for LammpsTopology '{}'".format(name))

  def __init__(self, bonds_per_atom, datanames, **kwargs):
    """
    This constructor ...
    [Arguments]
    * bonds_per_atom: <list>; not containing image flags
    * datanames: <tuple>
    """

    self._sequences = self._derive_sequences(bonds_per_atom)
    self._num = len(self._sequences)

    self._datanames = datanames
    self._lines = LammpsDataLines(
      self.__class__.__name__.replace("Lammps", ""),
      "{{{}}}".format("} {".join(self._datanames)))

    self._data = {k: None for k in datanames}

    if 0 < self._num:
      self._set_data(id=range(1, self._num+1))
      self._set_data(**dict(zip(
        self._datanames[2:], np.array(self._sequences, int).T)))

  def get_num(self):
    return self._num

  def get_num_type(self):
    try:
      return self._num_type
    except AttributeError:
      return None

  def set_types(self, seq_to_type, atom_types):
    """
    This method ...
    [Arguments]
    * seq_to_type: <dict>; (<tuple> => <int>)
    * atom_types: <list>
    """

    self._num_type = len(set(seq_to_type.values()))

    full_dict = self._make_full_seq_to_type(seq_to_type)

    atype_seqs = [
      tuple([atom_types[i-1] for i in seq]) for seq in self._sequences]

    self._set_data(type=[full_dict[seq] for seq in atype_seqs])

  def write_lines(self, path):
    """
    This method ...
    [Arguments]
    * path: <str>
    """
    self._lines.write(path, self._data)

  def _derive_sequences(self, bonds_per_atom):
    """
    This method ...
    [Arguments]
    * bonds_per_atom: <list>; not containing image flags
    [Return]
    <list>; <tuple>s representing topology
    """
    pass

  def _make_full_seq_to_type(self, seq_to_type):
    """
    This method ...
    [Arguments]
    * seq_to_type: <dict>; (<tuple> => <int>)
    [Return]
    <dict>; (<tuple> => <int>)
    """
    return {t: v for k, v in seq_to_type.items() for t in [k, k[::-1]]}

  def _set_data(self, **kwargs):
    """
    This method ...
    """
    for k, v in kwargs.items():
      data = list(v)  # ensure data is <list>
      if k in self._data:
        if self._data[k] is None:
          self._data[k] = data
          print("{}: '{}' have been set".format(
            self.__class__.__name__, k))
        else:
          self._data[k] = data
          print("{}: '{}' have been set again".format(
            self.__class__.__name__, k))

class LammpsBonds(LammpsTopology):
  """
  This class ...
  """

  def _derive_sequences(self, bonds_per_atom):

    bonds = set()

    for i, bs in enumerate(bonds_per_atom):
      bonds |= set([(i+1, j+1) for j in bs if (j+1, i+1) not in bonds])

    return list(bonds)

class LammpsAngles(LammpsTopology):
  """
  This class ...
  """

  def _derive_sequences(self, bonds_per_atom):

    angles = set()

    for j, bs in enumerate(bonds_per_atom):
      angles |= set([(i+1, j+1, k+1) for i, k in it.combinations(bs, 2)])

    return list(angles)

class LammpsDihedrals(LammpsTopology):
  """
  This class ...
  """

  def _derive_sequences(self, bonds_per_atom):

    bonds = set()
    dihedrals = set()

    for j, bs in enumerate(bonds_per_atom):

      for k in bs:

        if (k+1, j+1) in bonds:
          continue
        else:
          bonds.add((j+1, k+1))

        dihedrals |= set([
          (i+1, j+1, k+1, l+1)
          for i, l in it.product(set(bs)-{k}, set(bonds_per_atom[k])-{j})])

    return list(dihedrals)

class LammpsImpropers(LammpsTopology):
  """
  This class ...
  """

  def __init__(self, *args, class2=False, **kwargs):
    """
    This constructor ...
    """
    self._class2 = class2
    super().__init__(*args, **kwargs)

  def _derive_sequences(self, bonds_per_atom):

    impropers = set()

    c = 1 if self._class2 else 0

    for i, bs in enumerate(bonds_per_atom):
      impropers |= set([
        t[:c] + (i+1,) + t[c:]
        for t in it.combinations([n+1 for n in bs], 3)])

    return list(impropers)

  def _make_full_seq_to_type(self, seq_to_type):

    new_dict = {}

    c = 1 if self._class2 else 0

    for k, v in seq_to_type.items():
      new_dict.update({l: v for l in [
        t[:c] + (k[c],) + t[c:] for t in it.permutations(k[:c]+k[c+1:])]})

    return new_dict