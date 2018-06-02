"""
This file is for LammpsData class and its subclasses.

create: 2018/05/30 by Takayuki Kobayashi
"""

from .lammps_dataformats \
  import lmp_dataformats_atoms, lmp_dataformats_velocities

import itertools as it
import numpy as np

class LammpsAtoms:
  """
  This class ...
  """

  def __init__(self, atom_style, positions, types, velocities):
    """
    This constructor ...
    [Arguments]
    * atom_style: <str>
    * lammps_unit: <str>
    * positions: <numpy.ndarray>
    * types: <list>
    * velocities: <numpy.ndarray>
    """

    self._num = len(positions)

    #-------------------------------#
    # preparation for Atoms section #
    #-------------------------------#

    datanames_a = lmp_dataformats_atoms[atom_style]
    self._required_datanames_a = set(datanames_a)
    self._lines_a = LammpsDataLines(
      "Atoms", "{{{}}}".format("} {".join(datanames_a)))

    self._set_data_atom(id=range(1, self._num+1))
    self._set_data_atom(**dict(zip(["x","y","z"], positions.T)))

    if types is not None:
      self._num_type = len(set(types))
      self._set_data_atom(type=types)

    #------------------------------------#
    # preparation for Velocities section #
    #------------------------------------#

    if velocities is not None:

      try:
        datanames_v = lmp_dataformats_velocities[atom_style]
      except KeyError:
        datanames_v = lmp_dataformats_velocities["*"]

      self._required_datanames_v = set(datanames_v)
      self._lines_v = LammpsDataLines(
        "Velocities", "{{{}}}".format("} {".join(datanames_v)))

      self._set_data_velocity(id=range(1, self._num+1))
      self._set_data_velocity(**dict(zip(["vx","vy","vz"], velocities.T)))

  def get_num(self):
    return self._num

  def get_num_type(self):
    return self._num_type

  def get_required_datanames(self):
    if hasattr(self, "_required_datanames_v"):
      return self._required_datanames_a + self._required_datanames_v
    else:
      return self._required_datanames_a

  def set_data(self, **kwargs):
    """
    This method ...
    [Arguments]
    * kwargs:
      * (key) name of atom-data
      * (val) <list/tuple/numpy.ndarray>
    """

    if "type" in kwargs.keys():
      self._num_type = len(set(kwargs["type"]))

    self._set_data_atom(**{
      k: v for k, v in kwargs.items()
      if k in self._required_datanames_a})

    if hasattr(self, "_required_datanames_v"):
      self._set_data_velocity(**{
        k: v for k, v in kwargs.items()
        if k in self._required_datanames_v})

  def write_lines(self, path, velocity=False, **kwargs):
    """
    This method ...
    [Arguments]
    * path: <str>
    * velocity: <bool>
    """
    if self._required_datanames_a:
      raise RuntimeError("You have not set '{}'".format(
        "', '".join(self._required_datanames_a)))
    else:
      if not hasattr(self, "_lines_velocity") or not velocity:
        self._lines_a.write(path)
      else:
        if self._required_datanames_v:
          raise RuntimeError("You have not set '{}'".format(
            "', '".join(self._required_datanames_v)))
        else:
          self._lines_a.write(path)
          self._lines_v.write(path)

  def _set_data_atom(self, **kwargs):
    self._lines_a.set_data(**kwargs)
    self._required_datanames_a -= kwargs.keys()

  def _set_data_velocity(self, **kwargs):
    self._lines_v.set_data(**kwargs)
    self._required_datanames_v -= kwargs.keys()

class LammpsSequences:
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
        "Invalid name for LammpsSequences '{}'".format(name))

  def __init__(self, bonds_per_atom, datanames, **kwargs):
    """
    This constructor ...
    [Arguments]
    * bonds_per_atom: <list>; not containing image flags
    * datanames: <tuple>
    """
    self._datanames = datanames
    self._sequences = self._derive_sequences(bonds_per_atom)
    self._num = len(self._sequences)
    self._lines = LammpsDataLines(
      self.__class__.__name__.replace("Lammps", ""),
      "{{{}}}".format("} {".join(self._datanames)))
    self._lines.set_data(id=range(1, self._num+1))
    self._lines.set_data(**dict(zip(
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
    self._lines.set_data(type=[full_dict[seq] for seq in atype_seqs])

  def write_lines(self, path):
    """
    This method ...
    [Arguments]
    * path: <str>
    """
    self._lines.write(path)

  def _derive_sequences(self, bonds_per_atom):
    """
    This method ...
    [Arguments]
    * bonds_per_atom: <list>; not containing image flags
    [Return]
    <list>; <tuple>s representing sequences
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

class LammpsBonds(LammpsSequences):
  """
  This class ...
  """

  def _derive_sequences(self, bonds_per_atom):

    bonds = set()

    for i, bs in enumerate(bonds_per_atom):
      bonds.update(
        set([(i+1, j+1) for j in bs if (j+1, i+1) not in bonds]))

    return list(bonds)

class LammpsAngles(LammpsSequences):
  """
  This class ...
  """

  def _derive_sequences(self, bonds_per_atom):

    angles = set()

    for j, bs in enumerate(bonds_per_atom):
      angles.update(
        set([(i+1, j+1, k+1) for i, k in it.combinations(bs, 2)]))

    return list(angles)

class LammpsDihedrals(LammpsSequences):
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

        dihedrals.update(set([
          (i+1, j+1, k+1, l+1)
          for i, l in it.product(set(bs)-{k}, set(bonds_per_atom[k])-{j})]))

    return list(dihedrals)

class LammpsImpropers(LammpsSequences):
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

    for i, bs in enumerate(bonds_per_atom):
      impropers.update(set([
        (j+1, i+1, k+1, l+1) if self._class2 else (i+1, j+1, k+1, l+1)
        for j, k, l in it.combinations(bs, 3)]))

    return list(impropers)

  def _make_full_seq_to_type(self, seq_to_type):

    new_dict = {}

    c = 1 if self._class2 else 0

    for k, v in seq_to_type.items():
      others = k[:c] + k[c+1:]
      new_dict.update({tuple(l): v for l in [
        [t[i] if i < c else t[i-1] if c < i else k[c] for i in range(4)]
        for t in it.permutations(others)]})

    return new_dict

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
    self._line_format = line_format

  def set_data(self, **kwargs):
    """
    This method ...
    [Arguments]
    * kwargs:
      * (key) name of data
      * (val) <list/tuple/numpy.ndarray>
    """
    for k, vs in kwargs.items():
      if not hasattr(self, "_line_dicts"):
        self._line_dicts = [{} for i in range(len(vs))]
      elif len(vs) != len(self._line_dicts):
        raise RuntimeError(
          "Length of '{}' is different from the other(s)".format(k))
      for dic, v in zip(self._line_dicts, vs):
        dic[k] = v

  def write(self, path):
    """
    This method ...
    [Arguments]
    * path: <str>
    """

    with open(path, 'a') as f:
      f.write("\n{}\n\n".format(self._section_header))
      f.write("\n".join([
        self._line_format.format(**dic) for dic in self._line_dicts]))
      f.write("\n")

    print("'{}' section was written".format(self._section_header))
