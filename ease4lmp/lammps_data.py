"""
This file is for LammpsData class and its subclasses.

create: 2018/05/30 by Takayuki Kobayashi
"""

from .lammps_dataformats \
  import lmp_dataformats_atoms, lmp_dataformats_velocities

import itertools as it
import numpy as np

datanames_molecule = {"id", "type", "q", "x", "y", "z"}  # "diameter", "mass"

class LammpsAtoms:
  """
  This class corresponds to 'Atoms' (and 'Velocities') section in
  Lammps DATA (or MOLECULE) file.
  """

  def __init__(self, atom_style, positions, velocities=None, types=None):
    """
    [Arguments]
    * atom_style: <str>
    * positions: <numpy.ndarray>
    * velocities: <numpy.ndarray>
    * types: <list>
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
    This method stores data.
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

  def shift_positions(self, shift=(0.0, 0.0, 0.0)):
    """
    This method shifts positions.

    [Arguments]
    * dr: <list/tuple>; shift in x, y, z directions
    """
    for dim, d in zip(["x", "y", "z"], shift):
      self._data[dim] = list(np.array(self._data[dim]) + d)

  def write_lines(self, path, velocity=False, **kwargs):
    """
    This method makes LammpsDataLines instance to write data in Lammps
    DATA format to the specified path.

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

  def write_lines_for_molecule(self, path):
    """
    This method makes LammpsDataLines instance to write data in Lammps
    MOLECULE format to the specified path.

    [Arguments]
    * path: <str>
    """

    list_header = ["Coords", "Types"]
    list_dataformats = [
      ("id:4d", "x:15.8e", "y:15.8e", "z:15.8e"), ("id:4d", "type:4d")]

    if "q" in self._data:  # "diameter", "mass"
      list_header.append("Charges")
      list_header.append(("id:4d", "q:9.6f"))

    for h, dfs in zip(list_header, list_dataformats):
      LammpsDataLines(h, "{{{}}}".format("} {".join(dfs))).write(path, {
        k: v for k, v in self._data.items()
        if k in [s.split(":")[0] for s in dfs]})

class LammpsDataLines:
  """
  This class corresponds to data content in section in Lammps file.
  """

  def __init__(self, section_header, line_format):
    """
    [Arguments]
    * section_header: <str>
    * line_format: <str>
    """
    self._section_header = section_header
    self._line_format = line_format + "\n"

  def set_data(self, data):
    self._data = data
    return self

  def write(self, path, data=None):
    """
    This method writes data to a file in the specified path.

    [Arguments]
    * path: <str>
    * data: <dict>
    """

    if data is None:
      data = self._data

    keys, values = data.keys(), data.values()
    dicts = [dict(zip(keys, t)) for t in zip(*values)]

    with open(path, 'a') as f:
      f.write("\n{}\n\n".format(self._section_header))
      f.write("".join([
        self._line_format.format(**dic) for dic in dicts]))

    print("'{}' section was written".format(self._section_header))

class LammpsTopology:
  """
  This class corresponds to a section of topology ('Bonds', 'Angles',
  'Dihedrals', 'Impropers') in Lammps file.
  """

  @staticmethod
  def create(name, sequences, **kwargs):
    """
    This is a factory method of LammpsTopology.

    [Arguments]
    * name: <str>
    * data: <numpy.ndarray>
    * kwargs: contains ...
      1. 'class2' for improper
    [Return]
    <LammpsBonds/LammpsAngles/LammpsDihedrals/LammpsImpropers>
    """
    if name == "bond":
      return LammpsBonds(
        sequences, ("id", "type", "atom1", "atom2"), **kwargs)
    elif name == "angle":
      return LammpsAngles(
        sequences, ("id", "type", "atom1", "atom2", "atom3"), **kwargs)
    elif name == "dihedral":
      return LammpsDihedrals(
        sequences, ("id", "type", "atom1", "atom2", "atom3", "atom4"), **kwargs)
    elif name == "improper":
      return LammpsImpropers(
        sequences, ("id", "type", "atom1", "atom2", "atom3", "atom4"), **kwargs)
    else:
      raise RuntimeError(
        "Invalid name for LammpsTopology '{}'".format(name))

  def __init__(self, sequences, datanames, **kwargs):
    """
    [Arguments]
    * bonds_per_atom: <list>; not containing image flags
    * datanames: <tuple>
    """

    self._sequences = sequences + 1
    self._num = len(self._sequences)

    self._datanames = datanames
    self._lines = LammpsDataLines(
      self.__class__.__name__.replace("Lammps", ""),
      "{{{}}}".format("} {".join(self._datanames)))

    self._data = {k: None for k in datanames}

    if 0 < self._num:
      self._set_data(id=range(1, self._num+1))
      self._set_data(**dict(zip(
        self._datanames[2:], self._sequences.T)))

  def get_num(self):
    return self._num

  def get_num_type(self):
    try:
      return self._num_type
    except AttributeError:
      return 0

  def get_sequence_patterns(self, atom_types):
    """
    This method returns a set of sequences (tuple) of atom types.

    [Arguments]
    * atom_types: <list>
    """
    return set([
      tuple([atom_types[i-1] for i in seq]) for seq in self._sequences])

  def get_maximum_per_atom(self):
    """
    This method returns the maximum number of topologies per atom.
    """
    # see Lammps source code about angles/dihedrals/impopers per atom;
    # Atom::data_angles(), Atom::data_dihedrals(), Atom::data_impropers()
    unique, counts = np.unique(self._sequences.T[1], return_counts=True)
    return max(counts)

  def set_types(self, seq_to_type, atom_types):
    """
    This method set types to topologies.

    [Arguments]
    * seq_to_type: <dict>; (<tuple> => <int>)
    * atom_types: <list>
    """

    self._num_type = len(set(seq_to_type.values()))

    full_dict = self._make_full_seq_to_type(seq_to_type)
    typed_seqs = [
      tuple([atom_types[i-1] for i in seq]) for seq in self._sequences]

    self._set_data(type=[full_dict[seq] for seq in typed_seqs])

  def write_lines(self, path):
    """
    This method makes LammpsDataLines instance to write data to the
    specified Lammps file.

    [Arguments]
    * path: <str>
    """
    self._lines.write(path, self._data)

  def _make_full_seq_to_type(self, seq_to_type):
    """
    [Arguments]
    * seq_to_type: <dict>; (<tuple> => <int>)
    [Return]
    <dict>; (<tuple> => <int>)
    """
    return {t: v for k, v in seq_to_type.items() for t in [k, k[::-1]]}

  def _set_data(self, **kwargs):
    """
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
  This class corresponds to 'Bonds' section in Lammps file.
  """

  def get_maximum_per_atom(self):
    # see Lammps source code about bonds per atom; Atom::data_bonds()
    unique, counts = np.unique(self._sequences.T[0], return_counts=True)
    return max(counts)

class LammpsAngles(LammpsTopology):
  """
  This class corresponds to 'Angles' section in Lammps file.
  """

class LammpsDihedrals(LammpsTopology):
  """
  This class corresponds to 'Dihedrals' section in Lammps file.
  """

class LammpsImpropers(LammpsTopology):
  """
  This class corresponds to 'Impropers' section in Lammps file.
  """

  def __init__(self, sequences, datanames, class2=False, **kwargs):
    """
    LammpsImpropers takes arguments of LammpsTopology and an argument
    to determine whether forcefiled is class2 or not.
    """

    if class2:
      sequences[:, 0:2] = sequences[:, 0:2][:, ::-1]

    super().__init__(sequences, datanames, **kwargs)

    self._class2 = class2

  def _make_full_seq_to_type(self, seq_to_type):

    c = 1 if self._class2 else 0

    return {
      t[:c] + (k[c],) + t[c:]: v for k, v in seq_to_type.items()
      for t in it.permutations(k[:c]+k[c+1:])}

class LammpsSpecialBonds:
  """
  This class corresponds to 'Special Bond Counts' and 'Special Bonds'
  sections in Lammps MOLECULE file.
  """

  def __init__(self, bonds_per_atom):
    """
    [Arguments]
    * bonds_per_atom: <list>; not containing image flags
    """

    self._num = len(bonds_per_atom)
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
    return max([sum([len(v) for v in d.values()]) for d in self._data])

  def write_lines(self, path):
    """
    This method writes data to Lammps MOLECULE file.

    [Arguments]
    * path: <str>
    """

    with open(path, 'a') as f:
      f.write("\nSpecial Bond Counts\n\n")
      for i in range(self._num):
        data = self._data[i]
        f.write("{} {} {} {}\n".format(
          i+1, len(data["1-2"]), len(data["1-3"]), len(data["1-4"])))
      f.write("\nSpecial Bonds\n\n")
      for i in range(self._num):
        data = self._data[i]
        f.write("{} {}\n".format(
          i+1, " ".join(map(str, data["1-2"]+data["1-3"]+data["1-4"]))))

    print("'Special Bond Counts' and 'Special Bonds' section was written")
