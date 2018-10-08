"""
@file ease4lmp/lammps_data.py
@brief This file is for classes to write each section
in *Lammps data file* (or *molecule file*).
@author Takayuki Kobayashi
@date 2018/05/30
"""

from .lammps_dataformats \
  import lmp_dataformats_atoms, lmp_dataformats_velocities

import itertools as it
import numpy as np


class LammpsAtoms:
  """
  This class is an interface to write *Atoms* (and *Velocities*)
  section in *Lammps data file* (or *molecule file*).
  """

  def __init__(
    self, atom_style, positions, velocities=None, types=None):
    """
    @param atom_style
      A string denoting *Lammps atom style*.
    @param positions
      A `numpy.ndarray` of which shape is (*N*, 3), where *N* is
      the number of atoms. Elements of the array at (*i*, 0), (*i*, 1)
      and (*i*, 2) are *x*, *y* and *z* coordinates of the *i* th atom,
      respectively.
    @param velocities
      A `numpy.ndarray` of which shape is (*N*, 3), where *N* is
      the number of atoms. Elements of the array at (*i*, 0), (*i*, 1)
      and (*i*, 2) are velocities of *i* th atom in the *x*, *y*
      and *z* direction, respectively.
    @param types
      A list of which length is equal to the number of atoms.
      Each element denotes a type of each atom used in Lammps.
    """

    ## Number of atoms.
    self._num = len(positions)

    dataformats = lmp_dataformats_atoms[atom_style]
    datanames = [s.split(":")[0] for s in dataformats]

    ## A LammpsDataLines instance for *Atoms* section.
    self._lines = LammpsDataLines(
      "Atoms", "{{{}}}".format("} {".join(dataformats)))

    ## Dictionary from property names to lists (or other array-like
    # objects) containing property values for atoms.
    self._data = {k: None for k in datanames}

    self.set_data(id=range(1, self._num+1))  # `id` starts from 1 in Lammps
    self.set_data(**dict(zip(["x", "y", "z"], positions.T)))

    if types is not None:

      ## Number of atom types.
      self._num_type = len(set(types))

      self.set_data(type=types)

    if velocities is not None:

      try:
        dataformats_vel = lmp_dataformats_velocities[atom_style]
      except KeyError:
        dataformats_vel = lmp_dataformats_velocities["*"]

      datanames_vel = [s.split(":")[0] for s in dataformats_vel]

      ## A LammpsDataLines instance for *Velocities* section.
      self._lines_vel = LammpsDataLines(
        "Velocities", "{{{}}}".format("} {".join(datanames_vel)))

      ## Dictionary from property names to lists (or other array-like
      # objects) containing property values for velocities.
      self._data_vel = {k: None for k in datanames_vel}

      self._data_vel["id"] = self._data["id"]

      self.set_data(**dict(zip(["vx", "vy", "vz"], velocities.T)))


  def get_num(self):
    """
    @return
      The number of atoms.
    """

    return self._num


  def get_num_type(self):
    """
    @return
      The number of atom types.
    """

    return self._num_type


  def get_required_datanames(self, molecule=False):
    """
    @param molecule
      Whether returned data names are for *Lammps molecule file* or not.
      This method returns data names for *Lammps data file* by default.

    @return
      A set of required data names:
      *Lammps data file* (or *molecule file*) needs those data,
      but this instance does not have them yet.
    """

    names = set([k for k, v in self._data.items() if v is None])

    if molecule:
      names &= {"id", "type", "q", "x", "y", "z"}  # "diameter", "mass" are currently not supported
    elif hasattr(self, "_data_vel"):
      names |= set([k for k, v in self._data_vel.items() if v is None])

    return names


  def set_data(self, **kwargs):
    """
    @param **kwargs
      A variable number of named arguments.
      Keywords are property names, and arguments are property values.

    Properties given as `**kwargs` are stored in #_data.
    """

    if "type" in kwargs:
      self._num_type = len(set(kwargs["type"]))

    for k, v in kwargs.items():
      data = list(v)  # ensure data is `list`
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
    @param shift
      A tuple or list defining a three-dimension vector by which
      atom positions shift. The first, second and third element of
      the vector are for *x*, *y* and *z* direction, respectively.

    This method shifts atom positions stored in #_data.
    """

    for dim, d in zip(["x", "y", "z"], shift):
      self._data[dim] = list(np.array(self._data[dim]) + d)


  def write_lines(self, path, velocity=False, **kwargs):
    """
    @param path
      File path to *Lammps data file*.
    @param velocity
      Whether to write *Velocities* section or not. Default is `False`.
    @param **kwargs
      A variable number of named arguments (for excessive arguments).

    This method calls LammpsDataLines.write of #_lines to write
    *Atoms* section (and #_lines_vel to write *Velocities* section)
    in *Lammps data file*.
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
      self._lines.set_data(self._data).write(path)
      if velocity:
        self._lines_vel.set_data(self._data_vel).write(path)


  def write_lines_for_molecule(self, path):
    """
    @param path
      File path to *Lammps molecule file*.

    This method creates a new LammpsDataLines instance to write
    *Coords* and *Types* sections (and *Charges* section)
    in *Lammps molecule file*.
    """

    section_formats = {
      "Coords": ("id:4d", "x:15.8e", "y:15.8e", "z:15.8e"),
      "Types": ("id:4d", "type:4d")
    }

    if "q" in self._data:
      section_formats["Charges"] = ("id:4d", "q:9.6f")

    for h, dfs in section_formats.items():
      LammpsDataLines(h, "{{{}}}".format("} {".join(dfs))).set_data({
          k: v for k, v in self._data.items()
          if k in [s.split(":")[0] for s in dfs]
        }).write(path)



class LammpsDataLines:
  """
  This class writes a section in *Lammps data file*.

  Instance(s) of this class is created in,
  and used via @LammpsWriter, LammpsAtoms and LammpsTopology class.
  """

  def __init__(self, section_header, line_format):
    """
    @param section_header
      A string denoting a section header.
    @param line_format
      A string representing a format of each line in the section.
      The format should be described in the same manner
      with Python's `format` method.

    Consistency between a given header and format is not checked.
    """

    ## A string denoting a section header.
    self._section_header = section_header

    ## A formattable string used for writing each line.
    self._line_format = line_format + "\n"


  def set_data(self, data):
    """
    @param data
      Data (dictionary from strings for property name
      to lists containing property values) to be written.

    @return
      This instance.
    """

    ## Data to be written in the section.
    self._data = data

    return self


  def write(self, path):
    """
    @param path
      File path to *Lammps data file* (or *molecule file*),
      to which this instance writes a section.
    """

    keys, values = self._data.keys(), self._data.values()
    dicts = [dict(zip(keys, t)) for t in zip(*values)]

    with open(path, 'a') as f:
      f.write("\n{}\n\n".format(self._section_header))
      f.write("".join([
          self._line_format.format(**dic) for dic in dicts
        ]))

    print("'{}' section was written".format(self._section_header))



class LammpsTopology:
  """
  This class is an (abstract) interface to write a section of
  topology components in *Lammps data file* (or *molecule file*).
  """

  @staticmethod
  def create(name, sequences, **kwargs):
    """
    @param name
      Name of topology component:
      `bond`, `angle`, `dihedral` or `improper`.
    @param sequences
      A `numpy.ndarray` representing topology components.
      Shapes of the array for `bond`, `angle`, `dihedral`
      and `improper` are (*N*, 2), (*N*, 3), (*N*, 4) and (*N*, 4),
      respectively: where *N* is the number of topology components.
      Elements of the array are zero-based indices for atoms.
    @param **kwargs
      A variable number of named arguments.
      One can use `class2` keyword for impropers.

    @return
      A created instance of LammpsTopology's subclass.

    This static method is a factory method
    for LammpsTopology's subclass.
    """

    if name == "bond":
      return LammpsBonds(
        sequences,
        ("id", "type", "atom1", "atom2"),
        **kwargs)
    elif name == "angle":
      return LammpsAngles(
        sequences,
        ("id", "type", "atom1", "atom2", "atom3"),
        **kwargs)
    elif name == "dihedral":
      return LammpsDihedrals(
        sequences,
        ("id", "type", "atom1", "atom2", "atom3", "atom4"),
        **kwargs)
    elif name == "improper":
      return LammpsImpropers(
        sequences,
        ("id", "type", "atom1", "atom2", "atom3", "atom4"),
        **kwargs)
    else:
      raise RuntimeError(
        "Invalid name for LammpsTopology '{}'".format(name))


  def __init__(self, sequences, datanames, **kwargs):
    """
    @param sequences
      A `numpy.ndarray` representing topology components.
      Please see also the same name argument of #create.
    @param datanames
      Tuple of strings denoting data names written in a line.
      Please see also a source of #create.
    @param **kwargs
      A variable number of named arguments.
      Please see also the same name argument of #create.
    """

    ## A `numpy.ndarray` for topology components
    # of which elements are `atom-id`s used in Lammps
    self._sequences = sequences + 1

    ## Number of topology components.
    self._num = len(self._sequences)

    ## Tuple of strings denoting data names written in a line.
    self._datanames = datanames

    ## A LammpsDataLines instance to write one of topology section.
    self._lines = LammpsDataLines(
      self.__class__.__name__.replace("Lammps", ""),
      "{{{}}}".format("} {".join(self._datanames)))

    ## Dictionary from property names to lists (or other array-like
    # objects) containing property values for topology components.
    self._data = {k: None for k in datanames}

    if 0 < self._num:
      self._set_data(id=range(1, self._num+1))  # `id` starts from 1 in Lammps
      self._set_data(**dict(zip(
        self._datanames[2:], self._sequences.T)))


  def get_num(self):
    """
    @return
      The number of topology components.
    """

    return self._num


  def get_num_type(self):
    """
    @return
      The number of topology types.
    """

    try:
      return self._num_type
    except AttributeError:
      return 0


  def get_sequence_patterns(self, atom_types):
    """
    @param atom_types
      A one dimensional array-like object of which *i* th element is
      a type of the *i* th atom.

    @return
      A set of unique sequences of atom types
      appearing in the topology components.
      Note that the number of elements in the set is not necessarily
      equal to the number of types of topology component:
      for example, sequence `(1, 2, 3)` and `(3, 2, 1)` are
      reduced to one angle type.
    """

    return set([
        tuple([atom_types[i-1] for i in seq])
        for seq in self._sequences
      ])


  def get_maximum_per_atom(self):
    """
    @return
      The maximum number of topology components per atom.

    See also Lammps source code associated with
    angles/dihedrals/impropers per atom: `Atom::data_angles()`,
    `Atom::data_dihedrals()` and `Atom::data_impropers()`.
    """

    unique, counts = np.unique(
      self._sequences.T[1], return_counts=True)

    return max(counts)


  def set_types(self, seq_to_type, atom_types):
    """
    @param seq_to_type
      A dictionary from sequences of atom types to types of topology
      components consisting of those atoms. Note that one does not
      have to consider sequences in both directions: for example,
      using either `(1, 2, 3)` or `(3, 2, 1)` works fine.
    @param atom_types
      A one dimensional array-like object of which *i* th element is
      a type of the *i* th atom.
    """

    ## Number of types of the topology components.
    self._num_type = len(set(seq_to_type.values()))

    full_dict = self._make_full_seq_to_type(seq_to_type)
    typed_seqs = [
      tuple([atom_types[i-1] for i in seq]) for seq in self._sequences
    ]

    self._set_data(type=[full_dict[seq] for seq in typed_seqs])


  def write_lines(self, path):
    """
    @param path
      File path to *Lammps data file* (or *molecule file*).

    This method calls LammpsDataLines.write of #_lines to write
    a topology section in *Lammps data file* (or *molecule file*).
    """

    self._lines.set_data(self._data).write(path)


  def _make_full_seq_to_type(self, seq_to_type):
    """
    @param seq_to_type
      A dictionary from sequences of atom types to types of topology
      components consisting of those atoms.

    @return
      A dictionary from sequences of atom types to types of topology
      components consisting of those atoms. This dictionary is
      different from the `seq_to_type` by containing sequences
      in both directions as its keys.
    """

    return {t: v for k, v in seq_to_type.items() for t in [k, k[::-1]]}


  def _set_data(self, **kwargs):
    """
    @param **kwargs
      A variable number of named arguments.
      Keywords are property names, and arguments are property values.

    Properties given as `**kwargs` are stored in #_data.
    """

    for k, v in kwargs.items():
      data = list(v)  # ensure data is <list>
      if k in self._data:
        if self._data[k] is None:
          self._data[k] = data
          print(
            "{}: '{}' have been set"
            .format(self.__class__.__name__, k))
        else:
          self._data[k] = data
          print(
            "{}: '{}' have been set again"
          .format(self.__class__.__name__, k))



class LammpsBonds(LammpsTopology):
  """
  This class is an interface to write *Bonds* section
  in *Lammps data file* (or *molecule file*).
  """

  def get_maximum_per_atom(self):
    """
    @return
      The maximum number of topology components per atom.

    This method overrides LammpsTopology.get_maximum_per_atom,
    since data for a bond is linked to the first atom in the bond
    whereas data for a topology component except for bond is linked
    to the second atom in the component. See also Lammps source code
    associated with bonds per atom, `Atom::data_bonds()`.
    """

    unique, counts = np.unique(
      self._sequences.T[0], return_counts=True)

    return max(counts)



class LammpsAngles(LammpsTopology):
  """
  This class is an interface to write *Angles* section
  in *Lammps data file* (or *molecule file*).
  """



class LammpsDihedrals(LammpsTopology):
  """
  This class is an interface to write *Dihedrals* section
  in *Lammps data file* (or *molecule file*).
  """



class LammpsImpropers(LammpsTopology):
  """
  This class is an interface to write *Impropers* section
  in *Lammps data file* (or *molecule file*).
  """

  def __init__(self, sequences, datanames, class2=False, **kwargs):
    """
    @param sequences
      A `numpy.ndarray` representing topology components.
      Please see also the same name argument of #create.
    @param datanames
      Tuple of strings denoting data names written in a line.
      Please see also a source of #create.
    @param class2
      Whether forcefiled is *class2* or not.
      If forcefield is *class2*, a center atom of each improper must be
      the second atom in each sequence representing the improper.
      Default is `False`.
    @param **kwargs
      A variable number of named arguments (for excessive arguments).
    """

    if class2:
      sequences[:, 0:2] = sequences[:, 0:2][:, ::-1]

    super().__init__(sequences, datanames, **kwargs)

    ## Whether forcefiled is *class2* or not.
    self._class2 = class2


  def _make_full_seq_to_type(self, seq_to_type):
    """
    @param seq_to_type
      A dictionary from sequences of atom types to types of topology
      components consisting of those atoms.

    @return
      A dictionary from sequences of atom types to types of topology
      components consisting of those atoms. This dictionary is
      different from the `seq_to_type` by containing sequences
      in both directions as its keys.

    This method overrides LammpsTopology._make_full_seq_to_type
    to be compatible with impropers.
    Different from other topology components
    where atoms are connected *linearly*, impropers have a center atom
    to which all the other three atoms connect. For example,
    sequence `(1, 2, 3, 4)`, `(1, 2, 4, 3)`, `(1, 3, 2, 4)`,
    `(1, 3, 4, 2)`, `(1, 4, 2, 3)` and `(1, 4, 3, 2)` are reduced
    to one improper type, where the first atom in the sequences is
    the center atom of improper.
    """

    c = 1 if self._class2 else 0

    return {
        t[:c] + (k[c],) + t[c:]: v
        for k, v in seq_to_type.items()
        for t in it.permutations(k[:c]+k[c+1:])
      }



class LammpsSpecialBonds:
  """
  This class is an interface to write *Special Bond Counts*
  and *Special Bonds* section in *Lammps molecule file*.
  """

  def __init__(self, bonds_per_atom):
    """
    @param bonds_per_atom
      A nested list containing bond data for each atom.
      Each element of the outer list corresponds to each atom.
      The inner list consists of *absolute* indices for atoms bonded
      with the each atom.
    """

    ## Number of atoms.
    self._num = len(bonds_per_atom)

    ## Dictionary from property names to lists (or other array-like
    # objects) containing property values for special bonds.
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
    """
    @return
      The maximum number of special bonds per atom (sum over *1-2*,
      *1-3*, and *1-4* special bonds).
    """

    return max([sum([len(v) for v in d.values()]) for d in self._data])

  def write_lines(self, path):
    """
    @param path
      File path to *Lammps molecule file* to which this instance
      writes *Special Bond Counts* and *Special Bonds* section.
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
