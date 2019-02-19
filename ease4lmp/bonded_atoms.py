"""Submodule for a class extending ``ase.Atoms`` with bonds."""

import sys
import itertools

import ase
import numpy as np

from .lammps_reader import (
  read_box, read_bonds, read_atoms_from_data, read_atoms_from_molecule)

def create_atoms_from_data(path, atom_style, pbc=False):
  """Create a BondedAtoms instance from Lammps' data file.

  Parameters:

  path: str
    File path to Lammps' data file.

  atom_style: str
    Specifies an *atom style* used in Lammps.

  pbc: array-like object of bool
    For each axis in the unit cell decides whether the positions
    will be moved along this axis.

  """
  if isinstance(pbc, bool):
    pbc = (pbc,) * 3

  atoms = read_atoms_from_data(path, atom_style)
  bonds = read_bonds(path)

  obj = create_atoms_from_json(atoms, bonds)
  obj.set_cell(read_box(path))
  obj.set_pbc(pbc)
  obj.adjust_pbc_bonds()

  return obj

def create_atoms_from_molecule(path):
  """Create a BondedAtoms instance from Lammps' molecule file.

  Parameters:

  path: str
    File path to Lammps' molecule file.

  """
  atoms = read_atoms_from_molecule(path)
  bonds = read_bonds(path)

  return create_atoms_from_json(atoms, bonds)

def create_atoms_from_json(atoms, bonds=[]):
  """Create a BondedAtoms instance from JSON object(s).

  Parameters:

  atoms: list of dict
    JSON describing atoms.
    All the dictionary contained in the list must have the same keys.
    The following keys are required (some are optional):

    * ``xu`` (float) : In Lammps' unit.
    * ``yu`` (float) : In Lammps' unit.
    * ``zu`` (float) : In Lammps' unit.
    * ``type`` (int) : Optional.
    * ``mass`` (float) : In Lammps' unit. Optional.
    * ``vx`` (float) : In Lammps' unit. Optional.
    * ``vy`` (float) : In Lammps' unit. Optional.
    * ``vz`` (float) : In Lammps' unit. Optional.
    * ``id`` (int) : Required when ``bonds`` is supplied.

  bonds: list of dict
    JSON describing bonds.
    All the dictionary contained in the list must have the same keys.
    The following keys are required:

    * ``atom1-id`` (int) : The first atom id.
    * ``atom2-id`` (int) : The second atom id.

  """
  obj = BondedAtoms(positions=[
      [atom["xu"], atom["yu"], atom["zu"]]
      for atom in atoms if all(k in atom for k in ["xu", "yu", "zu"])
    ])

  types = [atom["type"] for atom in atoms if "type" in atom]
  if len(types) == len(obj):
    obj.set_types(types)

  masses = [atom["mass"] for atom in atoms if "mass" in atom]
  if len(masses) == len(obj):
    obj.set_masses(masses)

  velocities = [
    [atom["vx"], atom["vy"], atom["vz"]]
    for atom in atoms if all(k in atom for k in ["vx", "vy", "vz"])
  ]
  if len(velocities) == len(obj):
    obj.set_velocities(velocities)

  if len(bonds) > 0:

    bondarray = np.array([
        [bond["atom1-id"], bond["atom2-id"]]
        for bond in bonds if "atom1-id" in bond and "atom2-id" in bond
      ])

    _, counts = np.unique(bondarray, return_counts=True)

    max_bonds_per_atom = max(counts)

    if 4 < max_bonds_per_atom:
      obj.change_max_bonds(max_bonds_per_atom)

    id2idx = {
      atom["id"]: i for i, atom in enumerate(atoms) if "id" in atom
    }

    for id1, id2 in bondarray:
      obj.add_bond(id2idx[id1], id2idx[id2])

  return obj

def _compute_image_flags(vec, cell, pbc):
  """Compute image flags for bond data.

  Parameters:

  vec: numpy.ndarray of number
      Vector from an atom owning a bond to the other atom
      connected by the bond.

  cell: numpy.ndarray of number
      Unit cell vectors.

  pbc: array-like object of bool
      For each axis in the unit cell decides whether the positions
      will be moved along this axis.

  """
  coeffs = np.linalg.solve(cell.T, vec) + 0.5
  coeffs[np.logical_not(pbc)] = 0.0
  coeffs //= 1.0

  return -coeffs.astype(int)

#=======================================================================

class BondedAtoms(ase.Atoms):
  """Inheriting ``ase.Atoms`` class, and has functionalities
  to handle (covalent) bonds connecting a pair of atoms.

  An instance of this class has bond data as a ``numpy.ndarray``
  as with other atomic properties like positions and velocities.
  Shape of this array for bond data is (*N*, *M*, 4), where *N* is
  the number of atoms and *M* is the maximum number of bonds per atom.
  The first axis of the array corresponds to atoms, and the second axis
  corresponds to bonds connected to each atom.
  Four integers in the third axis describes each bond:
  the first integer is a relative index for the other atom,
  the second to fourth integers are relative image flags
  (used for resolving a periodic boundary condition)
  for the other atom in the *x*, *y* and *z* direction, respectively.

  In addition, an instance of this class has type of atoms
  as a one-dimensional ``numpy.ndarray``. Length of this array is
  equal to the number of atoms. Atom's type is typically required
  for making Lammps' data file and molecule file.

  Basically, you can use methods of ``ase.Atoms`` in the same manner.

  """

  def __init__(self, *args, **kwargs):
    """
    Parameters:

    args:
      A variable number of arguments passed to ``super().__init__()``.

    kwargs:
      A variable number of named arguments passed to ``super().__init__()``.

    """
    super().__init__(*args, **kwargs)

    # maximum number of bonds per atom.
    self._max_bonds = 4

    if "bonds" not in self.arrays:
      self.new_array(
        "bonds", np.zeros((len(self), self._max_bonds, 4)), int)

    if "types" not in self.arrays:
      self.new_array(
        "types", np.ones(len(self)), int)

  def add_bond(self, atom1, atom2, img1=(0,0,0), img2=(0,0,0)):
    """Add a bond connecting two atoms.

    Parameters:

    atom1: int
      Index for the first atom.

    atom2: int
      Index for the second atom.

    img1: tuple/list of int
      Three-element tuple or list of integers indicating which image of
      a periodic simulation box the first atom is in.

    img2: tuple/list of int
      Three-element tuple or list of integers indicating which image of
      a periodic simulation box the second atom is in.

    .. note::

      * Bond data stored in the first atom:
        ``[atom2-atom1, img2[0]-img1[0], img2[1]-img1[1], img2[2]-img1[2]]``.
      * Bond data stored in the second atom:
        ``[atom1-atom2, img1[0]-img2[0], img1[1]-img2[1], img1[2]-img2[2]]``.

    """
    if not isinstance(img1, (tuple, list)) or len(img1) != 3:
      raise RuntimeError(
        "Image flag 'img1' should be 3-membered tuple/list")
    elif not isinstance(img2, (tuple, list)) or len(img2) != 3:
      raise RuntimeError(
        "Image flag 'img2' should be 3-membered tuple/list")

    relative_data = np.asarray([  # not `np.array` due to pylint alarm
        atom2-atom1, img2[0]-img1[0], img2[1]-img1[1], img2[2]-img1[2]
      ], int)

    idx1 = self._get_available_bond_index(
      atom1, relative_data)
    idx2 = self._get_available_bond_index(
      atom2, -relative_data, 1 if atom1 == atom2 else 0)

    self.arrays["bonds"][atom1][idx1] = relative_data
    self.arrays["bonds"][atom2][idx2] = -relative_data

  def change_max_bonds(self, n=4):
    """Change the maximum number of bonds per atom.

    Like original properties of ``ase.Atoms``,
    bond data is also stored as a ``numpy.ndarray``,
    which is a homogeneous array of fixed-size items.
    This means that all the atoms have the same capacity for bond data,
    and the capacity is fixed until manually changed.
    This method changes the capacity, namely,
    the maximum number of bonds per atom.
    Note that existing data stay unchanged after changing the capacity.

    Parameters:

    n: int
      New maximum number of bonds per atom.

    """
    if n < self._max_bonds:
      self.arrays["bonds"] = np.array(
        [bs[:n] for bs in self.arrays["bonds"]])
    elif self._max_bonds < n:
      zeros = np.zeros((n - self._max_bonds, 4), int)
      self.arrays["bonds"] = np.array(
        [np.append(bs, zeros, axis=0) for bs in self.arrays["bonds"]])

    self._max_bonds = n

  def extend(self, other):
    """Extend this instance by appending the given instance.

    This method returns this extended instance,
    which is the same behavior with ``ase.Atoms.extend()``.

    Parameters:

    other: instance of ase.Atom/ase.Atoms/BondedAtoms
      For ase.Atom or ase.Atoms instance,
      the instance will be down-casted to BondedAtoms before extension.

    """
    if type(other) == ase.Atom:
      other = self.__class__([other])

    if type(other) == ase.Atoms:
      other = self.__class__(other)

    return super().extend(other)

  def get_bonds(self):
    """Return an array of integers containing bond data.

    An instance of ``ase.Atoms`` stores data in ``self.arrays``;
    for example, atoms' positions are accessed by key 'positions'.
    For BondedAtoms instance, bond data is also stored
    in ``self.arrays`` with key of 'bonds'.
    This method returns a copy of ``self.arrays["bonds"]``.

    """
    return self.arrays["bonds"].copy()

  def get_bonds_per_atom(self):
    """Return a nested list of bond data.

    The returned list is two-dimensional:
    each element of the first list corresponds to each atom,
    and the second list consists of *absolute* indices for atoms
    bonded with the each atom.

    .. note::

      This method is different from ``self.get_bonds()`` in terms of
      the following two points:

      * No data about image flags.
      * Contains data of *existing* bonds only.

    """
    return [
      [i + b[0] for b in bs if b[0] != 0]
      for i, bs in enumerate(self.arrays["bonds"])
    ]

  def get_types(self):
    """Return an array of integers containing atom's type.

    An instance of ``ase.Atoms`` stores data in ``self.arrays``;
    for example, atoms' positions are accessed by key 'positions'.
    For BondedAtoms instance, atom's type data is also stored
    in ``self.arrays`` with key of 'types'.
    This method returns a copy of ``self.arrays["types"]``.

    """
    return self.arrays["types"].copy()

  def get_bonded_bonds(self):
    """Return a two-dimensional ``numpy.ndarray`` describing all bonds.

    Shape of the returned array is (*N_b*, 2),
    where *N_b* is the number of bonds.
    Each row (the first axis) corresponds to each bond,
    and each column (the second axis) consists of *absolute* index
    of the first and second atom connected by the each bond.

    """
    bonds = set()

    bonds_per_atom = self.get_bonds_per_atom()

    for i, bs in enumerate(bonds_per_atom):
      bonds |= set((i, j) for j in bs if (j, i) not in bonds)

    return np.array(sorted(list(bonds)), int)

  def get_bonded_angles(self):
    """Return a two-dimensional ``numpy.ndarray``
    describing all angles defined by two consecutive bonds.

    Shape of the returned array is (*N_a*, 3),
    where *N_a* is the number of angles
    defined by two consecutive bonds.
    Each row (the first axis) corresponds to each angle,
    and each column (the second axis) consists of *absolute* index
    of the first, second and third atom in the each angle.

    """
    bonds_per_atom = self.get_bonds_per_atom()

    return np.array([
        (i, j, k)
        for j, bs in enumerate(bonds_per_atom)
        for i, k in itertools.combinations(bs, 2)
      ], int)

  def get_bonded_dihedrals(self):
    """Return a two-dimensional ``numpy.ndarray``
    describing all dihedrals defined by three consecutive bonds.

    Shape of the returned array is (*N_d*, 4),
    where *N_d* is the number of dihedrals
    defined by three consecutive bonds.
    Each row (the first axis) corresponds to each dihedral,
    and each column (the second axis) consists of *absolute* index
    of the first, second, third and fourth atom in the each dihedral.

    """
    bonds = set()
    dihedrals = []

    bonds_per_atom = self.get_bonds_per_atom()

    for j, bs in enumerate(bonds_per_atom):
      ks = [k for k in bs if (k, j) not in bonds]
      for k in ks:
        bonds.add((j, k))
        dihedrals.extend(
          (i, j, k, l) if n%2 else (l, k, j, i)  # equally distribute dihedrals to the second atom
          for n, (i, l) in enumerate(itertools.product(
            set(bs)-{k}, set(bonds_per_atom[k])-{j})))

    dihedrals.sort(key=lambda x: x[1])

    return np.array(dihedrals, int)

  def get_bonded_impropers(self):
    """Return a two-dimensional ``numpy.ndarray``
    describing all impropers defined by three bonds
    connected to the same atom.

    Shape of the returned array is (*N_i*, 4),
    where *N_i* is the number of impropers
    defined by three bonds connected to the same atom.
    Each row (the first axis) corresponds to each improper,
    and each column (the second axis) consists of *absolute* index
    of the first, second, third and fourth atom in the each improper;
    the first atom is the center one.

    """
    bonds_per_atom = self.get_bonds_per_atom()

    return np.array([
        (i,) + t
        for i, bs in enumerate(bonds_per_atom)
        for t in itertools.combinations(bs, 3)
      ], int)

  def remove_bond(self, atom1, atom2, img1=(0,0,0), img2=(0,0,0)):
    """Remove a bond connecting the given two atoms.

    Parameters:

    atom1: int
      Index for the first atom.

    atom2: int
      Index for the second atom.

    img1: tuple/list of int
      Three-element tuple or list of integers indicating which image of
      a periodic simulation box the first atom is in.

    img2: tuple/list of int
      Three-element tuple or list of integers indicating which image of
      a periodic simulation box the second atom is in.

    """
    if not isinstance(img1, (tuple, list)) or len(img1) != 3:
      raise RuntimeError(
        "Image flag 'img1' should be 3-membered tuple/list")
    elif not isinstance(img2, (tuple, list)) or len(img2) != 3:
      raise RuntimeError(
        "Image flag 'img2' should be 3-membered tuple/list")

    relative_data = np.asarray([  # not `np.array` due to pylint alarm
        atom2-atom1, img2[0]-img1[0], img2[1]-img1[1], img2[2]-img1[2]
      ], int)

    self._remove_bond(
      atom1, self._get_matched_bond_index(atom1, relative_data))
    self._remove_bond(
      atom2, self._get_matched_bond_index(atom2, -relative_data))

  def set_bonds(self, bonds):
    """Store the given bond data.

    This method sets bond data to ``self.array``
    using ``self.set_array()``.

    Parameters:

    bonds: array-like object of int
      Array describing bonds. Shape of the array must be (*N*, *M*, 4),
      where *N* is the number of atoms and *M* is the maximum number
      of bonds per atom. Four integers of innermost-axis are as follows:
      the first integer is a relative index for the other atom,
      the second to fourth integers are relative image flags
      for the other atom in the *x*, *y* and *z* direction, respectively.

    """
    self.set_array("bonds", bonds, int, ())

  def set_types(self, types):
    """Store the given type of atoms.

    This method sets type of atoms to ``self.array``
    using ``self.set_array()``.

    Parameters:

    types: array-like object of int
      Per-atom array containing atom's type.
      Length of the array must be equal to the number of atoms.

    """
    self.set_array("types", types, int, ())

  def sort_bonds(self):
    """Sort bonds in each atom.

    Note that the order of atoms stays unchanged.

    """
    def key_for_sorting_bonds(x):
      return (
        0 if x[0] > 0 else 1 if x[0] < 0 else 2,
        0 if all(i == 0 for i in x[1:]) else 1) + tuple(x)

    self.set_array("bonds", [
        sorted(bs, key=key_for_sorting_bonds)
        for bs in self.arrays["bonds"]
      ], int, ())

  def adjust_pbc_bonds(self):
    """Adjust bonds under periodic boundary conditions.

    This method corrects image flags of each bond if distance
    between two atoms connected with the bond is longer than
    a half length of the unit cell side.

    .. note::

      If the *i* th atom connects to one of images of the *j* th atom,
      it **must not** connect to another image of the the *j* th atom.
      In addition, all the bond lengths must be less than
      a half length of the unit cell side.

    """
    cell = self.get_cell()
    pbc = self.get_pbc()

    if pbc.any() and np.asarray(cell)[pbc].any(axis=1).all() and (cell, pbc):

      positions = self.get_positions()

      for i, bs in enumerate(self.arrays["bonds"]):
        for b in bs:
          if not np.all(b == 0):
            b[1:] = _compute_image_flags(
              positions[i+b[0]] - positions[i], cell, pbc)

  def wrap(self, center=(0.5, 0.5, 0.5), pbc=None, eps=1e-7):
      """Wrap positions to unit cell.

      Parameters:

      center: array-like object of float
          The positons in fractional coordinates that the new positions
          will be nearest possible to.

      pbc: bool or list/tuple of bool
        For each axis in the unit cell decides whether the positions
        will be moved along this axis.  By default, the boundary
        conditions of the Atoms object will be used.

      eps: float
        Small number to prevent slightly negative coordinates from being
        wrapped.

      """
      super().wrap(center, pbc, eps)
      self.adjust_pbc_bonds()

  def __delitem__(self, idx):
    """Delete a selected atom.

    This method overrides ``ase.Atoms.__delitem__``.
    Before calling the parent method to remove a selected atom,
    this method deletes bond data associated with the atom.

    Parameters:

    idx: int or list/tuple of int
      Index for an atom to be deleted.
      A list or tuple of integers is also acceptable for multiple atoms.

    Examples:

    >>> from ease4lmp import BondedAtoms
    >>> atoms = BondedAtoms('CO', positions=[(0, 0, 0), (0, 0, 1.1)])
    >>> atoms.set_bonds([
    ...   [[1, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0]],
    ...   [[-1, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0]]
    ... ])
    >>> atoms.get_bonds()
    array([[[ 1,  0,  0,  0],
            [ 0,  0,  0,  0],
            [ 0,  0,  0,  0],
            [ 0,  0,  0,  0]],
           [[-1,  0,  0,  0],
            [ 0,  0,  0,  0],
            [ 0,  0,  0,  0],
            [ 0,  0,  0,  0]]])
    >>> del atoms[0]
    >>> atoms.get_bonds()
    array([[[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]]])

    """
    if isinstance(idx, int):
      idx = np.array([idx])
    elif isinstance(idx, (tuple, list)) and len(idx) > 0:
      idx = np.array(idx)

    mask = np.ones(len(self), bool)
    mask[idx] = False

    for i, bs in enumerate(self.arrays["bonds"]):

      if not mask[i]:
        continue

      ib = 0
      while ib < self._max_bonds:
        j = i + bs[ib][0]
        if not mask[j]:
          self._remove_bond(i, ib)
        else:
          bs[ib][0] += (-np.sum(mask[i:j] == False)
            if i < j else np.sum(mask[j:i] == False))
          ib += 1

    super().__delitem__(idx)

  def __imul__(self, m):
    """Create a super cell and returns it.

    This method overrides ``ase.Atoms.__imul__``.
    After calling the parent method to expand the simulation box,
    this method recalculates image flags over the expanded bond data.
    Assuming that two atoms connected with a bond are
    in different images of the periodic simulation box before expansion.
    If the two atoms are in the same image after expansion,
    their image falgs for the opponent must be updated to (0, 0, 0).

    Parameters:

    m: int or list/tuple of int
      Three-element tuple or list defining how many times
      a periodic simulation box is repeated in each direction.
      The first, second and third element corresponds to the *x*, *y*
      and *z* direction, respectively.
      If an integer is given, a tuple ``(m, m, m)`` is used.

    Examples:

    >>> from ease4lmp import BondedAtoms
    >>> wire = BondedAtoms(
    ...   'Au',
    ...   positions=[[0.0, 5.0, 5.0]],
    ...   cell=[2.9, 5.0, 5.0],
    ...   pbc=[1, 0, 0])
    >>> wire.add_bond(0, 0, img2=(1, 0, 0))
    >>> wire.get_bonds()
    array([[[ 0,  1,  0,  0],
            [ 0, -1,  0,  0],
            [ 0,  0,  0,  0],
            [ 0,  0,  0,  0]]])
    >>> wire *= (3,1,1)
    >>> wire.get_bonds()
    array([[[ 1,  0,  0,  0],
            [ 2, -1,  0,  0],
            [ 0,  0,  0,  0],
            [ 0,  0,  0,  0]],
           [[ 1,  0,  0,  0],
            [-1,  0,  0,  0],
            [ 0,  0,  0,  0],
            [ 0,  0,  0,  0]],
           [[-2,  1,  0,  0],
            [-1,  0,  0,  0],
            [ 0,  0,  0,  0],
            [ 0,  0,  0,  0]]])

    """
    if isinstance(m, int):
      m = (m, m, m)

    n = len(self)

    super().__imul__(m)

    bonds = self.arrays["bonds"]
    i0 = 0

    for m0 in range(m[0]):
      for m1 in range(m[1]):
        for m2 in range(m[2]):
          i1 = i0 + n
          for i in range(i0, i1):
            bs = [b for b in bonds[i] if np.any(b[1:] != 0)]
            for b in bs:
              img_old = b[1:]
              img_new = np.array([
                (m0 + img_old[0]) // m[0],
                (m1 + img_old[1]) // m[1],
                (m2 + img_old[2]) // m[2]], int)
              tmp = img_old - img_new * m  # m will be converted to <numpy.ndarray>
              b[0] += n * (tmp[0]*m[1]*m[2] + tmp[1]*m[2] + tmp[2])
              b[1:] = img_new
          i0 = i1

    return self

  def _get_available_bond_index(self, atom, relative_data, offset=0):
    """Return the minimum index available for bond data.

    Parameters:

    atom: int
      Index for an atom which a new bond is connected to.

    relative_data: numpy.ndarray of int
      Four-element one-dimensional ``numpy.ndarray``
      describing a new bond data to be added the given atom.

    offset: int
      Returned index is offset by this value.

    """
    for i, b in enumerate(self.arrays["bonds"][atom]):
      if (b == relative_data).all():
        raise RuntimeError(
          "Bond between atom '{}' and '{}' already exists"
          .format(atom, atom+relative_data[0]))
      elif np.all(b == 0):
        if offset == 0:
          return i
        else:
          offset -= 1

    raise RuntimeError(
      "Too many bonds attache to atom '{}'".format(atom))

  def _get_matched_bond_index(self, atom, relative_data):
    """Return an index for bond in the given atom.

    Parameters:

    atom: int
      Index for an atom.

    relative_data: numpy.ndarray of int
      Four-element one-dimensional ``numpy.ndarray``.
      This method returns an index for a bond described by this parameter.

    """
    for i, b in enumerate(self.arrays["bonds"][atom]):
      if (b == relative_data).all():
        return i

    raise RuntimeError(
      "No such a bond attaches to atom '{}'".format(atom))

  def _remove_bond(self, atom, bond_idx):
    """Remove a specified bond from the given atom.

    Parameters:

    atom: int
      Index for an atom.

    bond_idx: int
      Index for a bond (in the given atom) to be removed.

    """
    bonds = self.arrays["bonds"][atom]
    bonds[bond_idx:-1] = bonds[bond_idx+1:]
    bonds[-1] = np.zeros(4, int)
