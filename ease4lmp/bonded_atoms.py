"""
@file ease4lmp/bonded_atoms.py
@brief This file is for `BondedAtoms` class.
@author Takayuki Kobayashi
@date 2018/05/24
"""

import sys

import ase

import itertools as it
import numpy as np


class BondedAtoms(ase.Atoms):
  """
  This class inherits from `ase.Atoms` class, and has functionalities
  to deal with (covalent) bonds connecting a pair of atoms.

  An instance of this class has bond data as a `numpy.ndarray`
  as with other atomic properties like positions and velocities.
  Shape of the array for bond data is (*N*, #_max_bonds, 4),
  where *N* is the number of atoms. The first axis of the array
  corresponds to atoms, and the second axis corresponds to
  bonds connected to each atom. The third axis contains data for
  each bond: an integer at 0 is relative index (of first axis)
  for the other atom, integers from 1 to 3 are relative image flags
  (to resolve a periodic boundary condition) for the other atom
  in the *x*, *y* and *z* direction, respectively.

  Basically, you can use methods of `ase.Atoms` in the same manner.
  """

  @staticmethod
  def inherit(atoms):
    """
    @param atoms
      An `ase.Atoms` instance to be down-casted to BondedAtoms.

    @return
      A BondedAtoms instance.

    This static method takes an `ase.Atoms` instance, converts it
    to a BondedAtoms instance, and returns the instance.
    This method does not adds any bond data, it only reserves array
    in which bond data will be stored.
    """

    atoms.__class__ = BondedAtoms
    atoms._max_bonds = 4

    atoms.new_array(
      "bonds", np.zeros((len(atoms), atoms._max_bonds, 4)), int)

    return atoms


  def __init__(self, *args, **kwargs):
    """
    @param *args
      A variable number of arguments.
    @param **kwargs
      A variable number of named arguments.

    Constructor of BondedAtoms takes the same arguments
    as that of `ase.Atoms`.
    """

    super().__init__(*args, **kwargs)

    ## Maximum number of bonds per atom.
    self._max_bonds = 4

    self.new_array(
      "bonds", np.zeros((len(self), self._max_bonds, 4)), int)


  def add_bond(self, atom1, atom2, img1=(0,0,0), img2=(0,0,0)):
    """
    @param atom1
      Index for the first atom.
    @param atom2
      Index for the second atom.
    @param img1
      Three-element tuple or list of integers indicating which image of
      a periodic simulation box the first atom is in.
    @param img2
      Three-element tuple or list of integers indicating which image of
      a periodic simulation box the second atom is in.

    This method adds a bond connecting two atoms. The bond data is
    stored as following.

    * Bond data in the first atom: `[atom2-atom1, *(img2-img1)]`
    * Bond data in the second atom: `[atom1-atom2, *(img1-img2)]`
    """

    if not isinstance(img1, (tuple, list)) or len(img1) != 3:
      raise RuntimeError(
        "Image flag 'img1' should be 3-membered tuple/list")
    elif not isinstance(img2, (tuple, list)) or len(img2) != 3:
      raise RuntimeError(
        "Image flag 'img2' should be 3-membered tuple/list")

    relative_data = np.insert(
      np.array(img2, int) - np.array(img1, int), 0, atom2 - atom1)

    idx1 = self._get_available_bond_index(
      atom1, relative_data)
    idx2 = self._get_available_bond_index(
      atom2, -relative_data, 1 if atom1 == atom2 else 0)

    self.arrays["bonds"][atom1][idx1] = relative_data
    self.arrays["bonds"][atom2][idx2] = -relative_data


  def change_max_bonds(self, n=4):
    """
    @param n
      New maximum number of bonds per atom.

    This method changes the maximum number of bonds per atom.
    Default is four bonds per atom.
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
    """
    @param other
      An instance of BondedAtoms, `ase.Atoms` or `ase.Atom` class.

    @return
      An extended BondedAtoms instance.

    This method extends this instance by appending another instance of
    BondedAtoms, `ase.Atoms` or `ase.Atom` class.
    """

    if isinstance(other, ase.Atom):
      other = self.__class__([other])

    if isinstance(other, ase.Atoms):
      other = self.__class__(other)

    return super().extend(other)


  def get_bonds(self):
    """
    @return
      A `numpy.ndarray` containing bond data.

    This method returns an array of integers containing bond data.
    For more details about the array, please see documentation of
    BondedAtoms class.
    """

    return self.arrays["bonds"].copy()


  def get_bonds_per_atom(self):
    """
    @return
      An atom-wise list containing *absolute* indices for bonded atoms.

    This method returns a nested list of bond data. Unlike #get_bonds,
    this method returns data of existing bonds only.
    Each element of the outer list corresponds to each atom.
    The inner list consists of *absolute* indices for atoms bonded with
    the each atom.
    """

    return [
      [i + b[0] for b in bs if b[0] != 0]
      for i, bs in enumerate(self.arrays["bonds"])
    ]


  def get_bonded_bonds(self):
    """
    @return
      A `numpy.ndarray` of *absolute* atom-index describing all bonds.

    This method returns a two-dimensional array describing bonds.
    Each row corresponds to each bond, and each column corresponds to
    the first and second atom connected by each bond.
    """

    bonds = set()

    bonds_per_atom = self.get_bonds_per_atom()

    for i, bs in enumerate(bonds_per_atom):
      bonds |= set([(i, j) for j in bs if (j, i) not in bonds])

    return np.array(list(bonds), int)


  def get_bonded_angles(self):
    """
    @return
      A `numpy.ndarray` of *absolute* atom-index describing
      all angles defined by bonds.

    This method returns a two-dimensional array describing angles.
    Each row corresponds to each angle, and each column corresponds to
    the first, second and third atom forming each angle.
    The second atom is the center one.
    """

    bonds_per_atom = self.get_bonds_per_atom()

    return np.array([
        (i, j, k)
        for j, bs in enumerate(bonds_per_atom)
        for i, k in it.combinations(bs, 2)
      ], int)


  def get_bonded_dihedrals(self):
    """
    @return
      A `numpy.ndarray` of *absolute* atom-index describing
      all dihedrals defined by bonds.

    This method returns a two-dimensional array describing dihedrals.
    Each row corresponds to each dihedral, and each column corresponds
    to the first, second, third and fourth atom forming each dihedral.
    The four atoms are connected linearly by three bonds.
    """

    bonds = set()
    dihedrals = []

    bonds_per_atom = self.get_bonds_per_atom()

    for j, bs in enumerate(bonds_per_atom):
      ks = [k for k in bs if (k, j) not in bonds]
      for k in ks:
        bonds.add((j, k))
        dihedrals.extend([
            (i, j, k, l)
            for i, l in it.product(
              set(bs)-{k}, set(bonds_per_atom[k])-{j})
          ])

    return np.array(dihedrals, int)


  def get_bonded_impropers(self):
    """
    @return
      A `numpy.ndarray` of *absolute* atom-index describing
      all impropers defined by bonds.

    This method returns a two-dimensional array describing impropers.
    Each row corresponds to each improper, and each column corresponds
    to the first, second, third and fourth atom forming each improper.
    The first atom is the center one.
    """

    bonds_per_atom = self.get_bonds_per_atom()

    return np.array([
        (i,) + t
        for i, bs in enumerate(bonds_per_atom)
        for t in it.combinations(bs, 3)
      ], int)


  def remove_bond(self, atom1, atom2, img1=(0,0,0), img2=(0,0,0)):
    """
    @param atom1
      Index for the first atom.
    @param atom2
      Index for the second atom.
    @param img1
      Three-element tuple or list of integers indicating which image of
      a periodic simulation box the first atom is in.
    @param img2
      Three-element tuple or list of integers indicating which image of
      a periodic simulation box the second atom is in.

    This method removes a bond connecting two atoms.
    """

    relative_data = np.insert(
      np.array(img2, int) - np.array(img1, int), 0, atom2 - atom1)

    self._remove_bond(
      atom1, self._get_matched_bond_index(atom1, relative_data))
    self._remove_bond(
      atom2, self._get_matched_bond_index(atom2, -relative_data))


  def set_bonds(self, bonds):
    """
    @param bonds
      A `numpy.ndarray` containing bond data.

    This method sets bond data by directly assigning an array.
    For more details about the array, please see documentation of
    BondedAtoms class.
    """

    self.set_array("bonds", bonds, int, ())


  def __delitem__(self, idx):
    """
    @param idx
      Index for an atom to be deleted.
      A list or tuple is also acceptable for multiple atoms.

    This method overrides `ase.Atoms.__delitem__`.
    Before calling the parent method to remove a selected atom,
    this method deletes bond data associated with the atom.
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
          bs[ib][0] += -np.sum(mask[i:j] == False) \
            if i < j else np.sum(mask[j:i] == False)
          ib += 1

    super().__delitem__(idx)


  def __imul__(self, m):
    """
    @param m
      Three-element tuple or list defining how many times
      a periodic simulation box is repeated in each direction.
      The first, second and third element corresponds to the *x*, *y*
      and *z* direction, respectively.
      If an integer is given as `m`, a tuple `(m, m, m)` is used.

    @return
      A BondedAtoms instance containing all the repeated atoms
      with the expanded simulation box.

    This method overrides `ase.Atoms.__imul__`.
    After calling the parent method to expand the simulation box,
    this method resolves bonds originally connecting two atoms in
    different images of the periodic simulation box.
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
    """
    @param atom
      Index for an atom.
    @param relative_data
      A `numpy.ndarray` of which shape is (4): new bond data
      to be added by #add_bond.
    @param offset
      Returned index is offset by this value.

    @return
      The minimum available index in the second axis of an array
      for bond data.

    For more details about the array for bond data,
    please see documentation of BondedAtoms class.
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
    """
    @param atom
      Index for an atom.
    @param relative_data
      A `numpy.ndarray` of which shape is (4): existing bond data
      to be removed by #remove_bond.

    @return
      Index in the second axis of an array for bond data,
      where data matches `relative_data`.

    For more details about the array for bond data, please see
    documentation of BondedAtoms class.
    """

    for i, b in enumerate(self.arrays["bonds"][atom]):
      if (b == relative_data).all():
        return i

    raise RuntimeError(
      "No such a bond attaches to atom '{}'".format(atom))


  def _remove_bond(self, atom, bond_idx):
    """
    @param atom
      Index for an atom to which a bond to be removed belongs.
    @param bond_idx
      Index for the bond to be removed (index in the second axis
      of an array for bond data).

    For more details about the array for bond data, please see
    documentation of BondedAtoms class.
    """

    bonds = self.arrays["bonds"][atom]
    bonds[bond_idx:-1] = bonds[bond_idx+1:]
    bonds[-1] = np.zeros(4, int)
