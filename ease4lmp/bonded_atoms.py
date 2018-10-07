"""
This file is for BondedAtoms class.

create: 2018/05/24 by Takayuki Kobayashi
"""

import sys

import ase

import itertools as it
import numpy as np

class BondedAtoms(ase.Atoms):
  """
  This class inherits from ase.Atoms, and has functionalities to deal
  with bonds connetcting a pair of atoms. Basically, you can use methods
  of ase.Atoms in the same manner.
  """

  @staticmethod
  def inherit(atoms):
    """
    This method is a static method which takes an ase.Atoms instance,
    converts it to a BondedAtoms instance, and returns the instance.

    **Arguments**
      * TETE
      * atoms: <ase.Atoms>
    """
    atoms.__class__ = BondedAtoms
    atoms._max_bonds = 4
    atoms.new_array(
      "bonds", np.zeros((len(atoms), atoms._max_bonds, 4)), int)
    return atoms

  def __init__(self, *args, **kwargs):
    """
    This constructor takes the same arguments with ase.Atoms' one.
    """
    super().__init__(*args, **kwargs)

    self._max_bonds = 4

    # For each bond, index 0 is for relative index of the other atom,
    # indices from 1 to 3 are (relative) image flags of the other atom.
    self.new_array(
      "bonds", np.zeros((len(self), self._max_bonds, 4)), int)

  def add_bond(self, atom1, atom2, img1=(0,0,0), img2=(0,0,0)):
    """
    This method add a bond connecting two atoms (specified by index).
    The bond data is stored as a property of the first atom and the
    second atoms as following.

    [
      [rel-id, rel-img_x, rel-img_y, rel-img_z],  # 1st bond
      [rel-id, rel-img_x, rel-img_y, rel-img_z],  # 2nd bond
      ...
    ]

    'rel-id' is a relative id; id of the second atom minus id of the
    first atom, for the first atom. 'rel-img's is relative image flags
    used in the case of periodic boundary conditions.

    [Arguments]
    * atom1: <int>
    * atom2: <int>
    * img1: <tuple/list>
    * img2: <tuple/list>
    """

    if not isinstance(img1, (tuple, list)) or len(img1) != 3:
      raise RuntimeError(
        "Image flag 'img1' should be 3-membered tuple/list")
    elif not isinstance(img2, (tuple, list)) or len(img2) != 3:
      raise RuntimeError(
        "Image flag 'img2' should be 3-membered tuple/list")

    relative_data = np.insert(
      np.array(img2, int) - np.array(img1, int), 0, atom2 - atom1)

    idx1 = self._get_available_bond_index(atom1, relative_data)
    idx2 = self._get_available_bond_index(
      atom2, -relative_data, 1 if atom1 == atom2 else 0)

    self.arrays["bonds"][atom1][idx1] = relative_data
    self.arrays["bonds"][atom2][idx2] = -relative_data

  def change_max_bonds(self, n=4):
    """
    This method change the maximum number of bonds per atom. Default is
    four bonds per atom.

    [Arguments]
    * n: <int>; new maximum number of bonds per atom
    """

    if n < self._max_bonds:
      self.arrays["bonds"] = np.array([
        bs[:n] for bs in self.arrays["bonds"]])
    elif self._max_bonds < n:
      zeros = np.zeros((n - self._max_bonds, 4), int)
      self.arrays["bonds"] = np.array([
        np.append(bs, zeros, axis=0) for bs in self.arrays["bonds"]])

    self._max_bonds = n

  def extend(self, other):
    """
    Extend atoms object by appending atoms from *other*.

    [Arguments]
    * other: <ase.Atom/ase.Atoms/BondedAtoms>
    [Return]
    <BondedAtoms>
    """

    if isinstance(other, ase.Atom):
      other = self.__class__([other])

    if isinstance(other, ase.Atoms):
      other = self.__class__(other)

    return super().extend(other)

  def get_bonds(self):
    """
    This method returns an integer array of bond data (an element in the
    first axis per atom).

    [Return]
    <numpy.ndarray>
    """
    return self.arrays["bonds"].copy()

  def get_bonds_per_atom(self):
    """
    This method returns a nested list of bond data. Unlike get_bonds(),
    this method returns data of existing bonds only.

    [Return]
    <list>
    """
    return [
      [i + b[0] for b in bs if b[0] != 0]
      for i, bs in enumerate(self.arrays["bonds"])]

  def get_bonded_bonds(self):
    """
    This method gets an array of 2-membered sequence of atom-index
    representing bonds.

    [Return]
    <numpy.ndarray>
    """

    bonds = set()

    bonds_per_atom = self.get_bonds_per_atom()

    for i, bs in enumerate(bonds_per_atom):
      bonds |= set([(i, j) for j in bs if (j, i) not in bonds])

    return np.array(list(bonds), int)

  def get_bonded_angles(self):
    """
    This method gets an array of 3-membered sequence of atom-index
    representing angles.

    [Return]
    <numpy.ndarray>
    """

    bonds_per_atom = self.get_bonds_per_atom()

    return np.array([
      (i, j, k)
      for j, bs in enumerate(bonds_per_atom)
      for i, k in it.combinations(bs, 2)], int)

  def get_bonded_dihedrals(self):
    """
    This method gets an array of 4-membered sequence of atom-index
    representing dihedrals.

    [Return]
    <numpy.ndarray>
    """

    bonds = set()
    dihedrals = []

    bonds_per_atom = self.get_bonds_per_atom()

    for j, bs in enumerate(bonds_per_atom):
      ks = [k for k in bs if (k, j) not in bonds]
      for k in ks:
        bonds.add((j, k))
        dihedrals.extend([
          (i, j, k, l) for i, l in it.product(
            set(bs)-{k}, set(bonds_per_atom[k])-{j})])

    return np.array(dihedrals, int)

  def get_bonded_impropers(self):
    """
    This method gets an array of 4-membered sequence of atom-index
    representing impropers.

    [Return]
    <numpy.ndarray>
    """

    bonds_per_atom = self.get_bonds_per_atom()

    return np.array([
      (i,) + t
      for i, bs in enumerate(bonds_per_atom)
      for t in it.combinations(bs, 3)], int)

  def remove_bond(self, atom1, atom2, img1=(0,0,0), img2=(0,0,0)):
    """
    This method removes a bond connecting tow atoms.

    [Arguments]
    * atom: <int>
    * bond: <int>
    """

    relative_data = np.insert(
      np.array(img2, int) - np.array(img1, int), 0, atom2 - atom1)

    self._remove_bond(
      atom1, self._get_matched_bond_index(atom1, relative_data))
    self._remove_bond(
      atom2, self._get_matched_bond_index(atom2, -relative_data))

  def set_bonds(self, bonds):
    """
    This method updates 'bonds' array.

    [Arguments]
    * bonds: <numpy.ndarray>
    """
    self.set_array("bonds", bonds, int, ())

  def __delitem__(self, idx):
    """
    [Arguments]
    * idx: <tuple/list/int>
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
    In-place repeat of atoms.

    [Arguments]
    * m: <tuple/list/int>
    [Return]
    <BondedAtoms>
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
    [Arguments]
    * atom: <int>
    * relative_data: <numpy.ndarray>
    * offset: <int>; postpone returning index until this becomes zero
    [Return]
    <int>
    """

    for i, b in enumerate(self.arrays["bonds"][atom]):
      if (b == relative_data).all():
        raise RuntimeError(
          "Bond between atom '{}' and '{}' already exists".format(
            atom, atom+relative_data[0]))
      elif np.all(b == 0):
        if offset == 0:
          return i
        else:
          offset -= 1

    raise RuntimeError(
      "Too many bonds attache to atom '{}'".format(atom))

  def _get_matched_bond_index(self, atom, relative_data):
    """
    [Arguments]
    * atom: <int>
    * relative_data: <numpy.ndarray>
    [Return]
    <int>
    """

    for i, b in enumerate(self.arrays["bonds"][atom]):
      if (b == relative_data).all():
        return i

    raise RuntimeError(
      "No such a bond attaches to atom '{}'".format(atom))

  def _remove_bond(self, atom, bond_idx):
    """
    [Arguments]
    * atom: <int>; bonds belonging to an atom
    * bond_idx: <int>; index of a bond to be removed
    """
    bonds = self.arrays["bonds"][atom]
    bonds[bond_idx:-1] = bonds[bond_idx+1:]
    bonds[-1] = np.zeros(4, int)
