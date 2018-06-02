"""
This file is for BondedAtoms class.

create: 2018/05/24 by Takayuki Kobayashi
"""

import sys

import ase
import numpy as np

class BondedAtoms(ase.Atoms):
  """
  This class ...
  """

  @staticmethod
  def inherit(atoms):
    """
    This method ...
    [Arguments]
    * atoms: <ase.Atoms>
    """
    atoms.__class__ = BondedAtoms
    atoms._max_bonds = 4
    atoms.new_array(
      "bonds", np.zeros((len(atoms), atoms._max_bonds, 4)), int)
    return atoms

  def __init__(self, *args, **kwargs):
    """
    This constructor ...
    """
    super().__init__(*args, **kwargs)

    self._max_bonds = 4

    # For each bond, index 0 is for relative index of the other atom,
    # indices from 1 to 3 are (relative) image flags of the other atom.
    self.new_array(
      "bonds", np.zeros((len(self), self._max_bonds, 4)), int)

  def add_bond(self, atom1, atom2, img1=(0,0,0), img2=(0,0,0)):
    """
    This method ...
    [Arguments]
    * atom1: <int>
    * atom2: <int>
    * img1: <tuple/list>
    * img2: <tuple/list>
    """

    if not all([
      isinstance(i, (tuple, list)) and len(i) == 3 for i in (img1, img2)
    ]):
      raise RuntimeError("Image flags should be 3-membered tuple/list")

    # get empty bonds
    bs1 = [b for b in self.arrays["bonds"][atom1] if np.all(b == 0)]
    bs2 = [b for b in self.arrays["bonds"][atom2] if np.all(b == 0)]

    i2 = 1 if atom1 == atom2 else 0

    if bs1 and bs2[i2:]:

      rel_idx = atom2 - atom1
      rel_imgs = np.array(img2, int) - np.array(img1, int)

      for b, sign in [(bs1[0], 1), (bs2[i2], -1)]:
        # apply data to the first empty bond
        b[0] = sign * rel_idx
        b[1:] = sign * rel_imgs

    else:
      raise RuntimeError("Too many bonds")

  def change_max_bonds(self, n=4):
    """
    This method ...
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
    This method gets integer array of bonds.
    [Return]
    <numpy.ndarray>
    """
    return self.arrays["bonds"].copy()

  def remove_bond(self, atom, bond):
    """
    This method ...
    [Arguments]
    * atom: <int>
    * bond: <int>
    """

    bs1 = self.arrays["bonds"][atom]
    ib1 = bond

    rel_idx = bs1[ib1][0]
    rel_imgs = bs1[ib1][1:]

    bs2 = self.arrays["bonds"][atom+rel_idx]
    ib2 = [
      i for i in range(self._max_bonds)
      if bs2[i][0] == -rel_idx and (bs2[i][1:] == -rel_imgs).all()][0]

    self._remove_bond(bs1, ib1)
    self._remove_bond(bs2, ib2)

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
          self._remove_bond(bs, ib)
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

  def _remove_bond(self, bonds, bond_idx):
    """
    This method ...
    [Arguments]
    * bonds: <numpy.ndarray>; bonds belonging to an atom
    * bond_idx: <int>; index of a bond to be removed
    """
    bonds[bond_idx:-1] = bonds[bond_idx+1:]
    bonds[-1] = np.zeros(4, int)
