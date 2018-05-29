"""
This file is for BondedAtoms class.

create: 2018/05/24 by Takayuki Kobayashi
"""

import sys

import ase
import numpy as np

tl = (tuple, list)

class BondedAtoms(ase.Atoms):
  """
  This class ...
  """

  def __init__(self, *args, **kwargs):
    """
    This constructor ...
    """

    super().__init__(*args, **kwargs)

    self._max_bonds = 4

    # For each bond, index 0 is for relative index of the other atom,
    # indices from 1 to 3 are (relative) image flags of the other atom.
    self.new_array(
      "bonds", np.zeros((len(self), self._max_bonds, 4)), int
    )

  def add_bond(self, atom1, atom2, img1=(0,0,0), img2=(0,0,0)):
    """
    This method ...
    [Arguments]
    * atom1: <int>
    * atom2: <int>
    * img1: <tuple/list>
    * img2: <tuple/list>
    """

    if not all([a < len(self) for a in (atom1, atom2)]):
      raise IndexError("Atom index(s) is out of range")

    if not all([isinstance(i, tl) and len(i) == 3 for i in (img1, img2)]):
      raise RuntimeError("Image flags should be 3-membered tuple/list")

    # get empty bonds
    bs1, bs2 = [[
      b for b in self.arrays["bonds"][a] if np.all(b == 0)
    ] for a in (atom1, atom2)]

    i2 = 1 if atom1 == atom2 else 0

    if bs1 and bs2[i2:]:

      rel_idx = atom2 - atom1
      rel_imgs = np.array(img2) - np.array(img1)

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
    * n: <int> new maximum number of bonds per atom
    """

    if n < self._max_bonds:
      self.arrays["bonds"] = np.array([
        bs[:n] for bs in self.arrays["bonds"]
      ])
    elif self._max_bonds < n:
      zeros = np.zeros((n - self._max_bonds, 4), int)
      self.arrays["bonds"] = np.array([
        np.append(bs, zeros, axis=0) for bs in self.arrays["bonds"]
      ])

    self._max_bonds = n

  def extend(self, other):
    """
    Extend atoms object by appending atoms from *other*.
    [Arguments]
    * other: <ase.Atom> or <ase.Atoms> or <BondedAtoms>
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
    * idx: <tuple/list> or <int>
    """

    n = len(self)

    for c in self._constraints:
      if not isinstance(c, ase.constraints.FixAtoms):
        raise RuntimeError("Remove constraint using set_constraint() before deleting atoms.")

    if isinstance(idx, int):
      idx = np.array([idx])
    elif isinstance(idx, tl) and len(idx) > 0:
      idx = np.array(idx)

    if len(self._constraints) > 0:
      idx = np.arange(n)[idx]
      constraints = []
      for c in self._constraints:
        c = c.delete_atoms(idx, n)
        if c is not None:
          constraints.append(c)
      self.constraints = constraints

    mask = np.ones(n, bool)
    mask[idx] = False

    for i, bs in enumerate(self.arrays["bonds"]):

      if not mask[i]:
        continue

      ib = 0
      while ib < self._max_bonds:
        j = i + bs[ib][0]
        if not mask[j]:
          bs[ib:-1] = bs[ib+1:]
          bs[-1] = np.zeros(4, int)
        else:
          bs[ib][0] += -np.sum(mask[i:j] == False) \
            if i < j else np.sum(mask[j:i] == False)
          ib += 1

    for name, a in self.arrays.items():
      self.arrays[name] = a[mask]

  def __imul__(self, m):
    """
    In-place repeat of atoms.
    [Arguments]
    * m: <tuple/list> or <int>
    [Return]
    <BondedAtoms>
    """

    if isinstance(m, int):
      m = (m, m, m)

    for x, vec in zip(m, self._cell):
      if x != 1 and not vec.any():
        raise ValueError("Cannot repeat along undefined lattice vector")

    M = np.product(m)
    n = len(self)

    for name, a in self.arrays.items():
      self.arrays[name] = np.tile(a, (M,) + (1,) * (len(a.shape) - 1))

    positions = self.arrays["positions"]
    i0 = 0
    for m0 in range(m[0]):
      for m1 in range(m[1]):
        for m2 in range(m[2]):
          i1 = i0 + n
          positions[i0:i1] += np.dot((m0, m1, m2), self._cell)
          i0 = i1

    if self.constraints is not None:
      self.constraints = [c.repeat(m, n) for c in self.constraints]

    self._cell = np.array([m[c] * self._cell[c] for c in range(3)])

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
                (m2 + img_old[2]) // m[2],
              ])
              tmp = img_old - img_new * m  # m will be converted to <numpy.ndarray>
              b[0] += n * (tmp[0]*m[1]*m[2] + tmp[1]*m[2] + tmp[2])
              b[1:] = img_new
          i0 = i1

    return self
