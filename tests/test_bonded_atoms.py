import unittest

from ease4lmp import BondedAtoms

from ase import Atoms
from ase.build import molecule

import numpy as np


class TestBondedAtoms(unittest.TestCase):

  def test_casting(self):
    """Test for down-casting from ``ase.Atoms`` to BondedAtoms."""
    atoms = molecule("CH3CH2OH")
    bonded_atoms = BondedAtoms(atoms)

    # the instance is BondedAtoms?
    self.assertTrue(isinstance(bonded_atoms, BondedAtoms))
    # 'bonds' property has been created?
    self.assertTrue("bonds" in bonded_atoms.arrays)
    # 'types' property has been created?
    self.assertTrue("types" in bonded_atoms.arrays)
    # all data except for 'bonds' stays unchanged?
    self.assertTrue(all(
      np.allclose(v, atoms.arrays[k])
      for k, v in bonded_atoms.arrays.items()
      if k not in ["bonds", "types"]))

  def test_gold(self):
    """Test for basic operations."""
    a = 4.05  # gold lattice constant
    b = a / 2
    au = BondedAtoms(
      "Au", cell=[(0, b, b), (b, 0, b), (b, b, 0)], pbc=True)

    # initial types are 1?
    self.assertTrue(np.allclose(au.get_types(), np.ones(len(au))))

    au.add_bond(0, 0, img2=(1,0,0))

    # a bond is added?
    self.assertTrue(np.allclose(
      au.get_bonds(), [[[ 0,  1,  0,  0],
                        [ 0, -1,  0,  0],
                        [ 0,  0,  0,  0],
                        [ 0,  0,  0,  0]]]))

    # multiplying and repeating give the same result?
    self.assertTrue(np.allclose(
      (au*(2, 2, 2)).get_bonds(), au.repeat((2, 2, 2)).get_bonds()))

    chain = au * (3, 1, 1)

    # types are correctly copied?
    self.assertTrue(np.allclose(chain.get_types(), np.ones(len(chain))))

    chain.change_max_bonds(2)

    # the maximum number of bonds per atoms can be reduced?
    self.assertTrue(np.allclose(
      chain.get_bonds(), [[[ 1,  0,  0,  0],
                           [ 2, -1,  0,  0]],
                          [[ 1,  0,  0,  0],
                           [-1,  0,  0,  0]],
                          [[-2,  1,  0,  0],
                           [-1,  0,  0,  0]]]))

    chain.remove_bond(0, 2, img2=(-1,0,0))

    # a specified bond is removed?
    self.assertTrue(np.allclose(
      chain.get_bonds(), [[[ 1,  0,  0,  0],
                           [ 0,  0,  0,  0]],
                          [[ 1,  0,  0,  0],
                           [-1,  0,  0,  0]],
                          [[-1,  0,  0,  0],
                           [ 0,  0,  0,  0]]]))

    chain.set_types([1, 2, 3])

    # types are correctly set?
    self.assertTrue(np.allclose(chain.get_types(), [1, 2, 3]))

    del chain[2]

    # a specified atom is deleted?
    self.assertTrue(np.allclose(chain.get_types(), [1, 2]))

    chain.change_max_bonds(4)

    # the maximum number of bonds per atoms can be increased?
    self.assertTrue(np.allclose(
      chain.get_bonds(), [[[ 1,  0,  0,  0],
                           [ 0,  0,  0,  0],
                           [ 0,  0,  0,  0],
                           [ 0,  0,  0,  0]],
                          [[-1,  0,  0,  0],
                           [ 0,  0,  0,  0],
                           [ 0,  0,  0,  0],
                           [ 0,  0,  0,  0]]]))

  def test_methanol(self):
    """Test for topology computation."""
    methanol = BondedAtoms(molecule("CH3OH"))

    # confirm atomic numbers
    self.assertTrue(np.allclose(
      methanol.get_atomic_numbers(), np.array([6, 8, 1, 1, 1, 1])))

    # confirm O-H distance
    self.assertTrue(np.allclose(methanol.get_distance(1, 3), 0.97))

    bond_list = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 3)]
    for t in bond_list:
      methanol.add_bond(*t)

    # bonds are calculated correctly?
    self.assertEqual(
      set(tuple(b) for b in methanol.get_bonded_bonds()),
      {(0, 1), (1, 3), (0, 5), (0, 4), (0, 2)})

    # angles are calculated correctly?
    self.assertEqual(
      set(tuple(a) for a in methanol.get_bonded_angles()),
      {(1, 0, 2), (1, 0, 4), (1, 0, 5), (2, 0, 4), (2, 0, 5), (4, 0, 5), (0, 1, 3)})

    # dihedrals are calculated correctly?
    dihedrals_ref = {(2, 0, 1, 3), (4, 0, 1, 3), (5, 0, 1, 3)}
    for d in methanol.get_bonded_dihedrals():
      if tuple(d) in dihedrals_ref:
        dihedrals_ref.remove(tuple(d))
      elif tuple(d[::-1]) in dihedrals_ref:
        dihedrals_ref.remove(tuple(d[::-1]))
    self.assertTrue(len(dihedrals_ref) == 0)

    # impropers are calculated correctly?
    self.assertEqual(
      set(tuple(i) for i in methanol.get_bonded_impropers()),
      {(0, 1, 2, 4), (0, 1, 2, 5), (0, 1, 4, 5), (0, 2, 4, 5)})


def suite():
  suite = unittest.TestSuite()
  suite.addTest(TestBondedAtoms("test_casting"))
  suite.addTest(TestBondedAtoms("test_gold"))
  suite.addTest(TestBondedAtoms("test_methanol"))
  return suite
