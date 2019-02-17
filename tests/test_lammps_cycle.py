import unittest

from ease4lmp import (
  BondedAtoms, LammpsWriter,
  create_atoms_from_data, create_atoms_from_molecule)

from ase import Atoms
from ase.build import molecule

import numpy as np

import os
import sys
import itertools

from io import StringIO

io = StringIO()
sys.stdout = io


def write_files(atoms):
  writer = LammpsWriter(atoms, atom_style="molecular")

  writer.set_atom_data(mol=[0]*len(atoms))

  writer.set_bond_types({
      seq: i+1 for i, seq in enumerate(writer.get_bond_patterns())
    })

  writer.set_angle_types({
      seq: i+1 for i, seq in enumerate(writer.get_angle_patterns())
    })

  writer.set_dihedral_types({
      seq: i+1 for i, seq in enumerate(writer.get_dihedral_patterns())
    })

  writer.set_improper_types({
      seq: i+1 for i, seq in enumerate(writer.get_improper_patterns())
    })

  writer.write_lammps_data("data.tmp", mass=True)
  writer.write_lammps_molecule("molecule.tmp", mass=True)

def remove_files():
  os.remove("data.tmp")
  os.remove("molecule.tmp")


class TestLammpsCycle(unittest.TestCase):

  def test_methanol(self):
    """Test for equivalence between original and written/read data."""
    atoms = BondedAtoms(molecule("CH3OH"))

    # confirm atomic numbers
    self.assertTrue(np.allclose(
      atoms.get_atomic_numbers(), np.array([6, 8, 1, 1, 1, 1])))

    # confirm O-H distance
    self.assertTrue(np.allclose(atoms.get_distance(1, 3), 0.97))

    atoms.set_types([1, 2, 3, 4, 3, 3])

    positions = atoms.get_positions()

    bonded_pairs = (
      (i, j) for i, j in itertools.combinations(range(len(atoms)), 2)
      if np.linalg.norm(positions[i] - positions[j]) < 1.5)

    for pair in bonded_pairs:
      atoms.add_bond(*pair)

    atoms.sort_bonds()

    atoms.set_cell([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]])
    atoms.center()

    write_files(atoms)

    atoms_from_data = create_atoms_from_data("data.tmp", "molecular")
    atoms_from_molecule = create_atoms_from_molecule("molecule.tmp")

    # atoms from Lammps' data and molecule file must be eaqual.
    self.assertTrue(np.allclose(
      atoms_from_data.get_positions(), atoms_from_molecule.get_positions()))
    self.assertTrue(np.allclose(
      atoms_from_data.get_masses(), atoms_from_molecule.get_masses()))
    self.assertTrue(np.allclose(
      atoms_from_data.get_types(), atoms_from_molecule.get_types()))
    self.assertTrue(np.allclose(
      atoms_from_data.get_bonds(), atoms_from_molecule.get_bonds()))

    # comparison with original atoms
    self.assertTrue(np.allclose(
      atoms_from_data.get_positions(), atoms.get_positions()))
    self.assertTrue(np.allclose(
      atoms_from_data.get_masses(), atoms.get_masses()))
    self.assertTrue(np.allclose(
      atoms_from_data.get_types(), atoms.get_types()))

    # storing order of bonds might be changed
    atoms_from_data.sort_bonds()
    self.assertTrue(np.allclose(
      atoms_from_data.get_bonds(), atoms.get_bonds()))

    remove_files()

  def test_acetic(self):
    """Test for equivalence between original and written/read data."""
    atoms = BondedAtoms(molecule("CH3COOH"))

    # confirm atomic numbers
    self.assertTrue(np.allclose(
      atoms.get_atomic_numbers(), np.array([6, 8, 8, 1, 6, 1, 1, 1])))

    # confirm O-H distance < C=H distance
    self.assertTrue(all(
      atoms.get_distance(2, 3) < atoms.get_distance(i, j)
      for i, j in [(4, 5), (4, 6), (4, 7)]))

    atoms.set_types([1, 2, 3, 4, 5, 6, 6, 6])

    positions = atoms.get_positions()

    bonded_pairs = (
      (i, j) for i, j in itertools.combinations(range(len(atoms)), 2)
      if np.linalg.norm(positions[i] - positions[j]) < 1.5)

    for pair in bonded_pairs:
      atoms.add_bond(*pair)

    atoms.sort_bonds()

    atoms.set_cell([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]])
    atoms.center()

    write_files(atoms)

    atoms_from_data = create_atoms_from_data("data.tmp", "molecular")
    atoms_from_molecule = create_atoms_from_molecule("molecule.tmp")

    # atoms from Lammps' data and molecule file must be eaqual.
    self.assertTrue(np.allclose(
      atoms_from_data.get_positions(), atoms_from_molecule.get_positions()))
    self.assertTrue(np.allclose(
      atoms_from_data.get_masses(), atoms_from_molecule.get_masses()))
    self.assertTrue(np.allclose(
      atoms_from_data.get_types(), atoms_from_molecule.get_types()))
    self.assertTrue(np.allclose(
      atoms_from_data.get_bonds(), atoms_from_molecule.get_bonds()))

    # comparison with original atoms
    self.assertTrue(np.allclose(
      atoms_from_data.get_positions(), atoms.get_positions()))
    self.assertTrue(np.allclose(
      atoms_from_data.get_masses(), atoms.get_masses()))
    self.assertTrue(np.allclose(
      atoms_from_data.get_types(), atoms.get_types()))

    # storing order of bonds might be changed
    atoms_from_data.sort_bonds()
    self.assertTrue(np.allclose(
      atoms_from_data.get_bonds(), atoms.get_bonds()))

    remove_files()

def suite():
  suite = unittest.TestSuite()
  suite.addTest(TestLammpsCycle("test_methanol"))
  suite.addTest(TestLammpsCycle("test_acetic"))
  return suite