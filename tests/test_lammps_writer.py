import unittest

from ease4lmp import BondedAtoms, LammpsWriter

from ase import Atoms
from ase.build import molecule

import numpy as np


class TestBondedAtoms(unittest.TestCase):

  def test_acetic(self):
    """Test for down-casting from ``ase.Atoms`` to BondedAtoms."""
    atoms = BondedAtoms(molecule("CH3COOH"))