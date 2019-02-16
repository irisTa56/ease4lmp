import unittest

from ease4lmp import (
  read_bonds, read_angles, read_dihedrals, read_impropers)


class TestLammpsReader(unittest.TestCase):

  def test_data_and_molecule(self):
    """Two data objects from Lammps' data and molecule file must be eaqual."""
    for name in ["acetic", "methanol"]:
      path1 = "lammps_files/data.{}".format(name)
      path2 = "lammps_files/molecule.{}".format(name)
      self.assertEqual(read_bonds(path1), read_bonds(path2))
      self.assertEqual(read_angles(path1), read_angles(path2))
      self.assertEqual(read_dihedrals(path1), read_dihedrals(path2))
      self.assertEqual(read_impropers(path1), read_impropers(path2))


def suite():
  suite = unittest.TestSuite()
  suite.addTest(TestLammpsReader("test_data_and_molecule"))
  return suite