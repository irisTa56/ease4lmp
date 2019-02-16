import unittest

from ease4lmp import (
  read_bonds, read_angles, read_dihedrals, read_impropers,
  read_atoms_from_data, read_atoms_from_molecule)


class TestLammpsReader(unittest.TestCase):

  def test_data_and_molecule_atom(self):
    """Two data objects from Lammps' data and molecule file must be eaqual."""
    for name in ["acetic", "methanol"]:

      path1 = "lammps_files/data.{}".format(name)
      path2 = "lammps_files/molecule.{}".format(name)

      self.assertEqual(read_bonds(path1), read_bonds(path2))
      self.assertEqual(read_angles(path1), read_angles(path2))
      self.assertEqual(read_dihedrals(path1), read_dihedrals(path2))
      self.assertEqual(read_impropers(path1), read_impropers(path2))

  def test_data_and_molecule_topology(self):
    """Two data objects from Lammps' data and molecule file must be eaqual."""
    for name in ["acetic", "methanol"]:

      path1 = "lammps_files/data.{}".format(name)
      path2 = "lammps_files/molecule.{}".format(name)

      for atom1, atom2 in zip(
        read_atoms_from_data(path1, "molecular", mass=True),
        read_atoms_from_molecule(path2)):

        keys = set(atom1.keys()) & set(atom2.keys())

        self.assertEqual(keys, {"id", "type", "mass", "xu", "yu", "zu"})
        self.assertTrue(all(atom1[k] == atom2[k] for k in keys))

def suite():
  suite = unittest.TestSuite()
  suite.addTest(TestLammpsReader("test_data_and_molecule_atom"))
  suite.addTest(TestLammpsReader("test_data_and_molecule_topology"))
  return suite