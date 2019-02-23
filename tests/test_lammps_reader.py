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

      bonds = read_bonds(path1)
      angles = read_angles(path1)
      dihedrals = read_dihedrals(path1)
      impropers = read_impropers(path1)

      self.assertEqual(bonds, read_bonds(path2))
      self.assertEqual(angles, read_angles(path2))
      self.assertEqual(dihedrals, read_dihedrals(path2))
      self.assertEqual(impropers, read_impropers(path2))

      self.assertTrue(all(
        set(b.keys())
          == {"id", "type"} | set("atom{}-id".format(i) for i in range(1, 3))
        for b in bonds))
      self.assertTrue(all(
        set(a.keys())
          == {"id", "type"} | set("atom{}-id".format(i) for i in range(1, 4))
        for a in angles))
      self.assertTrue(all(
        set(d.keys())
          == {"id", "type"} | set("atom{}-id".format(i) for i in range(1, 5))
        for d in dihedrals))
      self.assertTrue(all(
        set(i.keys())
          == {"id", "type"} | set("atom{}-id".format(i) for i in range(1, 5))
        for i in impropers))

  def test_data_and_molecule_topology(self):
    """Two data objects from Lammps' data and molecule file must be eaqual."""
    for name in ["acetic", "methanol"]:

      path1 = "lammps_files/data.{}".format(name)
      path2 = "lammps_files/molecule.{}".format(name)

      for atom1, atom2 in zip(
        read_atoms_from_data(path1, "molecular"),
        read_atoms_from_molecule(path2)):

        keys = set(atom1.keys()) & set(atom2.keys())

        self.assertEqual(keys, {"id", "type", "mass", "xu", "yu", "zu"})
        self.assertTrue(all(atom1[k] == atom2[k] for k in keys))

def suite():
  suite = unittest.TestSuite()
  suite.addTest(TestLammpsReader("test_data_and_molecule_atom"))
  suite.addTest(TestLammpsReader("test_data_and_molecule_topology"))
  return suite