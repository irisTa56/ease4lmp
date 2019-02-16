"""Submodule for functions to read Lammps' data file."""

from .lammps_dataformats import lmp_datanames


def _read_section(path, section):
  """Reads lines of a specified section in a specified file.

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  section: str
    Header of a section to be read.

  """
  with open(path, "r") as f:
    lines = (line.lstrip().rstrip() for line in f.readlines())

  blank_counter = 0
  in_section = False
  splitted_lines = []

  for line in lines:
    if in_section:
      if line == '':
        blank_counter += 1
        if 1 < blank_counter:
          break
      else:
        splitted_lines.append(line.split())
    else:
      if line.startswith(section):
        in_section = True

  return splitted_lines

def _read_topology_components(path, name, header):
  """Reads data of specified topology components
  from a specified Lammps' data (or molecule) file.

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  name: str
    Name of topology component.

  header: str
    Section header of the topology components.

  """
  return [
    dict(zip(lmp_datanames[name], map(int, line)))
    for line in _read_section(path, header)
  ]

#=======================================================================

def read_bonds(path):
  """Reads bonds data from a specified Lammps' data (or molecule) file.

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  """
  return _read_topology_components(path, "bond", "Bonds")

def read_angles(path):
  """Reads angles data from a specified Lammps' data (or molecule) file.

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  """
  return _read_topology_components(path, "angle", "Angles")

def read_dihedrals(path):
  """Reads dihedrals data from a specified Lammps' data (or molecule) file.

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  """
  return _read_topology_components(path, "dihedral", "Dihedrals")

def read_impropers(path):
  """Reads impropers data from a specified Lammps' data (or molecule) file.

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  """
  return _read_topology_components(path, "improper", "Impropers")