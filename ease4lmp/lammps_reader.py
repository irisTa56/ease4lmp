"""Submodule for functions to read Lammps' data file."""

from .lammps_dataformats import lmp_datanames


def _read_section(path, section):
  """Read lines of a specified section in a specified file.

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  section: str
    Header of a section to be read.

  """
  lines = (line.lstrip().rstrip() for line in open(path))

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
  """Read data of specified topology components
  from a specified Lammps' data (or molecule) file.

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  name: str
    Name of topology component.

  header: str
    Section header of the topology components.

  """
  datanames = ["{}-id".format(s) for s in lmp_datanames[name]]

  return [
    dict(zip(datanames, map(int, line)))
    for line in _read_section(path, header)
  ]

def _str2num(s):
  """Converts a string to a number checking whether it is int or not"""
  return int(s) if s.lstrip('-').isdigit() else float(s)

#-----------------------------------------------------------------------

def read_box(path):
  """Read side lengths of the simulation box.

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  """
  lines = (line.lstrip().rstrip() for line in open(path))

  lx, ly, lz = (None,) * 3
  done = False

  for line in lines:
    if not done:
      if line.endswith("xlo xhi"):
        tmp = line.split()
        lx = float(tmp[1]) - float(tmp[0])
      elif line.endswith("ylo yhi"):
        tmp = line.split()
        ly = float(tmp[1]) - float(tmp[0])
      elif line.endswith("zlo zhi"):
        tmp = line.split()
        lz = float(tmp[1]) - float(tmp[0])

      if all(l is not None for l in [lx, ly, lz]):
        done = True

    if line.endswith("xy xz yz"):
      raise RuntimeError("Skewed system is not supported")

  return lx, ly, lz

def read_bonds(path):
  """Read bonds data from a Lammps' data (or molecule) file.

  Returned value is a JSON object (list of dict).

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  """
  return _read_topology_components(path, "bond", "Bonds")

def read_angles(path):
  """Read angles data from a Lammps' data (or molecule) file.

  Returned value is a JSON object (list of dict).

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  """
  return _read_topology_components(path, "angle", "Angles")

def read_dihedrals(path):
  """Read dihedrals data from a Lammps' data (or molecule) file.

  Returned value is a JSON object (list of dict).

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  """
  return _read_topology_components(path, "dihedral", "Dihedrals")

def read_impropers(path):
  """Read impropers data from a Lammps' data (or molecule) file.

  Returned value is a JSON object (list of dict).

  Parameters:

  path: str
    File path to Lammps' data file (or molecule file).

  """
  return _read_topology_components(path, "improper", "Impropers")

def read_atoms_from_data(path, atom_style, mass=True, velocity=True):
  """Read atoms data from a Lammps' data file.

  Returned value is a JSON object (list of dict).

  Parameters:

  path: str
    File path to Lammps' data file.

  atom_style: str
    Specifies an *atom style* used in Lammps.

  mass: bool
    Whether to include mass data if *Masses* section extis.

  velocity: bool
    Whether to include velocity data if *Velocities* section extis.

  """
  atoms = []
  lx, ly, lz = read_box(path)

  datanames = tuple(
    n if n not in {"x", "y", "z"} else n + "u"
    for n in lmp_datanames["atom"][atom_style])

  lines = _read_section(path, "Atoms")

  for line in lines:
    if len(datanames) == len(line):
      atoms.append(dict(zip(datanames, map(_str2num, line))))
    elif len(datanames) + 3 == len(line):
      tmp = dict(zip(datanames, map(_str2num, line[:-3])))
      ix, iy, iz = map(int, line[-3:])
      tmp["xu"] += ix * lx
      tmp["yu"] += iy * ly
      tmp["zu"] += iz * lz
      atoms.append(tmp)
    else:
      raise RuntimeError("Invalid number of values in a line")

  if mass:

    type2mass = {
      int(s1): float(s2) for s1, s2 in _read_section(path, "Masses")
    }

    for atom in atoms:
      atom["mass"] = type2mass[atom["type"]]

  if velocity:

    datanames_vel = lmp_datanames["velocity"].get(
      atom_style, lmp_datanames["velocity"]["*"])

    lines_vel = _read_section(path, "Velocities")

    if len(lines_vel) == len(atoms):
      for atom, line in zip(atoms, lines_vel):
        atom.update(dict(zip(datanames_vel, map(_str2num, line))))

  return atoms

def read_atoms_from_molecule(path):
  """Read atoms data from a Lammps' molecule file.

  Returned value is a JSON object (list of dict).

  Parameters:

  path: str
    File path to Lammps' molecule file.

  """
  atoms = [
    dict(zip(("id", "xu", "yu", "zu"), map(_str2num, line)))
    for line in _read_section(path, "Coords")
  ]

  lines_type = _read_section(path, "Types")

  if len(lines_type) == len(atoms):
    for atom, line in zip(atoms, lines_type):
      atom["type"] = int(line[1])

  lines_q = _read_section(path, "Charges")

  if len(lines_q) == len(atoms):
    for atom, line in zip(atoms, lines_q):
      atom["q"] = float(line[1])

  lines_mass = _read_section(path, "Masses")

  if len(lines_mass) == len(atoms):
    for atom, line in zip(atoms, lines_mass):
      atom["mass"] = float(line[1])

  return atoms
