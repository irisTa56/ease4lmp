"""Submodule containing format data for each *atom style*
used in *Atoms* and *Velocities* section of Lammps' data file.

Pease see also the following for details:

* http://lammps.sandia.gov/doc/read_data.html.

"""


# Dictionary from *atom styles* to tuples representing formats
# of each line in *Atoms* section.
lmp_dataformats_atoms = {
  "angle": ("id:8d", "mol:6d", "type:4d", "x:15.8e", "y:15.8e", "z:15.8e"),
  "atomic": ("id:8d", "type:4d", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"body": ("id:8d", "type:4d", "bodyflag", "mass", "x:15.8e", "y:15.8e", "z:15.8e"),
  "bond": ("id:8d", "mol:6d", "type:4d", "x:15.8e", "y:15.8e", "z:15.8e"),
  "charge": ("id:8d", "type:4d", "q:9.6f", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"dipole": ("id:8d", "type:4d", "q:9.6f", "x:15.8e", "y:15.8e", "z:15.8e", "mux", "muy", "muz"),
  #"dpd": ("id:8d", "type:4d", "theta", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"edpd": ("id:8d", "type:4d", "edpd_temp", "edpd_cv", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"mdpd": ("id:8d", "type:4d", "rho", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"tdpd": ("id:8d", "type:4d", "x:15.8e", "y:15.8e", "z:15.8e", "cc1", "cc2", "...", "ccNspecies"),
  #"electron": ("id:8d", "type:4d", "q:9.6f", "spin", "eradius", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"ellipsoid": ("id:8d", "type:4d", "ellipsoidflag", "density", "x:15.8e", "y:15.8e", "z:15.8e"),
  "full": ("id:8d", "mol:6d", "type:4d", "q:9.6f", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"line": ("id:8d", "mol:6d", "type:4d", "lineflag", "density", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"meso": ("id:8d", "type:4d", "rho", "e", "cv", "x:15.8e", "y:15.8e", "z:15.8e"),
  "molecular": ("id:8d", "mol:6d", "type:4d", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"peri": ("id:8d", "type:4d", "volume", "density", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"smd": ("id:8d", "type:4d", "molecule", "volume", "mass", "kernel_radius", "contact_radius", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"sphere": ("id:8d", "type:4d", "diameter", "density", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"template": ("id:8d", "mol:6d", "template_index", "template_atom", "type:4d", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"tri": ("id:8d", "mol:6d", "type:4d", "triangleflag", "density", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"wavepacket": ("id:8d", "type:4d", "charge", "spin", "eradius", "etag", "cs_re", "cs_im", "x:15.8e", "y:15.8e", "z:15.8e"),
  #"hybrid": ("id:8d", "type:4d", "x:15.8e", "y:15.8e", "z:15.8e", "sub_style1", "sub_style2", "...")
}

# dictionary from *atom styles* to tuples representing formats
# of each line in *Velocities* section.
lmp_dataformats_velocities = {
  "*": ("id:8d", "vx:15.8e", "vy:15.8e", "vz:15.8e"),
  #"electron": ("id:8d", "vx:15.8e", "vy:15.8e", "vz:15.8e", "ervel"),
  #"ellipsoid": ("id:8d", "vx:15.8e", "vy:15.8e", "vz:15.8e", "lx", "ly", "lz"),
  #"sphere": ("id:8d", "vx:15.8e", "vy:15.8e", "vz:15.8e", "wx", "wy", "wz"),
  #"hybrid": ("id:8d", "vx:15.8e", "vy:15.8e", "vz:15.8e", "sub-style1", "sub-style2", "...")
}
