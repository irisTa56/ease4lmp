"""
This script is for Lammps" Atoms section styles used in DATA file.
Please see the following for details.
* http://lammps.sandia.gov/doc/read_data.html

create: 2018/05/30 by Takayuki Kobayashi
"""

lmp_dataformats_atoms = {
  "angle": ("id", "mol", "type", "x", "y", "z"),
  "atomic": ("id", "type", "x", "y", "z"),
  "body": ("id", "type", "bodyflag", "mass", "x", "y", "z"),
  "bond": ("id", "mol", "type", "x", "y", "z"),
  "charge": ("id", "type", "q", "x", "y", "z"),
  "dipole": ("id", "type", "q", "x", "y", "z", "mux", "muy", "muz"),
  "dpd": ("id", "type", "theta", "x", "y", "z"),
  "edpd": ("id", "type", "edpd_temp", "edpd_cv", "x", "y", "z"),
  "mdpd": ("id", "type", "rho", "x", "y", "z"),
  #"tdpd": ("id", "type", "x", "y", "z", "cc1", "cc2", "…", "ccNspecies"),
  "electron": ("id", "type", "q", "spin", "eradius", "x", "y", "z"),
  "ellipsoid": ("id", "type", "ellipsoidflag", "density", "x", "y", "z"),
  "full": ("id", "mol", "type", "q", "x", "y", "z"),
  "line": ("id", "mol", "type", "lineflag", "density", "x", "y", "z"),
  "meso": ("id", "type", "rho", "e", "cv", "x", "y", "z"),
  "molecular": ("id", "mol", "type", "x", "y", "z"),
  "peri": ("id", "type", "volume", "density", "x", "y", "z"),
  "smd": ("id", "type", "molecule", "volume", "mass", "kernel_radius", "contact_radius", "x", "y", "z"),
  "sphere": ("id", "type", "diameter", "density", "x", "y", "z"),
  "template": ("id", "mol", "template_index", "template_atom", "type", "x", "y", "z"),
  "tri": ("id", "mol", "type", "triangleflag", "density", "x", "y", "z"),
  "wavepacket": ("id", "type", "charge", "spin", "eradius", "etag", "cs_re", "cs_im", "x", "y", "z"),
  #"hybrid": ("id", "type", "x", "y", "z", "sub_style1", "sub_style2", "…}"
}

lmp_dataformats_velocities = {
  "*": ("id", "vx", "vy", "vz"),
  "electron": ("id", "vx", "vy", "vz", "ervel"),
  "ellipsoid": ("id", "vx", "vy", "vz", "lx", "ly", "lz"),
  "sphere": ("id", "vx", "vy", "vz", "wx", "wy", "wz"),
  #"hybrid": ("id", "vx", "vy", "vz", "sub-style1", "sub-style2", "…}"
}
