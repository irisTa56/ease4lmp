import unittest

import os
import sys
sys.path.insert(0, os.path.abspath("."))

from io import StringIO
io = StringIO()
sys.stdout = io

from test_bonded_atoms import suite as suite_bonded_atoms
from test_lammps_reader import suite as suite_lammps_reader
from test_lammps_cycle import suite as suite_lammps_cycle

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite_bonded_atoms())
    runner.run(suite_lammps_reader())
    runner.run(suite_lammps_cycle())
