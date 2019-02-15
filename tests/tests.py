import unittest

import os
import sys
sys.path.insert(0, os.path.abspath("."))

from test_bonded_atoms import suite as suite_bonded_atoms

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite_bonded_atoms())
