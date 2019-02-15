import os

from setuptools import setup, find_packages

# read version & description

version_ns = {}
with open(os.path.join("ease4lmp", "_version.py")) as f:
  exec(f.read(), {}, version_ns)

with open("README.md", "r") as f:
  long_description = f.read()

# setup

setup(
  name="ease4lmp",
  version=version_ns["__version__"],
  description="Extension of Atomic Simulation Environment for LAMMPS",
  author="Takayuki Kobayashi",
  author_email="iris.takayuki@gmail.com",
  license="MIT",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/irisTa56/ease4lmp.git",
  packages=find_packages(),
  classifiers=(
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License"))
