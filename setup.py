import setuptools as st

with open("README.md", "r") as fh:
    long_description = fh.read()

st.setup(
    name="ease4lmp",
    version="0.0.2a",
    author="Takayuki Kobayashi",
    author_email="iris.takayuki@gmail.com",
    description="Extension of Atomic Simulation Environment for LAMMPS",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/irisTa56/ease4lmp",
    packages=st.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    install_requires=["ase>=3.16.0", "numpy>=1.14.3"],
)