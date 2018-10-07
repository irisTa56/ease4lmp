#!/usr/bin/env python

import subprocess

from glob import glob

# pre-process

paths = glob("ease4lmp/*.py")
original_sources = {}

for path in paths:

  with open(path, "r") as f:
    data = f.read()

  original_sources[path] = data

  data_parts = data.split("\"\"\"")

  new_data = data_parts[0]

  for i, part in enumerate(data_parts[1:]):
    new_data += "\"\"\"" if i%2 else "\"\"\"!"
    new_data += part

  with open(path, "w") as f:
    f.write(new_data)

# main-process

p = subprocess.run(
  ['doxygen', 'Doxyfile'],
  stdout=subprocess.PIPE, stderr=subprocess.PIPE)

for line in p.stderr.decode().split("\n"):
  print(line)

# post-process

for path in paths:
  with open(path, "w") as f:
    f.write(original_sources[path])
