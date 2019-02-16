#!/bin/bash

cd tests
python tests.py
cd ..

for ipynb in $(ls ./examples/*.ipynb)
do
  jupyter nbconvert --execute $ipynb
done
