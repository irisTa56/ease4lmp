#!/bin/bash

cd tests
python test.py
cd ..

for ipynb in $(ls ./examples/*.ipynb)
do
  jupyter nbconvert --execute $ipynb
done
