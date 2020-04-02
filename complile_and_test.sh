#!/bin/bash

rm -r build
rm -r build-python
rm -r tmp
clear
python setup.py build_ext --inplace
cd build-python
cd CGAL

cp ../../test_module.py ./

python test_module.py