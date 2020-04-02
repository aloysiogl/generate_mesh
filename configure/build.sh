#! /bin/bash

cd ../

#Clearing old build files
rm -rf build
rm -rf build-python
rm -rf tmp

#Building project
python setup.py build_ext --inplace

#Clearing build directories and leaving only the important files inside build
rm -rf build
rm -rf tmp
mkdir build
cp ./build-python/CGAL/* ./build
rm -rf build-python

