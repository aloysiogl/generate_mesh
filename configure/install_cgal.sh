#!/bin/bash
sudo apt-get install libgmp3-dev libmpfr-dev swig

cd ~
FILE=~/CGAL-4.14.3.tar.xz
if [ -f "$FILE" ]; then
    echo "$FILE exist, installing directly..."
else 
    echo "$FILE does not exist, downloading..."
    wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.14.3/CGAL-4.14.3.tar.xz
fi

tar -C ~/ -xvf CGAL-4.14.3.tar.xz

cd CGAL-4.14.3

#Installing
mkdir build
cd build
cmake ..
sudo make install

#Removing install files
cd ~
rm -rf CGAL-4.14.3
rm -f CGAL-4.14.3.tar.xz
