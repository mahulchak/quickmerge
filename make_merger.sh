#!/bin/bash
cd merger
make
mergerpath=`pwd`
cd ../MUMmer3.23
make check
make install
mummerpath=`pwd`
cd ..
export PATH=$mergerpath:$mummerpath:$PATH


