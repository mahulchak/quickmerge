#!/bin/bash
cd merger
g++ -Wall work_in_prog_temp.cpp exp_testlib.cpp -o merger
mergerpath=`pwd`
cd ../MUMmer3.23
make check
make install
mummerpath=`pwd`
cd ..
export PATH=$mergerpath:$mummerpath:$PATH


