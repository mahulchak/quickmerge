#!/bin/bash
#make if statement for 32 or 64
cd merger
make clean
make
cd ..
mergerpath=`pwd`
cd MUMmer3.23
mkdir aux_bin
make clean
if [ $1 = "64" ] || [ $# -eq 0 ]; then
   make CPPFLAGS="-O3 -DSIXTYFOURBITS"
elif [ $1 = "32" ]; then
   make
else
   echo "please state either 32 or 64 bit compilation. i.e. bash make_merger.sh 64"
   exit 1
fi
make check
make install
mummerpath=`pwd`
cd ..

echo ""
echo ""
echo ""
echo ""
echo ""
echo ""
echo "##############################################"
echo "Please add to your path by:"
echo "export PATH=$mergerpath:$mummerpath:\$PATH"
echo ""
echo "or by:"
echo "source $mergerpath/.quickmergerc"
echo ""
echo "Please paste the follwing in your cluster submission script prior to running the quickmerge wrapper"
echo "export PATH=$mergerpath:$mummerpath:\$PATH" > .quickmergerc
