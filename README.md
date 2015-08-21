# quickmerge
This package contains all necessary components to run the assembly merger program quickmerge.

######################

INSTALL:

UNIX:

to install on a unix-based system, enter the following into the command line from the directory that this readme originated from:
bash make_merger.sh

This will compile 'quickmerge' and 'MUMmer', its primary dependency.

NON-UNIX:

On a non-unix system, you will have to manually compile these two programs, like so:

first, enter the 'merger' directory and enter the following command to make the merger program:
make

then, enter the 'MUMmer3.23' directory and enter the following commands, as specified in the MUMmer readme:
make check
make install

#####################

RUNNING MERGER:

WRAPPER:

The simplest way to run 'merger' is to use the python wrapper 'merge_wrapper_v2.py':
merge_wrapper_v2.py hybrid_assembly.fasta self_assembly.fasta

try the command 'merge_wrapper_v2.py -h' for detail on options available with this wrapper.

MANUAL:

To manually run 'merger', first make a call to 'nucmer'.  Nucmer aligns the two assemblies so that the merger can find the correct splice sites:

nucmer -l  100 -prefix out  self_assembly.fasta hybrid_assembly.fasta

Then, use delta-filter to filter out alignments due to repeats and duplicates:

delta-filter -i 95 -r -q out.delta > out.rq.delta

Finally, use 'quickmerge' to merge the two assemblies (note: the order of the self and hybrid assembly is important:

quickmerge -d out.rq.delta -q hybrid_assembly.fasta -r self_assembly.fasta -hco 5 -c 1.5

