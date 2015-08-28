# quickmerge

This package contains all necessary components to run the assembly merger program quickmerge. Please send questions and comments to mchakrab@uci.edu

######################

INSTALL:

UNIX:

to install on a unix-based system, enter the following into the command line from the directory that this readme originated from:

	bash make_merger.sh

This will compile 'quickmerge' and 'MUMmer', its primary dependency. Also requires GNU c++ compiler.

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

	quickmerge -d out.rq.delta -q hybrid_assembly.fasta -r self_assembly.fasta -hco 5 -c 1.5 -l n

-hco: controls the overlap cutoff for anchor contigs. Bigger the number, more stringent the criteria, fewer overlaps are selected.

-c: controls the overlap cutoff for contigs used for extension of the anchor contig.

For both "hco" and "c", bigger the number, more stringent is the criteria for contig selection (which will lead to fewer contig merging). If they are too small (<1), spurious overlaps will be used for contig merging.

-l: controls the length cutoff for anchor contigs. A good rule of thumb is to start with the N50 of the self_assembly.fasta. E.g. if the N50 of your self_assembly.fasta is 2Mb. Then use 2000000 as your cutoff. Lowering this value will lead to more merging but may increase the probability of mis-joins. 

Although this program was written to merge a hybrid assembly and a PB-only assembly, it can also be used to two different PB-only assemblies (e.g. one generated with FALCON and another generated with PBcR).
####################
If you find this program useful, please cite our paper.
