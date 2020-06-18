# quickmerge


What is quickmerge?

The program uses complementary information from genomes assembled with long reads in order to improve contiguity, and works with assemblies derived from both Pacific Biosciences or Oxford Nanopore. Quickmerge will even work with hybrid assemblies made by combining long reads and Illumina short reads. This can be useful when long read coverage is limiting. For more details, please see <a href="http://nar.oxfordjournals.org/content/early/2016/07/25/nar.gkw654.full">the paper</a> that describes it. Citation for ONP assembly merging is <a href="https://www.g3journal.org/content/8/10/3143.abstract">here</a>     

Why use quickmerge?

 * Saves money. Illumina sequences are cheaper than PacBio or ONP long reads. So quickmerge allows you to cut your long molecule requirement by as much as half by replacing the same with Illumina short reads. E.g. if you think you would get a N50 of 8Mb from 75X long reads (ONP or PacBio), try sequencing 45X long and 70X Illumina reads instead of 75X long reads. You may not need that extra 35X long reads.
 * It is fast. Takes less than a minute to run on most genomes. You run nucmer once (nucmer is the most time consuming step) and then you can run quickmerge over a large number of parameters in a very short time.
 * Requires only fasta files and does not depend on any special data or computational resources.
 
The package contains all necessary components to run quickmerge. We also provide a set of test data (currrently available on request) so that you can check that the program is working correctly in your computer. Please send questions and comments to mchakrab@uci.edu


1. DOWNLOAD

   To download the latest version of quickmerge and <a href = "http://mummer.sourceforge.net/">MUMmer</a>, its primary dependency, you can clone the repository using 
   ```
    git clone
   ```
   
2. INSTALL:

CONDA

[![Anaconda-Server Badge](https://anaconda.org/bioconda/quickmerge/badges/latest_release_date.svg)](https://anaconda.org/bioconda/quickmerge)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/quickmerge/badges/platforms.svg)](https://anaconda.org/bioconda/quickmerge)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/quickmerge/badges/downloads.svg)](https://anaconda.org/bioconda/quickmerge)

To install quickmerge, first install [conda](https://conda.io/en/latest/miniconda.html) then run:

```
conda install -c conda-forge -c bioconda quickmerge
```


   UNIX:
   To install on a unix-based system, enter the following into the command line from the directory that this readme originated from:
   ```
	bash make_merger.sh
   ```
   This will compile 'quickmerge' and MUMMer. Requires GNU c++ compiler.

  NON-UNIX:

   On a non-unix system, you will have to manually compile these two programs, like so:

   first, enter the 'merger' directory and enter the following command to make the merger program:
   ```
 	make
   ```
   then, enter the 'MUMmer3.23' directory and enter the following commands, as specified in the MUMmer readme:
   ```
	make check

	make install
   ```
3. RUNNING QUICKMERGE:
   WRAPPER:

   The simplest way to run 'merger' is to use the python wrapper 'merge_wrapper.py':
   ```
	merge_wrapper.py hybrid_assembly.fasta self_assembly.fasta
   ```
   try the command 'merge_wrapper.py -h' for detail on options available with this wrapper.

   Note that, if you use the wrapper with MUMmer version 4, the `-v` flag is
   required. Multithreading of MUMmer is available using the `-t` tag if using
   MUMmer v4 or above.

   MANUAL:

   To manually run 'merger', first make a call to 'nucmer'.  Nucmer aligns the two assemblies so that the merger can find the correct splice sites:
   ```
	nucmer -l 100 -prefix out  self_assembly.fasta hybrid_assembly.fasta
   ```
   Then, use delta-filter to filter out alignments due to repeats and duplicates. Use any length for filtering out alignments between non-syntenic repeats, but higher is better (here we show 10kb) :
   ```   
	delta-filter -r -q -l 10000 out.delta > out.rq.delta
   ```
   Finally, use 'quickmerge' to merge the two assemblies:
   ```
	quickmerge -d out.rq.delta -q hybrid_assembly.fasta -r self_assembly.fasta -hco 5.0 -c 1.5 -l n -ml m -p prefix
   ```
   Description of the parameters:
   
   -q: Hybrid assembly. (this can also be a PacBio or a ONP only assembly). see <a href ="https://github.com/mahulchak/quickmerge/wiki">quickmerge wiki</a> for details
   
   -r: Self assembly. (can also be a hybrid assembly).
   
   -hco: controls the overlap cutoff used in selection of anchor contigs. Default is 5.0. 

   -c: controls the overlap cutoff for contigs used for extension of the anchor contig. Default is 1.5.

   For both "hco" and "c", bigger the number, more stringent is the criteria for contig selection (which will lead to fewer contigs being merged). If they are too small (<1), chances of spurious merging will increase. It is better to be conservative while merging contigs!

   -l: controls the length cutoff for anchor contigs. A good rule of thumb is to start with the N50 of the self assembly. E.g. if the N50 of your self assembly is 2Mb then use 2000000 as your cutoff. Lowering this value may lead to more merging but may increase the probability of mis-joins.
   
   -ml: controls the minimum alignment length to be considered for merging. This is especially helpful for repeat-rich genomes. Default is 5000 but higher values (>5000) are recommended.
   
   -p: A prefix that is added to the output from the run.

4. SOME HELPFUL TIPS:

  * Although this program was written to merge a hybrid assembly (e.g. as generated by <a href="https://sites.google.com/site/dbg2olc/">DBG2OLC</a>) and a PacBio or ONP only assembly, it can also be used to merge two different long molecule only assemblies (e.g. one generated with <a href="https://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/">PBcR</a> or <a href="https://github.com/marbl/canu">canu</a> and another generated with <a href="https://github.com/PacificBiosciences/FALCON-integrate">FALCON</a>).
 
  * For optimal merging results, identify the major misassemblies (especially translocations and inversions) in the component assemblies and break the contigs at such misassembly boundaries. Alignment of the component assemblies to the merged assembly may help to identify such assembly errors because a specific error typically occurs in only one of the assemblies.  
  
  * The fasta sequence headers should not have white spaces in them. In case they do, as might happen for assemblies obtained from  <a href="https://github.com/PacificBiosciences/FALCON-integrate">FALCON</a> assembler, the white space needs to be removed before launching the merging python script or before running mummer. Our merging script now takes care of this issue.  

  * You can run Ka-kit's <a href="https://github.com/kakitone/finishingTool">finisherSC</a> after running quickmerge to improve the contiguity even further.

  * Assembly polishing with <a href="https://github.com/PacificBiosciences/GenomicConsensus">Quiver</a> and <a href="https://github.com/broadinstitute/pilon/wiki">pilon</a> before and after assembly merging is strongly recommended. However, if you are running finisherSC, you may perform the quiver polishing after the finisher step.

  * Check the merged assembly by aligning the hybrid and/or non-hybrid assembly to the merged assembly (you can use nucmer -mumreference and mummerplot for alignment and dot plot visualization). 

5. KNOWN ISSUES:

 * All sequence names (that is, fasta headers) in each of the two input assemblies must be unique, i.e., each assembly must have all unique names, and the two assemblies must not share any names.  The quickmerge python wrapper now automatically fix this.

