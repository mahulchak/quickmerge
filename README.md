# quickmerge
What is quickmerge?

quickmerge uses a simple concept to improve contiguity of genome assemblies based on long molecule sequences, often with dramatic outcomes. The program uses information from assemblies made with illumina short reads and PacBio long reads to improve contiguities of an assembly generated with PacBio long reads alone. This is counterintuitive because illumina short reads are not typically considered to cover genomic regions which PacBio long reads cannot. Although we have not evaluated this program for assemblies generated with Oxford nanopore sequences, the program should work with ONP-assemblies too. For more details, please see <a href="http://biorxiv.org/content/early/2015/10/16/029306">the paper</a> that describes it.    

Why use quickmerge?

 * Saves money. Illumina sequences are much cheaper than PacBio or ONP long reads. So quickmerge allows you to cut your long molecule requirement by half (or more) by replacing the same with Illumina short reads. E.g. if you think you would get a N50 of 8Mb from 75X PacBio reads, try sequencing 40X PacBio and 70X Illumina reads instead of 75X PacBio reads. You may not need that extra 35X PacBio reads.
 * It superfast. Takes less than a minute to run on most genomes. You run nucmer once (nucmer is the most time consuming step) and then you can run quickmerge over a large number of parameters in a very short time.
 * Requires only fasta files and does not depend on any special data or computational resources.
 
The package contains all necessary components to run quickmerge. We also provide a set of test data (currrently available on request) so that you can check that the program is working correctly in your computer. Please send questions and comments to mchakrab@uci.edu


1. DOWNLOAD

   To download the latest version of quickmerge and <a href "https://sourceforge.net/projects/mummer/files/">MUMmer</a>, its primary dependency, you can clone the repository using 
   ```
    git clone
   ```
   Alternatively, you can download it as a zip file.

2. INSTALL:

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

   The simplest way to run 'merger' is to use the python wrapper 'merge_wrapper_v2.py':
   ```
	merge_wrapper_v2.py hybrid_assembly.fasta self_assembly.fasta
   ```
   try the command 'merge_wrapper_v2.py -h' for detail on options available with this wrapper.

   MANUAL:

   To manually run 'merger', first make a call to 'nucmer'.  Nucmer aligns the two assemblies so that the merger can find the correct splice sites:
   ```
	nucmer -l  100 -prefix out  self_assembly.fasta hybrid_assembly.fasta
   ```
   Then, use delta-filter to filter out alignments due to repeats and duplicates:
   ```   
	delta-filter -i 95 -r -q out.delta > out.rq.delta
   ```
   Finally, use 'quickmerge' to merge the two assemblies (note: the order of the self and hybrid assembly is important:
   ```
	quickmerge -d out.rq.delta -q hybrid_assembly.fasta -r self_assembly.fasta -hco 5 -c 1.5 -l n
   ```
   Description of the parameters:

   -hco: controls the overlap cutoff used in selection of anchor contigs. Default is 5.0. 

   -c: controls the overlap cutoff for contigs used for extension of the anchor contig. Default is 1.5.

   For both "hco" and "c", bigger the number, more stringent is the criteria for contig selection (which will lead to fewer contigs being merged). If they are too small (<1), chances of spurious merging will increase.

   -l: controls the length cutoff for anchor contigs. A good rule of thumb is to start with the N50 of the 'self_assembly.fasta'. E.g. if the N50 of your self_assembly.fasta is 2Mb. Then use 2000000 as your cutoff. Lowering this value will lead to more merging but may increase the probability of mis-joins. 

4. SOME HELPFUL TIPS:

 * Although this program was written to merge a hybrid assembly and a PB-only assembly, it can also be used to merge two different PB-only assemblies (e.g. one generated with <a href="https://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/">PBcR</a> and another generated with <a href="https://github.com/PacificBiosciences/FALCON-integrate">FALCON</a>).

  * If using a self-assembly as the reference did not improve your contiguity much (unlikely for low coverage PacBio sequence but possible for high coverage PacBio sequences), use the hybrid assembly as your reference and the self-assembly as your query.

  * The fasta sequence headers should not have white spaces in them. In case they do, as might happen for assemblies obtained from  <a href="https://github.com/PacificBiosciences/FALCON-integrate">FALCON</a> assembler, the white space needs to be removed before launcing the merging python script or before running mummer. Also, there cannot be any line breaks in the fasta sequence. We are currently in the process of adding additional tools to our merging pipeline which will take care of these format issues.  

  * You can run Ka-kit's <a href="https://github.com/kakitone/finishingTool">finisherSC</a> after running quickmerge to improve the contiguity even further.

  * Assembly polishing with <a href="https://github.com/PacificBiosciences/GenomicConsensus">Quiver</a> before and after assembly merging is strongly recommended. However, if you are running finisherSC, you may perform the quiver polishing after the finisher step.

  * Check the merged assembly by aligning the hybrid and/or PB only assembly to the merged assembly (you can use nucmer -mumreference and mummerplot for alignment and dot plot visualization).


5. UPDATES

  i) If you use illumina sequences to improve contiguity of your assembly, and used <a href="https://sites.google.com/site/dbg2olc/">DBG2OLC</a> (recommended) to generate the hybrid assembly, you no longer have to run the time consuming consensus calling step for DBG2OLC. Instead, obtain the backbone_raw.fasta file after running DBG2OLC and use that as your input for quickmerge. If you decide to go this route, please run nucmer as 
  ```
  	nucmer -prefix hyb2pb pb.fasta backbone_raw.fasta
  	delta-filter -r -q hyb2pb.delta > hyb2pb.rq.delta
  ```
  then run quickmerge as mentioned above.
