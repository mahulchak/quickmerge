#!/bin/python
import subprocess
import sys
import argparse

#Required:
#the 'nucmer','delta-filter', and 'merger' executables must be in the user's PATH

parser = argparse.ArgumentParser(description = 'run mummer and the merge program.')
parser.add_argument("hybrid_assembly_fasta", help="the output of a hybrid assembly program such as DBG2OLC")
parser.add_argument("self_assembly_fasta",help="the output of a self assembly program such as PBcR")
parser.add_argument("-pre","--prefix", help="the prefix for all output files")
parser.add_argument("-hco","--hco", help="the quickmerge hco parameter (default=5.0)")
parser.add_argument("-c","--c", help="the quickmerge c parameter (default=1.5)")
parser.add_argument("-l","--length_cutoff", help="minimum contig length to be merged (default=0)")
parser.add_argument("--no_nucmer", help="skip the nucmer step",action="store_true")
parser.add_argument("--no_delta", help="skip the nucmer and delta-filter steps",action="store_true")
parser.add_argument("--stop_after_nucmer", help="do not perform the delta-filter and merger steps",action="store_true")
parser.add_argument("--stop_after_df", help="do not perform the merger step",action="store_true")

args=parser.parse_args()
if args.prefix:
  prefix = args.prefix
else:
  prefix = "out"
hypath = args.hybrid_assembly_fasta
selfpath = args.self_assembly_fasta
if args.hco:
  hco = args.hco
else:
  hco = 5.0
if args.c:
  c = args.c
else:
  c = 1.5
if args.length_cutoff:
  length_cutoff = args.length_cutoff
else:
  length_cutoff = 0

#check for fasta compatibility:
testpaths = [str(hypath),str(selfpath)]
order = ["hybrid","self"]
for iteration in range(0,2):
  file = testpaths[iteration]
  with open(file,"r") as tfile:
    i = 0
    ok = True
    for linea in tfile:
      line = linea.rstrip('\n')
      if i % 2 == 0 and line[0] != ">": ok = False
      if " " in line: ok = False
      i += 1
    
    #if fasta isn't compatible, make oneline version of fasta:
    if not ok:
      tfile.seek(0)
      tempoutpath = order[iteration] + "_oneline.fa"
      with open(tempoutpath,"w") as tout:
        header = ""
        seq = ""
        for linea in tfile:
          line = linea.rstrip('\n')
          print line # debug line
          if line[0] == ">":
            print line # debug line
            if len(header) > 0 and len(seq) > 0:
              tout.write(header + "\n" + seq + "\n")
            header = line
            seq = ""
          else:
            seq = seq + line
        tout.write(header + "\n" + seq + "\n")
      #change hypath and selfpath to use the temporary (oneline) fastas:
      if iteration == 0:
        hypath = tempoutpath
      else:
        selfpath = tempoutpath
      

#run nucmer:
if not args.no_nucmer and not args.no_delta:
  subprocess.call(['nucmer','-l','100','-prefix',str(prefix),str(selfpath),str(hypath)])

#run the delta filter on the nucmer alignment:
if not args.no_delta and not args.stop_after_nucmer:
  with open(str(str(prefix)+'.rq.delta'), "w") as outfile:
    subprocess.call(['delta-filter','-i','95','-r','-q',str(str(prefix)+'.delta')],stdout=outfile)

#append the correct options to the merger call:
mergercall = ['quickmerge','-d',str(str(prefix)+'.rq.delta'),'-q',str(hypath),'-r',str(selfpath)]
#if args.hco:
mergercall.append('-hco')
mergercall.append(str(hco))
#if args.c:
mergercall.append('-c')
mergercall.append(str(c))
#if args.length_cutoff:
mergercall.append('-l')
mergercall.append(str(length_cutoff))

#run the merging program
if not args.stop_after_nucmer and not args.stop_after_df:
  subprocess.call(mergercall)

