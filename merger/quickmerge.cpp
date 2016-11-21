//This program takes two different assemblies in fasta format and uses the delta output from mummer to generate a merged asembly
//The program does not correct any assembly errors. It only joins contigs. So if a contig was misassembled  by the main assembly program,
//it will most likely remain misassembled. Currently, the program is intended to use when 30X or above pacbio data is available, such that at least a 
//respectable self-assembly is available. 
//Copyright (C) 2015  Mahul Chakraborty 
//Please send questions and bug reports to mchakrab@uci.edu

#include<iostream>
#include<fstream>
#include<cstdlib>
#include "qmerge.h"

using namespace std;


int main(int argc, char * argv[])
{
        if((argc<14) || (argc==1))
        {
		cerr<<"All options were not supplied :("<<endl;
		cerr<<"Usage: "<<argv[0]<<" -d delta_file.out -q hybrid.fasta -r self.fasta -hco (default=5.0) -c (default=1.5) -l seed_length_cutoff -ml merging_length_cutoff "<<endl;
       		exit(EXIT_FAILURE);
        }

	ifstream fin,hyb,pb;
	ofstream fout;

	fout.open("merged.fasta");

	asmMerge merge,merge1; 
	fastaSeq hybrid,pbOnly,merged;

	string header,ref_name,qu_name,tempname,name;
	double hco,cutoff;
	hco = 5.0;
	cutoff = 1.5;
	int qu_st = 0;
	int qu_end = 0;
	int r_st = 0;
	int r_end = 0;
	const int length = atoi(argv[12]);
	const int absLenCutoff = atoi(argv[14]);
	if(*argv[8])
	{
		hco = strtod(argv[8],NULL);
	}
	if(*argv[10])
	{
		cutoff = strtod(argv[10],NULL);	
	}

	fin.open(argv[2]);

	while(getline(fin,header))
	{
		if(header[0] =='>')
             	   {
		
            		ref_name = xtractcol(header,' ',1);	
			ref_name = ref_name.substr(1);
			merge.r_name.push_back(ref_name); 
			qu_name = xtractcol(header,' ',2); 
			merge.q_name.push_back(qu_name);
			tempname = ref_name.append(qu_name); // tempname is the index for the map. they describe alignment pairs( e.g. BackboneX ctgY). should be unique. 
			merge.ref_len[tempname] = atoi(xtractcol(header,' ',3).c_str()); 
			merge.q_len[tempname] = atoi(xtractcol(header,' ',4).c_str()); 
		}
		
		if(header[0] != '>' && header.size()>10)
		{
			r_st = atoi(xtractcol(header,' ',1).c_str());
			merge.ref_st[tempname].push_back(r_st); // storing the coordinates for each alignment hit
			r_end = atoi(xtractcol(header,' ',2).c_str());
			merge.ref_end[tempname].push_back(r_end);
			qu_st = atoi(xtractcol(header,' ',3).c_str());
			merge.q_st[tempname].push_back(qu_st);
			qu_end = atoi(xtractcol(header,' ',4).c_str());
			merge.q_end[tempname].push_back(qu_end);
		}
	}
	fin.close();

	writeToFile(merge);
	ovlStoreCalculator(merge);
	innieChecker(merge);
	sideChecker(merge);
	sideCheckerR(merge);
	sideCheckerQ(merge);
	assignStrand(merge);
	//ovlStoreCalculator(merge);
	nOvlStoreCalculator(merge);
	ovrHngCal(merge);
	overHangSideR(merge);
	writeSummary(merge);
	
	if(*argv[4])
	{
		hyb.open(argv[4]);
	}
	if(*argv[6])
	{
	pb.open(argv[6]);
	}
	fillSeq(hybrid,hyb,' ');
	fillSeq(pbOnly,pb);
	fillAnchor(merge,merge1,hco,cutoff,length,absLenCutoff);
	writeAnchorSummary(merge);


	findChain(merge,merge1,pbOnly,merged,cutoff);
	createMseq(merge,merge1);

	fillOri(merge,merge1);	
	for(map<string,vector<string> >::iterator it=merge1.lseq.begin();it != merge1.lseq.end();it++)
	{
		name = it->first;
		cout<<name<<"\t";
		for(unsigned int j = 0; j<merge1.lseq[name].size();j++)
		{
			cout<<merge1.lseq[name][j]<<"\t"<<merge.Ori[name][j]<<"\t";
		}
		cout<<endl;
	}

	removeSeq(merge,hybrid,merged); // copy unaligned hybrid contigs to merged 
	ctgJoiner(merge,merge1,hybrid,pbOnly,merged);
	writeMerged(merged,fout);
	fout.close();
	return 0;
}
