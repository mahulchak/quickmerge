
#ifndef MASM_H_
#define MASM_H_

#include<vector>
#include<string>
#include<map>
#include<fstream>
#include<algorithm>

using namespace std;


class asmMerge {

public:

vector<string> r_name;
vector<string> q_name;
map<string,int> ref_len;
map<string,int> q_len;
map<string,vector<int> > ref_st;
map<string,vector<int> > ref_end;
map<string,vector<int> > q_st;
map<string,vector<int> > q_end;
map<string,bool> innie;
map<string,char> sideInfo;
map<string,char> sideInfoQ;
map<string,char> sideInfoR;
map<string,int> ovlStore;
map<string,int> nOvlStore;
map<string,int> ovrHangQ;
map<string,int> ovrHangR;
map<string,vector<int> >ovrHangQR;
map<string,vector<string> > anchor;
map<string,vector<string> > cAnchor;
map<string,vector<string> > lseq; 
map<string,vector<string> > rseq;
map<string,vector<int> > Ori;
map<string,bool> qToRemove; 
map<string,bool> strandInfo;
map<string,vector<int> > newEnd;
};

class fastaSeq{
public:
vector<string> seqName;
map<string,string> seq;
};

string xtractcol(string str,char c, int n);
void writeToFile(asmMerge & merge);
void innieChecker(asmMerge & merge);
void sideChecker(asmMerge & merge);
void sideCheckerQ(asmMerge & merge);
void sideCheckerR(asmMerge & merge);
void ovlStoreCalculator(asmMerge & merge);
int ovlCalculator(vector<int>& q_st, vector<int>& q_end);
void nOvlStoreCalculator(asmMerge & merge);
void writeSummary(asmMerge & merge);
string revCom(string & str);
void ovrHngCal(asmMerge & merge);
void fillAnchor(asmMerge & merge, asmMerge & merge1, double propAnchor, double propCutoff,const int & length,const int & absLenCutoff);
void fillSeq(fastaSeq & fasta, ifstream& fin, char c);
void fillSeq(fastaSeq & fasta, ifstream& fin);
void writeAnchorSummary(asmMerge & merge);
void findChain(asmMerge & merge, asmMerge & merge1,fastaSeq & pbOnly, fastaSeq & merged,double propCutoff);
vector<string> vfind(string tempname,vector<string>& temp_rname, vector<string>& temp_qname,asmMerge & merge,string & guruRef,double propCutoff);
string longestLeft(string tempname, vector<string>& seq,asmMerge & merge, char RorQ, char sideQ, string & prevElem);
string longestRt(string tempname, vector<string>& seq,asmMerge & merge, char RorQ, char sideQ, string & prevElem);
void createMseq(asmMerge & merge, asmMerge & merge1);
void fillOri(asmMerge & merge, asmMerge & merge1);
void fillToRemove(asmMerge & merge, vector<string> & myvector);
void fillToRemove(asmMerge & merge, string & str);
void removeSeq(asmMerge & merge, fastaSeq & hyb, fastaSeq & merged);
void writeMerged(fastaSeq & merged, ofstream & fout);
void ctgJoiner(asmMerge & merge,asmMerge & merge1,fastaSeq & hybrid, fastaSeq& pbObly, fastaSeq & merged);
int maxIntD(int & qf1,int & qe1, int & qf2, int & qe2);
int mapQonRef(string & tempname1,string & tempname2,asmMerge & merge);
bool chkOvl(int f1,int e1,int f2,int e2);
int max(int & a, int & b);
int min(int & a, int & b);
vector<int> minD(int & qf1,int & qe1, int & qf2, int & qe2);
vector<int> maxD(int & qf1,int & qe1, int & qf2, int & qe2);
string reversed(string & str);
int returnIndex (vector<string> & myvector, string & str);
void assignStrand(asmMerge & merge);
void discAnchor(string & guruQ, asmMerge & merge,string & guruRef,double propCutoff);
void lisCalculator(asmMerge & merge,string & tempname,vector<int> &ref_st,vector<int> & ref_end,vector<int> & q_st,vector<int> & q_end);
#endif
