#include<iostream>
#include "qmerge.h"

string xtractcol(string str, char c, int n)
{
int count =0;
int j =0; // j is the index for the string that'll be returned
string elem;
for(unsigned int i=0;i<str.size() && count<n;i++)
        {
        if(str[i] == c)
                {
                count++; //keep a count of the field separator
                }
        if(count == n-1)
                {
                elem.push_back(str[i]);
                j++;
                }
        }
return elem;
}
//////////////////////////////////////////////////////////////////////////////////

void writeToFile(asmMerge & merge)
	{
	ofstream fout;
	fout.open("aln_summary.tsv");
	fout<<"REF"<<'\t'<<"QUERY"<<'\t'<<"REF-LEN"<<'\t'<<"Q-LEN"<<'\t'<<"REF-ST"<<'\t'<<"REF-END"<<'\t'<<"Q-ST"<<'\t'<<"Q-END"<<endl;
		string tempname;

		for(unsigned int i =0;i<merge.r_name.size();i++)
		{
			tempname = merge.r_name[i] + merge.q_name[i];
			for(unsigned int j =0; j<merge.ref_st[tempname].size();j++)
			{
				fout<<merge.r_name[i]<<'\t'<<merge.q_name[i]<<'\t'<<merge.ref_len[tempname]<<'\t'<<merge.q_len[tempname]<<'\t'<<merge.ref_st[tempname][j]<<'\t'<<merge.ref_end[tempname][j]<<'\t'<<merge.q_st[tempname][j]<<'\t'<<merge.q_end[tempname][j]<<'\t'<<endl;
			}
		}
	fout.close();
	} 
		
			
/////////////////////////////////////////////////////////////////////////this function checks whether a query is completely contained within a reference

void  innieChecker(asmMerge & merge)
{
bool flag =0;
int ref_last = 0;
int q_last = 0;
string tempname;
	for(unsigned int i =0; i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i]+merge.q_name[i];
		ref_last = merge.newEnd[tempname][1];
		q_last = merge.newEnd[tempname][3];
		flag =0; //reset flag 

		if(merge.newEnd[tempname][2] < merge.newEnd[tempname][3])
		{ 	
			if( ((merge.newEnd[tempname][0] -1) >= (merge.newEnd[tempname][2] - 1)) && ((merge.ref_len[tempname] - ref_last)>= (merge.q_len[tempname] - q_last)))
			{ 
				flag = 1;
			}
		}
		if(merge.newEnd[tempname][2] > merge.newEnd[tempname][3])//^ 
		{
			if(((merge.q_len[tempname] - merge.newEnd[tempname][2]) < (merge.newEnd[tempname][0] -1)) && ((merge.ref_len[tempname]-merge.newEnd[tempname][1])>(merge.newEnd[tempname][3] - 1)))
			{									
				flag = 1;
			}
		}
		merge.innie[tempname] = flag;
	}

}			
////////////////////////////////////////////////////////////////////////////////////////////
int abs(int n)
{
int m = n;
	if(n<0)
	{
	m = n* (-1);
	}
return m;
}
/////////////////////////////////////////////////////////////////////////////////////////
void sideCheckerR(asmMerge & merge)  //checks whether the midpoint of alignment is on the right or left side of the reference
{
string tempname;
int rMid;
int r_f,r_last;
int rAlnMid;
int Dist =0;
        for(unsigned int i = 0; i<merge.r_name.size();i++)
        {
                tempname = merge.r_name[i]+merge.q_name[i];
                rMid = int(merge.ref_len[tempname] * 0.5);
               	r_f = merge.newEnd[tempname][0];
                r_last = merge.newEnd[tempname][1];
                Dist = abs(r_last-r_f); // total alignment length

                rAlnMid = r_f + int(Dist *0.5); //midpoint of alignment
                if(rAlnMid > rMid)
                {
                	merge.sideInfoR[tempname] = 'R';
	        }
                if(rAlnMid < rMid)
                {
                	merge.sideInfoR[tempname] = 'L';
                }
                if(rAlnMid == rMid) //when ref is innie
                {
                	merge.sideInfoR[tempname] = 'L';
                }
        }

}
////////////////////////////////////////////////////////////////////////////////////////
void sideChecker(asmMerge & merge) //checks whether midpoint of query is 5' or 3' end of the reference
{
string tempname;
int qMid;
int rMid;
int qMidonR;
	for(unsigned int i = 0; i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i]+merge.q_name[i];
		rMid = merge.ref_len[tempname] * 0.5;
		qMid = merge.q_len[tempname] * 0.5;
		if(merge.newEnd[tempname][2] < merge.newEnd[tempname][3]) //same strand
		{
			if(merge.newEnd[tempname][2] > qMid)
			{	
				qMidonR = merge.newEnd[tempname][0] - (merge.newEnd[tempname][2] -qMid);
				if(qMidonR > rMid)
				{
					merge.sideInfo[tempname] = 'R';
				}
				if(qMidonR < rMid)
				{
					merge.sideInfo[tempname] = 'L';
				}
			}
			if(merge.newEnd[tempname][2] < qMid)
			{
				qMidonR = merge.newEnd[tempname][0] + (qMid - merge.newEnd[tempname][2]);
				if(qMidonR > rMid)
				{
					merge.sideInfo[tempname] = 'R';
				}
				if(qMidonR < rMid)
				{
					merge.sideInfo[tempname] = 'L';
				}
			}
		}
		if(merge.newEnd[tempname][2] > merge.newEnd[tempname][3]) //if opposite strands
		{
			if(qMid < merge.newEnd[tempname][2])
			{
				qMidonR = merge.newEnd[tempname][0] + (merge.newEnd[tempname][2] - qMid);
				if(qMidonR > rMid)
				{
					merge.sideInfo[tempname] = 'R';
				}
				if(qMidonR < rMid)
				{
					merge.sideInfo[tempname] = 'L';
				}
			}
			if(qMid > merge.newEnd[tempname][2])
			{
				qMidonR = merge.newEnd[tempname][0] - (qMid - merge.newEnd[tempname][2]);
				if(qMidonR > rMid)
				{
					merge.sideInfo[tempname] = 'R';
				}
				if(qMidonR < rMid)
				{
					merge.sideInfo[tempname] = 'L';
				}
			}
		}
	}
	
}
////////////////////////////////////////////////////////////////////////////////////////
int ovlCalculator(vector<int>& q_st, vector<int>& q_end)
{
int ovl = 0;

	for(unsigned int j =0;j<q_st.size();j++)
	{
	ovl = ovl + abs(q_st[j] - q_end[j]);
	}

	
return ovl;	
}
/////////////////////////////////////////////////////////////////////////////////////
void ovlStoreCalculator(asmMerge & merge)
{
	string tempname;
	for(unsigned int i =0;i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i]+merge.q_name[i];
		merge.ovlStore[tempname] = ovlCalculator(merge.q_st[tempname],merge.q_end[tempname]);
		lisCalculator(merge,tempname,merge.ref_st[tempname],merge.ref_end[tempname],merge.q_st[tempname],merge.q_end[tempname]);
	}
}
/////////////////////////////////////////////////////////////////////////////////////
void lisCalculator(asmMerge & merge,string & tempname,vector<int>& ref_st,vector<int>& ref_end,vector<int>& q_st, vector<int>& q_end)
{
int ovl = 0;
long int storedQend = 0,storedRefEnd = 0,storedRefSt = 0, storedQst = 0;
vector<int> v;
map<int,int> ovl2refSt;
map<int,int> ovl2qSt;
map<int,int> ovl2refEnd;
map<int,int> ovl2qEnd;

	if(q_st.size() >1)
	{
	        for(unsigned int j =0;j<q_st.size();j++)
	        {
			if(q_st[j]<q_end[j]) //if this MUM is forward oriented
			{
				storedQend = 0;
				storedRefEnd = 0;
				for(unsigned k = j;k<q_st.size();k++)
				{
					if((q_st[k] > storedQst) && (abs(q_st[k]-storedQst)<merge.q_len[tempname]) && (ref_st[k] > storedRefSt))
					{
						ovl = ovl + abs(q_st[k] - q_end[k]);
						storedQend = q_end[k]; //confused? will change this later
						storedQst = q_st[k];
						storedRefEnd = ref_end[k];
						storedRefSt = ref_st[k];
					}
				}
			}
			if(q_st[j]>q_end[j]) //if this MUM is reverse oriented
			{
				storedQend = 100000000000;
				storedQst = 100000000000;
				storedRefEnd = 0;
				for(unsigned k = j;k<q_st.size();k++)
				{
					if(j == k)
					{
						if((q_st[k] < storedQst) && (ref_st[k] > storedRefSt))
					 	{
							ovl = ovl + abs(q_st[k] - q_end[k]);
							storedQend = q_end[k];
							storedQst = q_st[k];
							storedRefEnd = ref_end[k];
							storedRefSt = ref_st[k];
						}
					}
					if(k>j)
					{
						if((q_st[k] < storedQst) && (abs(q_st[k]-storedQst)<merge.q_len[tempname]) && (ref_st[k] > storedRefSt))
						{
							ovl = ovl + abs(q_st[k] - q_end[k]);
                                                        storedQend = q_end[k];
							storedQst = q_st[k];
							storedRefEnd = ref_end[k];
							storedRefSt = ref_st[k];
						}
					}
				}
			}
			ovl2refSt[ovl] = ref_st[j];
			ovl2qSt[ovl] = q_st[j];
			ovl2refEnd[ovl] = storedRefEnd;
			ovl2qEnd[ovl] = storedQend;
			v.push_back(ovl);
			ovl = 0;

		}
	}
	if(q_st.size() == 1)
	{
		ovl = abs(q_st[0] - q_end[0]);
		ovl2refSt[ovl] = ref_st[0];
		ovl2qSt[ovl] = q_st[0];
		ovl2refEnd[ovl] = ref_end[0];
		ovl2qEnd[ovl] = q_end[0];
		v.push_back(abs(q_st[0] - q_end[0]));
	}
		sort(v.begin(),v.end());
	        ovl = v[v.size()-1];
	merge.newEnd[tempname].push_back(ovl2refSt[ovl]); //stores the LIS start and end for each alignment.elements:refSt,refEnd,qSt,qEnd
	merge.newEnd[tempname].push_back(ovl2refEnd[ovl]);
	merge.newEnd[tempname].push_back(ovl2qSt[ovl]);
	merge.newEnd[tempname].push_back(ovl2qEnd[ovl]);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void nOvlStoreCalculator(asmMerge & merge)
{
string tempname;
int ref_last;
int q_last,noLovl,noRovl;
	for(unsigned int i=0; i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i] + merge.q_name[i];
		q_last = merge.newEnd[tempname][3];
		ref_last = merge.newEnd[tempname][1];
		
		if(merge.newEnd[tempname][2] < merge.newEnd[tempname][3])
		{
			if(((merge.newEnd[tempname][2]-1)>(merge.newEnd[tempname][0]-1)) && !((merge.q_len[tempname]-q_last) > (merge.ref_len[tempname] - ref_last)))// if left side of the query end is hanging but right side not
 			{
				noLovl = merge.newEnd[tempname][0];
				noRovl = merge.q_len[tempname] - q_last;
				
			}
			if(((merge.q_len[tempname]-q_last) > (merge.ref_len[tempname] - ref_last)) && !((merge.newEnd[tempname][2]-1)>(merge.newEnd[tempname][0]-1)))// if right side of the query is hanging but left side is not
			{
				noRovl = merge.ref_len[tempname] - ref_last;
				noLovl = merge.newEnd[tempname][2];
				
			}
			if(merge.innie[tempname] == 1)
			{
				noLovl = merge.newEnd[tempname][2];
				noRovl = merge.q_len[tempname] - q_last;
			}
			if(((merge.newEnd[tempname][2]-1)>(merge.newEnd[tempname][0]-1)) && ((merge.q_len[tempname]-q_last) > (merge.ref_len[tempname] - ref_last)))
			{
				noLovl = merge.newEnd[tempname][0];
				noRovl = merge.ref_len[tempname] - ref_last ; 
			}
		}

		if(merge.newEnd[tempname][2] > merge.newEnd[tempname][3])
		{
		
			if(((merge.q_len[tempname] - merge.newEnd[tempname][2]) > (merge.newEnd[tempname][0] - 1)) && !((q_last -1) > (merge.ref_len[tempname] - ref_last)))
			{
				noLovl = merge.newEnd[tempname][0];
				noRovl = q_last -1 ;
				
			}
			
			if(((q_last -1) > (merge.ref_len[tempname] - ref_last)) && !((merge.q_len[tempname] - merge.newEnd[tempname][2]) > (merge.newEnd[tempname][0] - 1)))
			{ 
				noRovl = merge.ref_len[tempname] - ref_last;
				noLovl = merge.q_len[tempname] - merge.newEnd[tempname][2];
				
			}
			if(merge.innie[tempname] == 1)
			{
				noLovl = merge.q_len[tempname] - merge.newEnd[tempname][2];
				noRovl = merge.newEnd[tempname][3] - 1;
			}
			
			if(((merge.q_len[tempname] - merge.newEnd[tempname][2]) > (merge.newEnd[tempname][0] - 1)) && ((q_last -1) > (merge.ref_len[tempname] - ref_last)))
			{
				noLovl = merge.newEnd[tempname][0];
				noRovl = merge.ref_len[tempname] - ref_last;
			}
					
		}
		
		merge.nOvlStore[tempname] = noLovl+noRovl;
		
	}
}
///////////////////////////////////////////////////////////////////////////
void writeSummary(asmMerge & merge)
{
	string tempname;
	double ovl;
	int lastElemRef,lastElemQ;
	ofstream fout;
	fout.open("summaryOut.txt");
	fout<<"REF"<<'\t'<<"QUERY"<<'\t'<<"REF_START"<<'\t'<<"REF_END"<<'\t'<<"Q_START"<<'\t'<<"Q_END"<<'\t'<<"ORIENTATION"<<'\t'<<"INNIE(1/0)"<<'\t'<<"OVERLAP_LEN"<<'\t'<<"OVERLAP_PROP"<<'\t'<<"NO_OVERLAP_AT_ENDS"<<'\t'<<"OVERHANG"<<endl;
	for(unsigned int i =0;i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i]+merge.q_name[i];
		ovl = (merge.ovlStore[tempname])/(double)(merge.nOvlStore[tempname]);
		lastElemRef = merge.ref_end[tempname].size()-1;
		lastElemQ = merge.q_end[tempname].size()-1;
		fout<<merge.r_name[i]<<'\t'<<merge.q_name[i]<<'\t'<<merge.ref_st[tempname][0]<<'\t'<<merge.ref_end[tempname][lastElemRef]<<'\t'<<merge.q_st[tempname][0]<<'\t'<<merge.q_end[tempname][lastElemQ]<<'\t'<<merge.sideInfo[tempname]<<'\t'<<merge.innie[tempname]<<'\t'<<merge.ovlStore[tempname]<<'\t'<<ovl<<'\t'<<merge.nOvlStore[tempname]<<'\t'<<merge.ovrHangQ[tempname]<<endl;
		
	}
	fout.close();
}
////////////////////////////////////////////////////////////////////
string revCom(string & str)
{
string revcom;
unsigned int k = str.size();
	for(int i=k-1; i>-1;--i)
	{
		if(str[i] == 'A')
		{
		revcom.push_back('T');
		}
		if(str[i] == 'T')
		{
		revcom.push_back('A');
		}
		if(str[i] == 'G')
		{
		revcom.push_back('C');
		}
		if(str[i] == 'C')
		{
		revcom.push_back('G');
		}
		if(str[i] == 'a')
                {
                revcom.push_back('T');
                }
                if(str[i] == 't')
                {
                revcom.push_back('A');
                }
                if(str[i] == 'g')
                {
                revcom.push_back('C');
                }
                if(str[i] == 'c')
                {
                revcom.push_back('G');
                }
		if(str[i] == 'N')
		{
		revcom.push_back('N');
		}
	}
return revcom;
}
///////////////////////////////////////////////////////////////////////
void ovrHngCal(asmMerge & merge)
{
int Dist = 0;
int ovrhangQ = 0,ovrhangR = 0;
int q_f,q_last,ref_f,ref_last;	
string tempname;
int ovrhangRTq,ovrhangLTq,Dist1,ovrhangRTr,ovrhangLTr;
	for(unsigned int i=0;i<merge.r_name.size();i++) 
	{
		tempname = merge.r_name[i] + merge.q_name[i];		
		q_f = merge.newEnd[tempname][2];
		q_last = merge.newEnd[tempname][3];
		ref_f = merge.newEnd[tempname][0];
		ref_last = merge.newEnd[tempname][1];

		if(merge.newEnd[tempname][2] < merge.newEnd[tempname][3])
		{
			if((!(((ref_f-1) < (q_f-1)) && ((merge.ref_len[tempname] - ref_last)<(merge.q_len[tempname] - q_last)))) && (merge.innie[tempname]!=1)) //when ref or query not innie
			{
				Dist = abs((q_f - q_last)) + merge.nOvlStore[tempname]; 
				ovrhangQ = merge.q_len[tempname] - Dist; // this would be zero for innie (overhangRef = merge.ref_len[tempname] - Dist) 
				ovrhangR = merge.ref_len[tempname] - Dist;
				if(ovrhangQ == -1) //this is to fix edge cases when alignment edges coincide with sequence ends
                                     {
                                             ovrhangQ = 0 ;
                                     }
                                if(ovrhangR == -1)
                                     {
                                             ovrhangR = 0;
                                     }

			}
			if(((ref_f-1) < (q_f-1)) && ((merge.ref_len[tempname] - ref_last)<(merge.q_len[tempname] - q_last))) // when ref is innie but query is not
			{
				ovrhangQ = -1;
				ovrhangR = -2; //changed from 0 to -2
				Dist = q_last + (merge.ref_len[tempname] - ref_last);
				ovrhangRTq = merge.q_len[tempname] - Dist;
				Dist1 = ref_f;
				ovrhangLTq = q_f - Dist1;
				merge.ovrHangQR[tempname].push_back(ovrhangLTq); // these maps are populated only when these alignments occur
				merge.ovrHangQR[tempname].push_back(ovrhangRTq);
			}
			if(merge.innie[tempname] == 1) // query is innie, ref is not
			{
				ovrhangQ = 0;

				ovrhangR = -1;
				Dist = ref_last + (merge.q_len[tempname] - q_last);
				ovrhangRTr = merge.ref_len[tempname] - Dist;
				ovrhangLTr = ref_f - q_f;
				merge.ovrHangQR[tempname].push_back(ovrhangLTr);
				merge.ovrHangQR[tempname].push_back(ovrhangRTr);
			}
				
		}
		if(merge.newEnd[tempname][2] > merge.newEnd[tempname][3])
		{
			if((!(((ref_f-1) < (merge.q_len[tempname] - q_f)) && ((merge.ref_len[tempname] - ref_last) < (q_last-1)))) && (merge.innie[tempname]!=1))
			{
				Dist = abs((q_f - q_last)) + merge.nOvlStore[tempname]; 
				ovrhangQ = merge.q_len[tempname] - Dist; // this would be zero for innie (ovrhangRef = merge.ref_len[tempname] -Dist)
				ovrhangR = merge.ref_len[tempname] -Dist;
				if(ovrhangQ == -1) //this is to fix edge cases when alignment edges coincides with sequence ends
				{
					ovrhangQ = 0 ;
				}
				if(ovrhangR == -1)
				{
					ovrhangR = 0;
				}
			}
			if(((ref_f-1) < (merge.q_len[tempname] - q_f)) && ((merge.ref_len[tempname] - ref_last) < (q_last-1))) // when ref is innie but query is not
			{
				ovrhangQ = -1;
				ovrhangR = -2; //changed from 0 to -2
				Dist = merge.ref_len[tempname] - ref_last;
				ovrhangRTq = q_last - Dist;
				Dist1 = q_f + merge.newEnd[tempname][0];
				ovrhangLTq = merge.q_len[tempname] - Dist1; 
				merge.ovrHangQR[tempname].push_back(ovrhangLTq); // these maps get populated only when these alignments occur
				merge.ovrHangQR[tempname].push_back(ovrhangRTq);
				
			}
			if(merge.innie[tempname] == 1)
			{
				ovrhangQ = 0;
				ovrhangR = -1;
				Dist = ref_last + (merge.q_len[tempname] - q_last);
				ovrhangRTr = merge.ref_len[tempname] - ref_last;
				Dist1 = merge.q_len[tempname] - q_f;
				ovrhangLTr = ref_f - Dist1;
				merge.ovrHangQR[tempname].push_back(ovrhangLTr);
				merge.ovrHangQR[tempname].push_back(ovrhangRTr);
			}
		}
	merge.ovrHangQ[tempname] = ovrhangQ; // this should be ovrHangQ and add ovrHangRef
	merge.ovrHangR[tempname] = ovrhangR;
	}
}
////////////////////////////////////////////////////////////////////
void fillAnchor(asmMerge & merge,asmMerge & merge1, double propAnchor, double propCutoff, const int & length, const int & absLenCutoff,fastaSeq & hyb)
{
	string tempname;
	double prop;
	for(unsigned int i=0;i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i]+merge.q_name[i];
		
		prop = merge.ovlStore[tempname]/(double)merge.nOvlStore[tempname];
		if((prop>propAnchor) && (merge.ref_len[tempname]>merge.q_len[tempname]) && (merge.ref_len[tempname]>length))
		{
			merge.anchor[tempname].push_back(merge.r_name[i]);
			merge.anchor[tempname].push_back(merge.q_name[i]);
			if(merge.innie[tempname] == 1)
			{
				fillToRemove(merge,merge.q_name[i]);//here add only those for which query is an innie
				
			}
		}
		if((prop > propCutoff) && (merge.ovlStore[tempname]>absLenCutoff))
		{
			merge1.r_name.push_back(merge.r_name[i]); //populating the new asmMerge names vector
			merge1.q_name.push_back(merge.q_name[i]);

		}
	}
}
///////////////////////////////////////////////////////////////////////
void trimSeq(asmMerge & merge, fastaSeq & hyb)
{
	vector<string> qToTrim;
	string tempname,tempSeq;
	for(unsigned int i=0;i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i] + merge.q_name[i];
		if((merge.innie[tempname] == 1) && (find(qToTrim.begin(),qToTrim.end(),merge.q_name[i]) == qToTrim.end()))
                        {
                                qToTrim.push_back(merge.q_name[i]);
                                if(merge.newEnd[tempname][2] > merge.newEnd[tempname][3])
                                {
                                        tempSeq = hyb.seq[merge.q_name[i]].substr(0,merge.newEnd[tempname][3]);
                                        hyb.seq[merge.q_name[i]] = tempSeq;
                                        tempSeq.clear();
                                }
                                if(merge.newEnd[tempname][2] < merge.newEnd[tempname][3])
                                {
                                        tempSeq = hyb.seq[merge.q_name[i]].substr(merge.newEnd[tempname][3]);
                                        hyb.seq[merge.q_name[i]] = tempSeq;
                                        tempSeq.clear();
                                }
                        }
	}
}
//////////////////////////////////////////////////////////////////////////
void fillSeq(fastaSeq & fasta, ifstream& fin, char c)
{
	string str,str1;
	
	while(getline(fin,str))
	{
		if(str[0] == '>')
		{	
			fasta.seqName.push_back(str);
			str1 = str;
			str1[0] = c;
		}
		if(str[0] != '>')
		{
			fasta.seq[str1].append(str);
		}
	}
}
//////////////////////////////////////////////////////////////////////////
void fillSeq(fastaSeq & fasta, ifstream& fin) // overloaded function. will be templated
{
	string str,str1;
	
	while(getline(fin,str))
	{
		if(str[0] == '>')
		{	
			fasta.seqName.push_back(str);
			str1 = str.substr(1);
		}
		if(str[0] != '>')
		{
			fasta.seq[str1].append(str);
		}
	}
}

///////////////////////////////////////////////////////////////////////
void writeAnchorSummary(asmMerge & merge)
{
	ofstream fout;
	fout.open("anchor_summary.txt");
	string tempname,str;
	fout<<"REF_NAME"<<'\t'<<"Q_NAME"<<'\t'<<"REF_LENGTH"<<'\t'<<"Q_LENGTH"<<'\t'<<"REF-ST"<<'\t'<<"REF-END"<<'\t'<<"Q-ST"<<'\t'<<"Q-END"<<endl;
	for( map<string,vector<string> >::iterator it = merge.anchor.begin();it != merge.anchor.end();it++)
	{
		tempname = it->first;
		merge.cAnchor[merge.anchor[tempname][0]].push_back(merge.anchor[tempname][1]);
		fout<<merge.anchor[tempname][0]<<'\t'<<merge.anchor[tempname][1]<<'\t'<<merge.ref_len[tempname]<<'\t'<<merge.q_len[tempname]<<'\t'<<merge.ref_st[tempname][0]<<'\t'<<merge.ref_end[tempname][merge.ref_end[tempname].size()-1]<<'\t'<<merge.q_st[tempname][0]<<'\t'<<merge.q_end[tempname][merge.q_end[tempname].size()-1]<<endl;
		
	}
	fout.close();
}
///////////////////////////////////////////////////////////////////
vector<string> vfind(string tempname,vector<string>& temp_rname, vector<string>& temp_qname,asmMerge & merge,string & guruRef,double propCutoff) // use tempvectors for names for searching. temp_rname is where tempname is. temp_qname is where the alignment partner of r_name is.
{
unsigned int i =0;
vector<string> v1;
string storeRef,storeQ;
char C;
        while(i<temp_rname.size())
        {
               	if(temp_rname[i] == tempname)
        	{
			if (find(merge.q_name.begin(),merge.q_name.end(),tempname) != merge.q_name.end()) //if tempname is a query
			{
				discAnchor(tempname,merge,guruRef,propCutoff);
			}

                	v1.push_back(temp_qname[i]);//may be add the condition that only overlaps above a certain point makes it into the list
			
			storeRef = temp_rname[i]; //storing them before erasing, in case I have to restore them
			storeQ = temp_qname[i];
			temp_rname[i] = ""; //akin to erasing the element. actual erasing is more time consuming.
			temp_qname[i] = "";
			if(merge.cAnchor.find(storeRef) != merge.cAnchor.end()) //if the reference is a contig anchor
			{
				C = 'Y';
			}
			if((find(temp_rname.begin(),temp_rname.end(),storeRef) == temp_rname.end())&& (merge.innie[storeRef+storeQ] ==1) && (C == 'Y')) // if the ref occurs only once and query is an innie, C means whether it is the anchor contig  
			{
				if(find(temp_qname.begin(), temp_qname.end(),storeQ) != temp_qname.end()) // if query occurs more than once
				{
					temp_qname[i] = storeQ;
					temp_rname[i] = storeRef;
				}
			} 
		}
        i++;
	}
        return v1;
}
////////////////////////////////////////////////////////////////
string longestLeft(string tempname, vector<string>& seq,asmMerge & merge, char RorQ, char sideQ,string & prevElem)
{
string str,lStr;
int size,longest;
longest =0;
size = 0;
	for(unsigned int i=0;i<seq.size();i++)
	{
		if(RorQ == 'Q')
		{
			str = tempname+seq[i];
		}
		if(RorQ == 'R')
		{
			str = seq[i]+tempname;
		}
		if(RorQ == 'Q') //if looking for query overhang length 
		{
			if (merge.ovrHangQ[str] == -1) //if ref is innie			
			{
				size = 	merge.ovrHangQR[str][0]; //left overhang is the first element
			}
			if((merge.overHangSideQ[str] == 'L') && (merge.ovrHangQ[str] != -1))
			{
				size = merge.ovrHangQ[str];
			}
		}
               if(RorQ == 'R') //if looking for reference overhang length
		{
			if(merge.ovrHangR[str] == -1) // if query is innie
			{
				size = merge.ovrHangQR[str][0];
			}
			if(merge.ovrHangR[str] == -2)
			{
				merge.innie[seq[i]] = 1; //ref is innie
			}
			if((merge.strandInfo[str] == merge.strandInfo[prevElem]) && (merge.ovrHangR[str] != -1)) //if the references for the query are on the same strand
			{
				if((merge.sideInfoQ[str] != sideQ) || (merge.sideInfo[str] != merge.sideInfo[prevElem]))
				{
					size = merge.ovrHangR[str];
				}
			}
			if((merge.strandInfo[str] != merge.strandInfo[prevElem]) && (merge.ovrHangR[str] != -1))
			{
				if((merge.sideInfoQ[str] != sideQ) ||(merge.sideInfo[str] == merge.sideInfo[prevElem]))
				{
					size = merge.ovrHangR[str];
				}
			}
		}

		if(size>longest)
		{
			longest = size;
			lStr = seq[i];
		}
	}
return lStr;
}
////////////////////////////////////////////////////////////////////////////
string longestRt(string tempname, vector<string>& seq,asmMerge & merge, char RorQ,char sideQ, string & prevElem)
{
	string str,lStr; // l is longest not left 
	int size,longest;
	longest =0;
	size =0;
	for(unsigned int i=0;i<seq.size();i++)
	{
		if(RorQ == 'Q')
		{
			str = tempname+seq[i];
		}
		if(RorQ == 'R')
		{
			str = seq[i]+tempname;
		}

		if(RorQ == 'Q') // if looking for query overhang length
		{
			if(merge.ovrHangQ[str] == -1)
			{
				size = merge.ovrHangQR[str][1]; // second element is the right element
			}
			if((merge.overHangSideQ[str] == 'R') && (merge.ovrHangQ[str] != -1))
			{
				size = merge.ovrHangQ[str];
			}
		
		}
		if(RorQ == 'R') //if looking for reference overhang length
		{
			if(merge.ovrHangR[str] == -1) // if query is innie
			{
				if(merge.strandInfo[str] == merge.strandInfo[prevElem])
				{
					size = merge.ovrHangQR[str][1]; //right overhang
				}
				if((merge.strandInfo[str] != merge.strandInfo[prevElem]) && (merge.sideInfo[str] == 'L'))
				{
					size = merge.ovrHangQR[str][0];//if the previous was right, and this is left, take left overhang
				}
				if((merge.strandInfo[str] != merge.strandInfo[prevElem]) && (merge.sideInfo[str] == 'R'))
				{
					size = merge.ovrHangQR[str][1];
				} 
			}
			if(merge.ovrHangR[str] == -2)
			{
				merge.innie[seq[i]] = 1; //ref is innie
			}
			if((merge.strandInfo[str] == merge.strandInfo[prevElem]) && (merge.ovrHangR[str] != -1)) //if the references for the query are on the same strand
			{
				if((merge.sideInfoQ[str] != sideQ) || (merge.sideInfo[str] != merge.sideInfo[prevElem]))
				{
					size = merge.ovrHangR[str];
				}
			}
			if((merge.strandInfo[str] != merge.strandInfo[prevElem]) && (merge.ovrHangR[str] != -1))
			{
				if((merge.sideInfoQ[str] != sideQ) ||(merge.sideInfo[str] == merge.sideInfo[prevElem]))
				{
					size = merge.ovrHangR[str];
				}
			}
		}
		if(size>longest)
		{
			longest = size;
			lStr = seq[i];
		}
	}
return lStr;
}
///////////////////////////////////////////////////////////////////////////////////
void findChain(asmMerge & merge, asmMerge & merge1,fastaSeq & pbOnly, fastaSeq & merged,double propCutoff)
{
	vector<string> temp_rname;
	vector<string> temp_qname;
	
	temp_rname = merge1.r_name;
	temp_qname = merge1.q_name;


	vector<string> seqs;
	vector<string> refSeqToAdd;
	vector<string> refSeqToRemove;
	string tempname,lRseq,lQseq,rRseq,rQseq,tlRseq,trRseq,tlQseq,trQseq;
	string ref,refForSideL,refForSideR,prevElem,guruR; //guruR stores the name of the anchor contig
	for(map<string,vector<string> >::iterator it= merge.cAnchor.begin();it!=merge.cAnchor.end();it++)
	{
		lRseq = ""; //left reference	
		lQseq = ""; //left query
		rRseq = ""; //right reference
		rQseq = ""; //right query
		tlRseq = "";
		trRseq = "";
		tlQseq = "";
		trQseq = "";
		refForSideL = "";
		refForSideR = "";
		prevElem = "";
		tempname = it->first; // tempname is the reference name


		ref = tempname;
		guruR = tempname;
		seqs = vfind(tempname,temp_rname,temp_qname,merge,guruR,propCutoff);
		fillToRemove(merge,seqs);
		lQseq = longestLeft(tempname,seqs,merge,'Q','N',prevElem);
		rQseq = longestRt(tempname,seqs,merge,'Q','N',prevElem);
		if((lQseq == "") && (rQseq == "")) //for those which are innie
		{
		refSeqToAdd.push_back(tempname);
		}
		if(merge.cAnchor[tempname].size() == 0) //anchors which are partners of previously used query chain elements
		{
			lQseq = "";
			rQseq = "";
		}
		while((rQseq != "")|(rRseq != "")|(lQseq != "")|(lRseq != ""))	
		{
	 
			if(lQseq != "")
			{
				merge.lseq[tempname].push_back(lQseq);
				seqs = vfind(lQseq,temp_qname,temp_rname,merge,guruR,propCutoff); // temp_qname is where lQseq should be present. temp_rname is where the references will be searched from
				
				if(!seqs.empty()) //check that empty returns null when vfind does not resturn an argument
				{
					if(ref != "") //if tempname is currently the reference
					{
						prevElem = ref+lQseq;
						tlRseq = longestLeft(lQseq,seqs,merge,'R',merge.sideInfoQ[tempname+lQseq],prevElem);
						prevElem = "";
					}
					if(ref == "") //if lRseq was the previous reference
					{
						prevElem = refForSideL+lQseq;
						tlRseq = longestLeft(lQseq,seqs,merge,'R',merge.sideInfoQ[refForSideL+lQseq],prevElem); //tlRseq inside the parentheses are from previous alignment
						prevElem = "";
					}
					lRseq = tlRseq;
					refSeqToRemove.push_back(lRseq);
				}
				tlQseq = lQseq; //store the info to get info about the alignment NEW
				lQseq = ""; //reset lQseq
			}
		
			if(lRseq != "")
			{
				merge.lseq[tempname].push_back(lRseq);
				seqs = vfind(lRseq,temp_rname,temp_qname,merge,guruR,propCutoff);
				fillToRemove(merge,seqs);
				if(!seqs.empty())
				{
					if(merge.sideInfo[lRseq+tlQseq] == 'L') //if the previous alignment was on the left side of the reference
					{
						lQseq = longestRt(lRseq,seqs,merge,'Q','N',prevElem); //'N' because this one does not require to check for query side
					}
					if(merge.sideInfo[lRseq+tlQseq] == 'R')
					{
						lQseq = longestLeft(lRseq,seqs,merge,'Q','N',prevElem); //'N' because this one does not require to check for query side
					}				
				}
				refForSideL = lRseq;
				lRseq = "";
				tlQseq = "";//reset tlQseq
			}
			if(rQseq != "")
			{
				merge.rseq[tempname].push_back(rQseq);
				seqs = vfind(rQseq,temp_qname,temp_rname,merge,guruR,propCutoff);
				if(!seqs.empty())
				{
					if(ref != "") //because rRseq is empty when tempname was the previous reference
					{
						prevElem = ref + rQseq;
						trRseq = longestRt(rQseq,seqs,merge,'R',merge.sideInfoQ[tempname+rQseq],prevElem);
						prevElem = "";
					}
					if(ref == "")
					{
						prevElem = refForSideR + rQseq;
						trRseq = longestRt(rQseq,seqs,merge,'R',merge.sideInfoQ[refForSideR+rQseq],prevElem);//refForSideR is the previous reference for this query
						prevElem = "";
					}
					rRseq = trRseq;
					refSeqToRemove.push_back(rRseq);
				}
				ref = "";
				trQseq = rQseq;
				rQseq = "";
			}
			if(rRseq != "")
			{
				merge.rseq[tempname].push_back(rRseq);
				seqs = vfind(rRseq,temp_rname,temp_qname,merge,guruR,propCutoff);
				fillToRemove(merge,seqs);
				if(!seqs.empty())
				{
					if(merge.sideInfo[rRseq+trQseq] == 'R') //if previous alignment was on the right side of the reference
					{
						rQseq = longestLeft(rRseq,seqs,merge,'Q','N',prevElem); // again N because no need to check for query side 
					}
					if(merge.sideInfo[rRseq+trQseq] == 'L') //if previous alignment was on the right side of the reference
					{
						rQseq = longestRt(rRseq,seqs,merge,'Q','N',prevElem); // again N because no need to check for query side NEW
					}
				}
				refForSideR = rRseq;
				rRseq = "";
				trQseq = ""; //reset trQseq
			}
		}
	}

	for(unsigned int i =0 ;i<refSeqToAdd.size();i++)
	{
		if((find(refSeqToRemove.begin(),refSeqToRemove.end(),refSeqToAdd[i]) == refSeqToRemove.end()) && merge.innie[refSeqToAdd[i]] != 1)//ref is not innie
		{
			merged.seq[refSeqToAdd[i]] = pbOnly.seq[refSeqToAdd[i]];
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
void createMseq(asmMerge & merge,asmMerge & merge1) // this function creates the list of contig names (lseq) and the list of their corresponding index
{

	string name,tempname;

	for(map<string,vector<string> >::iterator it = merge.rseq.begin();it!=merge.rseq.end();it++)
	{
		name = it->first;
		merge1.lseq[name].push_back(name); // add the anchor reference as the first element
		for(unsigned int i =0;i<merge.rseq[name].size();i++)
		{	
			if(i==0)
			{
				tempname = name + merge.rseq[name][i];
			}
			if(i%2 ==1 && i>0)
			{
				tempname = merge.rseq[name][i] + merge.rseq[name][i-1];
			}
			if(i%2 == 0 && i>0)
			{
				tempname = merge.rseq[name][i-1] + merge.rseq[name][i]; //ith element is query
			}
			merge1.rseq[name].push_back(tempname);		//fill up rseq with index names
			merge1.lseq[name].push_back(merge.rseq[name][i]); // fill up lseq with contig names
		}
	}

	for(map<string,vector<string> >::iterator it = merge.lseq.begin();it != merge.lseq.end();it++)
	{
		name = it->first;
		for(unsigned int j = 0;j<merge.lseq[name].size();j++)
		{
			if(j == 0)
			{
				tempname = name + merge.lseq[name][j];
			}
			if(j%2 == 1 && j>0)
			{
				tempname = merge.lseq[name][j] + merge.lseq[name][j-1];//jth element is a reference
			}
			if(j%2 == 0 && j>0)
			{
				tempname = merge.lseq[name][j-1] + merge.lseq[name][j];//j-1 th element is a reference
			}
			if(merge.rseq[name].size()!= 0)
			{
				merge1.rseq[name].insert(merge1.rseq[name].begin(),tempname); // insert the left contig indexes
				merge1.lseq[name].insert(merge1.lseq[name].begin(),merge.lseq[name][j]); // insert the left contigs
			}
			if(merge.rseq[name].size() == 0)
			{
				merge1.rseq[name].insert(merge1.rseq[name].begin(),tempname);
				if( j == 0)
				{
					merge1.lseq[name].insert(merge1.lseq[name].begin(),name);
					merge1.lseq[name].insert(merge1.lseq[name].begin(),merge.lseq[name][j]);
				}
				if(j>0)
				{
					merge1.lseq[name].insert(merge1.lseq[name].begin(),merge.lseq[name][j]);
				}
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void fillOri(asmMerge & merge, asmMerge & merge1)
{
	string name;
	int prev;
	for(map<string,vector<string> >::iterator it=merge1.lseq.begin();it != merge1.lseq.end();it++) 
	{
		name = it->first;
		prev =1; //reset prev
		for(unsigned int i =0; i<merge1.lseq[name].size();i++)
		{
			

			if(i==0 && (find(merge.q_name.begin(),merge.q_name.end(),merge1.lseq[name][i])!=merge.q_name.end())) //if first element is query
			{
				if(merge.newEnd[merge1.rseq[name][i]][2] > merge.newEnd[merge1.rseq[name][i]][3])
				{
					merge.Ori[name].push_back(1); // starting query will not be changed, the next reference will be reverse complemented
					prev = prev * (-1);
				}
				if(merge.newEnd[merge1.rseq[name][i]][2] < merge.newEnd[merge1.rseq[name][i]][3])
				{
					merge.Ori[name].push_back(prev);
				}

			}
			if(find(merge.r_name.begin(),merge.r_name.end(),merge1.lseq[name][i]) != merge.r_name.end()) // if it is reference
			{
				if(i == 0)
				{
					merge.Ori[name].push_back(prev); //prev = prev *1 because reference follows orientation of the query
				}
				
				if(i>0 && (merge.newEnd[merge1.rseq[name][i-1]][2] < merge.newEnd[merge1.rseq[name][i-1]][3]))
				{				
					merge.Ori[name].push_back(prev);	
				}
				if(i>0 && (merge.newEnd[merge1.rseq[name][i-1]][2] > merge.newEnd[merge1.rseq[name][i-1]][3]))
				{	
					if((i==1) && (merge.Ori[name][i-1] != -1)) // previous alignment was reverse but was recorded as 1
					{
						merge.Ori[name].push_back(prev);
					}
					else
					{
						prev = prev * (-1);
						merge.Ori[name].push_back(prev);
					}
						
				}
				
			}
			if((i>0) && (find(merge.q_name.begin(),merge.q_name.end(),merge1.lseq[name][i])!=merge.q_name.end()))
			{
				if(merge.newEnd[merge1.rseq[name][i-1]][2] > merge.newEnd[merge1.rseq[name][i-1]][3])
				{	
					prev = prev * (-1);
					merge.Ori[name].push_back(prev);
				}
				if(merge.newEnd[merge1.rseq[name][i-1]][2] < merge.newEnd[merge1.rseq[name][i-1]][3])
				{
					merge.Ori[name].push_back(prev);
				}
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////
void fillToRemove(asmMerge & merge, vector<string> & myvector)
{
	for(unsigned int i=0;i<myvector.size();i++)
	{
		if(merge.qToRemove[myvector[i]] == 0)
		{
			merge.qToRemove[myvector[i]] = 1;
		}

	}
}
void fillToRemove(asmMerge & merge, string & str) //function overloaded
{
	if(merge.qToRemove[str] == 0)
	{
		merge.qToRemove[str] = 1;
	}
}			
////////////////////////////////////////////////////////////////////////////
void removeSeq(asmMerge & merge, fastaSeq & hyb, fastaSeq & merged)
{
	string seqname;
	for(map<string,string>::iterator it = hyb.seq.begin();it!=hyb.seq.end();it++)
	{
		seqname = it->first;
		if(merge.qToRemove[seqname]!= 1)
		{
			merged.seq[seqname] = hyb.seq[seqname];
		}
		
	}
}
////////////////////////////////////////////////////////////////////////////
void writeMerged(fastaSeq & merged, ofstream & fout)
{
	string seqname,index;
	for(map<string,string>::iterator it = merged.seq.begin();it != merged.seq.end();it++)
	{
		seqname = it->first;
		index = seqname;
		if(index[0] == ' ') // for hybrid names
		{
			index[0] = '>';
			fout<<index<<endl;
		}
		else //for reference names
		{
			fout<<">"<<index<<endl;
		}
		fout<<merged.seq[seqname]<<endl;
	}
}	
///////////////////////////////////////////////////////////////////////////
void ctgJoiner(asmMerge & merge,asmMerge & merge1,fastaSeq & hybrid, fastaSeq & pbOnly, fastaSeq & merged)
{
	string name,subseq,subseqR,indexL1,indexL2,seqHolder;
	vector<int> v1;
	int q1_f,q1_last,q2_f,q2_last,r1_f,r1_last,r2_f,r2_last,tempRef_st,cheatF,cheatR;
	int begin_insrt = 1;
	bool skip = false;
	for(map<string,vector<string> >::iterator it = merge1.lseq.begin(); it!=merge1.lseq.end();it++) // merge1.lseq is container for the contig names
	{
		name = it->first;
		begin_insrt = 1;
		tempRef_st = 0;
		for(unsigned int i=0;i<merge1.lseq[name].size();i++)
		{
			q1_f = 0;			
			q1_last = 0;
			q2_f = 0;
			q2_last = 0;
			r1_last =0;
			r1_f = 0;
			r2_f = 0;
			r2_last = 0;
			if( i ==0 && (find(merge.q_name.begin(),merge.q_name.end(),merge1.lseq[name][i])!=merge.q_name.end()))//if first element is query
			{
				indexL2 = merge1.rseq[name][i]; // index name coresponding to the query element
				q2_f = merge.newEnd[indexL2][2];
				q2_last = merge.newEnd[indexL2][3];
			
				if(q2_f>q2_last)
				{
					if(merge.sideInfoR[indexL2] == 'L')
					{
						subseq = hybrid.seq[merge1.lseq[name][i]].substr(q2_f);
						begin_insrt = -1;
					}
					if(merge.sideInfoR[indexL2] == 'R')
					{
						subseq = hybrid.seq[merge1.lseq[name][i]].substr(0,q2_last);
					}	
				}
				if(q2_f < q2_last)
				{
					if(merge.sideInfoR[indexL2] == 'L')
					{					
						subseq = hybrid.seq[merge1.lseq[name][i]].substr(0,q2_f);
					}
					if(merge.sideInfoR[indexL2] == 'R') 
					{	
						subseq = hybrid.seq[merge1.lseq[name][i]].substr(q2_last);
						if(merge.sideInfoQ[indexL2] == 'L')
						{
							begin_insrt = -1;
						}
					}
				}
				if(merge1.lseq[name][0] == merge1.lseq[name][merge1.lseq[name].size()-1]) //if the first and last are same
				{
					subseq = hybrid.seq[merge1.lseq[name][i]];
					skip = true;
				}
			}

			if(i == 0 && (find(merge.r_name.begin(),merge.r_name.end(),merge1.lseq[name][i]) != merge.r_name.end())) // if first element is reference
			{
				indexL2 = merge1.rseq[name][i];
				r2_last = merge.newEnd[indexL2][1];
				r2_f = merge.newEnd[indexL2][0];
				if(merge1.lseq[name][i] == name)
				{
					subseq = pbOnly.seq[merge1.lseq[name][i]].substr(0,r2_last); 
				}
				if(merge1.lseq[name][i] != name)
				{
			
					if((merge.Ori[name][i] == 1) && (merge.Ori[name][i+1] == -1) && (merge.sideInfoQ[indexL2]=='R'))
					
					{
						subseq = pbOnly.seq[merge1.lseq[name][i]].substr(0,r2_last);
					}
					if((merge.Ori[name][i] == 1) && (merge.Ori[name][i+1] == -1) && (merge.sideInfoQ[indexL2]=='U') && (merge.sideInfo[indexL2]=='R'))
					{
						subseq = pbOnly.seq[merge1.lseq[name][i]].substr(0,r2_last);
					}
					
					if(merge.Ori[name][i] == merge.Ori[name][i+1])
					{
						subseq = pbOnly.seq[merge1.lseq[name][i]].substr(0,r2_last);
					}
					if((merge.Ori[name][i] == 1) && (merge.Ori[name][i+1] == -1) && (merge.sideInfoQ[indexL2]=='L'))
                                        {
                                                subseq = pbOnly.seq[merge1.lseq[name][i]].substr(r2_f,(merge.ref_len[indexL2] - r2_f));
                                                begin_insrt = -1;
                                        }
					if((merge.Ori[name][i] == 1) && (merge.Ori[name][i+1] == -1) && (merge.sideInfoQ[indexL2]=='U') && (merge.sideInfo[indexL2]=='L'))
					{
						subseq = pbOnly.seq[merge1.lseq[name][i]].substr(r2_f,(merge.ref_len[indexL2] - r2_f));
						begin_insrt = -1;
					}

				}
			}

			if(i > 0 && (find(merge.q_name.begin(),merge.q_name.end(),merge1.lseq[name][i])!=merge.q_name.end())) // if the element is query
			{
				indexL1 = merge1.rseq[name][i-1];
				q1_f = merge.newEnd[indexL1][2];
				q1_last = merge.newEnd[indexL1][3];
			if((merge.Ori[name][i-1] ==1) && merge.sideInfo[merge1.rseq[name][i-1]] == 'L')
			{
				begin_insrt = -1;
			}
			if(i != (merge1.lseq[name].size()-1)) //if it isn't the last element
			{
				indexL2 = merge1.rseq[name][i];
				q2_f = merge.newEnd[indexL2][2];
				q2_last = merge.newEnd[indexL2][3];
				if(chkOvl(q1_f,q1_last,q2_f,q2_last) == 0) //alignments are non-overlapping
				{
					v1 = minD(q1_f,q1_last,q2_f,q2_last);
					cheatF = min(v1[0],v1[1]);
					cheatR = max(v1[0],v1[1]);
					subseq = hybrid.seq[merge1.lseq[name][i]].substr(cheatF,(cheatR-cheatF));
				}
				if(chkOvl(q1_f,q1_last,q2_f,q2_last) == 1) // if alignments are overlapping
				{
					tempRef_st = mapQonRef(indexL1,indexL2,merge);
					subseq = ""; // may be add a single base if this causes segfault for revcom
				}
			}
			if( i == (merge1.lseq[name].size()-1)) // if this is the last element
			{
			
				if(merge.newEnd[merge1.rseq[name][i-1]][2] > merge.newEnd[merge1.rseq[name][i-1]][3]) //if reverse complement
				{
					subseq = hybrid.seq[merge1.lseq[name][i]].substr(0,q1_last);
					
					if((i>1) && (merge1.lseq[name][i] == merge1.lseq[name][i-2])) //if two queries for the reference are same
					{
						if(merge.sideInfoQ[merge1.rseq[name][i-3]] == 'L')
						{
							subseq = hybrid.seq[merge1.lseq[name][i]].substr(q1_f,(merge.q_len[merge1.rseq[name][i-1]] - q1_f));
						}
						if(merge.sideInfoQ[merge1.rseq[name][i-3]] == 'R')
						{
							subseq = hybrid.seq[merge1.lseq[name][i]].substr(0,q1_last);
						}

					}
					if(merge1.lseq[name][i] == merge1.lseq[name][0]) //if first element and last elements are same
					{
						subseq = hybrid.seq[merge1.lseq[name][i]].substr(0,q1_last);
						
					}
					
					if((merge.overHangSideQ[merge1.rseq[name][i-1]] == 'L') && (i>1) && (!(merge1.lseq[name][i] == merge1.lseq[name][i-2])))
					{
						subseq = hybrid.seq[merge1.lseq[name][i]].substr(q1_f);
					}
					if((merge.overHangSideQ[merge1.rseq[name][i-1]] == 'R') && (i>1) && (!(merge1.lseq[name][i] == merge1.lseq[name][i-2])))
					{
						subseq = hybrid.seq[merge1.lseq[name][i]].substr(0,q1_last);
					}
					
					
				}
				if(merge.newEnd[merge1.rseq[name][i-1]][2] < merge.newEnd[merge1.rseq[name][i-1]][3]) //same strand
				{
					if(merge.sideInfo[merge1.rseq[name][i-1]] == 'L')
					{
						subseq = hybrid.seq[merge1.lseq[name][i]].substr(0,q1_f);
					}
					if(merge.sideInfo[merge1.rseq[name][i-1]] == 'R')
					{
						subseq = hybrid.seq[merge1.lseq[name][i]].substr(q1_last);
					}
					if(merge1.lseq[name][i] == merge1.lseq[name][0]) //if first element and last elements are same
					{
						subseq = hybrid.seq[merge1.lseq[name][i]].substr(q1_last);
					}
					if((i>1) && (merge1.lseq[name][i] == merge1.lseq[name][i-2]))//if two queries for the reference are same
					{
						if(merge.sideInfoQ[merge1.rseq[name][i-3]] == 'L')
						{
							subseq = hybrid.seq[merge1.lseq[name][i]].substr(q1_last);
						}
						if(merge.sideInfoQ[merge1.rseq[name][i-3]] == 'R')
						{
							subseq = hybrid.seq[merge1.lseq[name][i]].substr(0,q1_f);
						}

					}
				}

			}

			if(merge.Ori[name][i] == -1)
			{
				subseqR = revCom(subseq); // if needs to be reverse complemented
				subseq = subseqR;
				
			}
		}


		if(i>0 && (find(merge.r_name.begin(),merge.r_name.end(),merge1.lseq[name][i]) != merge.r_name.end())) // if the element is reference
			{
				indexL1 = merge1.rseq[name][i-1];
				r1_f = merge.newEnd[indexL1][0];
				r1_last = merge.newEnd[indexL1][1];
				

				if(i != (merge1.lseq[name].size()-1)) // if this is not the last element
				{
					indexL2 = merge1.rseq[name][i];
					r2_f = merge.newEnd[indexL2][0];			
					r2_last = merge.newEnd[indexL2][1];
					
					if(tempRef_st != 0)
					{
						r1_f = tempRef_st;
					}
					v1 = maxD(r1_f,r1_last,r2_f,r2_last); // maximum distance between coords from the two alignments
					r1_f = min(v1[0],v1[1]);
					r2_last = max(v1[0],v1[1]);

					subseq = pbOnly.seq[merge1.lseq[name][i]].substr(r1_f,(r2_last - r1_f)); 
					tempRef_st = 0; //reset tempRef_st
					if((merge.sideInfoQ[indexL1] == 'L') && (merge.sideInfoR[indexL1] == 'R')&& (merge.Ori[name][i-1] == 1) && (merge.Ori[name][i] == 1))
					{
						begin_insrt = -1;
					}
					if((merge.sideInfoQ[indexL1] == 'L') && (merge.sideInfoR[indexL1] == 'L')&& (merge.Ori[name][i-1] == 1) && (merge.Ori[name][i] == -1))
					{
						begin_insrt = -1;
					}
				
				}
				if(i == (merge1.lseq[name].size()-1))
				{
					if((merge.Ori[name][i] == -1) && (name == merge1.lseq[name][i]))
					{
						begin_insrt = -1;
					}
				//	if(tempRef_st != 0)
				//	{
		
					//	if((merge.Ori[name][i] == -1) && (merge.Ori[name][i-1] == 1))
					//	{
							if(tempRef_st != 0)
							{
								r1_last = tempRef_st;
								r1_f = tempRef_st; // to use for the second condition
								tempRef_st = 0;
							}
							//if(name != merge1.lseq[name][i])
							if(merge.sideInfoR[indexL1] == 'R')
							{
								subseq = pbOnly.seq[merge1.lseq[name][i]].substr(0,r1_last);

							}
							//if(name == merge1.lseq[name][i])
							if(merge.sideInfoR[indexL1] == 'L')
							{
							
								subseq = pbOnly.seq[merge1.lseq[name][i]].substr(r1_f,(merge.ref_len[indexL1] - r1_f));//r_last is actually r1_f
							}
					//	}
					//	else
					//	{
					//		if(tempRef_st != 0)
					//		{
					//			r1_f = tempRef_st; //transfer the number so that it is used in case of an overlap;
					//			tempRef_st = 0; //reset tempRef
//cout<<name<<"\t"<<r1_f<<endl;
					//		}
					//	}
					//}	
					if(subseq.size() == 0) //if subseq has not been modified by condition 1 in the preceding section
					{
					subseq = pbOnly.seq[merge1.lseq[name][i]].substr(r1_f,(merge.ref_len[indexL1] - r1_f));
					}
				}
				if(merge.Ori[name][i] == -1)
				{
					subseqR = revCom(subseq);
					subseq = subseqR;
				}
			}
			if(begin_insrt == -1)
			{
					subseq.append(seqHolder);
					seqHolder = subseq;
			}

			if(begin_insrt == 1)
			{
				seqHolder.append(subseq);
			}
			subseq.clear(); //experimental. to use the zero size of subseq as a condition 
			if(skip == true)
			{
				i = merge1.lseq[name].size();
				skip = false;
			}
		}
		merged.seq[name] = seqHolder;
cout<<"Merged Contig and length: "<<name<<"\t"<<seqHolder.size()<<endl;
		seqHolder.clear();
		
	}
}

///////////////////////////////////////////////////////////////////////////////
vector<int> maxD(int & qf1,int & qe1, int & qf2, int & qe2)
{
	map<int, vector<int> > Dist;
	int d1,d2,d3,d4;
	vector<int> d;

	d1 = abs(qf1 - qe2);
	Dist[d1].push_back(qf1);
	Dist[d1].push_back(qe2);
	d.push_back(d1);
	
	d2 = abs(qf1 - qf2);
	Dist[d2].push_back(qf1);
	Dist[d2].push_back(qf2);
	d.push_back(d2);

	d3 = abs(qe1 - qf2);
	Dist[d3].push_back(qe1);
	Dist[d3].push_back(qf2);
	d.push_back(d3);

	d4 = abs(qe1 - qe2);
	Dist[d4].push_back(qe1);
	Dist[d4].push_back(qe2);
	d.push_back(d4);
	sort(d.begin(),d.end());
	return Dist[d[3]];
}

/////////////////////////////////////////////////////////////////////////////////
vector<int> minD(int & qf1,int & qe1, int & qf2, int & qe2)
{
	map<int, vector<int> > Dist;
	int d1,d2,d3,d4;
	vector<int> d;

	d1 = abs(qf1 - qe2);
	Dist[d1].push_back(qf1);
	Dist[d1].push_back(qe2);
	d.push_back(d1);
	
	d2 = abs(qf1 - qf2);
	Dist[d2].push_back(qf1);
	Dist[d2].push_back(qf2);
	d.push_back(d2);

	d3 = abs(qe1 - qf2);
	Dist[d3].push_back(qe1);
	Dist[d3].push_back(qf2);
	d.push_back(d3);

	d4 = abs(qe1 - qe2);
	Dist[d4].push_back(qe1);
	Dist[d4].push_back(qe2);
	d.push_back(d4);
	sort(d.begin(),d.end());
	return Dist[d[0]];
}
////////////////////////////////////////////////////////////////////////////////

int maxIntD(int & qf1,int & qe1, int & qf2, int & qe2) 
{
	map<int, vector<int> > Dist;
	int d1,d2,d3,d4;
	vector<int> d;
	d1 = abs(qf1 - qe2);
	Dist[d1].push_back(qf1);
	Dist[d1].push_back(qe2);
	d.push_back(d1);
	
	d2 = abs(qf1 - qf2);
	Dist[d2].push_back(qf1);
	Dist[d2].push_back(qf2);
	d.push_back(d2);

	d3 = abs(qe1 - qf2);
	Dist[d3].push_back(qe1);
	Dist[d3].push_back(qf2);
	d.push_back(d3);

	d4 = abs(qe1 - qe2);
	Dist[d4].push_back(qe1);
	Dist[d4].push_back(qe2);
	d.push_back(d4);
	sort(d.begin(),d.end());
	return abs(Dist[d[3]][0]-Dist[d[3]][1]);
}
///////////////////////////////////////////////////////////////////////////////
int mapQonRef(string & tempname1,string & tempname2,asmMerge & merge)
{
int q1_f,q1_end,q2_f,q2_end,r2_f,r2_end;

int cord =0;

q1_f = merge.newEnd[tempname1][2];
q1_end = merge.newEnd[tempname1][3];
q2_f = merge.newEnd[tempname2][2];
q2_end = merge.newEnd[tempname2][3];
r2_f = merge.newEnd[tempname2][0];
r2_end = merge.newEnd[tempname2][1];

	if(q2_f>q2_end)
	{
		if((q1_end<q2_f) && (q1_end>q2_end)) // q1_end is inside
		{
			cord = r2_f + abs(q2_f-q1_end);
			if(cord > r2_end)
			{
				cord = r2_end - (cord - r2_end); // when cord exceeds r2_end
			}
		}
		if((q1_f<q2_f) && (q1_f>q2_end)) //q1_f is inside
		{
			cord = r2_f + abs(q2_f - q1_f);
			if(cord > r2_end)
			{
				cord = r2_end - (cord-r2_end);
			}
		}
	}
	if(q2_f<q2_end)
	{
		if((q1_end<q2_end) && (q1_end > q2_f)) //q1_end is inside
		{
			cord = r2_end - abs(q2_end - q1_end);
			if(cord <0) // in cases where mummer alignment is noisy
			{
				cord = r2_f + (q1_end - q2_f);          
			}
		}
		if((q1_f<q2_end) && (q1_f > q2_f)) // q1_f is inside		
		{
			cord = r2_end - abs(q2_end - q1_f);
			if(cord < 0) //in cases where mummer alignment is noisy
			{
				cord = r2_f + abs(q2_f - q1_f);
			}
		}
	}

return cord;
}
///////////////////////////////////////////////////////////////////////////////

bool chkOvl(int f1,int e1,int f2,int e2)
{
bool flag = 0;
	if (maxIntD(f1,e1,f2,e2) < (abs(f1-e1) + abs(f2-e2)))
	{
		flag = 1;
	}
return flag;
}
//////////////////////////////////////////////////////////////////////////////
int max(int & a, int & b)
{
 int max;
	if(a>b)
	{
		max = a;
	}
	else
	{
		max =b;
	}
return max;
}
////////////////////////////////////////////////////////////////
int min(int & a, int & b)
{
int min;
	if(a<b)
	{
		min = a;
	}
	else
	{
		min = b;
	}
return min;
}	
//////////////////////////////////////////////////////////////
string reversed(string &str)
{
	string str1;
	for(unsigned int i = 0;i<str.size();i++)
	{
		str1.push_back(str[(str.size()-1)-i]);
	}
	return str1;
}
///////////////////////////////////////////////////////////
int returnIndex (vector<string> & myvector, string & str)
{
	int find = 0;
	for(int i =0;myvector[i]!= str;i++)
	{
		find = i+1;
	}
return find;
}
////////////////////////////////////////////////////////////
void sideCheckerQ(asmMerge & merge)  //checks whether the midpoint of alignment is on the right or left side of the query 
{
string tempname;
int qMid;
int q_f,q_last;
int qAlnMid;
	for(unsigned int i = 0; i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i]+merge.q_name[i];
		qMid = int(merge.q_len[tempname] * 0.5);
		q_f = merge.newEnd[tempname][2];
		q_last = merge.newEnd[tempname][3];
		
		if(q_f < q_last)
		{
			qAlnMid = (q_last * 0.5);
			if(qAlnMid == qMid) //when query is innie
			{
				merge.sideInfoQ[tempname] = 'U';//having 'U' was leading to missing segments in pac to pac merging
			}
			if((merge.q_len[tempname] -q_last) < (q_f -1))
			{
				merge.sideInfoQ[tempname] = 'R';
			}
			if((merge.q_len[tempname] -q_last) > (q_f -1))
			{
				merge.sideInfoQ[tempname] = 'L';
			}
			
		}
		if(q_f > q_last) //reverse strand
		{
			qAlnMid = (q_f *0.5);
			if(qAlnMid == qMid)
			{
				merge.sideInfoQ[tempname] = 'U';
			}
			if((merge.q_len[tempname] - q_f) < (q_last -1))
			{
			 	merge.sideInfoQ[tempname] = 'R';
			}
			if((merge.q_len[tempname] - q_f) > (q_last -1))
			{
				merge.sideInfoQ[tempname] = 'L';
			}
			
		}
	}

}
////////////////////////////////////////////////////////
void assignStrand(asmMerge & merge)
{
	string tempname;
	for(unsigned int i =0;i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i]+merge.q_name[i];
		if(merge.newEnd[tempname][2] > merge.newEnd[tempname][3])
		{
			merge.strandInfo[tempname] = 1;
		}
		if(merge.newEnd[tempname][2] < merge.newEnd[tempname][3])
		{
			merge.strandInfo[tempname] = 0;
		}
	}
}
///////////////////////////////////////////////////////
void discAnchor(string & guruQ, asmMerge & merge, string & guruRef,double propCutoff) //discard anchor contig that is partner of a query contig being used for chaining
{
	string temp,pastTemp;
	double cutoff;
	int count =0;
	for(map<string,vector<string> >::iterator it = merge.cAnchor.begin(); it!= merge.cAnchor.end(); it++)
	{
		count ++;
		temp = it->first;
		pastTemp = temp; //storing temp for use in the next iteration
		temp = temp.append(guruQ);
		
		if(merge.nOvlStore[temp] >0)
		{
			cutoff = merge.ovlStore[temp]/merge.nOvlStore[temp];
		}

		if((merge.ovlStore[temp]>0) && (guruRef != pastTemp) && (cutoff>propCutoff) && (guruRef != it->first))
		{
			merge.cAnchor[it->first].clear();
		}
	
	}
}
////////////////////////////////////////////////////////
void overHangSideR(asmMerge & merge) // whether overhanging query is on the 5' or 3' end of the reference 
{
	string tempname;
	int LnOverlap = 0; // part of query at the 5'side of the reference that does not align to reference
	int RnOverlap = 0;// part of the query at the 3' side of the reference that does not align to the reference
	int LqMaponRef = 0;
	int RqMaponRef = 0;
	int q_end = 0, q_st =0, ref_st = 0, ref_end =0;
	for(unsigned int i = 0; i<merge.r_name.size();i++)
	{
		tempname = merge.r_name[i] + merge.q_name[i];
		q_st = merge.newEnd[tempname][2];
		q_end = merge.newEnd[tempname][3];
		ref_st = merge.newEnd[tempname][0];
		ref_end = merge.newEnd[tempname][1];
		if(q_end < q_st) //if the query is on the other stand
		{
			LnOverlap = merge.q_len[tempname] - q_st;
			RnOverlap = q_end - 1;
			LqMaponRef = ref_st - LnOverlap;
			RqMaponRef = ref_end + RnOverlap;
			if(LqMaponRef < 0)
			{
				merge.overHangSideQ[tempname] = 'L'; //overhang is on the left side of the reference
			}
			if(RqMaponRef > merge.ref_len[tempname])
			{
				merge.overHangSideQ[tempname] = 'R'; // overhang is on the right side of the reference
			}
		}
		if(q_end > q_st)
		{
			LnOverlap = q_st - 1;
			RnOverlap = merge.q_len[tempname] - q_end;
			LqMaponRef = ref_st - LnOverlap;
                        RqMaponRef = ref_end + RnOverlap;
			if(LqMaponRef < 0)
                        {
                                merge.overHangSideQ[tempname] = 'L'; //overhang is on the left side of the reference
                        }
                        if(RqMaponRef > merge.ref_len[tempname])
                        {
                                merge.overHangSideQ[tempname] = 'R'; // overhang is on the right side of the reference
                        }
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////
void splitHaplo(asmMerge & merge,fastaSeq & hybrid)
{
	string hybname;
	ofstream fout;
	fout.open("residuals.fasta");
	for(map<string,string>::iterator it=hybrid.seq.begin();it!=hybrid.seq.end();it++)
	{
		hybname = it->first;
		if(find(merge.q_name.begin(),merge.q_name.end(),hybname) == merge.q_name.end())//if the fasta seq is not present in the delta file
		{
			fout<<">"<<hybname<<endl;
			fout<<it->second<<endl;
			hybrid.seq.erase(it);
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////
				
