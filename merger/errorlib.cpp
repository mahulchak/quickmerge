#include "qmerge.h"
#include<iostream>
void checkAln(asmMerge & merge, asmMerge & merge1)
{
	size_t pos;
	string name,tempname,refName,qName, temp;
	int rL = 0,rR=0, rs =0, re =0;
	int rL1 =0, rL2=0, rR1 =0, rR2=0, qL1=0, qL2=0,qR1 =0, qR2=0;
	int qs =0,qe =0;
	int ovlLen = 0;
	vector<string> vs;
	double co =0.0;
	for(map<string,vector<string> >::iterator it=merge1.rseq.begin();it != merge1.rseq.end();it++)
	{
		name = it->first;
		for(unsigned int i=0;i<merge1.rseq[name].size();i++)
		{
			//if((rL>200) || (rR>200) || (qL>200) || (qR>200))
			tempname = merge1.rseq[name][i];
			rL = merge.nOvlCord[tempname][1] - merge.nOvlCord[tempname][0];
			rR = merge.nOvlCord[tempname][3] - merge.nOvlCord[tempname][2];
			//qL = abs(merge.nOvlCord[tempname][5] - merge.nOvlCord[tempname][4]);
			//qR = abs(merge.nOvlCord[tempname][7] - merge.nOvlCord[tempname][6]);
			rL1 = merge.nOvlCord[tempname][0]; 
			rL2 = merge.nOvlCord[tempname][1];
			rR1 = merge.nOvlCord[tempname][2];
			rR2 = merge.nOvlCord[tempname][3];
			qL1 = merge.nOvlCord[tempname][4];
			qL2 = merge.nOvlCord[tempname][5];
			qR1 = merge.nOvlCord[tempname][6];
			qR2 = merge.nOvlCord[tempname][7];
			pos = tempname.find(' ');
			refName = tempname.substr(0,pos);
			qName = tempname.substr(pos);
			map<string,vector<string> >::iterator it1 = merge.rseq.begin();
//cout<<name<<" "<<refName<<" "<<qName<<"\t"<<rL<<"\t"<<rR<<"\t"<<"merge r seq has "<<merge.rseq[it->first].size()<<" "<<it->first<<endl;
			if((rL>200) || (rR>200))
			{
			//	cout<<name<<" "<<refName<<" "<<qName<<"\t"<<rL<<"\t"<<rR<<"\t"<<"merge r seq has "<<merge.rseq[it->first].size()<<" "<<it->first<<endl;
				for(unsigned int k=0;k<merge1.q_name.size();k++)
				{
					temp = refName + merge1.q_name[k];
					if(merge.newEnd[temp].size() != 0) //alignment exists
					{
						rs = merge.newEnd[temp][0];
						re = merge.newEnd[temp][1];
						qs = min(merge.newEnd[temp][2],merge.newEnd[temp][3]);
						qe = max(merge.newEnd[temp][2],merge.newEnd[temp][3]);
//cout<<name<<"\t"<<temp<<endl;
						if((!(re <rL1) && !(re>rL2)) && (find(vs.begin(),vs.end(),temp) == vs.end()))
						{
							ovlLen = min(rL2,re) - max(rL1,rs);
							co = double(ovlLen)/double(rL);
							if(co >0.9)
							{
//								cout<<"R L "<<temp<<" "<<merge.newEnd[temp][0]<<" "<<merge.newEnd[temp][1]<<" "<<rL1<<" "<<rL2<<" "<<rL<<" "<<ovlLen<<" "<<merge.innie[temp]<<" "<<co<<endl;
								vs.push_back(temp);
							}
						}
						if((!(rs<rR1) && !(rs>rR2)) && (find(vs.begin(),vs.end(),temp) == vs.end()))
						{
							ovlLen = min(rR2,re) - max(rR1,rs);
							co = double(ovlLen)/double(rR);
							if(co >0.9)
							{
//								cout<<"R R "<<temp<<" "<<merge.newEnd[temp][0]<<" "<<merge.newEnd[temp][1]<<" "<<rR1<<" "<<rR2<<" "<<rR<<" "<<ovlLen<<" "<<merge.innie[temp]<<" "<<co<<endl;
								vs.push_back(temp);
							}
						}
						
					}
					
				}
				vs.clear();
				for(unsigned int k= 0; k<merge.r_name.size();k++)
				{
					temp = merge.r_name[k] +qName;
					if(merge.newEnd[temp].size() != 0)
					{
//cout<<temp<<endl;
						qs = min(merge.newEnd[temp][2],merge.newEnd[temp][3]);
						qe = max(merge.newEnd[temp][2],merge.newEnd[temp][3]);
					//	qs = merge.newEnd[temp][2];
					//	qe = merge.newEnd[temp][3];
						if(qL1 <qL2)  //ovl is forward
						{
//cout<<temp<<endl;
							if((!(qe <qL1) && !(qe>qL2)) && (find(vs.begin(),vs.end(),temp) == vs.end()))
							{
								ovlLen = min(qL2,qe) - max(qL1,qs);
								co = double(ovlLen)/double(rL);
								if(co>0.9)
								{
//									cout<<"qF L "<<temp<<" "<<merge.newEnd[temp][2]<<" "<<merge.newEnd[temp][3]<<" "<<qL1<<" "<<qL2<<" "<<rL<<" "<<ovlLen<<" "<< co<<endl;
									vs.push_back(temp);
								}
							}
							if((!(qs<qR1) && !(qs>qR2)) && (find(vs.begin(),vs.end(),temp) == vs.end()))
							{
								ovlLen = min(qe,qR2) - max(qs,qR1);
								co = double(ovlLen)/double(rR);
								if(co >0.9)
								{
  //                                                              	cout<<"qF R "<<temp<<" "<<merge.newEnd[temp][2]<<" "<<merge.newEnd[temp][3]<<" "<<qR1<<" "<<qR2<<" "<<rR<<" "<<ovlLen<<" "<<co<<endl;
									vs.push_back(temp);	
								}
							}
								
						}
						if(qL1 > qL2) //novl is reverse but ovl is forward
						{
							if((!(qs>qL1) && !(qs<qL2)) && (find(vs.begin(),vs.end(),temp) == vs.end()))
							{
								ovlLen = min(qe,qL1) - max(qs,qL2);
								co = double(ovlLen)/double(rL);
								if(co>0.9)
								{
//									cout<<"qR L "<<temp<<" "<<merge.newEnd[temp][2]<<" "<<merge.newEnd[temp][3]<<" "<<qL1<<" "<<qL2<<" "<<rL<<" "<<ovlLen<<" "<<co<<endl;
									vs.push_back(temp);	
								}
							}
							if((!(qe>qR1) && !(qe<qR2)) && (find(vs.begin(),vs.end(),temp) == vs.end()))
							{
								ovlLen = min(qe,qR1) - max(qs,qR2);
								co = double(ovlLen)/double(rR);
								if(co >0.9)
								{
//									cout<<"qR R "<<temp<<" "<<merge.newEnd[temp][2]<<" "<<merge.newEnd[temp][3]<<" "<<qR1<<" "<<qR2<<" "<<rR<<" "<<ovlLen<<" "<<co<<endl;
									vs.push_back(temp);
								}
							}
						}
					}
				}
				
			}
		}
	}
}


		
