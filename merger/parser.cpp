#include<iostream>
#include "qmerge.h"
vector<string> argParser(char *argv[])
{
	vector<string> vs(8);//7 arguments total

	for(unsigned int i = 0; i<17;i++)
	{
		cout<<i<<"\t"<<string(argv[i])<<endl;
		if(string(argv[i]) == "-d")
		{
			vs[0] =string(argv[i+1]); //i+1th element is delta file
		}
		if(string(argv[i]) == "-q")
		{
			vs[1] = string(argv[i+1]);
		}
		if(string(argv[i]) == "-r")
		{
			vs[2] = string(argv[i+1]);
		}
		if(string(argv[i]) == "-hco")
		{
			vs[3] = string(argv[i+1]);
		}
		if(string(argv[i]) == "-c")
		{
			vs[4] = string(argv[i+1]);
		}
		if(string(argv[i]) == "-l")
		{
			vs[5] = string(argv[i+1]);
		}
		if(string(argv[i]) == "-ml")
		{
			vs[6] = string(argv[i+1]);
		}
		if(string(argv[i]) == "-p") //index
		{
			vs[7] = string(argv[i+1]);
		}
	}
	return vs;
}			
